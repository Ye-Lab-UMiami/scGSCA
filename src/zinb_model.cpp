// ============================================================================
// zinb_model.cpp
//
// Rewritten ZINB fitting core for the ctfix branch.
//
// Goals of this rewrite:
//   1. Keep the model interface used by the rest of MetaGSCA unchanged.
//   2. Preserve the DC formula structure used in R:
//        y ~ ct_factor + (1 | subj_factor), zi ~ ct_factor
//   3. Use a cleaner scalar Laplace solver for each subject random effect.
//   4. Make optimization failures easier to diagnose without model fallback.
// ============================================================================
#include "MetaGSCA.h"

#include <nlopt.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace gsca {

namespace {

constexpr double kTiny = 1e-12;
constexpr double kProbEps = 1e-10;
constexpr double kEtaBound = 20.0;
constexpr double kGradTol = 1e-6;
constexpr int kModeMaxIter = 80;
constexpr double kMaxModeAbs = 40.0;

static inline double clamp_prob(double p) {
    return std::clamp(p, kProbEps, 1.0 - kProbEps);
}

static inline double expit(double x) {
    if (x >= 35.0) return 1.0 - kProbEps;
    if (x <= -35.0) return kProbEps;
    double ex = std::exp(x);
    return clamp_prob(ex / (1.0 + ex));
}

static inline double logit(double p) {
    p = clamp_prob(p);
    return std::log(p / (1.0 - p));
}

static inline double safe_mu_from_eta(double eta) {
    return std::max(std::exp(std::clamp(eta, -kEtaBound, kEtaBound)), kTiny);
}

static double dnbinom_log(int y, double mu, double theta) {
    mu = std::max(mu, kTiny);
    theta = std::max(theta, kTiny);
    return std::lgamma(y + theta) - std::lgamma(theta) - std::lgamma(y + 1.0)
         + theta * std::log(theta / (mu + theta))
         + y * std::log(mu / (mu + theta));
}

static int rnbinom(double mu, double theta, RNG& rng) {
    mu = std::max(mu, kTiny);
    theta = std::max(theta, kTiny);
    std::gamma_distribution<double> gamma(theta, mu / theta);
    double lambda = gamma(rng.engine());
    if (lambda <= kTiny) return 0;
    std::poisson_distribution<int> pois(lambda);
    return pois(rng.engine());
}

struct ZINBDataGrouped {
    int n = 0;
    int n_subj = 0;
    int n_ct = 0;
    bool dc_mode = false;
    std::vector<std::vector<int>> subj_obs;
    std::vector<int> y;
    std::vector<int> ct_ids;
};

static int param_count(bool dc_mode, int n_ct) {
    return dc_mode ? (2 * n_ct + 2) : 4;
}

static inline double cond_eta_for_ct(const Vec& beta_cond, int ct_id) {
    double eta = beta_cond[0];
    if (ct_id > 0 && ct_id < beta_cond.size()) eta += beta_cond[ct_id];
    return eta;
}

static inline double zi_eta_for_ct(const Vec& beta_zi, int ct_id) {
    double eta = beta_zi[0];
    if (ct_id > 0 && ct_id < beta_zi.size()) eta += beta_zi[ct_id];
    return eta;
}

struct ObsDerivs {
    double ll = 0.0;
    double dll = 0.0;
    double d2ll = 0.0;
};

static ObsDerivs obs_loglik_derivs(int yi, double mu, double theta, double pi_i) {
    ObsDerivs out;
    mu = std::max(mu, kTiny);
    theta = std::max(theta, kTiny);
    pi_i = clamp_prob(pi_i);

    if (yi == 0) {
        double nb0 = std::exp(theta * std::log(theta / (mu + theta)));
        double L = pi_i + (1.0 - pi_i) * nb0;
        L = std::max(L, 1e-300);
        out.ll = std::log(L);

        double dnb0_dmu = -theta * nb0 / (mu + theta);
        double dnb0_db = dnb0_dmu * mu;
        double dL_db = (1.0 - pi_i) * dnb0_db;
        out.dll = dL_db / L;

        double d2nb0_dmu2 = nb0 * theta / ((mu + theta) * (mu + theta))
                          * (theta * theta / (mu + theta) - theta + 1.0);
        double d2nb0_db2 = (d2nb0_dmu2 * mu + dnb0_dmu) * mu;
        double d2L_db2 = (1.0 - pi_i) * d2nb0_db2;
        out.d2ll = (d2L_db2 * L - dL_db * dL_db) / (L * L);
    } else {
        out.ll = std::log(std::max(1.0 - pi_i, 1e-300)) + dnbinom_log(yi, mu, theta);

        double dll_dmu = static_cast<double>(yi) / mu
                       - (static_cast<double>(yi) + theta) / (mu + theta);
        out.dll = dll_dmu * mu;

        double d2ll_dmu2 = -static_cast<double>(yi) / (mu * mu)
                         + (static_cast<double>(yi) + theta) / ((mu + theta) * (mu + theta));
        out.d2ll = (d2ll_dmu2 * mu + dll_dmu) * mu;
    }

    return out;
}

struct SubjectPosteriorEval {
    double hval = 0.0;
    double grad = 0.0;
    double hess = 0.0;
};

struct SubjectLaplaceDiag {
    bool newton_warn = false;
    bool hessian_clamped = false;
    double abs_grad = 0.0;
    double neg_h = 0.0;
};

struct LaplaceGrouped {
    const ZINBDataGrouped* dat = nullptr;
    const double* params = nullptr;

    void unpack(Vec& beta_cond,
                double& log_theta,
                double& log_sigma,
                Vec& beta_zi) const {
        if (dat->dc_mode) {
            beta_cond.resize(dat->n_ct);
            beta_zi.resize(dat->n_ct);
            beta_cond.setZero();
            beta_zi.setZero();

            beta_cond[0] = params[0];
            for (int ct = 1; ct < dat->n_ct; ++ct) beta_cond[ct] = params[ct];

            log_theta = params[dat->n_ct];
            log_sigma = params[dat->n_ct + 1];

            beta_zi[0] = params[dat->n_ct + 2];
            for (int ct = 1; ct < dat->n_ct; ++ct) {
                beta_zi[ct] = params[dat->n_ct + 2 + ct];
            }
        } else {
            beta_cond.resize(1);
            beta_zi.resize(1);
            beta_cond[0] = params[0];
            log_theta = params[1];
            log_sigma = params[2];
            beta_zi[0] = params[3];
        }
    }

    SubjectPosteriorEval eval_subject(int j, double b,
                                      const Vec& beta_cond,
                                      double theta,
                                      double sigma,
                                      const Vec& beta_zi) const {
        const auto& obs = dat->subj_obs[j];
        double inv_s2 = 1.0 / std::max(sigma * sigma, kTiny);

        SubjectPosteriorEval out;
        out.hval = -0.5 * b * b * inv_s2
                 - std::log(std::max(sigma, kTiny))
                 - 0.5 * std::log(2.0 * M_PI);
        out.grad = -b * inv_s2;
        out.hess = -inv_s2;

        for (int idx : obs) {
            double eta = cond_eta_for_ct(beta_cond, dat->ct_ids[idx]) + b;
            double mu = safe_mu_from_eta(eta);
            double pi_i = expit(zi_eta_for_ct(beta_zi, dat->ct_ids[idx]));
            ObsDerivs d = obs_loglik_derivs(dat->y[idx], mu, theta, pi_i);
            out.hval += d.ll;
            out.grad += d.dll;
            out.hess += d.d2ll;
        }
        return out;
    }

    bool bracket_mode(int j,
                      const Vec& beta_cond,
                      double theta,
                      double sigma,
                      const Vec& beta_zi,
                      double& lo,
                      double& hi,
                      SubjectPosteriorEval& lo_eval,
                      SubjectPosteriorEval& hi_eval) const {
        SubjectPosteriorEval center = eval_subject(j, 0.0, beta_cond, theta, sigma, beta_zi);
        if (!std::isfinite(center.grad)) return false;

        if (std::abs(center.grad) < kGradTol) {
            lo = hi = 0.0;
            lo_eval = hi_eval = center;
            return true;
        }

        double step = 1.0;
        if (center.grad > 0.0) {
            lo = 0.0;
            lo_eval = center;
            hi = step;
            hi_eval = eval_subject(j, hi, beta_cond, theta, sigma, beta_zi);
            while ((hi_eval.grad > 0.0 || !std::isfinite(hi_eval.grad)) && hi < kMaxModeAbs) {
                step *= 2.0;
                hi = std::min(hi + step, kMaxModeAbs);
                hi_eval = eval_subject(j, hi, beta_cond, theta, sigma, beta_zi);
            }
            return std::isfinite(hi_eval.grad) && hi_eval.grad <= 0.0;
        }

        hi = 0.0;
        hi_eval = center;
        lo = -step;
        lo_eval = eval_subject(j, lo, beta_cond, theta, sigma, beta_zi);
        while ((lo_eval.grad < 0.0 || !std::isfinite(lo_eval.grad)) && std::abs(lo) < kMaxModeAbs) {
            step *= 2.0;
            lo = std::max(lo - step, -kMaxModeAbs);
            lo_eval = eval_subject(j, lo, beta_cond, theta, sigma, beta_zi);
        }
        return std::isfinite(lo_eval.grad) && lo_eval.grad >= 0.0;
    }

    SubjectPosteriorEval solve_mode(int j,
                                    const Vec& beta_cond,
                                    double theta,
                                    double sigma,
                                    const Vec& beta_zi,
                                    SubjectLaplaceDiag* diag) const {
        double lo = 0.0, hi = 0.0;
        SubjectPosteriorEval lo_eval, hi_eval;
        if (!bracket_mode(j, beta_cond, theta, sigma, beta_zi, lo, hi, lo_eval, hi_eval)) {
            if (diag) diag->newton_warn = true;
            SubjectPosteriorEval edge = eval_subject(j, 0.0, beta_cond, theta, sigma, beta_zi);
            if (diag) {
                diag->abs_grad = std::abs(edge.grad);
                diag->neg_h = -edge.hess;
            }
            return edge;
        }

        if (lo == hi) {
            if (diag) {
                diag->abs_grad = std::abs(lo_eval.grad);
                diag->neg_h = -lo_eval.hess;
            }
            return lo_eval;
        }

        double b = 0.5 * (lo + hi);
        SubjectPosteriorEval cur = eval_subject(j, b, beta_cond, theta, sigma, beta_zi);

        for (int iter = 0; iter < kModeMaxIter; ++iter) {
            if (!std::isfinite(cur.grad) || !std::isfinite(cur.hess) || !std::isfinite(cur.hval)) {
                if (diag) diag->newton_warn = true;
                break;
            }
            if (std::abs(cur.grad) < kGradTol) break;

            if (cur.grad > 0.0) {
                lo = b;
                lo_eval = cur;
            } else {
                hi = b;
                hi_eval = cur;
            }

            double candidate = std::numeric_limits<double>::quiet_NaN();
            if (cur.hess < -1e-8) candidate = b - cur.grad / cur.hess;
            if (!std::isfinite(candidate) || candidate <= lo || candidate >= hi) {
                candidate = 0.5 * (lo + hi);
            }

            if (std::abs(candidate - b) < 1e-8 * (1.0 + std::abs(b))) {
                b = candidate;
                cur = eval_subject(j, b, beta_cond, theta, sigma, beta_zi);
                break;
            }

            b = candidate;
            cur = eval_subject(j, b, beta_cond, theta, sigma, beta_zi);
        }

        if (diag) {
            diag->abs_grad = std::abs(cur.grad);
            diag->neg_h = -cur.hess;
            if (std::abs(cur.grad) >= 1e-4) diag->newton_warn = true;
        }
        return cur;
    }

    double laplace_subject(int j,
                           const Vec& beta_cond,
                           double theta,
                           double sigma,
                           const Vec& beta_zi,
                           SubjectLaplaceDiag* diag = nullptr) const {
        const auto& obs = dat->subj_obs[j];
        if (obs.empty()) return 0.0;

        SubjectPosteriorEval mode_eval = solve_mode(j, beta_cond, theta, sigma, beta_zi, diag);
        double neg_h = -mode_eval.hess;
        if (!std::isfinite(neg_h) || neg_h < 1e-10) {
            neg_h = 1e-10;
            if (diag) diag->hessian_clamped = true;
        }
        if (diag) diag->neg_h = neg_h;
        return mode_eval.hval + 0.5 * std::log(2.0 * M_PI) - 0.5 * std::log(neg_h);
    }

    double marginal_loglik() const {
        Vec beta_cond, beta_zi;
        double log_theta = 0.0, log_sigma = 0.0;
        unpack(beta_cond, log_theta, log_sigma, beta_zi);

        if (!std::isfinite(log_theta) || !std::isfinite(log_sigma)) return -std::numeric_limits<double>::infinity();

        double theta = std::exp(log_theta);
        double sigma = std::exp(log_sigma);
        if (!std::isfinite(theta) || !std::isfinite(sigma) || theta <= 0.0 || sigma <= 0.0) {
            return -std::numeric_limits<double>::infinity();
        }

        double ll = 0.0;
        for (int j = 0; j < dat->n_subj; ++j) {
            double lj = laplace_subject(j, beta_cond, theta, sigma, beta_zi, nullptr);
            if (!std::isfinite(lj)) return -std::numeric_limits<double>::infinity();
            ll += lj;
        }
        return ll;
    }
};

static double nlopt_obj_grouped(unsigned /* n */, const double* x,
                                double* /* grad */, void* data) {
    auto* obj = static_cast<LaplaceGrouped*>(data);
    obj->params = x;
    double ll = obj->marginal_loglik();
    if (!std::isfinite(ll)) return std::numeric_limits<double>::infinity();
    return -ll;
}

static std::vector<double> make_initial_values(const ZINBDataGrouped& dat,
                                               const std::vector<int>& subj_ids,
                                               int np) {
    std::vector<double> x0(np, 0.0);

    double y_mean = 0.0;
    double y_var = 0.0;
    int zero_count = 0;
    for (int yi : dat.y) y_mean += yi;
    y_mean /= std::max(1, dat.n);
    y_mean = std::max(y_mean, 0.05);

    for (int yi : dat.y) {
        double d = static_cast<double>(yi) - y_mean;
        y_var += d * d;
        if (yi == 0) ++zero_count;
    }
    y_var /= std::max(1, dat.n - 1);

    double theta_init = 1.0;
    if (y_var > y_mean + 1e-8) theta_init = (y_mean * y_mean) / (y_var - y_mean);
    theta_init = std::clamp(theta_init, 0.05, 100.0);

    std::vector<double> subj_sum(dat.n_subj, 0.0);
    std::vector<int> subj_n(dat.n_subj, 0);
    for (int i = 0; i < dat.n; ++i) {
        subj_sum[subj_ids[i]] += static_cast<double>(dat.y[i]);
        subj_n[subj_ids[i]] += 1;
    }
    double subj_log_mean = 0.0;
    int subj_used = 0;
    std::vector<double> subj_logs(dat.n_subj, 0.0);
    for (int j = 0; j < dat.n_subj; ++j) {
        if (subj_n[j] == 0) continue;
        subj_logs[j] = std::log(std::max(subj_sum[j] / subj_n[j], 0.05));
        subj_log_mean += subj_logs[j];
        ++subj_used;
    }
    subj_log_mean /= std::max(1, subj_used);
    double subj_log_var = 0.0;
    for (int j = 0; j < dat.n_subj; ++j) {
        if (subj_n[j] == 0) continue;
        double d = subj_logs[j] - subj_log_mean;
        subj_log_var += d * d;
    }
    subj_log_var /= std::max(1, subj_used - 1);
    double sigma_init = std::clamp(std::sqrt(std::max(subj_log_var, 1e-4)), 0.05, 2.0);
    double zero_frac_all = static_cast<double>(zero_count) / std::max(1, dat.n);

    if (dat.dc_mode) {
        std::vector<double> ct_sum(dat.n_ct, 0.0);
        std::vector<int> ct_n(dat.n_ct, 0);
        std::vector<int> ct_zero(dat.n_ct, 0);

        for (int i = 0; i < dat.n; ++i) {
            int ct = dat.ct_ids[i];
            ct_sum[ct] += dat.y[i];
            ct_n[ct] += 1;
            if (dat.y[i] == 0) ++ct_zero[ct];
        }

        auto safe_log_mean = [&](int ct) {
            double mean = ct_sum[ct] / std::max(1, ct_n[ct]);
            double zero_frac = static_cast<double>(ct_zero[ct]) / std::max(1, ct_n[ct]);
            double nonzero_frac = std::max(1.0 - zero_frac, 0.02);
            double latent_mean = mean / nonzero_frac;
            return std::log(std::max(latent_mean, 0.05));
        };
        auto safe_zero_logit = [&](int ct) {
            double frac = static_cast<double>(ct_zero[ct]) / std::max(1, ct_n[ct]);
            return logit(std::clamp(frac, 0.02, 0.98));
        };

        x0[0] = safe_log_mean(0);
        for (int ct = 1; ct < dat.n_ct; ++ct) x0[ct] = safe_log_mean(ct) - x0[0];

        x0[dat.n_ct] = std::log(theta_init);
        x0[dat.n_ct + 1] = std::log(sigma_init);

        double base_zi = safe_zero_logit(0);
        x0[dat.n_ct + 2] = base_zi;
        for (int ct = 1; ct < dat.n_ct; ++ct) {
            x0[dat.n_ct + 2 + ct] = safe_zero_logit(ct) - base_zi;
        }
    } else {
        x0[0] = std::log(y_mean);
        x0[1] = std::log(theta_init);
        x0[2] = std::log(sigma_init);
        x0[3] = logit(std::clamp(zero_frac_all, 0.02, 0.98));
    }

    return x0;
}

}  // namespace

ZINBFit zinb_fit(const std::string& model_type,
                 const std::vector<int>& y,
                 const std::vector<int>& subj_ids,
                 const std::vector<int>& ct_ids,
                 int n_subj,
                 int n_ct) {
    ZINBFit result;
    bool dc_mode = (model_type == "dc");

    ZINBDataGrouped dat;
    dat.n = static_cast<int>(y.size());
    dat.n_subj = n_subj;
    dat.n_ct = dc_mode ? n_ct : 1;
    dat.dc_mode = dc_mode;
    dat.y = y;
    dat.ct_ids = dc_mode ? ct_ids : std::vector<int>(y.size(), 0);
    dat.subj_obs.resize(n_subj);
    for (int i = 0; i < dat.n; ++i) dat.subj_obs[subj_ids[i]].push_back(i);

    int np = param_count(dc_mode, dat.n_ct);
    std::vector<double> x0 = make_initial_values(dat, subj_ids, np);

    std::vector<double> lb(np, -12.0), ub(np, 12.0), step(np, dc_mode ? 0.35 : 0.25);
    if (dc_mode) {
        lb[dat.n_ct] = -6.0;
        ub[dat.n_ct] = 6.0;
        lb[dat.n_ct + 1] = -6.0;
        ub[dat.n_ct + 1] = 4.0;
        step[dat.n_ct] = 0.3;
        step[dat.n_ct + 1] = 0.2;
        for (int i = dat.n_ct + 2; i < np; ++i) step[i] = 0.45;
    } else {
        lb[1] = -6.0;
        ub[1] = 6.0;
        lb[2] = -6.0;
        ub[2] = 4.0;
        step[1] = 0.3;
        step[2] = 0.2;
        step[3] = 0.45;
    }

    for (int i = 0; i < np; ++i) {
        if (!std::isfinite(x0[i])) {
            x0[i] = 0.5 * (lb[i] + ub[i]);
            ++result.initial_clamped;
        } else if (x0[i] < lb[i]) {
            x0[i] = lb[i];
            ++result.initial_clamped;
        } else if (x0[i] > ub[i]) {
            x0[i] = ub[i];
            ++result.initial_clamped;
        }
    }

    LaplaceGrouped obj;
    obj.dat = &dat;

    nlopt_opt opt = nlopt_create(NLOPT_LN_BOBYQA, np);
    nlopt_set_lower_bounds(opt, lb.data());
    nlopt_set_upper_bounds(opt, ub.data());
    nlopt_set_initial_step(opt, step.data());
    nlopt_set_min_objective(opt, nlopt_obj_grouped, &obj);
    nlopt_set_xtol_rel(opt, dc_mode ? 1e-6 : 1e-5);
    nlopt_set_ftol_rel(opt, dc_mode ? 1e-8 : 1e-6);
    nlopt_set_maxeval(opt, dc_mode ? std::max(5000, 250 * np) : 1200);

    double minf = std::numeric_limits<double>::infinity();
    nlopt_result ret = nlopt_optimize(opt, x0.data(), &minf);
    nlopt_destroy(opt);

    result.converged = (ret > 0) && std::isfinite(minf);
    result.nlopt_status = static_cast<int>(ret);
    result.objective = minf;

    if (dc_mode) {
        result.beta_cond.resize(dat.n_ct);
        result.beta_zi.resize(dat.n_ct);
        result.beta_cond.setZero();
        result.beta_zi.setZero();
        result.beta_cond[0] = x0[0];
        for (int ct = 1; ct < dat.n_ct; ++ct) result.beta_cond[ct] = x0[ct];
        result.log_theta = x0[dat.n_ct];
        result.sigma_re = std::exp(x0[dat.n_ct + 1]);
        result.beta_zi[0] = x0[dat.n_ct + 2];
        for (int ct = 1; ct < dat.n_ct; ++ct) result.beta_zi[ct] = x0[dat.n_ct + 2 + ct];
    } else {
        result.beta_cond.resize(1);
        result.beta_zi.resize(1);
        result.beta_cond[0] = x0[0];
        result.log_theta = x0[1];
        result.sigma_re = std::exp(x0[2]);
        result.beta_zi[0] = x0[3];
    }

    for (int i = 0; i < np; ++i) {
        if (std::abs(x0[i] - lb[i]) < 1e-5 || std::abs(x0[i] - ub[i]) < 1e-5) {
            ++result.boundary_hits;
        }
    }

    Vec beta_cond, beta_zi;
    double log_theta = 0.0, log_sigma = 0.0;
    obj.params = x0.data();
    obj.unpack(beta_cond, log_theta, log_sigma, beta_zi);
    double theta = std::exp(log_theta);
    double sigma = std::max(std::exp(log_sigma), 1e-8);

    result.min_neg_h = std::numeric_limits<double>::infinity();
    result.max_abs_grad = 0.0;
    for (int j = 0; j < dat.n_subj; ++j) {
        SubjectLaplaceDiag diag;
        obj.laplace_subject(j, beta_cond, theta, sigma, beta_zi, &diag);
        if (diag.newton_warn) ++result.newton_warn_subjects;
        if (diag.hessian_clamped) ++result.hessian_clamp_subjects;
        result.max_abs_grad = std::max(result.max_abs_grad, diag.abs_grad);
        result.min_neg_h = std::min(result.min_neg_h, diag.neg_h);
    }
    if (!std::isfinite(result.min_neg_h)) result.min_neg_h = 0.0;

    return result;
}

std::vector<double> zinb_simulate_residuals(const ZINBFit& fit,
                                            const std::vector<int>& y,
                                            const std::vector<int>& subj_ids,
                                            const std::vector<int>& ct_ids,
                                            int n_subj,
                                            int n_sim,
                                            uint64_t seed) {
    int n = static_cast<int>(y.size());
    double theta = fit.theta();
    double sigma = fit.sigma_re;

    RNG rng(seed);
    std::vector<double> eta_fixed(n), zi_fixed(n);
    for (int i = 0; i < n; ++i) {
        int ct_id = (fit.beta_cond.size() > 1) ? ct_ids[i] : 0;
        eta_fixed[i] = cond_eta_for_ct(fit.beta_cond, ct_id);
        zi_fixed[i] = zi_eta_for_ct(fit.beta_zi, ct_id);
    }

    std::vector<int> sim_le(n, 0), sim_eq(n, 0);
    for (int s = 0; s < n_sim; ++s) {
        std::vector<double> b(n_subj);
        for (int j = 0; j < n_subj; ++j) b[j] = rng.normal() * sigma;

        for (int i = 0; i < n; ++i) {
            double mu_i = safe_mu_from_eta(eta_fixed[i] + b[subj_ids[i]]);
            double pi_i = expit(zi_fixed[i]);
            int y_sim = (rng.uniform() < pi_i) ? 0 : rnbinom(mu_i, theta, rng);

            if (y_sim <= y[i]) ++sim_le[i];
            if (y_sim == y[i]) ++sim_eq[i];
        }
    }

    std::vector<double> residuals(n);
    for (int i = 0; i < n; ++i) {
        double u = rng.uniform();
        residuals[i] = (static_cast<double>(sim_le[i]) - u * static_cast<double>(sim_eq[i]))
                     / static_cast<double>(n_sim);
        residuals[i] = std::clamp(residuals[i], 1e-10, 1.0 - 1e-10);
    }
    return residuals;
}

}  // namespace gsca
