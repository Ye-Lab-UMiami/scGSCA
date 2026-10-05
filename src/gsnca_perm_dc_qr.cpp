// ============================================================================
// gsnca_perm_dc_qr.cpp
//
// Permutation test for GSNCA with ZINB quantile residuals (DC mode).
//
// FLATTENED PARALLELISM: all ZINB fits across all permutations × groups ×
// genes are pooled into one big parallel-for loop.
// 100 perms × 2 groups × 57 genes = 11,400 independent tasks → N threads.
//
// Port of: gsnca_perm_dc_qr_parv7.R
// ============================================================================
#include "MetaGSCA.h"
#include <set>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <atomic>

#ifdef HAS_OPENMP
#include <omp.h>
#endif

static double qnorm_local(double p) {
    if (p <= 0.0) return -1e15;
    if (p >= 1.0) return  1e15;
    if (p == 0.5) return  0.0;
    static const double a[]={-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00};
    static const double b[]={-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01};
    static const double c[]={-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00};
    static const double d[]={7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,3.754408661907416e+00};
    double q,r;
    if(p<0.02425){q=std::sqrt(-2.0*std::log(p));return(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);}
    else if(p<=0.97575){q=p-0.5;r=q*q;return(((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);}
    else{q=std::sqrt(-2.0*std::log(1.0-p));return-(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);}
}

namespace gsca {

static bool check_perm_validity_dc(const DataFrame& group1, const DataFrame& group2,
                                    const DataFrame& objt, size_t gene_start,
                                    bool check_sd, double min_sd, double filter_count) {
    auto gene_cols = objt.gene_names(gene_start - 1);
    try {
        std::vector<std::string> removed;
        if (check_sd) {
            auto r1 = gsnca_sd(group1, group2, gene_cols, gene_start, min_sd);
            auto r2 = gsnca_filt(group1, group2, gene_cols, gene_start, filter_count);
            std::set<std::string> s(r1.begin(), r1.end());
            s.insert(r2.begin(), r2.end());
            removed.assign(s.begin(), s.end());
        } else {
            removed = gsnca_filt(group1, group2, gene_cols, gene_start, filter_count);
        }
        size_t surviving = 0;
        for (auto& g : gene_cols)
            if (std::find(removed.begin(), removed.end(), g) == removed.end()) ++surviving;
        return surviving >= 2;
    } catch (...) { return false; }
}

static void encode_factor_local(const std::vector<std::string>& vals,
                                std::vector<int>& ids, std::vector<std::string>& levels) {
    std::set<std::string> uq(vals.begin(), vals.end());
    levels.assign(uq.begin(), uq.end());
    std::sort(levels.begin(), levels.end());
    std::unordered_map<std::string, int> m;
    for (size_t i = 0; i < levels.size(); ++i) m[levels[i]] = static_cast<int>(i);
    ids.resize(vals.size());
    for (size_t i = 0; i < vals.size(); ++i) ids[i] = m[vals[i]];
}

// -------------------------------------------------------------------------
std::vector<double> gsnca_perm_dc_qr(
    const std::vector<std::string>& domain, int nperm,
    const DataFrame& object, size_t gene_start, int group_num_1,
    bool check_sd, double min_sd, int max_skip,
    double filter_count, bool na_counts,
    const std::string& cor_method, const std::string& subj_var,
    const std::string& ct_var,
    int num_cores, uint64_t perm_seed, uint64_t seed, int n_sim, bool paired)
{
    RNG rng(perm_seed);
    const auto& subj_col = object.str_col(subj_var);
    std::unordered_map<std::string, std::vector<size_t>> subj_rows;
    for (size_t i = 0; i < object.nrow(); ++i)
        subj_rows[subj_col[i]].push_back(i);
    size_t nv = object.nrow();

    // === Phase 1: generate valid permutations (sequential) ===
    std::set<std::string> perm_keys;
    std::vector<std::vector<std::string>> perm_indices;
    perm_indices.reserve(nperm);
    std::set<std::string> original_group1(domain.begin(), domain.begin() + group_num_1);
    int skip_total = 0, itr = 0;
    while (itr < nperm) {
        std::vector<std::string> group1_subjects;
        std::vector<std::string> group2_subjects;
        std::string key;

        if (paired) {
            // Paired ST permutation: each subject independently flips sample type
            // with probability 0.5, while all cells from the subject move together.
            for (size_t k = 0; k < domain.size(); ++k) {
                const auto& subj = domain[k];
                bool originally_group1 = static_cast<int>(k) < group_num_1;
                bool flip = rng.coin_flip();
                bool perm_group1 = flip ? !originally_group1 : originally_group1;
                if (perm_group1) group1_subjects.push_back(subj);
                else group2_subjects.push_back(subj);
                key += (perm_group1 ? '1' : '2');
            }
            if (group1_subjects.empty() || group2_subjects.empty()) continue;
        } else {
            auto randperm = domain;
            rng.shuffle(randperm);
            bool same_split = true;
            for (int k = 0; k < group_num_1; ++k) {
                if (!original_group1.count(randperm[k])) {
                    same_split = false;
                    break;
                }
            }
            if (same_split) continue;
            for (auto& s : randperm) key += s + "\r";
            group1_subjects.assign(randperm.begin(), randperm.begin() + group_num_1);
            group2_subjects.assign(randperm.begin() + group_num_1, randperm.end());
        }

        if (perm_keys.count(key)) continue;
        perm_keys.insert(key);

        std::vector<size_t> ord1, ord2;
        for (auto& subj : group1_subjects) { auto& r = subj_rows[subj]; ord1.insert(ord1.end(), r.begin(), r.end()); }
        for (auto& subj : group2_subjects) { auto& r = subj_rows[subj]; ord2.insert(ord2.end(), r.begin(), r.end()); }
        size_t nv1 = ord1.size();
        std::vector<size_t> ord = ord1;
        ord.insert(ord.end(), ord2.begin(), ord2.end());

        auto sub1 = object.row_subset(ord1);
        auto sub2 = object.row_subset(ord2);
        NumCol c1(nv1), c2(nv-nv1);
        std::iota(c1.begin(), c1.end(), 1.0);
        std::iota(c2.begin(), c2.end(), 1.0);
        auto g1 = DataFrame::cbind_num("cid", c1, sub1);
        auto g2 = DataFrame::cbind_num("cid", c2, sub2);

        if (!check_perm_validity_dc(g1, g2, object, gene_start, check_sd, min_sd, filter_count)) {
            if (++skip_total >= max_skip) { std::cerr << "max.skip reached\n"; break; }
            continue;
        }
        std::vector<std::string> perm_subjects = group1_subjects;
        perm_subjects.insert(perm_subjects.end(), group2_subjects.begin(), group2_subjects.end());
        perm_indices.push_back(std::move(perm_subjects));
        ++itr;
    }

    int n_valid = static_cast<int>(perm_indices.size());
    std::cout << "Number of valid DC permutations: " << n_valid << std::endl;
    if (n_valid == 0) return {};

    // === Phase 2: stream permutations one-by-one to keep memory bounded ===
    auto genes = object.gene_names(gene_start - 1);
    size_t ng = genes.size();
    size_t total = static_cast<size_t>(n_valid) * 2 * ng;

    std::atomic<size_t> done{0};
    std::cout << "  Flattened parallel: " << total << " ZINB tasks, "
              << num_cores << " threads" << std::endl;
    std::vector<double> D_perm(n_valid, std::numeric_limits<double>::quiet_NaN());

    for (int p = 0; p < n_valid; ++p) {
        auto& rp = perm_indices[p];
        std::vector<size_t> ord;
        for (auto& s : rp) { auto& r = subj_rows.at(s); ord.insert(ord.end(), r.begin(), r.end()); }
        size_t nv1 = 0;
        for (int k = 0; k < group_num_1; ++k) nv1 += subj_rows.at(rp[k]).size();

        DataFrame sub1 = object.row_subset(std::vector<size_t>(ord.begin(), ord.begin()+nv1));
        DataFrame sub2 = object.row_subset(std::vector<size_t>(ord.begin()+nv1, ord.end()));
        NumCol c1(nv1), c2(nv-nv1);
        std::iota(c1.begin(), c1.end(), 1.0);
        std::iota(c2.begin(), c2.end(), 1.0);
        DataFrame group1 = DataFrame::cbind_num("cid", c1, sub1);
        DataFrame group2 = DataFrame::cbind_num("cid", c2, sub2);

        Mat qr_g1 = Mat::Zero(static_cast<Eigen::Index>(sub1.nrow()), static_cast<Eigen::Index>(ng));
        Mat qr_g2 = Mat::Zero(static_cast<Eigen::Index>(sub2.nrow()), static_cast<Eigen::Index>(ng));

#ifdef HAS_OPENMP
        #pragma omp parallel for schedule(dynamic, 1) num_threads(num_cores)
#endif
        for (size_t rem = 0; rem < 2 * ng; ++rem) {
            size_t gi = rem / ng;   // 0=group1, 1=group2
            size_t ge = rem % ng;   // gene index

            const DataFrame& sub = (gi == 0) ? sub1 : sub2;
            size_t nc = sub.nrow();

            // Encode factors
            std::vector<int> sids, cids;
            std::vector<std::string> slvl, clvl;
            encode_factor_local(sub.str_col(subj_var), sids, slvl);
            encode_factor_local(sub.str_col(ct_var), cids, clvl);

            // Sort
            std::vector<size_t> so(nc);
            std::iota(so.begin(), so.end(), 0);
            std::sort(so.begin(), so.end(), [&](size_t a, size_t b) {
                return cids[a] != cids[b] ? cids[a] < cids[b] : sids[a] < sids[b];
            });
            std::vector<int> ss(nc), cs(nc);
            for (size_t i = 0; i < nc; ++i) { ss[i] = sids[so[i]]; cs[i] = cids[so[i]]; }

            // Gene counts
            const auto& gcol = sub.num_col(genes[ge]);
            std::vector<int> ys(nc);
            for (size_t i = 0; i < nc; ++i) {
                double v = gcol[so[i]];
                if (std::isnan(v) && na_counts) v = 0.0;
                ys[i] = std::max(0, static_cast<int>(std::round(v)));
            }

            uint64_t tseed = seed + 100000ULL * p + 50000ULL * gi + ge;
            ZINBFit fit = zinb_fit("dc", ys, ss, cs,
                                   static_cast<int>(slvl.size()),
                                   static_cast<int>(clvl.size()));
            if (fit.converged) {
                auto u = zinb_simulate_residuals(fit, ys, ss, cs,
                                                 static_cast<int>(slvl.size()), n_sim, tseed);
                Mat& tgt = (gi == 0) ? qr_g1 : qr_g2;
                for (size_t i = 0; i < nc; ++i)
                    tgt(static_cast<Eigen::Index>(so[i]), static_cast<Eigen::Index>(ge)) = qnorm_local(u[i]);
            } else {
#ifdef HAS_OPENMP
                #pragma omp critical
#endif
                std::cerr << "  [FAIL] ZINB fit failed: perm=" << p
                          << " grp=" << gi << " gene=" << genes[ge] << std::endl;
            }

            size_t d = ++done;
            if (d % 100 == 0 || d == total) {
#ifdef HAS_OPENMP
                #pragma omp critical
#endif
                std::cout << "    DC ZINB: " << d << "/" << total
                          << " (" << static_cast<int>(100.0*d/total) << "%)" << std::endl;
            }
        }

        QRMatrix qm1, qm2;
        qm1.data = std::move(qr_g1); qm1.gene_names = genes; qm1.row_names = group1.row_names;
        qm2.data = std::move(qr_g2); qm2.gene_names = genes; qm2.row_names = group2.row_names;
        auto qr = QRMatrix::rbind(qm1, qm2);
        try {
            D_perm[p] = gsnca_stat_qr(group1, group2, gene_start,
                                      na_counts, cor_method, subj_var, qr);
        } catch (...) {}
    }

    std::vector<double> result;
    for (double d : D_perm) if (!std::isnan(d)) result.push_back(d);
    std::cout << "  DC permutation done: " << result.size() << " valid stats." << std::endl;
    return result;
}

} // namespace gsca
