// ============================================================================
// gsnca_perm_ct_qr.cpp
//
// Permutation test for GSNCA with ZINB quantile residuals (CT mode).
//
// CT permutation scheme:
// perform a DC-style subject permutation across disease groups, then within
// each permuted disease group randomly relabel cell types within each subject
// while preserving that subject's original CT1 count.
//
// FLATTENED PARALLELISM: all ZINB fits across all permutations × cell types ×
// genes are pooled into one big parallel-for loop.
//
// Port of: gsnca_perm_ct_qr_parv8.R
// ============================================================================
#include "MetaGSCA.h"
#include <set>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <atomic>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>

#ifdef HAS_OPENMP
#include <omp.h>
#endif

static double qnorm_local_ct(double p) {
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

static void encode_factor_ct(const std::vector<std::string>& vals,
                              std::vector<int>& ids, std::vector<std::string>& levels) {
    std::set<std::string> uq(vals.begin(), vals.end());
    levels.assign(uq.begin(), uq.end());
    std::sort(levels.begin(), levels.end());
    std::unordered_map<std::string, int> m;
    for (size_t i = 0; i < levels.size(); ++i) m[levels[i]] = static_cast<int>(i);
    ids.resize(vals.size());
    for (size_t i = 0; i < vals.size(); ++i) ids[i] = m[vals[i]];
}

static bool check_perm_validity_ct(const DataFrame& ct1_df, const DataFrame& ct2_df,
                                   const DataFrame& object, size_t gene_start,
                                   bool check_sd, double min_sd, double filter_count) {
    auto gene_cols = object.gene_names(gene_start);
    try {
        std::vector<std::string> removed;
        if (check_sd) {
            auto r1 = gsnca_sd(ct1_df, ct2_df, gene_cols, gene_start, min_sd);
            auto r2 = gsnca_filt(ct1_df, ct2_df, gene_cols, gene_start, filter_count);
            std::set<std::string> s(r1.begin(), r1.end());
            s.insert(r2.begin(), r2.end());
            removed.assign(s.begin(), s.end());
        } else {
            removed = gsnca_filt(ct1_df, ct2_df, gene_cols, gene_start, filter_count);
        }
        size_t surv = 0;
        for (const auto& g : gene_cols)
            if (std::find(removed.begin(), removed.end(), g) == removed.end()) ++surv;
        return surv >= 2;
    } catch (...) {
        return false;
    }
}

static std::string sanitize_token_ct(const std::string& s) {
    std::string out;
    bool prev_us = false;
    for (unsigned char ch : s) {
        if (std::isalnum(ch)) {
            out.push_back(static_cast<char>(std::tolower(ch)));
            prev_us = false;
        } else if (!prev_us && !out.empty()) {
            out.push_back('_');
            prev_us = true;
        }
    }
    while (!out.empty() && out.back() == '_')
        out.pop_back();
    return out.empty() ? "level" : out;
}

static std::string perm_prefix_ct(int perm_index) {
    std::ostringstream oss;
    oss << "perm_" << std::setw(6) << std::setfill('0') << (perm_index + 1);
    return oss.str();
}

static void write_named_matrix_tsv_ct(const Mat& mat,
                                      const std::vector<std::string>& names,
                                      const std::string& path) {
    std::ofstream ofs(path);
    if (!ofs.is_open())
        throw std::runtime_error("write_named_matrix_tsv_ct: cannot open file: " + path);

    ofs << "gene";
    for (const auto& name : names)
        ofs << '\t' << name;
    ofs << '\n';

    for (Eigen::Index i = 0; i < mat.rows(); ++i) {
        ofs << names[static_cast<size_t>(i)];
        for (Eigen::Index j = 0; j < mat.cols(); ++j)
            ofs << '\t' << mat(i, j);
        ofs << '\n';
    }
}

// -------------------------------------------------------------------------
CTPermResult gsnca_perm_ct_qr(
    int nperm, const DataFrame& object, size_t gene_start,
    const std::vector<int>& group,
    const std::vector<std::string>& cell_type,
    const std::string& preserve_cell_type,
    bool check_sd, double min_sd, int max_skip,
    double filter_count, bool na_counts,
    const std::string& cor_method, const std::string& subj_var,
    int group_num_1, int num_cores, uint64_t perm_seed,
    uint64_t seed, int n_sim, bool paired,
    bool save_cor_mats,
    const std::string& cor_mat_dir)
{
    RNG rng(perm_seed);
    CTPermResult result;

    std::set<std::string> ct_set(cell_type.begin(), cell_type.end());
    std::vector<std::string> ct_levels(ct_set.begin(), ct_set.end());
    std::sort(ct_levels.begin(), ct_levels.end());
    if (ct_levels.size() != 2)
        throw std::runtime_error("gsnca_perm_ct_qr expects exactly 2 cell types.");
    if (!preserve_cell_type.empty()) {
        auto it = std::find(ct_levels.begin(), ct_levels.end(), preserve_cell_type);
        if (it == ct_levels.end()) {
            throw std::runtime_error("gsnca_perm_ct_qr: preserve_cell_type '" +
                                     preserve_cell_type + "' not found.");
        }
        std::rotate(ct_levels.begin(), it, it + 1);
    }

    std::set<int> grp_set(group.begin(), group.end());
    std::vector<int> grp_levels(grp_set.begin(), grp_set.end());
    std::sort(grp_levels.begin(), grp_levels.end());
    if (grp_levels.size() != 2)
        throw std::runtime_error("gsnca_perm_ct_qr expects exactly 2 disease groups.");

    const auto& subj_col = object.str_col(subj_var);
    std::set<std::string> subj_set(subj_col.begin(), subj_col.end());
    std::vector<std::string> subjects(subj_set.begin(), subj_set.end());
    std::sort(subjects.begin(), subjects.end());
    if (subjects.size() < 2)
        throw std::runtime_error("gsnca_perm_ct_qr: need >= 2 subjects.");
    size_t n_subj = subjects.size();

    std::unordered_map<std::string, std::vector<size_t>> subj_row_map;
    for (size_t i = 0; i < object.nrow(); ++i)
        subj_row_map[subj_col[i]].push_back(i);

    std::unordered_map<std::string, size_t> subj_ct1_count;
    std::unordered_map<std::string, int> subj_group;
    for (size_t i = 0; i < object.nrow(); ++i) {
        if (cell_type[i] == ct_levels[0]) ++subj_ct1_count[subj_col[i]];
        auto it = subj_group.find(subj_col[i]);
        if (it == subj_group.end()) {
            subj_group[subj_col[i]] = group[i];
        } else if (it->second != group[i]) {
            throw std::runtime_error("gsnca_perm_ct_qr: subject appears in multiple disease groups.");
        }
    }

    size_t n_group1_subjects = 0;
    for (const auto& subject : subjects)
        if (subj_group.at(subject) == grp_levels[0]) ++n_group1_subjects;
    if (static_cast<int>(n_group1_subjects) != group_num_1)
        throw std::runtime_error("gsnca_perm_ct_qr: group_num_1 does not match pooled subject counts.");

    // === Phase 1: generate valid joint DC+CT permutations ===
    struct PermSplit {
        std::vector<size_t> g1_ct1;
        std::vector<size_t> g1_ct2;
        std::vector<size_t> g2_ct1;
        std::vector<size_t> g2_ct2;
    };

    std::set<std::string> perm_keys;
    std::vector<PermSplit> perm_splits;
    perm_splits.reserve(nperm);
    int itr = 0, skip_total = 0;

    while (itr < nperm) {
        std::set<std::string> orig_group1;
        for (const auto& subject : subjects)
            if (subj_group.at(subject) == grp_levels[0]) orig_group1.insert(subject);

        std::vector<std::string> perm_subject_order;
        std::unordered_map<std::string, int> perm_subject_group;
        std::string subject_key;

        if (paired) {
            // Paired CT permutation: each subject independently flips sample type
            // with probability 0.5 before within-subject cell-type reassignment.
            for (const auto& subject : subjects) {
                int original_group = subj_group.at(subject);
                bool flip = rng.coin_flip();
                int perm_group = flip ? ((original_group == grp_levels[0]) ? grp_levels[1] : grp_levels[0])
                                      : original_group;
                perm_subject_group[subject] = perm_group;
                subject_key += (perm_group == grp_levels[0]) ? '1' : '2';
            }
        } else {
            auto randperm = subjects;
            rng.shuffle(randperm);

            std::set<std::string> perm_group1(randperm.begin(), randperm.begin() + group_num_1);
            if (perm_group1 == orig_group1) continue;

            for (size_t s = 0; s < n_subj; ++s) {
                const auto& subject = randperm[s];
                perm_subject_group[subject] = (static_cast<int>(s) < group_num_1) ? grp_levels[0] : grp_levels[1];
                subject_key += subject + "\r";
            }
        }

        std::vector<char> state(object.nrow(), 0); // 0=g1ct1,1=g1ct2,2=g2ct1,3=g2ct2
        std::vector<size_t> g1_ct1, g1_ct2, g2_ct1, g2_ct2;
        g1_ct1.reserve(object.nrow());
        g1_ct2.reserve(object.nrow());
        g2_ct1.reserve(object.nrow());
        g2_ct2.reserve(object.nrow());

        for (size_t s = 0; s < n_subj; ++s) {
            const auto& subject = subjects[s];
            const auto& rows = subj_row_map.at(subject);
            size_t n_ct1 = subj_ct1_count[subject];
            auto chosen = rng.sample_indices(rows.size(), n_ct1);
            std::vector<char> chosen_mask(rows.size(), 0);
            for (size_t idx : chosen) chosen_mask[idx] = 1;

            bool in_group1 = perm_subject_group.at(subject) == grp_levels[0];
            for (size_t local_idx = 0; local_idx < rows.size(); ++local_idx) {
                size_t row = rows[local_idx];
                if (in_group1) {
                    if (chosen_mask[local_idx]) {
                        state[row] = 0;
                        g1_ct1.push_back(row);
                    } else {
                        state[row] = 1;
                        g1_ct2.push_back(row);
                    }
                } else {
                    if (chosen_mask[local_idx]) {
                        state[row] = 2;
                        g2_ct1.push_back(row);
                    } else {
                        state[row] = 3;
                        g2_ct2.push_back(row);
                    }
                }
            }
        }

        std::string key = subject_key + "|";
        key.reserve(subject_key.size() + object.nrow() + 1);
        for (char b : state) key += static_cast<char>('A' + b);
        if (perm_keys.count(key)) continue;
        perm_keys.insert(key);

        if (g1_ct1.empty() || g1_ct2.empty() || g2_ct1.empty() || g2_ct2.empty()) {
            if (++skip_total >= max_skip) { std::cerr << "max.skip reached\n"; break; }
            continue;
        }

        bool valid_g1 = check_perm_validity_ct(
            object.row_subset(g1_ct1), object.row_subset(g1_ct2),
            object, gene_start, check_sd, min_sd, filter_count);
        bool valid_g2 = check_perm_validity_ct(
            object.row_subset(g2_ct1), object.row_subset(g2_ct2),
            object, gene_start, check_sd, min_sd, filter_count);

        if (!valid_g1 && !valid_g2) {
            if (++skip_total >= max_skip) { std::cerr << "max.skip reached\n"; break; }
            continue;
        }
        perm_splits.push_back({std::move(g1_ct1), std::move(g1_ct2),
                               std::move(g2_ct1), std::move(g2_ct2)});
        ++itr;
    }

    int n_valid = static_cast<int>(perm_splits.size());
    std::cout << "Number of valid joint CT permutations: " << n_valid << std::endl;
    if (n_valid == 0) return result;
    if (save_cor_mats) {
        if (cor_mat_dir.empty())
            throw std::runtime_error("gsnca_perm_ct_qr: save_cor_mats requested without cor_mat_dir.");
        std::filesystem::create_directories(cor_mat_dir);
    }

    // === Phase 2: build permuted DataFrames ===
    struct PD { DataFrame g1_ct1, g1_ct2, g2_ct1, g2_ct2; };
    std::vector<PD> pd(n_valid);

    for (int p = 0; p < n_valid; ++p) {
        pd[p].g1_ct1 = object.row_subset(perm_splits[p].g1_ct1);
        pd[p].g1_ct2 = object.row_subset(perm_splits[p].g1_ct2);
        pd[p].g2_ct1 = object.row_subset(perm_splits[p].g2_ct1);
        pd[p].g2_ct2 = object.row_subset(perm_splits[p].g2_ct2);
    }

    // === Phase 3: FLATTENED parallel ZINB ===
    auto genes = object.gene_names(gene_start);
    size_t ng = genes.size();
    size_t total = static_cast<size_t>(n_valid) * 4 * ng;

    struct PQR { Mat g1_ct1, g1_ct2, g2_ct1, g2_ct2; };
    std::vector<PQR> pqr(n_valid);
    for (int p = 0; p < n_valid; ++p) {
        pqr[p].g1_ct1 = Mat::Zero(static_cast<Eigen::Index>(pd[p].g1_ct1.nrow()), static_cast<Eigen::Index>(ng));
        pqr[p].g1_ct2 = Mat::Zero(static_cast<Eigen::Index>(pd[p].g1_ct2.nrow()), static_cast<Eigen::Index>(ng));
        pqr[p].g2_ct1 = Mat::Zero(static_cast<Eigen::Index>(pd[p].g2_ct1.nrow()), static_cast<Eigen::Index>(ng));
        pqr[p].g2_ct2 = Mat::Zero(static_cast<Eigen::Index>(pd[p].g2_ct2.nrow()), static_cast<Eigen::Index>(ng));
    }

    std::atomic<size_t> done{0};
    std::cout << "  Flattened parallel: " << total << " CT ZINB tasks, "
              << num_cores << " threads" << std::endl;

#ifdef HAS_OPENMP
    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_cores)
#endif
    for (size_t t = 0; t < total; ++t) {
        size_t pi = t / (4 * ng);
        size_t rem = t % (4 * ng);
        size_t ci = rem / ng;   // 0=g1ct1, 1=g1ct2, 2=g2ct1, 3=g2ct2
        size_t ge = rem % ng;

        const DataFrame& sub = (ci == 0) ? pd[pi].g1_ct1 :
                               (ci == 1) ? pd[pi].g1_ct2 :
                               (ci == 2) ? pd[pi].g2_ct1 :
                                           pd[pi].g2_ct2;
        size_t nc = sub.nrow();

        // Encode subject factor
        std::vector<int> sids;
        std::vector<std::string> slvl;
        encode_factor_ct(sub.str_col(subj_var), sids, slvl);

        // Sort by subject
        std::vector<size_t> so(nc);
        std::iota(so.begin(), so.end(), 0);
        std::sort(so.begin(), so.end(), [&](size_t a, size_t b) { return sids[a] < sids[b]; });
        std::vector<int> ss(nc);
        for (size_t i = 0; i < nc; ++i) ss[i] = sids[so[i]];

        // Gene counts
        const auto& gcol = sub.num_col(genes[ge]);
        std::vector<int> ys(nc);
        for (size_t i = 0; i < nc; ++i) {
            double v = gcol[so[i]];
            if (std::isnan(v) && na_counts) v = 0.0;
            ys[i] = std::max(0, static_cast<int>(std::round(v)));
        }

        std::vector<int> ct_dummy(nc, 0);
        uint64_t tseed = seed + 100000ULL * pi + 50000ULL * ci + ge;

        ZINBFit fit = zinb_fit("ct", ys, ss, ct_dummy,
                               static_cast<int>(slvl.size()), 1);
        if (fit.converged) {
            auto u = zinb_simulate_residuals(fit, ys, ss, ct_dummy,
                                             static_cast<int>(slvl.size()), n_sim, tseed);
            Mat& tgt = (ci == 0) ? pqr[pi].g1_ct1 :
                       (ci == 1) ? pqr[pi].g1_ct2 :
                       (ci == 2) ? pqr[pi].g2_ct1 :
                                   pqr[pi].g2_ct2;
            for (size_t i = 0; i < nc; ++i)
                tgt(static_cast<Eigen::Index>(so[i]), static_cast<Eigen::Index>(ge)) = qnorm_local_ct(u[i]);
        }else {
            int out_group = (ci < 2) ? 1 : 2;
            int out_ct = (ci % 2);
            #pragma omp critical
            std::cerr << "  [FAIL] ZINB fit failed: perm=" << pi
                      << " group=" << out_group << " ct=" << out_ct
                      << " gene=" << genes[ge] << std::endl;
        }

        size_t d = ++done;
        if (d % 100 == 0 || d == total) {
#ifdef HAS_OPENMP
            #pragma omp critical
#endif
            std::cout << "    CT ZINB: " << d << "/" << total
                      << " (" << static_cast<int>(100.0*d/total) << "%)" << std::endl;
        }
    }

    // === Phase 4: compute test statistics ===
    std::cout << "  Computing CT test statistics..." << std::endl;
    std::vector<double> D_perm_g1(n_valid, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> D_perm_g2(n_valid, std::numeric_limits<double>::quiet_NaN());

#ifdef HAS_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(num_cores)
#endif
    for (int p = 0; p < n_valid; ++p) {
        try {
            QRMatrix q1, q2;
            q1.data = pqr[p].g1_ct1; q1.gene_names = genes; q1.row_names = pd[p].g1_ct1.row_names;
            q2.data = pqr[p].g1_ct2; q2.gene_names = genes; q2.row_names = pd[p].g1_ct2.row_names;
            auto qr_g1 = QRMatrix::rbind(q1, q2);
            auto stat_g1 = gsnca_stat_qr_detail(pd[p].g1_ct1, pd[p].g1_ct2, gene_start,
                                                na_counts, cor_method, subj_var, qr_g1);
            D_perm_g1[p] = stat_g1.D;
            if (save_cor_mats) {
                const std::string prefix = cor_mat_dir + "/" + perm_prefix_ct(p);
                write_named_matrix_tsv_ct(
                    stat_g1.cor_mat1, genes,
                    prefix + ".group_" + std::to_string(grp_levels[0]) + "." +
                    sanitize_token_ct(ct_levels[0]) + "_corr.tsv");
                write_named_matrix_tsv_ct(
                    stat_g1.cor_mat2, genes,
                    prefix + ".group_" + std::to_string(grp_levels[0]) + "." +
                    sanitize_token_ct(ct_levels[1]) + "_corr.tsv");
            }
        } catch (...) {}
        try {
            QRMatrix q1, q2;
            q1.data = pqr[p].g2_ct1; q1.gene_names = genes; q1.row_names = pd[p].g2_ct1.row_names;
            q2.data = pqr[p].g2_ct2; q2.gene_names = genes; q2.row_names = pd[p].g2_ct2.row_names;
            auto qr_g2 = QRMatrix::rbind(q1, q2);
            auto stat_g2 = gsnca_stat_qr_detail(pd[p].g2_ct1, pd[p].g2_ct2, gene_start,
                                                na_counts, cor_method, subj_var, qr_g2);
            D_perm_g2[p] = stat_g2.D;
            if (save_cor_mats) {
                const std::string prefix = cor_mat_dir + "/" + perm_prefix_ct(p);
                write_named_matrix_tsv_ct(
                    stat_g2.cor_mat1, genes,
                    prefix + ".group_" + std::to_string(grp_levels[1]) + "." +
                    sanitize_token_ct(ct_levels[0]) + "_corr.tsv");
                write_named_matrix_tsv_ct(
                    stat_g2.cor_mat2, genes,
                    prefix + ".group_" + std::to_string(grp_levels[1]) + "." +
                    sanitize_token_ct(ct_levels[1]) + "_corr.tsv");
            }
        } catch (...) {}
    }

    for (double d : D_perm_g1) {
        if (!std::isnan(d))
            result.group1.push_back(d);
    }
    for (double d : D_perm_g2) {
        if (!std::isnan(d))
            result.group2.push_back(d);
    }
    for (int p = 0; p < n_valid; ++p) {
        if (std::isnan(D_perm_g1[p]) || std::isnan(D_perm_g2[p]))
            continue;
        const double signed_delta = D_perm_g2[p] - D_perm_g1[p];
        result.signed_delta.push_back(signed_delta);
        result.abs_delta.push_back(std::abs(signed_delta));
    }
    std::cout << "  CT permutation done: "
              << result.group1.size() << " valid stats for group 1, "
              << result.group2.size() << " valid stats for group 2." << std::endl;
    return result;
}

} // namespace gsca
