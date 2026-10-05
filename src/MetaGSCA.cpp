// ============================================================================
// MetaGSCA.cpp
//
// Main orchestration functions for the MetaGSCA pipeline.
//
// GSAR_qr_dc — Test differential co-expression across DISEASE CONDITIONS
//              (pooling cell types).
//              Permutes subject labels between conditions.
//
// GSAR_qr_ct — Test differential co-expression across CELL TYPES
//              within each disease condition.
//              Permutes cell-type labels at the subject level.
//
// Port of: MetaGSCA_qr_parv5.R (GSAR_qr_dc, GSAR_qr_ct)
// ============================================================================
#include "MetaGSCA.h"
#include "ct_tail.h"
#include <set>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <cctype>
#include <cmath>

#ifdef HAS_OPENMP
#include <omp.h>
#endif

namespace gsca {

// =========================================================================
// Helper: remove genes from a DataFrame by column name
// =========================================================================
static void remove_genes(DataFrame& df,
                         const std::vector<std::string>& to_remove) {
    for (auto& g : to_remove) {
        if (df.has_col(g)) df.remove_col(g);
    }
}

// =========================================================================
// Helper: collect unique values from a string column, ordered
// =========================================================================
static std::vector<std::string> unique_sorted(const StrCol& col) {
    std::set<std::string> s(col.begin(), col.end());
    return {s.begin(), s.end()};
}

// =========================================================================
// Helper: collect unique values from a string column in first-seen order
// =========================================================================
static std::vector<std::string> unique_in_order(const StrCol& col) {
    std::set<std::string> seen;
    std::vector<std::string> out;
    out.reserve(col.size());
    for (const auto& val : col) {
        if (seen.insert(val).second) out.push_back(val);
    }
    return out;
}

static std::string sanitize_token_meta(const std::string& s) {
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

static double mean_or_zero(const std::vector<double>& x) {
    if (x.empty())
        return 0.0;
    double s = 0.0;
    for (double v : x)
        s += v;
    return s / static_cast<double>(x.size());
}

static void write_named_matrix_tsv_meta(const Mat& mat,
                                        const std::vector<std::string>& names,
                                        const std::string& path) {
    std::ofstream ofs(path);
    if (!ofs.is_open())
        throw std::runtime_error("write_named_matrix_tsv_meta: cannot open file: " + path);

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

// =========================================================================
// GSAR_qr_dc
//
// Test differential co-expression across disease conditions.
// =========================================================================
GSCAOutput GSAR_qr_dc(const GSCAConfig& cfg,
                       DataFrame& gematrix,
                       const std::vector<int>& group,
                       int group_ref,
                       const std::vector<std::string>& genelist_in)
{
    auto t_start = std::chrono::steady_clock::now();
    GSCAOutput output;

#ifdef HAS_OPENMP
    omp_set_num_threads(cfg.num_cores);
    std::cout << "OpenMP enabled with " << cfg.num_cores << " threads." << std::endl;
#endif

    // Determine gene list
    std::vector<std::string> genelist = genelist_in;
    bool use_all = (genelist.empty() ||
                    (genelist.size() == 1 && genelist[0] == "NA"));

    if (use_all) {
        genelist = gematrix.gene_names(cfg.gene_start);
        output.genes_not_found = "-";
    } else {
        // Find genes not present in data
        auto all_cols = gematrix.col_names();
        std::set<std::string> col_set(all_cols.begin(), all_cols.end());
        std::string not_found;
        for (auto& g : genelist) {
            if (!col_set.count(g)) {
                if (!not_found.empty()) not_found += ",";
                not_found += g;
            }
        }
        output.genes_not_found = not_found.empty() ? "-" : not_found;
    }

    // Subset to metadata + gene columns
    std::vector<std::string> keep_cols;
    for (size_t c = 0; c < cfg.gene_start; ++c)
        keep_cols.push_back(gematrix.col_name(c));
    for (auto& g : genelist)
        if (gematrix.has_col(g)) keep_cols.push_back(g);

    DataFrame object = gematrix.col_subset(keep_cols);
    size_t gene_start = cfg.gene_start;

    // Round gene columns
    object.round_cols(gene_start);

    // Input validation
    gsnca_check_dc(cfg.projectname, object, group, group_ref, cfg.ct_var);

    // Split by group
    std::vector<bool> mask1(object.nrow()), mask2(object.nrow());
    for (size_t i = 0; i < object.nrow(); ++i) {
        mask1[i] = (group[i] == 1);
        mask2[i] = (group[i] == 2);
    }
    auto object1_1 = object.row_subset(mask1);
    auto object2_1 = object.row_subset(mask2);
    auto object_1  = DataFrame::rbind(object1_1, object2_1);

    // Add cid column
    NumCol cid1(object1_1.nrow()), cid2(object2_1.nrow());
    std::iota(cid1.begin(), cid1.end(), 1.0);
    std::iota(cid2.begin(), cid2.end(), 1.0);
    object1_1 = DataFrame::cbind_num("cid", cid1, object1_1);
    object2_1 = DataFrame::cbind_num("cid", cid2, object2_1);
    gene_start += 1;  // shifted by cid

    // SD check
    std::vector<std::string> genes_removed1;
    auto genes = object.gene_names(cfg.gene_start);
    if (cfg.check_sd) {
        std::cout << "Checking gene SD..." << std::endl;
        genes_removed1 = gsnca_sd(object1_1, object2_1, genes,
                                  gene_start, cfg.min_sd);
        remove_genes(object, genes_removed1);
        remove_genes(object1_1, genes_removed1);
        remove_genes(object2_1, genes_removed1);
        remove_genes(object_1, genes_removed1);
    }

    // Low-count filter
    std::cout << "Filtering low-count genes..." << std::endl;
    genes = object.gene_names(cfg.gene_start);
    auto genes_removed2 = gsnca_filt(object1_1, object2_1, genes,
                                     gene_start, cfg.filter_count);
    remove_genes(object, genes_removed2);
    remove_genes(object1_1, genes_removed2);
    remove_genes(object2_1, genes_removed2);
    remove_genes(object_1, genes_removed2);

    // Merge removed gene lists
    std::set<std::string> all_removed(genes_removed1.begin(), genes_removed1.end());
    all_removed.insert(genes_removed2.begin(), genes_removed2.end());
    output.genes_removed.assign(all_removed.begin(), all_removed.end());

    // --- Compute ZINB quantile residuals (observed data) ---
    std::cout << "Computing ZINB quantile residuals..." << std::endl;

    // Drop cid for ZINB fitting
    auto cols_no_cid1 = object1_1.col_names();
    cols_no_cid1.erase(cols_no_cid1.begin());  // remove "cid"
    auto objt_g1 = object1_1.col_subset(cols_no_cid1);

    auto cols_no_cid2 = object2_1.col_names();
    cols_no_cid2.erase(cols_no_cid2.begin());
    auto objt_g2 = object2_1.col_subset(cols_no_cid2);

    auto qr_g1 = zinb_qr_dc_sep(objt_g1, gene_start - 1,
                                 cfg.subj_var, cfg.ct_var,
                                 cfg.na_counts, cfg.n_sim, cfg.seed, false);
    qr_g1.row_names = objt_g1.row_names;

    auto qr_g2 = zinb_qr_dc_sep(objt_g2, gene_start - 1,
                                 cfg.subj_var, cfg.ct_var,
                                 cfg.na_counts, cfg.n_sim, cfg.seed, false);
    qr_g2.row_names = objt_g2.row_names;

    auto qr_mat_1 = QRMatrix::rbind(qr_g1, qr_g2);
    std::cout << "Quantile residuals computed." << std::endl;

    // --- Observed test statistic ---
    double D_obs = gsnca_stat_qr(object1_1, object2_1, gene_start,
                                 cfg.na_counts, cfg.cor_method,
                                 cfg.subj_var, qr_mat_1);
    std::cout << "Observed sample completed." << std::endl;

    // --- Permutation test ---
    // Build domain: unique subject IDs, group1 first
    auto subj1 = unique_in_order(objt_g1.str_col(cfg.subj_var));
    auto subj2 = unique_in_order(objt_g2.str_col(cfg.subj_var));
    std::vector<std::string> domain;
    domain.insert(domain.end(), subj1.begin(), subj1.end());
    for (auto& s : subj2)
        if (std::find(subj1.begin(), subj1.end(), s) == subj1.end())
            domain.push_back(s);

    int group_num_1 = static_cast<int>(subj1.size());

    std::vector<double> D_perm;
    double pval = 0.0;

    if (domain.size() <= 2) {
        std::cout << "We need more than one sample per disease condition "
                  << "to create permutation samples. "
                  << "We can only report the observed test statistic."
                  << std::endl;
    } else {
        D_perm = gsnca_perm_dc_qr(
            domain, cfg.nperm, object_1, gene_start, group_num_1,
            cfg.check_sd, cfg.min_sd, cfg.max_skip,
            cfg.filter_count, cfg.na_counts,
            cfg.cor_method, cfg.subj_var, cfg.ct_var,
            cfg.num_cores, 42, cfg.seed, cfg.n_sim, cfg.paired);

        if (!D_perm.empty()) {
            int count_ge = 0;
            for (double d : D_perm)
                if (d >= D_obs) ++count_ge;
            pval = static_cast<double>(count_ge + 1) /
                   static_cast<double>(D_perm.size() + 1);
        }
    }

    // --- Build result ---
    int n_tested = static_cast<int>(object.ncol() - cfg.gene_start);
    GSCAResult res;
    res.original_ts   = D_obs;
    res.perm_ts_mean  = 0.0;
    if (!D_perm.empty()) {
        double sum = 0.0;
        for (double d : D_perm) sum += d;
        res.perm_ts_mean = sum / D_perm.size();
    }
    res.original_pval = pval;
    res.n_gene        = n_tested + static_cast<int>(all_removed.size());
    res.n_gene_tested = n_tested;
    res.n_samp_g1     = static_cast<int>(object1_1.nrow());
    res.n_samp_g2     = static_cast<int>(object2_1.nrow());
    res.n_perm        = static_cast<int>(D_perm.size());

    output.results.push_back(res);
    output.n_samp_total = static_cast<int>(object.nrow());
    output.n_gene       = n_tested;

    // Write result TSV
    {
        DataFrame res_df;
        res_df.add_num_col("Original.TS",     {res.original_ts});
        res_df.add_num_col("Perm.TS.Mean",    {res.perm_ts_mean});
        res_df.add_num_col("Original.Pvalue", {res.original_pval});
        res_df.add_num_col("Ngene",           {static_cast<double>(res.n_gene)});
        res_df.add_num_col("Ngene.tested",    {static_cast<double>(res.n_gene_tested)});
        res_df.add_num_col("Nsamp.g1",        {static_cast<double>(res.n_samp_g1)});
        res_df.add_num_col("Nsamp.g2",        {static_cast<double>(res.n_samp_g2)});
        res_df.add_num_col("Nperm",           {static_cast<double>(res.n_perm)});
        write_tsv(res_df, cfg.resdir + "/" + cfg.projectname + "_result.tsv");
    }

    auto t_end = std::chrono::steady_clock::now();
    output.time_sec = std::chrono::duration<double>(t_end - t_start).count();
    std::cout << "Time taken: " << output.time_sec << " seconds" << std::endl;

    return output;
}

// =========================================================================
// GSAR_qr_ct
//
// Test differential co-expression across cell types within each
// disease condition.
// =========================================================================
GSCAOutput GSAR_qr_ct(const GSCAConfig& cfg,
                       DataFrame& gematrix,
                       const std::vector<int>& group,
                       int group_ref,
                       const std::vector<std::string>& genelist_in)
{
    auto t_start = std::chrono::steady_clock::now();
    GSCAOutput output;

#ifdef HAS_OPENMP
    omp_set_num_threads(cfg.num_cores);
#endif

    // Gene list handling
    std::vector<std::string> genelist = genelist_in;
    bool use_all = (genelist.empty() ||
                    (genelist.size() == 1 && genelist[0] == "NA"));
    if (use_all) {
        genelist = gematrix.gene_names(cfg.gene_start);
        output.genes_not_found = "-";
    } else {
        auto all_cols = gematrix.col_names();
        std::set<std::string> col_set(all_cols.begin(), all_cols.end());
        std::string not_found;
        for (auto& g : genelist)
            if (!col_set.count(g)) {
                if (!not_found.empty()) not_found += ",";
                not_found += g;
            }
        output.genes_not_found = not_found.empty() ? "-" : not_found;
    }

    // Subset columns
    std::vector<std::string> keep_cols;
    for (size_t c = 0; c < cfg.gene_start; ++c)
        keep_cols.push_back(gematrix.col_name(c));
    for (auto& g : genelist)
        if (gematrix.has_col(g)) keep_cols.push_back(g);

    DataFrame object = gematrix.col_subset(keep_cols);
    size_t gene_start = cfg.gene_start;
    object.round_cols(gene_start);

    gsnca_check_ct(cfg.projectname, object, group, group_ref, cfg.ct_var);

    // Check exactly 2 groups and 2 cell types
    std::set<int> grp_set(group.begin(), group.end());
    if (grp_set.size() != 2)
        throw std::runtime_error("GSAR_qr_ct requires exactly 2 disease groups.");

    auto ct_col = object.str_col(cfg.ct_var);
    auto ct_levels = unique_sorted(ct_col);
    if (ct_levels.size() != 2)
        throw std::runtime_error("GSAR_qr_ct supports exactly 2 cell types.");
    if (!cfg.ct_ref_level.empty()) {
        auto it = std::find(ct_levels.begin(), ct_levels.end(), cfg.ct_ref_level);
        if (it == ct_levels.end()) {
            throw std::runtime_error("Requested ct_ref_level '" + cfg.ct_ref_level +
                                     "' is not present in the CT subset.");
        }
        std::rotate(ct_levels.begin(), it, it + 1);
    }

    // Assign row names for alignment
    for (size_t i = 0; i < object.nrow(); ++i)
        object.row_names[i] = std::to_string(i + 1);

    // Determine group and cell-type levels
    std::vector<int> grp_levels(grp_set.begin(), grp_set.end());
    std::sort(grp_levels.begin(), grp_levels.end());
    int g1 = grp_levels[0], g2 = grp_levels[1];
    std::string ct1 = ct_levels[0], ct2 = ct_levels[1];

    // Split by group x cell type
    auto make_mask = [&](int grp, const std::string& ct) {
        std::vector<bool> mask(object.nrow(), false);
        for (size_t i = 0; i < object.nrow(); ++i)
            if (group[i] == grp && ct_col[i] == ct) mask[i] = true;
        return mask;
    };

    auto object1_1 = object.row_subset(make_mask(g1, ct1));
    auto object1_2 = object.row_subset(make_mask(g1, ct2));
    auto object2_1 = object.row_subset(make_mask(g2, ct1));
    auto object2_2 = object.row_subset(make_mask(g2, ct2));

    if (object1_1.nrow() == 0 || object1_2.nrow() == 0)
        throw std::runtime_error("Group " + std::to_string(g1) +
                                 " does not contain both cell types.");
    if (object2_1.nrow() == 0 || object2_2.nrow() == 0)
        throw std::runtime_error("Group " + std::to_string(g2) +
                                 " does not contain both cell types.");

    // Pool within disease group for permutation
    auto object_dc1 = DataFrame::rbind(object1_1, object1_2);
    auto object_dc2 = DataFrame::rbind(object2_1, object2_2);

    // Add cid
    auto add_cid = [](DataFrame& df) {
        NumCol cid(df.nrow());
        std::iota(cid.begin(), cid.end(), 1.0);
        df = DataFrame::cbind_num("cid", cid, df);
    };
    add_cid(object1_1);
    add_cid(object1_2);
    add_cid(object2_1);
    add_cid(object2_2);
    add_cid(object_dc1);
    add_cid(object_dc2);

    auto object_ct_perm = DataFrame::rbind(object_dc1, object_dc2);
    std::vector<int> group_ct_perm(object_ct_perm.nrow(), g2);
    for (size_t i = 0; i < object_dc1.nrow(); ++i)
        group_ct_perm[i] = g1;

    gene_start += 1;  // shifted by cid

    // Cell-type vectors for permutation objects
    auto ct_ct_perm_col = object_ct_perm.str_col(cfg.ct_var);

    // SD check
    auto genes = object.gene_names(cfg.gene_start);
    std::vector<std::string> genes_removed1;
    if (cfg.check_sd) {
        std::cout << "Checking gene SD..." << std::endl;
        auto r1a = gsnca_sd(object1_1, object1_2, genes, gene_start, cfg.min_sd);
        auto r1b = gsnca_sd(object2_1, object2_2, genes, gene_start, cfg.min_sd);
        std::set<std::string> s(r1a.begin(), r1a.end());
        s.insert(r1b.begin(), r1b.end());
        genes_removed1.assign(s.begin(), s.end());

        remove_genes(object, genes_removed1);
        remove_genes(object1_1, genes_removed1);
        remove_genes(object1_2, genes_removed1);
        remove_genes(object2_1, genes_removed1);
        remove_genes(object2_2, genes_removed1);
        remove_genes(object_dc1, genes_removed1);
        remove_genes(object_dc2, genes_removed1);
        remove_genes(object_ct_perm, genes_removed1);
    }

    // Low-count filter
    std::cout << "Filtering low-count genes..." << std::endl;
    genes = object.gene_names(cfg.gene_start);
    auto r2a = gsnca_filt(object1_1, object1_2, genes, gene_start, cfg.filter_count);
    auto r2b = gsnca_filt(object2_1, object2_2, genes, gene_start, cfg.filter_count);
    std::set<std::string> s2(r2a.begin(), r2a.end());
    s2.insert(r2b.begin(), r2b.end());
    std::vector<std::string> genes_removed2(s2.begin(), s2.end());

    remove_genes(object, genes_removed2);
    remove_genes(object1_1, genes_removed2);
    remove_genes(object1_2, genes_removed2);
    remove_genes(object2_1, genes_removed2);
    remove_genes(object2_2, genes_removed2);
    remove_genes(object_dc1, genes_removed2);
    remove_genes(object_dc2, genes_removed2);
    remove_genes(object_ct_perm, genes_removed2);

    // Merged removed genes
    std::set<std::string> all_removed(genes_removed1.begin(), genes_removed1.end());
    all_removed.insert(genes_removed2.begin(), genes_removed2.end());
    output.genes_removed.assign(all_removed.begin(), all_removed.end());

    // --- Compute ZINB quantile residuals per cell type per disease group ---
    std::cout << "Computing ZINB quantile residuals (observed data, per cell type)..."
              << std::endl;

    auto drop_cid = [](const DataFrame& df) {
        auto names = df.col_names();
        names.erase(names.begin()); // remove "cid"
        return df.col_subset(names);
    };

    // Disease group 1
    auto qr_ct1_g1 = zinb_qr_ct(drop_cid(object1_1), gene_start - 1,
                                  cfg.subj_var, cfg.na_counts, cfg.n_sim,
                                  cfg.seed, false);
    qr_ct1_g1.row_names = object1_1.row_names;

    auto qr_ct2_g1 = zinb_qr_ct(drop_cid(object1_2), gene_start - 1,
                                  cfg.subj_var, cfg.na_counts, cfg.n_sim,
                                  cfg.seed + 50000ULL, false);
    qr_ct2_g1.row_names = object1_2.row_names;

    auto qr_mat_dc1 = QRMatrix::rbind(qr_ct1_g1, qr_ct2_g1);

    // Disease group 2
    auto qr_ct1_g2 = zinb_qr_ct(drop_cid(object2_1), gene_start - 1,
                                  cfg.subj_var, cfg.na_counts, cfg.n_sim,
                                  cfg.seed, false);
    qr_ct1_g2.row_names = object2_1.row_names;

    auto qr_ct2_g2 = zinb_qr_ct(drop_cid(object2_2), gene_start - 1,
                                  cfg.subj_var, cfg.na_counts, cfg.n_sim,
                                  cfg.seed + 50000ULL, false);
    qr_ct2_g2.row_names = object2_2.row_names;

    auto qr_mat_dc2 = QRMatrix::rbind(qr_ct1_g2, qr_ct2_g2);

    std::cout << "Quantile residuals computed." << std::endl;

    // --- Observed test statistics per disease group ---
    std::string perm_corr_dir;
    if (cfg.save_ct_perm_corr) {
        perm_corr_dir = cfg.resdir + "/" + cfg.projectname + "_ct_perm_corr";
        std::filesystem::create_directories(perm_corr_dir);
    }

    auto obs_stat_g1 = gsnca_stat_qr_detail(object1_1, object1_2, gene_start,
                                            cfg.na_counts, cfg.cor_method,
                                            cfg.subj_var, qr_mat_dc1);
    double D_obs_dc1 = obs_stat_g1.D;
    std::cout << "Observed sample completed for group " << g1 << "." << std::endl;

    auto obs_stat_g2 = gsnca_stat_qr_detail(object2_1, object2_2, gene_start,
                                            cfg.na_counts, cfg.cor_method,
                                            cfg.subj_var, qr_mat_dc2);
    double D_obs_dc2 = obs_stat_g2.D;
    std::cout << "Observed sample completed for group " << g2 << "." << std::endl;

    if (cfg.save_ct_perm_corr) {
        const std::string ct1_token = sanitize_token_meta(ct1);
        const std::string ct2_token = sanitize_token_meta(ct2);
        const auto& genes_now = qr_mat_dc1.gene_names;

        write_named_matrix_tsv_meta(
            obs_stat_g1.cor_mat1, genes_now,
            perm_corr_dir + "/observed.group_" + std::to_string(g1) + "." +
            ct1_token + "_corr.tsv");
        write_named_matrix_tsv_meta(
            obs_stat_g1.cor_mat2, genes_now,
            perm_corr_dir + "/observed.group_" + std::to_string(g1) + "." +
            ct2_token + "_corr.tsv");
        write_named_matrix_tsv_meta(
            obs_stat_g2.cor_mat1, genes_now,
            perm_corr_dir + "/observed.group_" + std::to_string(g2) + "." +
            ct1_token + "_corr.tsv");
        write_named_matrix_tsv_meta(
            obs_stat_g2.cor_mat2, genes_now,
            perm_corr_dir + "/observed.group_" + std::to_string(g2) + "." +
            ct2_token + "_corr.tsv");
    }

    // --- Permutation tests ---
    int group_num_1 = 0;
    {
        const auto& subj_col = object_dc1.str_col(cfg.subj_var);
        std::set<std::string> subj_set(subj_col.begin(), subj_col.end());
        group_num_1 = static_cast<int>(subj_set.size());
    }

    auto D_perm = gsnca_perm_ct_qr(
        cfg.nperm, object_ct_perm, gene_start, group_ct_perm, ct_ct_perm_col, ct1,
        cfg.check_sd, cfg.min_sd, cfg.max_skip,
        cfg.filter_count, cfg.na_counts, cfg.cor_method, cfg.subj_var,
        group_num_1, cfg.num_cores, 42, cfg.seed, cfg.n_sim, cfg.paired,
        cfg.save_ct_perm_corr, perm_corr_dir);

    double pval_dc1 = 0.0;
    if (!D_perm.group1.empty()) {
        int cnt = 0;
        for (double d : D_perm.group1) if (d >= D_obs_dc1) ++cnt;
        pval_dc1 = static_cast<double>(cnt + 1) /
                   static_cast<double>(D_perm.group1.size() + 1);
    }

    double pval_dc2 = 0.0;
    if (!D_perm.group2.empty()) {
        int cnt = 0;
        for (double d : D_perm.group2) if (d >= D_obs_dc2) ++cnt;
        pval_dc2 = static_cast<double>(cnt + 1) /
                   static_cast<double>(D_perm.group2.size() + 1);
    }

    const double observed_signed_delta = D_obs_dc2 - D_obs_dc1;
    const double observed_abs_delta = std::abs(observed_signed_delta);
    const double perm_signed_delta_mean = mean_or_zero(D_perm.signed_delta);
    const double perm_abs_delta_mean = mean_or_zero(D_perm.abs_delta);

    const double pval_signed_delta = ct_right_tail_pvalue(
        observed_signed_delta, D_perm.signed_delta);

    double pval_abs_delta = 0.0;
    if (!D_perm.abs_delta.empty()) {
        int cnt = 0;
        for (double d : D_perm.abs_delta) {
            if (d >= observed_abs_delta)
                ++cnt;
        }
        pval_abs_delta = static_cast<double>(cnt + 1) /
                         static_cast<double>(D_perm.abs_delta.size() + 1);
    }

    // --- Build results ---
    int n_tested = static_cast<int>(object.ncol() - cfg.gene_start);
    int n_gene_total = n_tested + static_cast<int>(all_removed.size());

    auto make_result = [&](double D_obs, const std::vector<double>& D_perm,
                           double pval, int n_ct1, int n_ct2) {
        GSCAResult r;
        r.original_ts   = D_obs;
        r.perm_ts_mean  = 0.0;
        if (!D_perm.empty()) {
            double s = 0.0;
            for (double d : D_perm) s += d;
            r.perm_ts_mean = s / D_perm.size();
        }
        r.original_pval = pval;
        r.n_gene        = n_gene_total;
        r.n_gene_tested = n_tested;
        r.n_samp_g1     = n_ct1;
        r.n_samp_g2     = n_ct2;
        r.n_perm        = static_cast<int>(D_perm.size());
        return r;
    };

    output.results.push_back(
        make_result(D_obs_dc1, D_perm.group1, pval_dc1,
                    static_cast<int>(object1_1.nrow()),
                    static_cast<int>(object1_2.nrow())));
    output.results.push_back(
        make_result(D_obs_dc2, D_perm.group2, pval_dc2,
                    static_cast<int>(object2_1.nrow()),
                    static_cast<int>(object2_2.nrow())));

    // Preserve AbsDelta and record the fixed right-tailed signed test.
    {
        DataFrame delta_df;
        delta_df.add_num_col("Group1", {static_cast<double>(g1)});
        delta_df.add_num_col("Group2", {static_cast<double>(g2)});
        delta_df.add_num_col("Observed.TS.Group1", {D_obs_dc1});
        delta_df.add_num_col("Observed.TS.Group2", {D_obs_dc2});
        delta_df.add_num_col("Observed.Delta", {observed_signed_delta});
        delta_df.add_num_col("Observed.AbsDelta", {observed_abs_delta});
        delta_df.add_num_col("Perm.Delta.Mean", {perm_signed_delta_mean});
        delta_df.add_num_col("Perm.AbsDelta.Mean", {perm_abs_delta_mean});
        delta_df.add_num_col("Permutation.Pvalue.Delta", {pval_signed_delta});
        delta_df.add_num_col("Permutation.Pvalue.AbsDelta", {pval_abs_delta});
        delta_df.add_num_col("Nperm.Paired", {static_cast<double>(D_perm.abs_delta.size())});
        delta_df.add_str_col("CT.Tail", {"right"});
        delta_df.add_str_col("Delta.Definition", {"D(Group2)-D(Group1)"});
        write_tsv(delta_df, cfg.resdir + "/" + cfg.projectname + "_delta_result.tsv");
    }

    output.n_samp_total = static_cast<int>(object.nrow());
    output.n_gene       = n_tested;

    auto t_end = std::chrono::steady_clock::now();
    output.time_sec = std::chrono::duration<double>(t_end - t_start).count();
    std::cout << "Time taken: " << output.time_sec << " seconds" << std::endl;

    return output;
}

} // namespace gsca
