// ============================================================================
// MetaGSCA.h
//
// Master header for the MetaGSCA C++ pipeline.
//
// This is a C++ port of the R package originally developed by
// Zhuoli Jin and Rebecca Irlmeier (University of Miami).
//
// The pipeline tests differential co-expression of gene sets using
// ZINB (Zero-Inflated Negative Binomial) quantile residuals with
// permutation-based inference.
//
// Two analysis modes:
//   DC mode — differential co-expression across disease conditions
//             (pooling cell types)
//   CT mode — differential co-expression across cell types
//             within a disease condition
//
// Dependencies: Eigen (linear algebra), NLopt (optimization), OpenMP (optional)
// ============================================================================
#ifndef METAGSCA_H
#define METAGSCA_H

#include "DataFrame.h"
#include "rng.h"
#include "tsv_io.h"

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <functional>
#include <iostream>
#include <chrono>

namespace gsca {

// =========================================================================
// Common types
// =========================================================================
using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;

// =========================================================================
// Configuration struct for a single analysis run
// =========================================================================
struct GSCAConfig {
    std::string projectname;
    std::string pathname;
    std::string resdir = ".";

    // Gene columns start at this 0-based index in the DataFrame
    size_t gene_start = 3;

    // Group / factor variable names
    std::string ct_var  = "Celltype";
    std::string subj_var = "Samples";
    std::string ct_ref_level;

    // Filtering
    bool   check_sd     = true;
    double min_sd       = 0.001;
    double filter_count = 0.0;
    bool   na_counts    = true;

    // Correlation
    std::string cor_method = "pearson";  // only pearson is implemented

    // Permutation
    int  nperm    = 50;
    int  max_skip = 50;
    int  num_cores = 1;
    bool paired = false;

    // ZINB residual simulation
    int  n_sim = 100;

    // Significance level (for reporting)
    double level = 0.05;

    // Diagnostic plots (not implemented in C++; placeholder)
    bool diagnostic_plots = false;
    bool save_ct_perm_corr = false;

    // Random seed
    uint64_t seed = 42;
};

// =========================================================================
// Result struct for one group comparison
// =========================================================================
struct GSCAResult {
    double original_ts    = 0.0;
    double perm_ts_mean   = 0.0;
    double original_pval  = 0.0;
    int    n_gene         = 0;
    int    n_gene_tested  = 0;
    int    n_samp_g1      = 0;
    int    n_samp_g2      = 0;
    int    n_perm         = 0;
};

// =========================================================================
// Full output struct
// =========================================================================
struct GSCAOutput {
    std::vector<GSCAResult> results;   // one per group (DC) or two (CT)
    std::vector<std::string> genes_removed;
    std::string genes_not_found;
    int    n_samp_total = 0;
    int    n_gene       = 0;
    double time_sec     = 0.0;
};

// =========================================================================
// Input validation  (gsnca_check_dc.R / gsnca_check_ct.R)
// =========================================================================
void gsnca_check_dc(const std::string& projectname,
                    const DataFrame& object,
                    const std::vector<int>& group,
                    int group_ref,
                    const std::string& ct_var);

void gsnca_check_ct(const std::string& projectname,
                    const DataFrame& object,
                    const std::vector<int>& group,
                    int group_ref,
                    const std::string& ct_var);

// =========================================================================
// Gene filtering  (gsnca_sd.R / gsnca_filt.R)
// =========================================================================

// Return names of genes to remove because per-gene SD < min_sd in either group.
std::vector<std::string> gsnca_sd(const DataFrame& objt1,
                                  const DataFrame& objt2,
                                  const std::vector<std::string>& genes,
                                  size_t gene_start,
                                  double min_sd);

// Return names of genes to remove because column sum <= filter_count.
std::vector<std::string> gsnca_filt(const DataFrame& objt1,
                                    const DataFrame& objt2,
                                    const std::vector<std::string>& genes,
                                    size_t gene_start,
                                    double filter_count);

// =========================================================================
// ZINB quantile residuals  (zinb_qr_residuals_parv6.R)
// =========================================================================

// QR matrix: rows = cells, columns = genes.
// Stored as Eigen::MatrixXd with gene_names and row_names for alignment.
struct QRMatrix {
    Mat data;                              // n_cells x n_genes
    std::vector<std::string> gene_names;   // column labels
    std::vector<std::string> row_names;    // row labels for alignment

    size_t nrow() const { return static_cast<size_t>(data.rows()); }
    size_t ncol() const { return static_cast<size_t>(data.cols()); }

    // Subset rows by row_names matching the given names
    QRMatrix row_subset(const std::vector<std::string>& names) const;

    // Vertical concatenation
    static QRMatrix rbind(const QRMatrix& a, const QRMatrix& b);
};

// DC mode: y ~ ct_factor + (1 | subj_factor), zi ~ ct_factor
// Called once per disease group.
// parallel: if true, parallelize gene fitting with OpenMP (use for observed data).
//           if false, run serially (use when called inside permutation-level parallelism).
QRMatrix zinb_qr_dc_sep(const DataFrame& objt,
                         size_t gene_start,
                         const std::string& subj_var,
                         const std::string& ct_var,
                         bool na_counts = true,
                         int n_sim = 100,
                         uint64_t seed = 42,
                         bool verbose = false,
                         bool parallel = true);

// CT mode: y ~ 1 + (1 | subj_factor), zi ~ 1
// Called once per cell type within a disease group.
// parallel: if true, parallelize gene fitting with OpenMP.
QRMatrix zinb_qr_ct(const DataFrame& objt,
                     size_t gene_start,
                     const std::string& subj_var,
                     bool na_counts = true,
                     int n_sim = 100,
                     uint64_t seed = 42,
                     bool verbose = false,
                     bool parallel = true);

// =========================================================================
// ZINB model fitting  (zinb_model.h)
// =========================================================================

// Parameters for a fitted ZINB model (simplified for the pipeline needs).
struct ZINBFit {
    bool converged = false;
    int nlopt_status = 0;
    double objective = 0.0;
    int initial_clamped = 0;
    int boundary_hits = 0;
    int newton_warn_subjects = 0;
    int hessian_clamp_subjects = 0;
    double max_abs_grad = 0.0;
    double min_neg_h = 0.0;

    // Conditional model: log(mu_i) = X_cond * beta_cond + Z * b_i
    Vec beta_cond;        // fixed effects (conditional)
    double log_theta = 0; // log(dispersion)
    double sigma_re  = 0; // random-effect SD

    // Zero-inflation model: logit(pi_i) = X_zi * beta_zi
    Vec beta_zi;          // fixed effects (ZI)

    // Derived quantities needed for simulation
    double theta() const { return std::exp(log_theta); }
};

// Fit a ZINB GLMM via marginal likelihood (Laplace approximation + NLopt).
//
// model_type:
//   "dc" -> y ~ ct + (1|subj), zi ~ ct
//   "ct" -> y ~ 1  + (1|subj), zi ~ 1
//
// y            : response vector (counts)
// subj_ids     : integer-coded subject IDs (0-based)
// ct_ids       : integer-coded cell-type IDs (0-based), unused for "ct" mode
// n_subj       : number of unique subjects
// n_ct         : number of unique cell types (unused for "ct" mode)
ZINBFit zinb_fit(const std::string& model_type,
                 const std::vector<int>& y,
                 const std::vector<int>& subj_ids,
                 const std::vector<int>& ct_ids,
                 int n_subj,
                 int n_ct);

// Simulate quantile residuals from a fitted ZINB model (DHARMa-style).
// Returns a vector of scaled residuals in [0,1].
std::vector<double> zinb_simulate_residuals(const ZINBFit& fit,
                                            const std::vector<int>& y,
                                            const std::vector<int>& subj_ids,
                                            const std::vector<int>& ct_ids,
                                            int n_subj,
                                            int n_sim,
                                            uint64_t seed);

// =========================================================================
// GSNCA test statistic  (gsnca_stat_qr.R)
// =========================================================================

// Compute the GSNCA test statistic D from precomputed quantile residuals.
// objt1, objt2: DataFrames with cid column.
// qr_mat: precomputed quantile residual matrix (rows aligned by row_names).
struct GSNCAStatDetail {
    double D = 0.0;
    Mat cor_mat1;
    Mat cor_mat2;
};

GSNCAStatDetail gsnca_stat_qr_detail(const DataFrame& objt1,
                                     const DataFrame& objt2,
                                     size_t gene_start,
                                     bool na_counts,
                                     const std::string& cor_method,
                                     const std::string& subj_var,
                                     const QRMatrix& qr_mat);

double gsnca_stat_qr(const DataFrame& objt1,
                     const DataFrame& objt2,
                     size_t gene_start,
                     bool na_counts,
                     const std::string& cor_method,
                     const std::string& subj_var,
                     const QRMatrix& qr_mat);

// =========================================================================
// Permutation tests (gsnca_perm_dc_qr_parv7.R / gsnca_perm_ct_qr_parv8.R)
// =========================================================================

// DC mode permutation
std::vector<double> gsnca_perm_dc_qr(
    const std::vector<std::string>& domain,    // unique subject IDs
    int nperm,
    const DataFrame& object,                   // pooled data (no cid)
    size_t gene_start,
    int group_num_1,                           // # subjects in group 1
    bool check_sd, double min_sd, int max_skip,
    double filter_count, bool na_counts,
    const std::string& cor_method,
    const std::string& subj_var,
    const std::string& ct_var,
    int num_cores, uint64_t perm_seed,
    uint64_t seed, int n_sim, bool paired = false);

struct CTPermResult {
    std::vector<double> group1;
    std::vector<double> group2;
    std::vector<double> signed_delta;
    std::vector<double> abs_delta;
};

// CT mode permutation
CTPermResult gsnca_perm_ct_qr(
    int nperm,
    const DataFrame& object,                   // pooled data across both disease groups
    size_t gene_start,
    const std::vector<int>& group,             // disease group per row
    const std::vector<std::string>& cell_type, // cell type per row
    const std::string& preserve_cell_type,     // preserve this cell-type count within subject
    bool check_sd, double min_sd, int max_skip,
    double filter_count, bool na_counts,
    const std::string& cor_method,
    const std::string& subj_var,
    int group_num_1,                           // # subjects in permuted group 1
    int num_cores, uint64_t perm_seed,
    uint64_t seed, int n_sim, bool paired = false,
    bool save_cor_mats = false,
    const std::string& cor_mat_dir = "");

// =========================================================================
// Main orchestration functions (MetaGSCA_qr_parv5.R)
// =========================================================================

// DC mode: test differential co-expression across disease conditions
GSCAOutput GSAR_qr_dc(const GSCAConfig& cfg,
                       DataFrame& gematrix,
                       const std::vector<int>& group,
                       int group_ref,
                       const std::vector<std::string>& genelist);

// CT mode: test differential co-expression across cell types
GSCAOutput GSAR_qr_ct(const GSCAConfig& cfg,
                       DataFrame& gematrix,
                       const std::vector<int>& group,
                       int group_ref,
                       const std::vector<std::string>& genelist);

} // namespace gsca

#endif // METAGSCA_H
