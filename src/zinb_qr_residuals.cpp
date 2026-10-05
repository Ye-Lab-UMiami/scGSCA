// ============================================================================
// zinb_qr_residuals.cpp
//
// Wrapper functions that compute ZINB quantile residual matrices for
// DC and CT modes.  Each gene is fitted independently (serial within
// a single call; parallelism lives at the permutation layer).
//
// Port of: zinb_qr_residuals_parv6.R
//
// DC mode (zinb_qr_dc_sep):
//   Model per gene: y ~ ct_factor + (1 | subj_factor), zi ~ ct_factor
//   Called separately for each disease group.
//   NOT permutation-invariant -> must recompute per permutation.
//
// CT mode (zinb_qr_ct):
//   Model per gene: y ~ 1 + (1 | subj_factor), zi ~ 1
//   Called separately for each cell type within a disease group.
//   NOT permutation-invariant -> must recompute per permutation.
// ============================================================================
#include "MetaGSCA.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <set>
#include <atomic>
#include <iostream>

#ifdef HAS_OPENMP
#include <omp.h>
#endif

// Inverse normal CDF (probit).  Uses the Beasley-Springer-Moro algorithm.
static double qnorm(double p) {
    // Rational approximation for the inverse normal CDF
    if (p <= 0.0) return -1e15;
    if (p >= 1.0) return  1e15;
    if (p == 0.5) return  0.0;

    // Coefficients
    static const double a[] = {
        -3.969683028665376e+01, 2.209460984245205e+02,
        -2.759285104469687e+02, 1.383577518672690e+02,
        -3.066479806614716e+01, 2.506628277459239e+00
    };
    static const double b[] = {
        -5.447609879822406e+01, 1.615858368580409e+02,
        -1.556989798598866e+02, 6.680131188771972e+01,
        -1.328068155288572e+01
    };
    static const double c[] = {
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00,  2.938163982698783e+00
    };
    static const double d[] = {
         7.784695709041462e-03, 3.224671290700398e-01,
         2.445134137142996e+00, 3.754408661907416e+00
    };

    double q, r;
    if (p < 0.02425) {
        q = std::sqrt(-2.0 * std::log(p));
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    } else if (p <= 0.97575) {
        q = p - 0.5;
        r = q * q;
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5]) * q /
               (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
    } else {
        q = std::sqrt(-2.0 * std::log(1.0 - p));
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                 ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    }
}

namespace gsca {

// =========================================================================
// QRMatrix methods
// =========================================================================
QRMatrix QRMatrix::row_subset(const std::vector<std::string>& names) const {
    // Build row-name -> index lookup
    std::unordered_map<std::string, size_t> rn_map;
    for (size_t i = 0; i < row_names.size(); ++i)
        rn_map[row_names[i]] = i;

    QRMatrix out;
    out.gene_names = gene_names;
    out.row_names = names;
    out.data.resize(static_cast<Eigen::Index>(names.size()),
                    static_cast<Eigen::Index>(gene_names.size()));

    for (size_t i = 0; i < names.size(); ++i) {
        auto it = rn_map.find(names[i]);
        if (it == rn_map.end())
            throw std::runtime_error(
                "QRMatrix::row_subset: row name '" + names[i] + "' not found.");
        out.data.row(static_cast<Eigen::Index>(i)) =
            data.row(static_cast<Eigen::Index>(it->second));
    }
    return out;
}

QRMatrix QRMatrix::rbind(const QRMatrix& a, const QRMatrix& b) {
    if (a.gene_names != b.gene_names)
        throw std::runtime_error("QRMatrix::rbind: gene names mismatch.");

    QRMatrix out;
    out.gene_names = a.gene_names;
    out.row_names.reserve(a.row_names.size() + b.row_names.size());
    out.row_names.insert(out.row_names.end(),
                         a.row_names.begin(), a.row_names.end());
    out.row_names.insert(out.row_names.end(),
                         b.row_names.begin(), b.row_names.end());

    out.data.resize(a.data.rows() + b.data.rows(), a.data.cols());
    out.data.topRows(a.data.rows()) = a.data;
    out.data.bottomRows(b.data.rows()) = b.data;
    return out;
}

// =========================================================================
// Helper: encode a string factor column as 0-based integer IDs.
// Returns the mapping and the encoded vector.
// =========================================================================
static void encode_factor(const std::vector<std::string>& vals,
                          std::vector<int>& ids,
                          std::vector<std::string>& levels) {
    std::set<std::string> uq(vals.begin(), vals.end());
    levels.assign(uq.begin(), uq.end());
    std::sort(levels.begin(), levels.end());

    std::unordered_map<std::string, int> lvl_map;
    for (size_t i = 0; i < levels.size(); ++i)
        lvl_map[levels[i]] = static_cast<int>(i);

    ids.resize(vals.size());
    for (size_t i = 0; i < vals.size(); ++i)
        ids[i] = lvl_map[vals[i]];
}

// =========================================================================
// Helper: compute sort order (argsort) by factor levels.
// =========================================================================
static std::vector<size_t> argsort_factors(const std::vector<int>& primary,
                                           const std::vector<int>& secondary) {
    size_t n = primary.size();
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
        if (primary[a] != primary[b]) return primary[a] < primary[b];
        return secondary[a] < secondary[b];
    });
    return idx;
}

static std::vector<size_t> argsort_single(const std::vector<int>& v) {
    size_t n = v.size();
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
        return v[a] < v[b];
    });
    return idx;
}

// =========================================================================
// zinb_qr_dc_sep
//
// DC mode quantile residuals.
// For each gene:
//   1. Fit ZINB: y ~ ct_factor + (1 | subj_factor), zi ~ ct_factor
//   2. Simulate residuals (DHARMa-style)
//   3. Transform to normal quantile residuals: qnorm(u)
// =========================================================================
QRMatrix zinb_qr_dc_sep(const DataFrame& objt,
                         size_t gene_start,
                         const std::string& subj_var,
                         const std::string& ct_var,
                         bool na_counts,
                         int n_sim,
                         uint64_t seed,
                         bool verbose,
                         bool parallel)
{
    size_t n_cells = objt.nrow();
    auto genes = objt.gene_names(gene_start);
    size_t n_genes = genes.size();

    // Encode factors
    std::vector<int> subj_ids, ct_ids;
    std::vector<std::string> subj_levels, ct_levels;
    encode_factor(objt.str_col(subj_var), subj_ids, subj_levels);
    encode_factor(objt.str_col(ct_var), ct_ids, ct_levels);
    int n_subj = static_cast<int>(subj_levels.size());
    int n_ct = static_cast<int>(ct_levels.size());

    // Sort order: by ct then subj
    auto sort_ord = argsort_factors(ct_ids, subj_ids);

    // Sorted factor arrays
    std::vector<int> subj_sorted(n_cells), ct_sorted(n_cells);
    for (size_t i = 0; i < n_cells; ++i) {
        subj_sorted[i] = subj_ids[sort_ord[i]];
        ct_sorted[i]   = ct_ids[sort_ord[i]];
    }

    QRMatrix qr;
    qr.data.resize(static_cast<Eigen::Index>(n_cells),
                   static_cast<Eigen::Index>(n_genes));
    qr.gene_names = genes;
    qr.row_names = objt.row_names;

    std::atomic<int> fail_count{0};

    // Preload all gene columns as raw pointers for thread-safe access
    std::vector<const NumCol*> gene_cols(n_genes);
    for (size_t g = 0; g < n_genes; ++g)
        gene_cols[g] = &objt.num_col(genes[g]);

    std::cout << "  ZINB DC fitting " << n_genes << " genes";
#ifdef HAS_OPENMP
    if (parallel) std::cout << " (parallel, " << omp_get_max_threads() << " threads)";
    else std::cout << " (serial, inside permutation)";
#endif
    std::cout << "..." << std::endl;

#ifdef HAS_OPENMP
    #pragma omp parallel for schedule(dynamic, 1) if(parallel)
#endif
    for (size_t g = 0; g < n_genes; ++g) {
        // Get gene counts in sorted order
        const auto& gene_col = *gene_cols[g];
        std::vector<int> y_sorted(n_cells);
        for (size_t i = 0; i < n_cells; ++i) {
            double val = gene_col[sort_ord[i]];
            if (std::isnan(val) && na_counts) val = 0.0;
            y_sorted[i] = static_cast<int>(std::round(val));
            if (y_sorted[i] < 0) y_sorted[i] = 0;
        }

        // Fit ZINB model
        ZINBFit fit = zinb_fit("dc", y_sorted, subj_sorted, ct_sorted,
                               n_subj, n_ct);

        if (!fit.converged) {
            ++fail_count;
            #pragma omp critical
            std::cerr << "  [FAIL] Observed ZINB fit failed: gene=" << genes[g] << std::endl;
            qr.data.col(static_cast<Eigen::Index>(g)).setZero();
            continue;
        }

        // Simulate residuals
        auto u = zinb_simulate_residuals(fit, y_sorted, subj_sorted, ct_sorted,
                                         n_subj, n_sim, seed + g);

        // Transform to normal quantile residuals, un-sort
        for (size_t i = 0; i < n_cells; ++i) {
            double qr_val = qnorm(u[i]);
            qr.data(static_cast<Eigen::Index>(sort_ord[i]),
                    static_cast<Eigen::Index>(g)) = qr_val;
        }

        if (parallel) {
#ifdef HAS_OPENMP
            #pragma omp critical
#endif
            {
                std::cout << "  Gene " << (g + 1) << "/" << n_genes
                          << " (" << genes[g] << ") done" << std::endl;
            }
        }
    }

    if (verbose && fail_count.load() > 0) {
        std::cerr << "[zinb_qr_dc_sep] " << fail_count.load() << "/" << n_genes
                  << " genes returned zero residual vectors due to "
                  << "fit/simulation failure." << std::endl;
    }

    return qr;
}

// =========================================================================
// zinb_qr_ct
//
// CT mode quantile residuals.
// For each gene:
//   1. Fit ZINB: y ~ 1 + (1 | subj_factor), zi ~ 1
//   2. Simulate residuals (DHARMa-style)
//   3. Transform: qnorm(u)
// =========================================================================
QRMatrix zinb_qr_ct(const DataFrame& objt,
                     size_t gene_start,
                     const std::string& subj_var,
                     bool na_counts,
                     int n_sim,
                     uint64_t seed,
                     bool verbose,
                     bool parallel)
{
    size_t n_cells = objt.nrow();
    auto genes = objt.gene_names(gene_start);
    size_t n_genes = genes.size();

    // Encode subject factor
    std::vector<int> subj_ids;
    std::vector<std::string> subj_levels;
    encode_factor(objt.str_col(subj_var), subj_ids, subj_levels);
    int n_subj = static_cast<int>(subj_levels.size());

    // Sort order by subject
    auto sort_ord = argsort_single(subj_ids);

    std::vector<int> subj_sorted(n_cells);
    for (size_t i = 0; i < n_cells; ++i)
        subj_sorted[i] = subj_ids[sort_ord[i]];

    // Dummy ct_ids (unused in CT mode)
    std::vector<int> ct_dummy(n_cells, 0);

    QRMatrix qr;
    qr.data.resize(static_cast<Eigen::Index>(n_cells),
                   static_cast<Eigen::Index>(n_genes));
    qr.gene_names = genes;
    qr.row_names = objt.row_names;

    std::atomic<int> fail_count{0};

    // Preload gene column pointers for thread-safe access
    std::vector<const NumCol*> gene_cols_ptr(n_genes);
    for (size_t g = 0; g < n_genes; ++g)
        gene_cols_ptr[g] = &objt.num_col(genes[g]);

    std::cout << "  ZINB CT fitting " << n_genes << " genes";
#ifdef HAS_OPENMP
    if (parallel) std::cout << " (parallel, " << omp_get_max_threads() << " threads)";
    else std::cout << " (serial, inside permutation)";
#endif
    std::cout << "..." << std::endl;

#ifdef HAS_OPENMP
    #pragma omp parallel for schedule(dynamic, 1) if(parallel)
#endif
    for (size_t g = 0; g < n_genes; ++g) {
        const auto& gene_col = *gene_cols_ptr[g];
        std::vector<int> y_sorted(n_cells);
        for (size_t i = 0; i < n_cells; ++i) {
            double val = gene_col[sort_ord[i]];
            if (std::isnan(val) && na_counts) val = 0.0;
            y_sorted[i] = static_cast<int>(std::round(val));
            if (y_sorted[i] < 0) y_sorted[i] = 0;
        }

        ZINBFit fit = zinb_fit("ct", y_sorted, subj_sorted, ct_dummy,
                               n_subj, 1);

        if (!fit.converged) {
            ++fail_count;
            qr.data.col(static_cast<Eigen::Index>(g)).setZero();
            continue;
        }

        auto u = zinb_simulate_residuals(fit, y_sorted, subj_sorted, ct_dummy,
                                         n_subj, n_sim, seed + g);

        for (size_t i = 0; i < n_cells; ++i) {
            double qr_val = qnorm(u[i]);
            qr.data(static_cast<Eigen::Index>(sort_ord[i]),
                    static_cast<Eigen::Index>(g)) = qr_val;
        }

        if (parallel) {
#ifdef HAS_OPENMP
            #pragma omp critical
#endif
            {
                std::cout << "  Gene " << (g + 1) << "/" << n_genes
                          << " (" << genes[g] << ") done" << std::endl;
            }
        }
    }

    if (verbose && fail_count.load() > 0) {
        std::cerr << "[zinb_qr_ct] " << fail_count.load() << "/" << n_genes
                  << " genes returned zero residual vectors due to "
                  << "fit/simulation failure." << std::endl;
    }

    return qr;
}

} // namespace gsca
