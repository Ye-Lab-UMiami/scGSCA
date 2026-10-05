// ============================================================================
// gsnca_stat_qr.cpp
//
// Compute the GSNCA test statistic D from precomputed quantile residuals.
//
// Algorithm:
//   1. For each condition, extract quantile residual submatrix (aligned by
//      row names).
//   2. Compute the absolute Pearson correlation matrix (q x q) on residuals.
//   3. Compute leading eigenvectors of both correlation matrices.
//   4. D = sum(|p1_i * ||p1|| - p2_i * ||p2|| |)
//      where p1, p2 are the absolute leading eigenvectors.
//
// Port of: gsnca_stat_qr.R
// ============================================================================
#include "MetaGSCA.h"
#include <cmath>
#include <algorithm>
#include <Eigen/Eigenvalues>

namespace gsca {

// ---------------------------------------------------------------------------
// Build the absolute correlation matrix from a QRMatrix subset.
//
// Given a DataFrame objt (with row_names) and a QRMatrix, extract the
// rows corresponding to objt, compute pairwise Pearson correlation on
// columns, take absolute values, set diagonal to 1, NaN to 0.
//
// This mirrors the active .build_cor_mat() in gsnca_stat_qr.R (line 111-139).
// ---------------------------------------------------------------------------
static Mat build_cor_mat(const DataFrame& objt,
                         const QRMatrix& qr_mat,
                         size_t /* gene_start */,
                         const std::string& /* cor_method */,
                         const std::string& /* subj_var */)
{
    // Align qr_mat rows to objt row_names
    const auto& rn = objt.row_names;
    if (rn.empty())
        throw std::runtime_error("build_cor_mat: objt has no row_names.");
    if (qr_mat.row_names.empty())
        throw std::runtime_error("build_cor_mat: qr_mat has no row_names.");

    // Build row-name -> index lookup for qr_mat
    std::unordered_map<std::string, size_t> rn_map;
    for (size_t i = 0; i < qr_mat.row_names.size(); ++i)
        rn_map[qr_mat.row_names[i]] = i;

    size_t n = rn.size();
    Eigen::Index q = qr_mat.data.cols();

    // Extract subset
    Mat qr_sub(static_cast<Eigen::Index>(n), q);
    for (size_t i = 0; i < n; ++i) {
        auto it = rn_map.find(rn[i]);
        if (it == rn_map.end())
            throw std::runtime_error(
                "build_cor_mat: row name '" + rn[i] + "' not found in qr_mat.");
        qr_sub.row(static_cast<Eigen::Index>(i)) =
            qr_mat.data.row(static_cast<Eigen::Index>(it->second));
    }

    // Compute Pearson correlation matrix on columns
    // Center each column
    Vec col_mean = qr_sub.colwise().mean();
    Mat centered = qr_sub.rowwise() - col_mean.transpose();

    // Compute column norms
    Vec col_norms = centered.colwise().norm();
    // Avoid division by zero
    for (Eigen::Index j = 0; j < q; ++j) {
        if (col_norms(j) < 1e-15) col_norms(j) = 1e-15;
    }

    // Correlation = (centered^T * centered) / (n-1), then normalize
    Mat cor_mat(q, q);
    for (Eigen::Index a = 0; a < q; ++a) {
        cor_mat(a, a) = 1.0;
        for (Eigen::Index b = a + 1; b < q; ++b) {
            double r = centered.col(a).dot(centered.col(b)) /
                       (col_norms(a) * col_norms(b));
            if (std::isnan(r)) r = 0.0;
            double abs_r = std::abs(r);
            cor_mat(a, b) = abs_r;
            cor_mat(b, a) = abs_r;
        }
    }

    return cor_mat;
}

// ---------------------------------------------------------------------------
// gsnca_stat_qr_detail
//
// Compute the GSNCA test statistic D and retain the two correlation matrices.
// ---------------------------------------------------------------------------
GSNCAStatDetail gsnca_stat_qr_detail(const DataFrame& objt1,
                                     const DataFrame& objt2,
                                     size_t gene_start,
                                     bool na_counts,
                                     const std::string& cor_method,
                                     const std::string& subj_var,
                                     const QRMatrix& qr_mat)
{
    GSNCAStatDetail out;

    // Build correlation matrices for both conditions
    out.cor_mat1 = build_cor_mat(objt1, qr_mat, gene_start, cor_method, subj_var);
    out.cor_mat2 = build_cor_mat(objt2, qr_mat, gene_start, cor_method, subj_var);

    Eigen::Index q1 = out.cor_mat1.rows();
    Eigen::Index q2 = out.cor_mat2.rows();

    // Validation checks
    if (q1 != q2)
        throw std::runtime_error(
            "The correlation matrices do not have the same dimensions.");

    // Symmetry check (should be guaranteed by construction)
    double asym1 = (out.cor_mat1 - out.cor_mat1.transpose()).norm();
    double asym2 = (out.cor_mat2 - out.cor_mat2.transpose()).norm();
    if (asym1 > 1e-10 || asym2 > 1e-10)
        throw std::runtime_error("The matrices should be symmetric.");

    // Eigen decomposition
    // SelfAdjointEigenSolver returns eigenvalues in ascending order,
    // so the leading (largest) eigenvector is in the last column.
    Eigen::SelfAdjointEigenSolver<Mat> e1(out.cor_mat1);
    Eigen::SelfAdjointEigenSolver<Mat> e2(out.cor_mat2);

    // Extract leading eigenvector (last column = largest eigenvalue)
    Mat evecs1 = e1.eigenvectors();
    Mat evecs2 = e2.eigenvectors();

    // Take absolute values of leading eigenvector components
    Vec p1(q1), p2(q2);
    for (Eigen::Index i = 0; i < q1; ++i)
        p1(i) = std::abs(evecs1(i, q1 - 1));
    for (Eigen::Index i = 0; i < q2; ++i)
        p2(i) = std::abs(evecs2(i, q2 - 1));

    // Compute norms: R's norm(matrix(p)) with default type="O" is the
    // column-wise maximum absolute sum.  For a single-column matrix this
    // equals sum(abs(p)), i.e. the L1 norm of the vector.
    double norm_p1 = p1.sum();  // p1 is already abs, so sum = L1 norm
    double norm_p2 = p2.sum();

    // D = sum(|p1_i * ||p1|| - p2_i * ||p2||  |)
    out.D = 0.0;
    for (Eigen::Index i = 0; i < q1; ++i) {
        out.D += std::abs(p1(i) * norm_p1 - p2(i) * norm_p2);
    }

    return out;
}

// ---------------------------------------------------------------------------
// gsnca_stat_qr
//
// Compute the GSNCA test statistic D.
// ---------------------------------------------------------------------------
double gsnca_stat_qr(const DataFrame& objt1,
                     const DataFrame& objt2,
                     size_t gene_start,
                     bool na_counts,
                     const std::string& cor_method,
                     const std::string& subj_var,
                     const QRMatrix& qr_mat)
{
    return gsnca_stat_qr_detail(
        objt1, objt2, gene_start, na_counts, cor_method, subj_var, qr_mat).D;
}

} // namespace gsca
