// ============================================================================
// gsnca_sd.cpp
//
// Check per-gene standard deviation across groups.
// Genes whose SD is below min_sd in either group are flagged for removal.
//
// Port of: gsnca_sd.R
// ============================================================================
#include "MetaGSCA.h"
#include <cmath>
#include <numeric>

namespace gsca {

// ---------------------------------------------------------------------------
// Helper: compute column-wise SD for numeric columns starting at gene_start.
// Returns a vector of SDs, one per gene column.
// ---------------------------------------------------------------------------
static std::vector<double> col_sds(const DataFrame& df, size_t gene_start) {
    size_t n = df.nrow();
    size_t ncols = df.ncol();
    std::vector<double> sds;
    sds.reserve(ncols - gene_start);

    for (size_t c = gene_start; c < ncols; ++c) {
        const auto& v = std::get<NumCol>(df.col_by_idx(c));
        // Compute mean
        double sum = 0.0;
        size_t count = 0;
        for (double x : v) {
            if (!std::isnan(x)) { sum += x; ++count; }
        }
        if (count < 2) { sds.push_back(0.0); continue; }
        double mean = sum / static_cast<double>(count);

        // Compute variance (sample variance, N-1)
        double ss = 0.0;
        for (double x : v) {
            if (!std::isnan(x)) {
                double d = x - mean;
                ss += d * d;
            }
        }
        sds.push_back(std::sqrt(ss / static_cast<double>(count - 1)));
    }
    return sds;
}

// ---------------------------------------------------------------------------
// gsnca_sd
//
// Return names of genes that have SD < min_sd in either group.
// Throws if fewer than 2 genes survive.
// ---------------------------------------------------------------------------
std::vector<std::string> gsnca_sd(const DataFrame& objt1,
                                  const DataFrame& objt2,
                                  const std::vector<std::string>& genes,
                                  size_t gene_start,
                                  double min_sd)
{
    auto sd1 = col_sds(objt1, gene_start);
    auto sd2 = col_sds(objt2, gene_start);

    if (sd1.size() != genes.size() || sd2.size() != genes.size())
        throw std::runtime_error(
            "gsnca_sd: number of SDs does not match number of genes.");

    std::vector<std::string> removed;
    size_t kept = 0;
    for (size_t i = 0; i < genes.size(); ++i) {
        if (sd1[i] < min_sd || sd2[i] < min_sd) {
            removed.push_back(genes[i]);
        } else {
            ++kept;
        }
    }

    if (kept < 2) {
        throw std::runtime_error(
            "There must be at least 2 genes with standard deviation of "
            "gene expression data bigger than " + std::to_string(min_sd) +
            " in group 1 and group 2");
    }

    return removed;
}

} // namespace gsca
