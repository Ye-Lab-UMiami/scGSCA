// ============================================================================
// gsnca_filt.cpp
//
// Filter out genes whose column sum is <= filter_count in either group.
//
// Port of: gsnca_filt.R
// ============================================================================
#include "MetaGSCA.h"
#include <cmath>

namespace gsca {

// ---------------------------------------------------------------------------
// Helper: compute column sums for gene columns.
// ---------------------------------------------------------------------------
static std::vector<double> col_sums(const DataFrame& df, size_t gene_start) {
    size_t ncols = df.ncol();
    std::vector<double> sums;
    sums.reserve(ncols - gene_start);

    for (size_t c = gene_start; c < ncols; ++c) {
        const auto& v = std::get<NumCol>(df.col_by_idx(c));
        double s = 0.0;
        for (double x : v) {
            if (!std::isnan(x)) s += x;
        }
        sums.push_back(s);
    }
    return sums;
}

// ---------------------------------------------------------------------------
// gsnca_filt
//
// Return names of genes with column sum <= filter_count in either group.
// Throws if fewer than 2 genes survive.
// ---------------------------------------------------------------------------
std::vector<std::string> gsnca_filt(const DataFrame& objt1,
                                    const DataFrame& objt2,
                                    const std::vector<std::string>& genes,
                                    size_t gene_start,
                                    double filter_count)
{
    auto cs1 = col_sums(objt1, gene_start);
    auto cs2 = col_sums(objt2, gene_start);

    if (cs1.size() != genes.size() || cs2.size() != genes.size())
        throw std::runtime_error(
            "gsnca_filt: number of sums does not match number of genes.");

    std::vector<std::string> removed;
    size_t kept = 0;
    for (size_t i = 0; i < genes.size(); ++i) {
        if (cs1[i] <= filter_count || cs2[i] <= filter_count) {
            removed.push_back(genes[i]);
        } else {
            ++kept;
        }
    }

    if (kept < 2) {
        throw std::runtime_error(
            "There must be at least 2 genes with sum of counts of gene "
            "expression data bigger than " + std::to_string(filter_count) +
            " in group 1 and group 2");
    }

    return removed;
}

} // namespace gsca
