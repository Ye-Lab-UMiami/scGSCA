#pragma once
#include <cmath>
#include <stdexcept>
#include <vector>

namespace gsca {
inline double ct_right_tail_pvalue(double observed, const std::vector<double>& perm) {
    if (!std::isfinite(observed) || perm.empty())
        throw std::invalid_argument("CT P value requires a finite observed Delta and nonempty permutations.");
    size_t count = 0;
    for (double d : perm) {
        if (!std::isfinite(d))
            throw std::invalid_argument("Nonfinite CT permutation Delta.");
        if (d >= observed) ++count;
    }
    return static_cast<double>(count + 1) / static_cast<double>(perm.size() + 1);
}
}
