// ============================================================================
// rng.cpp — Random number generation utilities.
// ============================================================================
#include "rng.h"

namespace gsca {

std::vector<size_t> RNG::sample_indices(size_t n, size_t k) {
    std::vector<size_t> pool(n);
    std::iota(pool.begin(), pool.end(), 0);
    for (size_t i = 0; i < k; ++i) {
        std::uniform_int_distribution<size_t> dist(i, n - 1);
        std::swap(pool[i], pool[dist(gen_)]);
    }
    return std::vector<size_t>(pool.begin(), pool.begin() + k);
}

} // namespace gsca
