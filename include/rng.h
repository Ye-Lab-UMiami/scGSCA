// ============================================================================
// rng.h
//
// Reproducible random number generation utilities.
// Uses the Mersenne Twister (mt19937_64) seeded explicitly, mirroring R's
// set.seed() + L'Ecuyer-CMRG for reproducibility within a single stream.
// ============================================================================
#ifndef METAGSCA_RNG_H
#define METAGSCA_RNG_H

#include <random>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <numeric>

namespace gsca {

// -------------------------------------------------------------------------
// Thread-local RNG wrapper
// -------------------------------------------------------------------------
class RNG {
public:
    explicit RNG(uint64_t seed = 42) : gen_(seed) {}

    // Shuffle a vector of T in place (Fisher-Yates).
    template <typename T>
    void shuffle(std::vector<T>& v) {
        for (size_t i = v.size() - 1; i > 0; --i) {
            std::uniform_int_distribution<size_t> dist(0, i);
            std::swap(v[i], v[dist(gen_)]);
        }
    }

    // Sample k elements from [0, n) without replacement.
    std::vector<size_t> sample_indices(size_t n, size_t k);

    // Random boolean with probability 0.5.
    bool coin_flip() {
        return dist01_(gen_) < 0.5;
    }

    // Uniform [0,1).
    double uniform() { return dist01_(gen_); }

    // Normal(0,1).
    double normal() { return norm_(gen_); }

    // Access underlying generator.
    std::mt19937_64& engine() { return gen_; }

private:
    std::mt19937_64 gen_;
    std::uniform_real_distribution<double> dist01_{0.0, 1.0};
    std::normal_distribution<double> norm_{0.0, 1.0};
};

} // namespace gsca

#endif // METAGSCA_RNG_H
