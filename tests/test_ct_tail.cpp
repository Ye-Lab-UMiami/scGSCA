#include "ct_tail.h"
#include <iostream>
#include <limits>

int main() {
    const std::vector<double> perm{-2, -1, 0, 1, 2};
    auto check = [](double actual, double expected) {
        if (std::abs(actual - expected) > 1e-12) throw std::runtime_error("P value mismatch");
    };
    auto rejects = [](auto f) {
        try { f(); } catch (const std::invalid_argument&) { return; }
        throw std::runtime_error("Invalid input accepted");
    };
    check(gsca::ct_right_tail_pvalue(-1, perm), 5.0/6);
    check(gsca::ct_right_tail_pvalue(1, perm), 3.0/6);
    check(gsca::ct_right_tail_pvalue(0, perm), 4.0/6);
    check(gsca::ct_right_tail_pvalue(3, perm), 1.0/6);
    check(gsca::ct_right_tail_pvalue(-3, perm), 1.0);
    check(gsca::ct_right_tail_pvalue(0, {0,0}), 1.0);
    rejects([&] { gsca::ct_right_tail_pvalue(0, {}); });
    rejects([&] { gsca::ct_right_tail_pvalue(0, {std::numeric_limits<double>::quiet_NaN()}); });
    rejects([&] { gsca::ct_right_tail_pvalue(std::numeric_limits<double>::quiet_NaN(), perm); });
    std::cout << "Fixed-right-tail tests passed (both signs, ties, zero, extremes, validation).\n";
}
