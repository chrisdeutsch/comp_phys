#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include <cstddef>
#include <cmath>

double mc_sum(std::size_t n) {
    static std::mt19937 gen((std::random_device())());
    std::normal_distribution<> norm_dist(0, 1.0);

    std::size_t accept = 0;
    for (std::size_t i = 0; i != n; ++i) {
        if (norm_dist(gen) > 5.0) {
            ++accept;
        }
    }
    return static_cast<double>(accept) / n;
}

double instrumental_density(std::size_t n) {
    // N(0,1) pdf:
    auto f = [](double x) {
        return std::exp(-x * x / 2.0) / std::sqrt(2.0 * M_PI);
    };
    // Exp(1) pdf:
    auto g = [](double x) { return std::exp(-x); };

    static std::mt19937 gen((std::random_device())());

    // Sample from Exp(1) truncated at x = 5 using inverse transform
    // cdf: G(x) = 1 - exp(-lambda * x)
    std::uniform_real_distribution<> trunc_unif(1.0 - std::exp(-5.0), 1.0);
    std::uniform_real_distribution<> unif;

    // Get best envelope by setting f(5) == M g(5)
    const double M = f(5) / g(5);
    // Probability that x > 5 for Exp(1)
    const double area = M * std::exp(-5);

    std::size_t accept = 0;
    for (std::size_t i = 0; i != n; ++i) {
        // Inverse transform to generate Exp(1) truncated at x = 5
        auto x = -std::log(1 - trunc_unif(gen));
        auto y = unif(gen) * M * g(x);
        if (y < f(x)) {
            ++accept;
        }
    }

    return accept * area / n;
}

int main() {
    // exact value
    std::cout << 0.5 * std::erfc(5.0 / std::sqrt(2.0)) << std::endl;

    std::cout << instrumental_density(100000) << std::endl;

    //
    std::cout << mc_sum(10000000) << std::endl;
    return 0;
}
