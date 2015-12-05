#ifndef N_SPHERE
#define N_SPHERE

#include <cstddef>
#include <random>
#include <valarray>
#include <cmath>
#include <iostream>

// Random number generation on a sphere using Marsaglia algorithm

using n_vec = std::valarray<double>;

// Random vector on shell of unit-sphere
n_vec random_n_sphere_shell(std::size_t n) {
    static std::mt19937 gen((std::random_device())());
    static std::normal_distribution<> norm_dist;

    n_vec ret(n);
    for (auto &coord : ret) {
        coord = norm_dist(gen);
    }

    auto radius = std::sqrt(std::pow(ret, 2.0).sum());

    return ret / radius;
}

// Random vector in unit-sphere
n_vec random_n_sphere(std::size_t n) {
    static std::mt19937 gen((std::random_device())());
    static std::uniform_real_distribution<> r_dist;

    return std::pow(r_dist(gen), 1.0 / n) * random_n_sphere_shell(n);
}

// Print vector coordinates comma separated
std::ostream &print_vec(const n_vec &v, std::ostream &os = std::cout) {
    auto it = std::begin(v), end = std::end(v);
    while (os && it != end) {
        os << *it;
        if (++it != end) {
            os << ", ";
        }
    }
    os << "\n";
    return os;
}

#endif // N_SPHERE
