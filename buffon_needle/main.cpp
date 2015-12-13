#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>
#include <algorithm>

using needle = std::tuple<double, double>;

std::vector<needle> scatter_needles(double height, std::size_t n) {
    static std::mt19937 gen((std::random_device())());
    std::uniform_real_distribution<> x_dist(0.0, height);
    std::uniform_real_distribution<> phi_dist(0.0, 2 * M_PI);

    std::vector<needle> ret;
    ret.reserve(n);

    for (std::size_t i = 0; i != n; ++i) {
        ret.emplace_back(x_dist(gen), phi_dist(gen));
    }

    return ret;
}

int main() {
    const double height = 10.0;
    const double stripe_spacing = 1.0;
    const double needle_length = 0.75;
    const std::size_t needle_num = 1000000;

    auto needles = scatter_needles(height, needle_num);

    auto is_intersecting = [=](const needle &n) {
        double x, phi;
        std::tie(x, phi) = n;

        // project into [0.0, d];
        x -= floor(x / stripe_spacing) * stripe_spacing;
        return std::min(x, stripe_spacing - x) <
               0.5 * needle_length * std::abs(std::cos(phi));
    };

    auto intersections =
        std::count_if(needles.cbegin(), needles.cend(), is_intersecting);

    double pi = 2.0 * needle_num * needle_length / (intersections * stripe_spacing);
    std::cout << pi << std::endl;

    return 0;
}
