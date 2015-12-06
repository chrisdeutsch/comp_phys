#define _USE_MATH_DEFINES
#include <iostream>
#include <cstddef>
#include <random>
#include <vector>
#include <numeric>
#include <cmath>
#include "../rng_sphere/random_n_sphere.h"

double sinc(const double x) {
    if (x == 0.0) {
        return 1.0;
    }
    return std::sin(x) / x;
}

// double slit intensity (b: slit width, a: slit distance)
// I(alpha) = I0 * (sinc(k/2 b sin(alpha))^2 * cos(k/2 a * sin(alpha))^2
// normalized to I0 = 1
double double_slit_intensity(double alpha, double k, double slit_width,
                             double slit_distance) {
    const double arg = k * std::sin(alpha) / 2;
    return std::pow(sinc(arg * slit_width) * std::cos(arg * slit_distance),
                    2.0);
}

// Voltages in volts / volts^2
std::vector<double> simulate_double_slit(double mean_voltage,
                                         double stddev_voltage,
                                         std::size_t size) {
    // phys. constants
    constexpr double hbar = 1.0545718E-34;      // J s
    constexpr double e_charge = 1.60217662e-19; // C
    constexpr double e_mass = 9.109e-31;        // kg

    // wavenumber of electron accelerated by "voltage"
    auto k = [=](double voltage) {
        return std::sqrt(2.0 * e_mass * e_charge * voltage) / hbar;
    };

    // double slit parameters
    constexpr double slit_width = 62e-9;     // m
    constexpr double slit_distance = 272e-9; // m

    // position of the 2nd main minimum
    // sin(alpha) = 2pi / (k a) (1/2 + n)  -  n = 1 for 2nd minimum
    double alpha_min = std::asin(3 * M_PI / (k(mean_voltage) * slit_distance));

    // rng
    static std::mt19937 gen((std::random_device())());
    // distributions
    std::uniform_real_distribution<> alpha_dist(-alpha_min, alpha_min);
    std::uniform_real_distribution<> intensity_dist(0, 1.0);
    std::normal_distribution<> voltage_dist(mean_voltage, stddev_voltage);

    std::vector<double> angles;
    std::size_t accept = 0;
    angles.reserve(size);

    while (accept != size) {
        auto alpha = alpha_dist(gen);
        auto intensity = intensity_dist(gen);
        auto voltage = voltage_dist(gen);
        auto theo_intensity =
            double_slit_intensity(alpha, k(voltage), slit_width, slit_distance);

        if (intensity < theo_intensity) {
            angles.push_back(alpha);
            ++accept;
        }
    }
    return angles;
}

double viviani_body_vol(std::size_t n) {
    std::mt19937 gen((std::random_device())());

    // Random number generators for box around cylinder
    std::uniform_real_distribution<> x_dist(0.0, 2.0);
    std::uniform_real_distribution<> y_dist(-1.0, 1.0);
    std::uniform_real_distribution<> z_dist(-2.0, 2.0);

    auto in_sphere = [](double x, double y, double z) {
        return x * x + y * y + z * z < 4.0;
    };

    auto in_cylinder = [](double x, double y, double z) {
        const double dx = x - 1.0;
        return dx * dx + y * y < 1.0;
    };

    double x, y, z;
    std::size_t accept = 0;
    for (std::size_t i = 0; i < n; ++i) {
        x = x_dist(gen);
        y = y_dist(gen);
        z = z_dist(gen);

        if (in_sphere(x, y, z) && in_cylinder(x, y, z)) {
            ++accept;
        }
    }

    constexpr double vol = 16.0;
    return accept * vol / n;
}

double viviani_body_area(std::size_t n) {
    std::size_t accept = 0;

    auto in_cylinder = [](double x, double y, double z) {
        const double dx = x - 1.0;
        return dx * dx + y * y < 1.0;
    };

    for (std::size_t i = 0; i != n; ++i) {
        // Random vector on unit-sphere
        auto random_vector = random_n_sphere_shell(3);
        // Rescale to radius 2 and take modulus of x-coordinate (random vectors
        // with negative x-coordinate will not be accepted, since vectors in the
        // cylinder obey x >= 0)
        random_vector *= 2.0;
        random_vector[0] = std::abs(random_vector[0]);

        if (in_cylinder(random_vector[0], random_vector[1], random_vector[2])) {
            ++accept;
        }
    }

    // 1/2 * 4 pi r^2 (with r = 2)
    const double sphere_half_area = 8.0 * M_PI;
    return accept * sphere_half_area / n;
}

double n_simplex_vol(std::size_t n, std::size_t samp) {
    std::mt19937 gen((std::random_device())());
    std::uniform_real_distribution<> dist(0.0, 1.0);

    std::size_t accept = 0;
    for (std::size_t i = 0; i != samp; ++i) {
        std::vector<double> pts;
        for (std::size_t i = 0; i != n; ++i) {
            pts.push_back(dist(gen));
        }

        if (std::accumulate(pts.cbegin(), pts.cend(), 0.0) < 1.0) {
            ++accept;
        }
        pts.clear();
    }
    return 1.0 * accept / samp;
}

double n_simplex_vol_sphere(std::size_t n, std::size_t samp) {
    std::size_t accept = 0;
    for (std::size_t i = 0; i != samp; ++i) {
        auto x = std::abs(random_n_sphere(n));
        if (x.sum() < 1.0) {
            ++accept;
        }
    }
    // volume of unit n-sphere
    const double vol = std::pow(2.0, -1.0 * n) * std::pow(M_PI, n / 2.0) /
                       std::tgamma(n / 2.0 + 1.0);
    return accept * vol / samp;
}

#include <fstream>
#include <iterator>
#include <algorithm>
int main() {
    // Should be 8 (pi - 2) = 9.13274122871835
    // std::cout << viviani_body_area(100000000) << std::endl;

    // double slit simulation
    auto slit_stable_voltage = simulate_double_slit(600.0, 0.0, 10000000);
    auto slit_variable_voltage = simulate_double_slit(600.0, 10.0, 10000000);

    std::ofstream os("slit_stable.tsv");
    std::ostream_iterator<double> os_it(os, "\n");
    std::copy(slit_stable_voltage.cbegin(), slit_stable_voltage.cend(), os_it);
    os.close();

    os.open("slit_variable.tsv");
    std::copy(slit_variable_voltage.cbegin(), slit_variable_voltage.cend(),
              os_it);
    os.close();

    return 0;
}
