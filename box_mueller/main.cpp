#define _USE_MATH_DEFINES
#include <cstddef>
#include <random>
#include <cmath>

#include <iostream>
#include <chrono>

double box_mueller() {
    static std::mt19937 gen((std::random_device())());
    static std::uniform_real_distribution<> dist;

    static bool available = false;
    static double rnd_num;

    // if a random number is available return it
    if (available) {
        available = false;
        return rnd_num;
    }

    // else generate two new ones
    auto u = dist(gen), v = dist(gen);

    static const double two_pi = 2.0 * M_PI;
    auto sqrt_log = std::sqrt(-2.0 * std::log(u));

    // save the second random number for next call and return the first
    rnd_num = sqrt_log * std::sin(two_pi * v);
    available = true;

    return sqrt_log * std::cos(two_pi * v);
}

double box_mueller_polar() {
    static std::mt19937 gen((std::random_device())());
    static std::uniform_real_distribution<> dist(-1.0, 1.0);

    static bool available = false;
    static double rnd_num;

    // if a random number is available return it
    if (available) {
        available = false;
        return rnd_num;
    }

    auto u1 = dist(gen), u2 = dist(gen);
    double s;
    while ((s = u1 * u1 + u2 * u2) > 1.0) {
        u1 = dist(gen);
        u2 = dist(gen);
    }

    auto z = std::sqrt(-2.0 * std::log(s) / s);

    // save the second random number for next call and return the first
    rnd_num = z * u2;
    available = true;

    return z * u1;
}

int main() {
    const std::size_t count = 10000000;

    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i != count; ++i) {
        box_mueller();
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "box_mueller " << count << " generations: "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
              << " us" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i != count; ++i) {
        box_mueller();
    }
    end = std::chrono::high_resolution_clock::now();

    std::cout << "box_mueller_polar " << count << " generations: "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
              << " us" << std::endl;

    return 0;
}
