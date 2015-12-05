#include "random_n_sphere.h"
#include <iostream>
#include <vector>
#include <fstream>


int main() {
    std::vector<n_vec> pts;

    for (std::size_t i = 0; i < 1000; ++i) {
        pts.push_back(std::abs(random_n_sphere(2)));
    }

    std::ofstream os("sphere.csv");
    for (const auto &elem : pts) {
        print_vec(elem, os);
    }

    return 0;
}
