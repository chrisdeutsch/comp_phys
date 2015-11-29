#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

#include <cstddef>
#include <random>
#include <functional>
#include <vector>

// for std_dev
#include "binning_method.h"

template <typename InputIt>
double bootstrap_error(InputIt beg, InputIt end,
                       std::function<double(const std::vector<double> &)> observable,
                       std::size_t samples) {
    // size of sample to bootstrap
    std::size_t sample_sz = end - beg;

    // RNG to sample for bootstrap
    std::mt19937 gen((std::random_device())());
    std::uniform_int_distribution<std::size_t> dist(0, sample_sz - 1);

    // vector to hold the bootstrap sample and a list of observed values
    std::vector<double> b_sample(sample_sz);
    std::vector<double> observ_vals;

    for (std::size_t sample_no = 0; sample_no != samples; ++sample_no) {
        // create bootstrap sample
        for (std::size_t i = 0; i != sample_sz; ++i) {
            b_sample[i] = beg[dist(gen)];
        }
        // add observed value for this sample to list of observed values
        observ_vals.push_back(observable(b_sample));
    }
    // estimated observable error
    return std_dev(observ_vals.cbegin(), observ_vals.cend());
}

#endif // BOOTSTRAP_H
