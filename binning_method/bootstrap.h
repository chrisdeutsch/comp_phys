#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

#include <cstddef>
#include <random>
#include <vector>

// for std_dev
#include "binning_method.h"

// bootstrap error of estimator with "b" bootstrap samples
template <typename InputIt>
double bootstrap_error(InputIt beg, InputIt end,
                       std::function<double(const std::vector<double> &)> estimator,
                       std::size_t b) {
    // size of sample to bootstrap
    std::size_t sample_sz = end - beg;

    // RNG to sample for bootstrap
    static std::mt19937 gen((std::random_device())());
    std::uniform_int_distribution<std::size_t> dist(0, sample_sz - 1);

    // vector to hold the bootstrap sample and a list of observed values
    std::vector<double> bootstrap_sample(sample_sz);
    std::vector<double> estimated_vals;

    for (std::size_t sample_num = 0; sample_num != b; ++sample_num) {
        // create bootstrap sample
        for (std::size_t i = 0; i != sample_sz; ++i) {
            bootstrap_sample[i] = beg[dist(gen)];
        }
        // add estimated value for this sample to list
        estimated_vals.push_back(estimator(bootstrap_sample));
    }
    // estimator error
    return std_dev(estimated_vals.cbegin(), estimated_vals.cend());
}

#endif // BOOTSTRAP_H
