#ifndef AUTOCORRELATION_H
#define AUTOCORRELATION_H

#include <cstddef>
#include <numeric>
#include <vector>
#include <random>

template <typename InputIterator, typename OutputIterator>
void autocorrelation_function(InputIterator begin, InputIterator end,
                              OutputIterator dest, std::size_t t_end) {
    const std::size_t sz = end - begin;
    if (t_end >= sz) {
        t_end = sz - 1;
    }

    // mean of data
    double mean = std::accumulate(begin, end, 0.0) / sz;

    // variance of data
    auto add_squared = [](double x, double y) { return x + y * y; };
    double variance =
        std::accumulate(begin, end, 0.0, add_squared) / sz - mean * mean;

    // calculate autocorrelation for lags up to (excl.) t_end
    for (std::size_t t = 0; t < t_end; ++t) {
        double sum = 0.0;
        auto it = begin, lag_it = begin + t;
        while (lag_it != end) {
            sum += (*it++ - mean) * (*lag_it++ - mean);
        }
        sum /= sz - t;
        *dest++ = sum / variance;
    }
}

std::vector<double> time_series(unsigned n, double alpha, double x0) {
    static std::mt19937 gen((std::random_device())());
    static std::uniform_real_distribution<> dist;

    std::vector<double> series;
    series.reserve(n + 1);
    series.push_back(x0);

    for (unsigned i = 0; i != n; ++i) {
        x0 *= alpha;
        x0 += (1.0 - alpha) * dist(gen);
        series.push_back(x0);
    }
    return series;
}

#endif // AUTOCORRELATION_H
