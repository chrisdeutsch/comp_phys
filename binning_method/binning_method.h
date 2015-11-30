#ifndef BINNING_METHOD_H
#define BINNING_METHOD_H

#include <cstddef>
#include <numeric>
#include <functional>
#include <cmath>

template <typename Iterator>
double mean(Iterator beg, Iterator end) {
    return std::accumulate(beg, end, 0.0) / (end - beg);
}

template <typename Iterator>
double mean(Iterator beg, Iterator end, std::function<double(double)> f) {
    auto add_f = [&f](double x, double y) { return x + f(y); };
    return std::accumulate(beg, end, 0.0, add_f) / (end - beg);
}

// sample standard deviation
template <typename Iterator>
double std_dev(Iterator beg, Iterator end) {
    const double mu = mean(beg, end);
    auto add_sq_dev = [&mu](double x, double y) {
        double dev = (y - mu);
        return x + dev * dev;
    };
    return std::sqrt(std::accumulate(beg, end, 0.0, add_sq_dev) /
                     (end - beg - 1));
}

// standard error of sample mean
template <typename Iterator>
double std_err(Iterator beg, Iterator end) {
    return std_dev(beg, end) / std::sqrt(end - beg);
}

template <typename InputIterator, typename OutputIterator>
void block_mean(InputIterator beg, InputIterator end, OutputIterator dest,
                std::size_t b) {
    std::size_t no_blocks = (end - beg) / b;
    auto block_beg = beg, block_end = beg + b;

    for (std::size_t i = 0; i != no_blocks; ++i) {
        ++*dest = mean(block_beg, block_end);
        block_beg += b;
        block_end += b;
    }
}

#endif // BINNING_METHOD_H
