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

template <typename Iterator>
double std_dev(Iterator beg, Iterator end) {
    auto mu = mean(beg, end);
    auto square = [](double x) { return x * x; };
    return std::sqrt(mean(beg, end, square) - mu * mu);
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
