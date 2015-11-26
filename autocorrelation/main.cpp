#include <cstddef>
#include <numeric>

#include <vector>
#include <iostream>
#include <iterator>

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

template <typename Container>
void print_container(const Container &c) {
    for (const auto &elem : c) {
        std::cout << elem << std::endl;
    }
}

int main() {
    std::vector<double> vd = {
        0.727504, 0.718915, 0.449088,  0.797964, 0.0968756,  0.8185,   0.882323,
        0.164807, 0.607614, 0.455075,  0.416633, 0.209309,   0.595713, 0.459449,
        0.76053,  0.218027, 0.12488,   0.355834, 0.306831,   0.510878, 0.890206,
        0.137602, 0.305521, 0.0691475, 0.750782, 0.515048,   0.55674,  0.602874,
        0.630518, 0.617191, 0.884834,  0.216637, 0.398987,   0.075161, 0.223157,
        0.548678, 0.919536, 0.275558,  0.595751, 0.00478411, 0.181478, 0.214613,
        0.175442, 0.61162,  0.274793,  0.409027, 0.45021,    0.312003, 0.177765,
        0.375598};

    std::vector<double> corr;
    autocorrelation_function(vd.cbegin(), vd.cend(), std::back_inserter(corr), 10);
    print_container(corr);

    return 0;
}
