#include <fstream>
#include <iterator>
#include <algorithm>
#include "autocorrelation.h"

int main() {
    auto series = time_series(10000, 0.9, 0.0);

    std::vector<double> autocorrelation;
    autocorrelation_function(series.cbegin(), series.cend(),
                             std::back_inserter(autocorrelation), 100);

    std::ofstream os("autocorrelation.dat");
    std::ostream_iterator<double> os_it(os, "\n");
    std::copy(autocorrelation.cbegin(), autocorrelation.cend(), os_it);

    return 0;
}
