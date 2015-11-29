#include <iostream>
#include <iterator>
#include "binning_method.h"
#include "../autocorrelation/autocorrelation.h"

int main() {
    auto series = time_series(100000, 0.7, 0.1);
    std::vector<double> binned;
    block_mean(series.cbegin(), series.cend(), std::back_inserter(binned), 10);
    std::cout << "mean: " << mean(binned.cbegin(), binned.cend()) << std::endl;
    std::cout << "standard deviation: "
              << std_dev(binned.cbegin(), binned.cend()) << std::endl;

    return 0;
}
