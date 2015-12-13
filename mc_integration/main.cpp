#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include <cstddef>
#include <cmath>
#include <tuple>
#include <vector>
#include <fstream>

bool is_power_of_two(unsigned x) { return ((x != 0) && ((x & (~x + 1)) == x)); }

using result = std::tuple<unsigned, double, double>;

std::vector<result> monte_carlo_sum(unsigned n) {
    static std::mt19937 gen((std::random_device())());
    static std::normal_distribution<> norm_dist(0.0, 1.0);

    std::vector<result> ret;

    std::size_t sum = 0;
    for (std::size_t N = 0; N != n; ++N) {
        if (norm_dist(gen) > 5.0) {
            ++sum;
        }

        if (is_power_of_two(N) && (N != 1)) {
            auto mean = static_cast<double>(sum) / N;

            // X^2 = X, since x in {0, 1}
            // Var(X) = E(X^2) - E(X)^2 = E(X) - E(X)^2 = E(X) (1 - E(X))
            auto err = std::sqrt(mean * (1.0 - mean) / (N - 1));

            ret.emplace_back(N, mean, err);
        }
    }
    return ret;
}

std::vector<result> importance_sampling(unsigned n) {
    static std::mt19937 gen((std::random_device())());

    // Sample from Exp(1) truncated at x = 5 using inverse transform
    // cdf: G(x) = 1 - exp(-x)  ->  - log(1 - U[1 - e^-5, 1])
    static std::uniform_real_distribution<> trunc_unif(1.0 - std::exp(-5.0),
                                                       1.0);
    // N(0,1) pdf:
    auto f = [](double x) {
        return std::exp(-x * x / 2.0) / std::sqrt(2.0 * M_PI);
    };

    // Exp(1) pdf:
    auto g = [](double x) { return std::exp(-x + 5.0); };
    // + 5.0 due to normalization

    std::vector<result> ret;

    // for mean and variance
    double sum = 0.0;
    double sq_sum = 0.0;

    for (unsigned N = 0; N != n; ++N) {
        auto x = -std::log(1 - trunc_unif(gen));
        auto f_x = f(x);
        auto g_x = g(x);

        sum += f_x / g_x;
        sq_sum += f_x * f_x / (g_x * g_x);

        if (is_power_of_two(N) && (N != 1)) {
            auto mean = sum / N;
            auto err = std::sqrt((sq_sum / N - mean * mean) / (N - 1));

            ret.emplace_back(N, mean, err);
        }
    }
    return ret;
}

std::vector<result> result_mean(unsigned n, unsigned N,
                                std::vector<result> (*f)(unsigned)) {
    std::vector<std::vector<result>> data;
    for (std::size_t i = 0; i != N; ++i) {
        data.push_back(f(n));
    }

    std::vector<result> ret;

    // Iterates through rows
    for (std::size_t i = 0; i != data[0].size(); ++i) {
        double sum = 0.0;
        double sum_sq = 0.0;

        // Iterates through tables
        for (std::size_t j = 0; j != N; ++j) {
            auto x = std::get<1>(data[j][i]);
            sum += x;
            sum_sq += x * x;
        }
        auto mean = sum / N;
        auto mean_sq = sum_sq / N;

        // auto err = std::sqrt((mean_sq - mean*mean) / (N - 1));
        // ret.emplace_back(std::get<0>(data[0][i]), mean, err);
        auto err = std::sqrt((mean_sq - mean * mean) * N / (N - 1));
        ret.emplace_back(std::get<0>(data[0][i]), std::get<1>(data[0][i]), err);
    }
    return ret;
}

void write(std::ostream &os, std::vector<result> data) {
    unsigned N;
    double mean, err;

    auto it = data.cbegin(), end = data.cend();
    while (os && (it != end)) {
        std::tie(N, mean, err) = *it++;
        os << N << "\t" << mean << "\t" << err << "\n";
    }
}

int main() {
    // exact value
    std::cout << 0.5 * std::erfc(5.0 / std::sqrt(2.0)) << std::endl;

    // max. number of samples
    const unsigned n = 100000000;
    // number of independent runs
    const unsigned N = 10;

    std::ofstream os("mc_sum_running.tsv");
    auto data = monte_carlo_sum(n);
    write(os, data);
    os.close();

    os.open("mc_sum_indepedent_runs.tsv");
    data = result_mean(n, N, monte_carlo_sum);
    write(os, data);
    os.close();

    os.open("importance_sampling_running.tsv");
    data = importance_sampling(n);
    write(os, data);
    os.close();

    os.open("importance_sampling_independent_runs.tsv");
    data = result_mean(n, N, importance_sampling);
    write(os, data);
    os.close();

    return 0;
}
