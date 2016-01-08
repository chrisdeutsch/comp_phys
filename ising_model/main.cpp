#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>

/**
 * Spin lattice as an array of arrays of shorts
 * (-1 or 1 depending on orientation)
 */
template <std::size_t N>
using Lattice = std::array<std::array<short, N>, N>;

/**
  * Calculate energy at lattice site given by "x" and "y"
  */
template <std::size_t N>
double site_energy(const Lattice<N> &lattice, double J, std::size_t x,
                   std::size_t y) {
    assert(x < N && y < N);

    /**
     * Maps indices to fulfill the periodic boundary condition
     * e.g. maps -2 -> N - 2; -1 -> N - 1
     *           (0, ..., N-1) -> (0, ..., N-1); N -> 0 etc.
     */
    static auto pbc = [N_int = static_cast<int>(N)](int n)->std::size_t {
        assert(n >= -N_int && "n must be greater-equal than -N");
        return (N_int + n) % N_int;
    };

    return -0.5 * J * lattice[x][y] * (lattice[pbc(x - 1)][y]     // Up
                                       + lattice[pbc(x + 1)][y]   // Down
                                       + lattice[x][pbc(y - 1)]   // Left
                                       + lattice[x][pbc(y + 1)]); // Right
}

/**
 * Calculates the total energy of the spin lattice
 */
template <std::size_t N>
double energy(const Lattice<N> &lattice, double J) {
    double energy = 0.0;
    for (int x = 0; x != N; ++x) {
        for (int y = 0; y != N; ++y) {
            energy += site_energy(lattice, J, x, y);
        }
    }
    return energy;
}

/**
 * Prints the lattice using custom characters for spin-up and -down
 */
template <std::size_t N>
void print_lattice(Lattice<N> lattice, char up = 'x', char down = 'o') {
    for (const auto &row : lattice) {
        for (const auto &spin : row) {
            std::cout << (spin == 1 ? up : down);
        }
        std::cout << "\n";
    }
}

int main() {
    // Lattice size
    constexpr std::size_t lattice_size = 20;
    Lattice<lattice_size> lattice;

    // Hamiltonian
    double J = -1.0;

    // Simulated annealing settings
    const double temp_max = 1; // from linear cooling scheme
    unsigned kmax = 100000;    // maximum SA steps

    // RNG
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution bernoulli_dist;

    // Populate lattice with random spins
    for (auto &row : lattice) {
        for (auto &spin : row) {
            spin = bernoulli_dist(gen) ? 1 : -1;
        }
    }

    std::cout << "Start:" << std::endl;
    print_lattice(lattice);
    std::cout << "Energy: " << energy(lattice, J) << "\n" << std::endl;

    /**
      * Acceptance probability for a transition from energy "e0" to energy "e1"
      * (delta_e = e1 - e0) at temperature "temp"
      */
    auto acceptance_probability = [](double delta_e, double temp) {
        return std::exp(-delta_e / temp);
    };

    /**
      * Linear cooling scheme
      */
    auto annealing_schedule = [temp_max](double time_budget) {
        return temp_max * (1 - time_budget);
    };

    std::uniform_int_distribution<std::size_t> index_dist(0, lattice_size - 1);
    std::uniform_real_distribution<> unif_dist;

    for (unsigned k = 0; k != kmax; ++k) {
        // Choose random spin and flip it
        auto x_rand = index_dist(gen);
        auto y_rand = index_dist(gen);

        /**
         * Flipping the spin results in sign change of the energy at the lattice
         * site. Therefore the change in energy is -2 times the energy at the
         * lattice site.
         */
        auto delta_e = -2 * site_energy(lattice, J, x_rand, y_rand);

        auto temp = annealing_schedule(k / static_cast<double>(kmax));
        if (acceptance_probability(delta_e, temp) > unif_dist(gen)) {
            // Accept new state (flip the spin)
            lattice[x_rand][y_rand] *= -1;
        } // else: Reject new state
    }

    std::cout << "End:" << std::endl;
    print_lattice(lattice);
    std::cout << "Energy: " << energy(lattice, J) << "\n" << std::endl;

    return 0;
}
