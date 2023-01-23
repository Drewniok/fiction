//
// Created by Jan Drewniok on 11.01.23.
//

#ifndef FICTION_QUICKSIM_HPP
#define FICTION_QUICKSIM_HPP

#include "fiction/algorithms/simulation_sidb/get_energy_dist.hpp"
#include "fiction/algorithms/simulation_sidb/minimum_energy.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/traits.hpp"

#include <fmt/format.h>
#include <algorithm>
#include <vector>
#include <iostream>

namespace fiction
{
/**
 * This struct stores the simulation runtime and all physically valid charge layouts gained by the quicksim simulation (see quicksim.hpp).
 *
 * @paramt Lyt cell-level layout.
 */
template <typename Lyt>
struct quicksim_stats
{
    double                                        time_total{0};
    std::vector<charge_distribution_surface<Lyt>> valid_lyts{};

    void report(std::ostream& out = std::cout)
    {
        out << fmt::format("[i] total time  = {:.2f} millisecs\n", time_total / 1000);
        if (!get_statistics<Lyt>(valid_lyts).empty())
        {
            for (auto [energy, count] : get_statistics<Lyt>(valid_lyts))
            {
                out << fmt::format("[i] the lowest state energy is  = {:.4f} \n", minimum_energy(valid_lyts));
                out << fmt::format("energy: {} | occurance: {} \n", energy, count);
            }
        }
        else
        {
            std::cout << "no state found" << std::endl;
        }
        std::cout << "_____________________________________________________ \n";
    }
};

/**
 * quicksim determines physically valid charge configurations (with minimal energy) of a given (already initialized) charge distribution layout. Depending on the simulation paramaters, the ground state is found with a certain probability after one run.
 *
 * @tparam Lyt cell-level layout.
 * @param lyt charge distribution layout.
 * @param ps struct that stores the simulation results (simulation runtime, and all physically valid charge distribution layouts).
 * @param physical_params physical parameters, they are material-specific and may vary from experiment to experiment.
 */
template <typename Lyt>
void quicksim(charge_distribution_surface<Lyt>& lyt, quicksim_stats<Lyt>& ps,
              const physical_params& phys_params = physical_params{}, const uint64_t iteration_steps = 80,
              const double alpha = 0.7)
{

    // set the given physical parameters
    lyt.set_physical_parameters(phys_params);

    // measure run time
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<charge_distribution_surface<Lyt>> result{};

    lyt.set_charge_states(sidb_charge_state::NEUTRAL);
    lyt.local_potential();
    lyt.system_energy();
    lyt.validity_check();

    if (lyt.get_validity())
    {
        charge_distribution_surface<Lyt> lyt_new{lyt};
        ps.valid_lyts.push_back(lyt_new);
    }

    float best_energy = MAXFLOAT;
    auto  bound       = static_cast<uint64_t>(round(0.6 * static_cast<double>(lyt.num_cells())));
    for (uint64_t z = 0u; z < iteration_steps; z++)
    {
        for (uint64_t i = 0u; i < bound; i++)
        {
            std::vector<uint64_t> index_start = {i};
            lyt.set_charge_states(sidb_charge_state::NEUTRAL);
            lyt.assign_charge_state_index(i, sidb_charge_state::NEGATIVE);
            lyt.local_potential();
            lyt.system_energy();

            for (uint64_t num = 0; num < lyt.num_cells() / 1.5; num++)
            {
                lyt.adjacent_search(alpha, index_start);
                lyt.validity_check();

                if (lyt.get_validity() && (lyt.get_system_energy() <= best_energy))
                {
                    charge_distribution_surface<Lyt> lyt_new{lyt};
                    ps.valid_lyts.push_back(lyt_new);
                }
            }
        }
    }
    auto end      = std::chrono::high_resolution_clock::now();
    ps.time_total = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

}  // namespace fiction

#endif  // FICTION_QUICKSIM_HPP
