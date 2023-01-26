//
// Created by Jan Drewniok on 18.01.23.
//

#ifndef FICTION_CHECK_GROUNDSTATE_HPP
#define FICTION_CHECK_GROUNDSTATE_HPP

#include "fiction/algorithms/simulation_sidb/exhaustive_ground_state_simulation.hpp"
#include "fiction/algorithms/simulation_sidb/quicksim.hpp"
#include <cmath>

namespace fiction
{

/**
* This function checks if the ground state is found by the quicksim algorithm.
*
* @tparam Lyt cell-level layout.
* @param result_new_ap All found physically valid charge distribution surfaces obtained with the new quicksim algorithm (see quicksim.hpp).
* @param result_exact All valid charge distribution surfaces (ExGS, see ExGS.hpp).
* @return Returns true if the relative difference between the lowest energies of the two sets is less than
0.00001, false otherwise.
*/
template <typename Lyt>
bool check_groundstate(const quicksim_stats<Lyt>& quicksim_results, const exgs_stats<Lyt>& exhaustive__results)
{
    if (exhaustive__results.valid_lyts.empty())
    {
        return false;
    }
    const auto min_energy_exact  = minimum_energy(exhaustive__results.valid_lyts);
    const auto min_energy_new_ap = minimum_energy(quicksim_results.valid_lyts);

    return std::abs(min_energy_exact - min_energy_new_ap) / min_energy_exact < 0.00001;
}

}  // namespace fiction

#endif  // FICTION_CHECK_GROUNDSTATE_HPP
