//
// Created by Jan Drewniok on 18.01.23.
//

#ifndef FICTION_CHECK_GROUNDSTATE_HPP
#define FICTION_CHECK_GROUNDSTATE_HPP

#include "fiction/algorithms/simulation_sidb/ExGS.hpp"
#include "fiction/algorithms/simulation_sidb/quicksim.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"

namespace fiction
{

/**
* @brief This function checks if the ground state function is found by the quicksim algorithm.
* @tparam Lyt The type of the charge distribution surface layout.
* @param result_new_ap The set of valid charge distribution surfaces obtained from the new quicksim algorithm (see quicksim.hpp).
* @param result_exact The set of valid charge distribution surfaces obtained from the exact method (ExGS, see ExGS.hpp).
* @return Returns true if the relative difference between the ground state energies of the two sets is less than
0.00001, false otherwise.
*/
template <typename Lyt>
bool check_groundstate(const quicksim_stats<Lyt>& result_new_ap, const exgs_stats<Lyt>& result_exact)
{
    if (result_exact.valid_lyts.empty())
    {
        return false;
    }
    auto min_energy_exact  = minimum_energy(result_exact.valid_lyts);
    auto min_energy_new_ap = minimum_energy(result_new_ap.valid_lyts);

    return std::abs(min_energy_exact - min_energy_new_ap) / min_energy_exact < 0.00001;
}

}  // namespace fiction

#endif  // FICTION_CHECK_GROUNDSTATE_HPP
