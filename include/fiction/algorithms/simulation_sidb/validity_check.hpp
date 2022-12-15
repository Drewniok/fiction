//
// Created by Jan Drewniok on 13.12.22.
//

#ifndef FICTION_VALIDITY_CHECK_HPP
#define FICTION_VALIDITY_CHECK_HPP

#include <fiction/algorithms/simulation_sidb/local_potential.hpp>
#include <fiction/algorithms/simulation_sidb/potential_matrix.hpp>
#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/technology/charge_distribution_surface.hpp>

namespace fiction
{
/**
 * Check the physically validity (configuration and population stability).
 *
 * @tparam Lyt cell-level layout.
 * @tparam Potential Data type of the electrostatic potential.
 * @param lyt charge distribution layout.
 * @param loc_pot local electrostatic potential.
 * @param pot_mat electrostatic potential matrix.
 * @param sim_params simulation parameters.
 */
template <typename Lyt, typename Potential = double>
bool validity_check(const charge_distribution_surface<Lyt>& lyt, const local_pot<Lyt, Potential>& loc_pot,
                    const potential_matrix<charge_distribution_surface<Lyt>, Potential>& pot_mat,
                    const simulation_params& sim_params = simulation_params{})
{
    bool valid = false;
    for (auto& it : loc_pot)
    {
        valid = (((lyt.get_charge_state(it.first) == sidb_charge_state::NEGATIVE) &&
                  ((-it.second + sim_params.mu) < physical_sim_constants::POP_STABILITY_ERR)) ||
                 ((lyt.get_charge_state(it.first) == sidb_charge_state::POSITIVE) &&
                  ((-it.second + sim_params.mu_p) > physical_sim_constants::POP_STABILITY_ERR)) ||
                 ((lyt.get_charge_state(it.first) == sidb_charge_state::NEUTRAL) &&
                  ((-it.second + sim_params.mu) > physical_sim_constants::POP_STABILITY_ERR) &&
                  (-it.second + sim_params.mu_p) < physical_sim_constants::POP_STABILITY_ERR));

        if (!valid)
        {
            return false;
        }
    }

    auto hopDel = [&lyt, &loc_pot, &pot_mat](const cell<Lyt>& c1, const cell<Lyt>& c2)
    {
        int dn_i = (lyt.get_charge_state(c1) == sidb_charge_state::NEGATIVE) ? 1 : -1;
        int dn_j = -dn_i;

        return loc_pot.at(c1) * dn_i + loc_pot.at(c1) * dn_j - pot_mat.at(std::make_pair(c1, c2)) * 1;
    };

    for (auto& it : loc_pot)
    {
        if (lyt.get_charge_state(it.first) == sidb_charge_state::POSITIVE)
        {
            continue;
        }

        for (auto& it_second : loc_pot)
        {
            auto E_del = hopDel(it.first, it_second.first);
            if ((transform_to_sign(lyt.get_charge_state(it.first)) >
                 transform_to_sign(lyt.get_charge_state(it.first))) &&
                (E_del < -physical_sim_constants::POP_STABILITY_ERR))
            {
                return false;
            }
        }
    }
    return true;
}

}  // namespace fiction

#endif  // FICTION_VALIDITY_CHECK_HPP
