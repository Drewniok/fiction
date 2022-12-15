//
// Created by Jan Drewniok on 13.12.22.
//

#ifndef FICTION_SYSTEM_ENERGY_HPP
#define FICTION_SYSTEM_ENERGY_HPP

#include <fiction/algorithms/simulation_sidb/local_potential.hpp>
#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/technology/charge_distribution_surface.hpp>

namespace fiction
{
/**
 * Calculate the system's total electrostatic potential energy.
 *
 * @param lyt charge distribution layout
 * @param loc_pot local electrostatic potential
 */
template <typename Lyt, typename Energy = double, typename Potential = double>
Energy system_energy(const charge_distribution_surface<Lyt>& lyt, const local_pot<Lyt, Potential> loc_pot)

{
    Energy total_energy = 0;
    for (auto& it : loc_pot)
    {
        total_energy += 0.5 * it.second * transform_to_sign(lyt.get_charge_state(it.first));
    }
    return total_energy;
}

}  // namespace fiction

#endif  // FICTION_SYSTEM_ENERGY_HPP
