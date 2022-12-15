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
 * @tparam Lyt cell-level layout.
 * @tparam Energy data type of the system's total electrostatic potential energy.
 * @tparam Potential data type of the electrostatic potential.
 * @param lyt charge distribution layout
 * @param loc_pot local electrostatic potential
 */
template <typename Lyt, typename Energy = double, typename Potential = double>
Energy system_energy(const charge_distribution_surface<Lyt>& lyt, const local_pot<Lyt, Potential> loc_pot)

{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is cell-level layout");
    static_assert(std::is_floating_point_v<Energy>, "Potential is not a floating-point type");
    static_assert(std::is_floating_point_v<Potential>, "Potential is not a floating-point type");
    Energy total_energy = 0;
    for (auto& it : loc_pot)
    {
        total_energy += 0.5 * it.second * transform_to_sign(lyt.get_charge_state(it.first));
    }
    return total_energy;
}

}  // namespace fiction

#endif  // FICTION_SYSTEM_ENERGY_HPP
