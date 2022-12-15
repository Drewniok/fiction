//
// Created by Jan Drewniok on 13.12.22.
//

#ifndef FICTION_LOCAL_POTENTIAL_HPP
#define FICTION_LOCAL_POTENTIAL_HPP

#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/algorithms/simulation_sidb/potential_matrix.hpp>
#include <fiction/technology/sidb_charge_state.hpp>

namespace fiction
{
/**
 * Calculate local electrostatic potential at each SiDB position.
 *
 * @tparam Lyt cell-level layout.
 * @tparam Potential data type of the electrostatic potential.
 * @param lyt charge distribution layout.
 * @param pot_mat electrostatic potential matrix.
 */
template <typename Lyt, typename Potential = double>
using local_pot = std::unordered_map<cell<Lyt>, Potential>;

template <typename Lyt, typename Potential = double>
local_pot<Lyt, Potential> local_potential(const Lyt &lyt, const potential_matrix<Lyt,Potential> &pot_mat)
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is cell-level layout");
    static_assert(std::is_floating_point_v<Potential>, "Potential is not a floating-point type");

    local_pot<Lyt, Potential> loc_pot{};
    lyt.foreach_cell([&loc_pot, &pot_mat, &lyt](const auto& c){

        Potential collect = 0;
        for (auto& it:pot_mat)
        {

            if (it.first.first == c)
            {
                collect += it.second * transform_to_sign(lyt.get_charge_state(it.first.second));
            }
            else {
                continue;
}
        }
        loc_pot.insert(std::make_pair(c,collect));
    });
    return loc_pot;
}

} // namespace fiction

#endif  // FICTION_LOCAL_POTENTIAL_HPP
