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

template <typename Lyt, typename Potential = double>
using local_pot = std::unordered_map<cell<Lyt>, Potential>;

template <typename Lyt, typename Potential = double>
local_pot<Lyt, Potential> local_potential(Lyt &lyt, const potential_matrix<Lyt,Potential> &pot_mat)
{
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