//
// Created by Jan Drewniok on 12.12.22.
//

#ifndef FICTION_NM_POSITION_HPP
#define FICTION_NM_POSITION_HPP

#include "fiction/algorithms/simulation_sidb/simulation_parameters.hpp"
#include "fiction/traits.hpp"

namespace fiction
{
/**
 * SiQAD coordinates are converted to a nm location on the Si-substrate by taking Silicon's lattice constants into
 * account (see simulation_parameters.hpp).
 *
 * @tparam Lyt cell layout type (SiQAD coordinates are required).
 * @tparam Dist Data type for the distance.
 * @param c cell of the layout Lyt.
 * @return nm position
 */

template <typename Lyt, typename Dist = double>
std::pair<Dist, Dist> nm_position(const cell<Lyt>& c, const simulation_params &sim_params = simulation_params{})
{
    static_assert(std::is_same_v<cell<Lyt>, siqad::coord_t>, "Lyt is not based on SiQAD coordinates");
    static_assert(is_cell_level_layout_v<Lyt>,"Lyt is not a cell-level layout");
    static_assert(std::is_floating_point_v<Dist>, "Dist is not a floating-point type");
    const auto x = static_cast<Dist>(c.x * sim_params.lat_a);
    const auto y =
        static_cast<Dist>(c.y * sim_params.lat_b + c.z * sim_params.lat_c);
    return std::pair(x, y);
}

#endif  // FICTION_NM_POSITION_HPP

}  // namespace fiction
