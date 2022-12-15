//
// Created by Jan Drewniok on 12.12.22.
//

#ifndef FICTION_NM_POSITION_HPP
#define FICTION_NM_POSITION_HPP

#include "fiction/algorithms/simulation_sidb/simulation_parameter.hpp"
#include "fiction/traits.hpp"

namespace fiction
{
/**
 * SiQAD coordinates are converted to a nm location on the Si-substrate by taking Silicon's lattice constants into
 * account (lat_a, lat_b, lat_c).
 *
 * @tparam Lyt cell layout type (SiQAD coordinates are required).
 * @tparam Dist Floating-point type for the distance.
 * @param lyt Layout.
 * @return nm position
 */

template <typename Lyt, typename Dist = double>
std::pair<Dist, Dist> nm_position(const cell<Lyt>& c)
{
    const auto x = static_cast<Dist>(c.x * fiction::simulation_params{}.lat_a);
    const auto y =
        static_cast<Dist>(c.y * fiction::simulation_params{}.lat_b + c.z * fiction::simulation_params{}.lat_c);
    return std::pair(x, y);
}

#endif  // FICTION_NM_POSITION_HPP

}  // namespace fiction
