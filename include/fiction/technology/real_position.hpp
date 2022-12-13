//
// Created by Jan Drewniok on 12.12.22.
//

#ifndef FICTION_REAL_POSITION_HPP
#define FICTION_REAL_POSITION_HPP

#include "fiction/algorithms/simulation_sidb/simulation_parameter.hpp"
#include "fiction/traits.hpp"

namespace fiction
{
/**
 * SiQAD coordinates are converted into a position on the Si-H substrate
 *
 * @param c cell
 * @return pair (x,y)
 */
template <typename Lyt, typename Dist = double>
std::pair<Dist, Dist> real_position(const cell<Lyt>& c)
{
    const auto x = static_cast<Dist>(c.x * fiction::simulation_parameter{}.lat_a);
    const auto y =
        static_cast<Dist>(c.y * fiction::simulation_parameter{}.lat_b + c.z * fiction::simulation_parameter{}.lat_c);
    return std::pair(x, y);
}

}  // namespace fiction

#endif  // FICTION_REAL_POSITION_HPP
