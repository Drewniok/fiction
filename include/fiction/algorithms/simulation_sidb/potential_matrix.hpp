//
// Created by Jan Drewniok on 07.12.22.
//

#ifndef FICTION_POTENTIAL_MATRIX_HPP
#define FICTION_POTENTIAL_MATRIX_HPP

#include "fiction/algorithms/simulation_sidb/distance_matrix.hpp"
#include "fiction/technology/constants.hpp"
#include "fiction/technology/electrostatic_potential.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"

namespace fiction
{
/**
 * The potential matrix is an unordered map with pairs of cells as key and the corresponding electrostatic potential as
 * value.
 */
template <typename Lyt, typename Potential = double>
using potential_matrix = std::unordered_map<std::pair<cell<Lyt>, cell<Lyt>>, Potential>;

/**
 * All electrostatic inter-potentials between two cells are calculated. The Euclidean distance is provided by the
 * distance matrix. Electrostatic potential between identical cells is set to 0.
 *
 * @tparam Lyt Coordinate layout type (SiQAD coordinates are required).
 * @tparam Dist Floating-point type for the distance.
 * @tparam Potential Floating-point type for the electrostatic potential.
 * @param lyt Layout.
 * @return Potential matrix
 */
template <typename Lyt, typename Potential = double, typename Dist = double>
potential_matrix<Lyt, Potential> potential_sidbs(const distance_matrix<Lyt, Dist>& dist)
{
    potential_matrix<Lyt, Potential> potential_values{};
    for (const auto& it : dist)
    {
        potential_values.insert(std::make_pair(it.first, potential_sidb_pair<Potential>(it.second)));
    }
    return potential_values;
};

}  // namespace fiction

#endif  // FICTION_POTENTIAL_MATRIX_HPP
