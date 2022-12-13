#ifndef FICTION_DISTANCE_MATRIX_HPP
#define FICTION_DISTANCE_MATRIX_HPP
//
// Created by Jan Drewniok on 07.12.22.
//

#include "fiction/algorithms/path_finding/distance.hpp"
#include "fiction/algorithms/simulation_sidb/constants.hpp"
#include "fiction/layouts/coordinates.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"

namespace fiction
{
/**
 * The distance matrix is an unordered map with pairs of cells as key and the corresponding euclidean distance as value.
 */
template <typename Lyt, typename Dist = double>
using distance_matrix = std::unordered_map<std::pair<const cell<Lyt>, const cell<Lyt>>, Dist>;

/**
 * The distance matrix stores all euclidean inter-distances. The euclidean distance between
 * two identical SiDBs is set to 0. \f$ D = \sqrt{(x_1 - x_2)^2 + (y_1 - y_2)^2} \f$
 *
 * @tparam Lyt Coordinate layout type (SiQAD coordinates are required).
 * @tparam Dist Floating-point type for the distance.
 * @param lyt Layout.
 * @return Distance matrix
 */
template <typename Lyt, typename Dist = double>
distance_matrix<Lyt, Dist> distance_SiDBs(const Lyt& lyt)
{
    distance_matrix<Lyt, Dist> distance_values{};
    lyt.foreach_cell(
        [&distance_values, lyt](const auto& c1)
        {
            lyt.foreach_cell(
                [&distance_values, c1, lyt](const auto& c2) {
                    distance_values.insert(
                        std::make_pair(std::make_pair(c1, c2), distance_SiDB_pair<Lyt>(lyt, c1, c2)));
                });
        });

    return distance_values;
};

}  // namespace fiction

#endif  // FICTION_DISTANCE_MATRIX_HPP
