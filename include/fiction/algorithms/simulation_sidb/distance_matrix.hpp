//
// Created by Jan Drewniok on 07.12.22.
//

#ifndef FICTION_DISTANCE_MATRIX_HPP
#define FICTION_DISTANCE_MATRIX_HPP

#include "fiction/algorithms/path_finding/distance.hpp"
#include "fiction/layouts/coordinates.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/technology/constants.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"

namespace fiction
{
/**
 * The distance matrix is an unordered map with pairs of cells as key and the corresponding euclidean distance as value.
 */
template <typename Lyt, typename Dist = double>
using distance_matrix = std::unordered_map<std::pair<cell<Lyt>, cell<Lyt>>, Dist>;

/**
 * The Euclidean distance for every cell (needs to exhibit an assigned cell type) pair in the layout is calculated and
 * stored.
 *
 * @tparam Lyt Coordinate layout type (SiQAD coordinates are required).
 * @tparam Dist Floating-point type for the distance.
 * @param lyt Layout.
 * @return Distance matrix
 */
template <typename Lyt, typename Dist = double>
distance_matrix<Lyt, Dist> initialize_sidb_distance_matrix(const Lyt& lyt)
{
    static_assert(std::is_same_v<cell<Lyt>, siqad::coord_t>, "Cell level layout is not based on siqad coordinates");
    distance_matrix<Lyt, Dist> distance_values{};
    lyt.foreach_cell(
        [&distance_values, lyt](const auto& c1)
        {
            lyt.foreach_cell(
                [&distance_values, c1, lyt](const auto& c2) {
                    distance_values.insert(
                        std::make_pair(std::make_pair(c1, c2), distance_sidb_pair<Lyt>(lyt, c1, c2)));
                });
        });

    return distance_values;
};

}  // namespace fiction

#endif  // FICTION_DISTANCE_MATRIX_HPP
