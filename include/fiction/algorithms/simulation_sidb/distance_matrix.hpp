#ifndef FICTION_DISTANCE_MATRIX_HPP
#define FICTION_DISTANCE_MATRIX_HPP
//
// Created by Jan Drewniok on 07.12.22.
//

#include "fiction/traits.hpp"
#include "fiction/types.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/layouts/coordinates.hpp"
#include "fiction/algorithms/simulation_sidb/constants.hpp"
#include "fiction/algorithms/path_finding/distance.hpp"

namespace fiction
{

template <typename Lyt, typename Dist = double>
using distance_matrix = std::unordered_map<std::pair<const cell<Lyt>, const cell<Lyt>>, Dist>;

template <typename Lyt, typename Dist = double>
distance_matrix<Lyt,Dist> distance_SiDBs(const Lyt& lyt)
{
    distance_matrix<Lyt, Dist> distance_values{};
    lyt.foreach_cell(
        [&distance_values, lyt](const auto& c1)
        {
            lyt.foreach_cell([&distance_values, c1, lyt](const auto& c2)
                             { distance_values.insert(std::pair(c1, c2), distance_SiDB_pair<Lyt>(lyt, c1, c2)); });
        });
    return distance_values;
};

} // namespace fiction

#endif  // FICTION_DISTANCE_MATRIX_HPP
