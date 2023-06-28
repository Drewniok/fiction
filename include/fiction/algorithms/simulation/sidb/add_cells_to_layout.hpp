//
// Created by Jan Drewniok on 05.04.23.
//

#ifndef FICTION_RANDOM_LAYOUT_GENERATOR_HPP
#define FICTION_RANDOM_LAYOUT_GENERATOR_HPP

#include "fiction/algorithms/path_finding/distance.hpp"
#include "fiction/io/write_sqd_layout.hpp"
#include "fiction/layouts/cell_level_layout.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/technology/sidb_nm_position.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/hash.hpp"
#include "fiction/utils/layout_utils.hpp"

#include <cstdint>
#include <iostream>
#include <random>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace fiction
{

/**
 * Generates a random layout of SiDBs based on the provided parameters.
 *
 * @tparam Lyt The layout type.
 * @param params The parameters for generating the random layout.
 * @param lyt_skeleton A layout to which random cells are added (useful if you need to add random cells to a given
 * layout).
 * @return A randomly generated layout of SiDBs.
 */
template <typename Lyt>
Lyt add_cells_to_layout(const Lyt& lyt_skeleton, const std::vector<coordinate<Lyt>>& placed_cells)
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");

    Lyt lyt{};

    if (lyt_skeleton.num_cells() != 0)
    {
        lyt_skeleton.foreach_cell([&lyt, &lyt_skeleton](const auto& cell)
                                  { lyt.assign_cell_type(cell, lyt_skeleton.get_cell_type(cell)); });
    }

    bool     successful_generation = false;
    uint64_t attempt_counter       = 0;

    for (const auto& place_cell : placed_cells)
    {
        lyt.assign_cell_type(place_cell, Lyt::cell_type::NORMAL);
    }

    return lyt;
}

template <typename Lyt>
Lyt erase_cells_from_layout(const Lyt& lyt, const std::vector<coordinate<Lyt>>& placed_cells)
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");

    for (const auto& place_cell : placed_cells)
    {
        lyt.assign_cell_type(place_cell, Lyt::cell_type::EMPTY);
    }

    return lyt;
}

}  // namespace fiction

#endif  // FICTION_RANDOM_LAYOUT_GENERATOR_HPP
