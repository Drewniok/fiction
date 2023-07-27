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
 * This struct stores the parameters for the *generate_random_layout* algorithm.
 */
template <typename Lyt>
struct random_layout_params
{
    /**
     * Two coordinates that span the region where SiDBs may be placed (order is not important).
     */
    std::pair<typename Lyt::cell, typename Lyt::cell> coordinate_pair;
    /**
     * Number of SiDBs that are placed on the layout.
     */
    uint64_t number_of_sidbs = 0;
    /**
     * If positively charged SiDBs should be prevented, SiDBs are not placed closer than the minimal_spacing.
     */
    bool prevent_positive_charges = true;
    /**
     * If positively charged SiDBs should be prevented, SiDBs are not placed closer than this value (2 cells as
     * Euclidean distance by default).
     */
    double minimal_spacing = 3;
    /**
     * Maximal number of steps to place the given number of SiDBs.
     */
    uint64_t maximal_attempts = 10E6;
};

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
Lyt generate_random_layout(const random_layout_params<Lyt>& params, const Lyt& lyt_skeleton)
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");

    const uint64_t number_of_sidbs_of_final_layout = lyt_skeleton.num_cells() + params.number_of_sidbs;

    Lyt lyt{};

    if (lyt_skeleton.num_cells() != 0)
    {
        lyt_skeleton.foreach_cell([&lyt, &lyt_skeleton](const auto& cell)
                                  { lyt.assign_cell_type(cell, lyt_skeleton.get_cell_type(cell)); });
    }

    if (has_get_sidb_defect_v<Lyt>)
    {
        lyt_skeleton.foreach_sidb_defect(
            [&lyt, &lyt_skeleton](const auto& defect)
            { lyt.assign_sidb_defect(defect.first, lyt_skeleton.get_sidb_defect(defect.first)); });
    }

    bool     successful_generation = false;
    uint64_t attempt_counter       = 0;

    // Stops if either all SiDBs are placed or the maximum number of attempts were performed.
    while (lyt.num_cells() < number_of_sidbs_of_final_layout && attempt_counter < params.maximal_attempts)
    {
        const auto random_coord = random_coordinate(params.coordinate_pair.first, params.coordinate_pair.second);

        bool constraint_violation_positive_sidbs = false;

        if (params.prevent_positive_charges)
        {
            // Checks if the new coordinate is not closer than 2 cells (Euclidean distance) from an already
            // placed SiDB.
            lyt.foreach_cell(
                [&lyt, &random_coord, &constraint_violation_positive_sidbs, &params](const auto& c1)
                {
                    if (euclidean_distance<Lyt>(lyt, c1, random_coord) < params.minimal_spacing)
                    {
                        constraint_violation_positive_sidbs = true;
                    }
                });
            lyt.foreach_sidb_defect(
                [&lyt, &random_coord, &constraint_violation_positive_sidbs, &params](const auto& c1)
                {
                    if (euclidean_distance<Lyt>(lyt, c1.first, random_coord) < params.minimal_spacing)
                    {
                        constraint_violation_positive_sidbs = true;
                    }
                });
        }
        // If the constraint that no positive SiDBs occur is satisfied, the SiDB is added to the layout.
        if (!constraint_violation_positive_sidbs)
        {
            lyt.assign_cell_type(random_coord, Lyt::cell_type::NORMAL);
        }
        attempt_counter += 1;
    }

    return lyt;
}

/**
 * Generates multiple unique random layouts of SiDBs based on the provided parameters.
 *
 * @tparam Lyt The layout type.
 * @param params The parameters for generating the random layouts.
 * @param lyt_skeleton A layout to which random cells are added (useful if you need to add random cells to a given
 * layout).
 * @param number_of_unique_generated_layouts The desired number of unique layouts to be generated.
 * @param maximal_attempts The maximum number of attempts allowed to generate a unique layout (default: \f$ 10^{6} \f$).
 * @return A vector containing the unique randomly generated layouts.
 */
template <typename Lyt>
std::vector<Lyt> generate_multiple_random_layout(const random_layout_params<Lyt>& params, const Lyt& lyt_skeleton,
                                                 const uint64_t number_of_unique_generated_layouts,
                                                 const uint64_t maximal_attemps = 10E6)
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");

    std::vector<Lyt> unique_lyts{};
    unique_lyts.reserve(number_of_unique_generated_layouts);
    uint64_t counter = 0;
    while (unique_lyts.size() < number_of_unique_generated_layouts && counter < maximal_attemps)
    {
        const auto random_lyt = generate_random_layout(params, lyt_skeleton);

        uint64_t identical_layout_counter = 0;
        // Checks if two layouts have the same cell.

        for (const auto& old_lyt : unique_lyts)
        {
            uint64_t identical_cell_counter = 0;
            old_lyt.foreach_cell(
                [&identical_cell_counter, random_lyt](const auto& cell_old)
                {
                    random_lyt.foreach_cell(
                        [&identical_cell_counter, &cell_old](const auto& cell_new)
                        {
                            if (cell_new == cell_old)
                            {
                                identical_cell_counter += 1;
                            }
                        });
                });

            // all cells are identical, so the new layout is a duplicate.
            if (identical_cell_counter == random_lyt.num_cells())
            {
                identical_layout_counter += 1;
            }
        }

        if (identical_layout_counter == 0)
        {
            unique_lyts.push_back(random_lyt);
        }

        if (unique_lyts.size() % 1000 == 0)
        {
            std::cout << unique_lyts.size() << std::endl;
        }

        counter += 1;
    }

    return unique_lyts;
}

}  // namespace fiction

#endif  // FICTION_RANDOM_LAYOUT_GENERATOR_HPP
