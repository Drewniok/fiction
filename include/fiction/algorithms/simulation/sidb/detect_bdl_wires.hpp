//
// Created by Jan Drewniok on 15.07.24.
//

#ifndef FICTION_DETECT_BDL_WIRES_HPP
#define FICTION_DETECT_BDL_WIRES_HPP

#include "fiction/algorithms/simulation/sidb/detect_bdl_pairs.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/traits.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <optional>
#include <set>
#include <vector>

namespace fiction
{

/**
 * Possible directions of a BDL wire.
 */
enum class bdl_wire_direction
{
    /**
     * Defines the direction of the wire from north to south.
     * The starting point is defined by the position of input cells, and the ending point is defined by the
     * position of output cells.
     *
     * @note A wire running from west to east is interpreted as north-south direction.
     */
    NORTH_SOUTH,
    /**
     * Defines the direction of the wire from south to north.
     * The starting point is defined by the position of input cells, and the ending point is defined by the
     * position of output cells.
     *
     * @note A wire running from east to west is interpreted as south-north direction.
     */
    SOUTH_NORTH,
    /**
     * A wire without input or output cells does not exhibit a direction.
     */
    NO_DIRECTION
};

/**
 * An enumeration of the selection of different types of wires.
 */
enum class bdl_wire_selection
{
    /**
     * Select all BDL wires.
     */
    ALL,
    /**
     * Select only BDL wires that start with input cells.
     */
    INPUT,
    /**
     * Select only BDL wires that end with output cells.
     */
    OUTPUT
};

/**
 * This struct encapsulates parameters used for detecting BDL wires.
 */
struct detect_bdl_wires_params
{
    /**
     * A distance threshold, which is used to determine if two pairs of BDLs are part of the same wire. (unit: nm).
     */
    double threshold_bdl_interdistance = 2.0;
    /**
     * Parameters for the `detect_bdl_pairs` algorithm.
     */
    detect_bdl_pairs_params params_bdl_pairs{};
};

/**
 * This alias represents a BDL wire as a vector of BDL pairs.
 *
 * @tparam Lyt SiDB cell-level layout type.
 */
template <typename Lyt>
using bdl_wire = std::vector<bdl_pair<cell<Lyt>>>;

/**
 * Finds the first Binary-dot Logic (BDL) pair of a specified type in a BDL wire.
 *
 * @tparam Lyt SiDB cell-level layout type.
 * @param wire BDL wire containing BDL pairs to search through.
 * @param t Type of BDL pair to search for (`sidb_technology::cell_type::INPUT`, `sidb_technology::cell_type::OUTPUT`,
 * etc.).
 * @return Optional containing the first BDL pair with the specified type `t`,
 *         or `std::nullopt` if no such BDL pair is found.
 */
template <typename Lyt>
[[nodiscard]] std::optional<bdl_pair<cell<Lyt>>>
find_bdl_pair_in_wire_by_type(const bdl_wire<Lyt>& wire, const sidb_technology::cell_type t) noexcept
{
    auto it = std::find_if(wire.cbegin(), wire.cend(), [t](const auto& bdl) { return bdl.type == t; });

    if (it != wire.end())
    {
        return *it;  // Return the first BDL pair of the specified type found
    }

    return std::nullopt;  // Return std::nullopt if no BDL pair of the specified type is found
}
/**
 * Determines the direction of a Binary-dot Logic (BDL) wire based on the positions of input and output cells.
 *
 * @tparam Lyt SiDB cell-level layout type.
 * @param wire The BDL wire to analyze.
 * @return Direction of the wire: `NORTH_SOUTH`, `SOUTH_NORTH`, or `NO_DIRECTION`.
 */
template <typename Lyt>
[[nodiscard]] bdl_wire_direction determine_wire_direction(const bdl_wire<Lyt>& wire) noexcept
{
    // A BDL wire fewer than two BDL pairs does not have a direction.
    if (wire.size() < 2)
    {
        return bdl_wire_direction::NO_DIRECTION;
    }

    // check if the wire contains only normal cells
    if (std::all_of(wire.cbegin(), wire.cend(),
                    [](const auto& bdl) { return bdl.type == sidb_technology::cell_type::NORMAL; }))
    {
        return bdl_wire_direction::NO_DIRECTION;
    }

    // check if the wire exhibits an input BDL pair
    const auto input_exists = std::any_of(wire.cbegin(), wire.cend(), [](const auto& bdl)
                                          { return bdl.type == sidb_technology::cell_type::INPUT; });

    // check if the wire exhibits an output BDL pair
    const auto output_exists = std::any_of(wire.cbegin(), wire.cend(), [](const auto& bdl)
                                           { return bdl.type == sidb_technology::cell_type::OUTPUT; });

    if (input_exists && output_exists)
    {
        const auto input_bdl  = find_bdl_pair_in_wire_by_type<Lyt>(wire, sidb_technology::cell_type::INPUT).value();
        const auto output_bdl = find_bdl_pair_in_wire_by_type<Lyt>(wire, sidb_technology::cell_type::OUTPUT).value();

        // if the input BDL pair is above the output BDL pair, the wire runs from north to south
        if (input_bdl < output_bdl)
        {
            return bdl_wire_direction::NORTH_SOUTH;
        }

        // if the input BDL pair and the output BDL pair are at the same position, the wire does not have a direction
        if (input_bdl.equal_ignore_type(output_bdl))
        {
            assert(false && "input and output BDL pairs are at the same position");
            return bdl_wire_direction::NO_DIRECTION;
        }

        return bdl_wire_direction::SOUTH_NORTH;
    }

    // if the wire contains only input and normal BDL pairs
    if (input_exists)
    {
        const auto input_bdl = find_bdl_pair_in_wire_by_type<Lyt>(wire, sidb_technology::cell_type::INPUT);

        const auto bdl_above_exists =
            std::any_of(wire.cbegin(), wire.cend(), [&input_bdl](const auto& bdl) { return bdl > input_bdl; });

        if (bdl_above_exists)
        {
            return bdl_wire_direction::NORTH_SOUTH;
        }

        return bdl_wire_direction::SOUTH_NORTH;
    }

    const auto output_bdl = find_bdl_pair_in_wire_by_type<Lyt>(wire, sidb_technology::cell_type::OUTPUT);

    const auto bdl_above_exists =
        std::any_of(wire.cbegin(), wire.cend(), [&output_bdl](const auto& bdl) { return bdl > output_bdl; });

    if (bdl_above_exists)
    {
        return bdl_wire_direction::SOUTH_NORTH;
    }

    return bdl_wire_direction::NORTH_SOUTH;
}

/**
 * This function searches for the first Binary-dot Logic (BDL) pair in a given set of BDL pairs that is above a
 * specified BDL pair. The function returns the first BDL pair that meets the following criteria:
 * - The distance between the upper SiDB of the given BDL pair and the lower SiDB of the candidate BDL pair is less than
 * the specified inter-BDL distance.
 * - The candidate BDL pair is not equal to the given BDL pair (ignoring the type of SiDBs).
 * - The y-coordinate of the upper SiDB of the given BDL pair is less than the y-coordinate of the lower SiDB of the
 * candidate BDL pair.
 *
 * @tparam Lyt SiDB cell-level layout type.
 * @param given_bdl The BDL pair to find a neighbor for.
 * @param bdl_pairs A set of BDL pairs to search within.
 * @param inter_bdl_distance The maximum allowable inter distance between the lower SiDB of the given BDL pair and the
 * upper SiDB of the candidate BDL pair.
 * @return A std::optional containing the first BDL pair that meets the criteria, or std::nullopt if no such pair is
 * found.
 */
template <typename Lyt>
[[nodiscard]] std::optional<bdl_pair<cell<Lyt>>> find_bdl_neighbor_above(const bdl_pair<cell<Lyt>>&           given_bdl,
                                                                         const std::set<bdl_pair<cell<Lyt>>>& bdl_pairs,
                                                                         const double inter_bdl_distance) noexcept
{
    for (const auto& bdl : bdl_pairs)
    {
        if (sidb_nm_distance<Lyt>(Lyt{}, given_bdl.upper, bdl.lower) < inter_bdl_distance &&
            !given_bdl.equal_ignore_type(bdl) && given_bdl > bdl)
        {
            return std::optional<bdl_pair<cell<Lyt>>>(bdl);
        }
    }

    return std::nullopt;
};

/**
 * This function searches for the first Binary-dot Logic (BDL) pair in a given set of BDL pairs that is below a
 * specified BDL pair. The function returns the first BDL pair that meets the following criteria:
 * - The distance between the lower SiDB of the given BDL pair and the upper SiDB of the candidate BDL pair is less than
 * the specified inter-BDL distance.
 * - The candidate BDL pair is not equal to the given BDL pair (ignoring the type of SiDBs).
 * - The y-coordinate of the lower SiDB of the given BDL pair is less than the y-coordinate of the upper SiDB of the
 * candidate BDL pair.
 *
 * @tparam Lyt SiDB cell-level layout type.
 * @param given_bdl The BDL pair to find a neighbor for.
 * @param bdl_pairs A set of BDL pairs to search within.
 * @param inter_bdl_distance The maximum allowable inter distance between the lower SiDB of the given BDL pair and the
 * upper SiDB of the candidate BDL pair.
 * @return A std::optional containing the first BDL pair that meets the criteria, or std::nullopt if no such pair is
 * found.
 */
template <typename Lyt>
[[nodiscard]] std::optional<bdl_pair<cell<Lyt>>> find_bdl_neighbor_below(const bdl_pair<cell<Lyt>>&           given_bdl,
                                                                         const std::set<bdl_pair<cell<Lyt>>>& bdl_pairs,
                                                                         const double inter_bdl_distance) noexcept
{
    for (const auto& bdl : bdl_pairs)
    {
        if (sidb_nm_distance<Lyt>(Lyt{}, given_bdl.lower, bdl.upper) < inter_bdl_distance &&
            given_bdl.not_equal_ignore_type(bdl) && given_bdl < bdl)
        {
            return std::optional<bdl_pair<cell<Lyt>>>(bdl);
        }
    }

    return std::nullopt;
};

/**
 * This function identifies BDL wires in a given SiDB cell-level layout by detecting BDL pairs and linking them based on
 * their spatial relationships. The function supports selection of different types of wires, such as input wires, output
 * wires, or all wires.
 *
 * @tparam Lyt SiDB cell-level layout type.
 * @param lyt The SiDB layout to detect BDL wires in.
 * @param wire_selection The type of wires to detect, specified by the `bdl_wire_selection` enum. Default is
 * `bdl_wire_selection::ALL`.
 * @return A vector of BDL wires, where each wire is represented as a vector of BDL pairs.
 */
template <typename Lyt>
[[nodiscard]] std::vector<bdl_wire<Lyt>>
detect_bdl_wires(const Lyt& lyt, const detect_bdl_wires_params& params = {},
                 const bdl_wire_selection wire_selection = bdl_wire_selection::ALL) noexcept
{
    std::set<bdl_pair<cell<Lyt>>> bdl_pairs{};

    const auto all_input_bdls  = detect_bdl_pairs(lyt, Lyt::cell_type::INPUT, params.params_bdl_pairs);
    const auto all_output_bdls = detect_bdl_pairs(lyt, Lyt::cell_type::OUTPUT, params.params_bdl_pairs);
    const auto all_normal_bdls = detect_bdl_pairs(lyt, Lyt::cell_type::NORMAL, params.params_bdl_pairs);

    for (const auto& bdl : all_input_bdls)
    {
        bdl_pairs.insert(bdl);
    }

    for (const auto& bdl : all_output_bdls)
    {
        bdl_pairs.insert(bdl);
    }

    for (const auto& bdl : all_normal_bdls)
    {
        bdl_pairs.insert(bdl);
    }

    std::vector<bdl_wire<Lyt>> wires{};

    while (!bdl_pairs.empty())
    {
        bdl_wire<Lyt> wire{};

        bool neighbor_bdl_found = true;

        // choose the current bdl pair
        auto current_bdl_pair = *bdl_pairs.cbegin();
        // store the front bdl pair as an initial bdl pair
        const auto initial_bdl_pair = current_bdl_pair;

        // add current bdl pair to wire
        wire.push_back(current_bdl_pair);
        // delete bdl pair from a set
        bdl_pairs.erase(current_bdl_pair);

        while (neighbor_bdl_found)
        {
            // determine a neighbor bdl pair below
            const auto neighbor_below =
                find_bdl_neighbor_below<Lyt>(current_bdl_pair, bdl_pairs, params.threshold_bdl_interdistance);

            if (neighbor_below.has_value())
            {
                // add a neighbor bdl pair to wire
                wire.push_back(neighbor_below.value());
                // delete bdl pair from a set
                bdl_pairs.erase(neighbor_below.value());
                // update current bdl pair
                current_bdl_pair = neighbor_below.value();
            }
            // If no neighbor bdl pair is found below, the existence of a neighbor bdl pair above is checked
            else
            {
                current_bdl_pair = initial_bdl_pair;

                if (find_bdl_neighbor_above<Lyt>(current_bdl_pair, bdl_pairs, params.threshold_bdl_interdistance)
                        .has_value())
                {
                    wire.push_back(neighbor_below.value());
                    bdl_pairs.erase(neighbor_below.value());
                }
                else
                {
                    neighbor_bdl_found = false;
                    wires.push_back(wire);
                }
            }
        }
    }

    if (wire_selection == bdl_wire_selection::INPUT)
    {
        std::vector<bdl_wire<Lyt>> input_wires{};

        std::set<uint64_t> input_wire_length{};

        for (const auto& wire : wires)
        {
            for (const auto& bdl : wire)
            {
                if (bdl.type == sidb_technology::cell_type::INPUT)
                {
                    input_wires.push_back(wire);
                    input_wire_length.insert(wire.size());
                }
            }
        }

        if (input_wires.empty())
        {
            return {};
        }

        assert(input_wire_length.size() < 2 && "input wires have different length");

        return input_wires;
    }

    if (wire_selection == bdl_wire_selection::OUTPUT)
    {
        std::vector<bdl_wire<Lyt>> output_wires{};

        std::set<uint64_t> output_wire_length{};

        for (const auto& wire : wires)
        {
            for (const auto& bdl : wire)
            {
                if (bdl.type == sidb_technology::cell_type::OUTPUT)
                {
                    output_wires.push_back(wire);
                    output_wire_length.insert(wire.size());
                }
            }
        }

        assert(output_wire_length.size() < 2 && "output wires have different length");

        return output_wires;
    }

    return wires;
}

}  // namespace fiction

#endif  // FICTION_DETECT_BDL_WIRES_HPP