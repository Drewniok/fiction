//
// Created by Jan Drewniok on 28.12.23.
//

#ifndef FICTION_MAXIMUM_DEFECT_INFLUENCE_POSITION_AND_DISTANCE_OF_SIDB_GATE_HPP
#define FICTION_MAXIMUM_DEFECT_INFLUENCE_POSITION_AND_DISTANCE_OF_SIDB_GATE_HPP

#include "fiction/algorithms/iter/bdl_input_iterator.hpp"
#include "fiction/algorithms/simulation/sidb/maximum_defect_influence_position_and_distance.hpp"
#include "fiction/traits.hpp"

#include <cassert>
#include <utility>
#include <cmath>

namespace fiction
{

/**
 * Parameters for the `maximum_defect_influence_position_and_distance_of_sidb_gate` algorithm.
 */
struct maximum_defect_influence_position_and_distance_of_sidb_gate_params
{
    /**
     * Parameters for the defect influence simulation.
     */
    maximum_defect_influence_position_and_distance_params defect_influence_params{};
    /**
     * Parameters for the input BDL iterator.
     */
    bdl_input_iterator_params bdl_iterator_params{};
};

/**
 * This function calculates the maximum influence position and distance of a defect on the ground state
 * of an SiDB gate layout. It iterates over all input combinations and finds the defect position at maximum position
 * that affects the gate's ground state.
 *
 * @Note The `maximum defect influence distance` describes the maximum distance at which a defect influences the ground state. It does
 * not check when the successful operation starts to fail, since a change in the ground state can still lead to an
 * operational gate.
 *
 * @tparam Lyt SiDB cell-level layout type.
 * @param lyt Layout to compute the maximum defect influence position and distance for.
 * @param params Parameters for the defect influence simulation and BDL pair detection.
 * @return A pair containing the maximum influence defect position and its distance from the layout/gate.
 */
template <typename Lyt>
[[nodiscard]] std::pair<typename Lyt::cell, double> maximum_defect_influence_position_and_distance_of_sidb_gate(
    const Lyt& lyt,
    const maximum_defect_influence_position_and_distance_of_sidb_gate_params& params = {}) noexcept
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
    static_assert(!has_offset_ucoord_v<Lyt>, "Lyt should not be based on offset coordinates");
    static_assert(!is_charge_distribution_surface_v<Lyt>, "Lyt cannot be a charge distribution surface");

    assert(lyt.num_pis() > 0 && "skeleton needs input cells");
    assert(lyt.num_pos() > 0 && "skeleton needs output cells");

    bdl_input_iterator<Lyt> bii{lyt, params.bdl_iterator_params};
    double                  maximum_defect_influence_distance = 0.0;
    cell<Lyt>               defect_cell{};

    // number of different input combinations
    for (auto i = 0u; i < std::pow(2,lyt.num_pis()); ++i, ++bii)
    {
        maximum_defect_influence_position_and_distance_stats stats_defect{};
        const auto                                           influence_cell_distance =
            maximum_defect_influence_position_and_distance(lyt, params.defect_influence_params, &stats_defect);

        if (influence_cell_distance.second > maximum_defect_influence_distance)
        {
            maximum_defect_influence_distance = influence_cell_distance.second;
            defect_cell                       = influence_cell_distance.first;
        }
    }
    return {defect_cell, maximum_defect_influence_distance};
}

}  // namespace fiction

#endif  // FICTION_MAXIMUM_DEFECT_INFLUENCE_POSITION_AND_DISTANCE_OF_SIDB_GATE_HPP
