//
// Created by Jan Drewniok on 19.06.23.
//

#ifndef FICTION_BESTAGON_GATE_GENERATOR_HPP
#define FICTION_BESTAGON_GATE_GENERATOR_HPP

#include "fiction/algorithms/simulation/sidb/critical_temperature.hpp"
#include "fiction/algorithms/simulation/sidb/random_layout_generator.hpp"
#include "fiction/io/read_sqd_layout.hpp"
#include "fiction/technology/sidb_surface.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/truth_table_utils.hpp"

#include <kitty/dynamic_truth_table.hpp>

#include <cstdint>
#include <iostream>
#include <mutex>
#include <random>
#include <string_view>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fiction
{
enum class bestagon_direction
{
    /**
     * This simulation engine computes Critical Temperature values with 100 % accuracy.
     */
    STRAIGHT,
    /**
     * This simulation engine quickly calculates the Critical Temperature. However, there may be deviations from the
     * exact Critical Temperature. This mode is recommended for larger layouts (> 40 SiDBs).
     */
    DIAGONAL
};
/**
 * This struct stores the parameters for the *generate_random_layout* algorithm.
 */
template <typename Lyt>
struct bestagon_gate_generator_params
{
    /**
     * Canvas size.
     */
    std::pair<uint64_t, uint64_t> canvas_size{std::pair(10, 6)};
    /**
     * Parameter to generate the random layout.
     */
    random_layout_params<Lyt> random_lyt_params{};
    /**
     * Truth table of the given gate (if layout is simulated in `gate-based` mode).
     */
    tt truth_table{};

    Lyt defect_surface{};
    /**
     * Number of generated layouts.
     */
    uint64_t number_of_layouts = 0;
    /**
     * Number of layouts to find a working gate.
     */
    uint64_t number_of_maximal_iterations = 10E6;

    bestagon_direction gate_direction = bestagon_direction::STRAIGHT;
};

template <typename Lyt>
std::vector<Lyt> bestagon_gate_generator(bestagon_gate_generator_params<Lyt>& params)
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
    static_assert(has_siqad_coord_v<Lyt>, "Lyt is not based on SiQAD coordinates");

    Lyt                layout{};
    typename Lyt::cell top_left_cell{};
    typename Lyt::cell bottom_right_cell{};

    std::vector<std::vector<uint64_t>> deactivating_cell_indices{};

    if (params.truth_table.num_vars() == 3)      // truth table of (double wire, cx, etc.) consists of three variables.
    {
        if (params.truth_table.num_bits() == 8)  // number of bits of truth table.
        {
            layout = read_sqd_layout<Lyt>(
                "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton/skeleton_hex_inputsdbp_2i2o.sqd");
            top_left_cell             = {19 - (params.canvas_size.first - params.canvas_size.first % 2) / 2, 7, 0};
            bottom_right_cell         = {19 + (params.canvas_size.first + params.canvas_size.first % 2) / 2,
                                         7 + params.canvas_size.second, 0};
            deactivating_cell_indices = {{2, 3}, {1, 2}, {0, 3}, {0, 1}};
        }
    }

    else if (params.truth_table.num_vars() == 1 && params.truth_table.num_bits() == 2)
    {
        if (params.gate_direction == bestagon_direction::STRAIGHT)
        {
            layout            = read_sqd_layout<Lyt>("/Users/jandrewniok/CLionProjects/fiction_fork/experiments/"
                                                     "skeleton/skeleton_hex_inputsdbp_1i1o_straight.sqd");
            top_left_cell     = {14, 6, 0};
            bottom_right_cell = {14 + params.canvas_size.first, 6 + params.canvas_size.second, 0};
            deactivating_cell_indices = {{1}, {0}};
        }
        else
        {
            layout            = read_sqd_layout<Lyt>("/Users/jandrewniok/CLionProjects/fiction_fork/experiments/"
                                                     "skeleton/skeleton_hex_inputsdbp_1i1o_diagonal.sqd");
            top_left_cell     = {19 - (params.canvas_size.first - params.canvas_size.first % 2) / 2, 7, 0};
            bottom_right_cell = {19 + (params.canvas_size.first + params.canvas_size.first % 2) / 2,
                                 6 + params.canvas_size.second, 0};
            deactivating_cell_indices = {{1}, {0}};
        }
    }

    else if (params.truth_table.num_vars() == 2)
    {
        if (params.truth_table.num_bits() == 4 && params.truth_table != create_fan_out_tt())  // and, or, nand, etc.
        {
            layout = read_sqd_layout<Lyt>(
                "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton/skeleton_hex_inputsdbp_2i1o.sqd");
            top_left_cell             = {19 - (params.canvas_size.first - params.canvas_size.first % 2) / 2, 6, 0};
            bottom_right_cell         = {19 + (params.canvas_size.first + params.canvas_size.first % 2) / 2,
                                         6 + params.canvas_size.second, 0};
            deactivating_cell_indices = {{2, 3}, {1, 2}, {0, 3}, {0, 1}};
        }
        else
        {
            layout = read_sqd_layout<Lyt>(
                "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton/skeleton_hex_inputsdbp_1i2o.sqd");
            top_left_cell             = {19 - (params.canvas_size.first - params.canvas_size.first % 2) / 2, 8, 0};
            bottom_right_cell         = {19 + (params.canvas_size.first + params.canvas_size.first % 2) / 2,
                                         8 + params.canvas_size.second, 0};
            deactivating_cell_indices = {{1}, {0}};
        }
    }


    params.random_lyt_params.coordinate_pair = std::make_pair(top_left_cell, bottom_right_cell);

    std::cout << "layouts produced" << std::endl;

    std::vector<Lyt> found_gates{};

    uint64_t const           num_threads = 10;
    std::vector<std::thread> threads{};
    threads.reserve(num_threads);
    std::mutex mutex{};  // used to control access to shared resources

    std::atomic<bool> found(false);

    params.defect_surface.foreach_sidb_defect(
        [&params, &layout](const auto& defect)
        { layout.assign_sidb_defect(defect.first, params.defect_surface.get_sidb_defect(defect.first)); });

    std::vector<typename Lyt::cell> cells{};
    cells.reserve(layout.num_cells());
    layout.foreach_cell([&cells](const auto& cell) { cells.push_back(cell); });
    std::sort(cells.begin(), cells.end());

    for (uint64_t z = 0u; z < num_threads; z++)
    {
        threads.emplace_back(
            [&, z, params]
            {
                critical_temperature_params temp_params{
                    simulation_engine::EXACT,
                    critical_temperature_mode::GATE_BASED_SIMULATION,
                    quicksim_params{sidb_simulation_parameters{2, -0.32}, 100, 0.65},
                    0.99,
                    350,
                    params.truth_table};

                while (!found)
                {

                    auto lyt = generate_random_layout(params.random_lyt_params, layout);

                    double temp = 1000;

                    for (auto i = 0u; i < deactivating_cell_indices.size(); i++)
                    {
                        critical_temperature_stats<Lyt> criticalstats{};
                        for (const auto& deactive_cell : deactivating_cell_indices[i])
                        {
                            lyt.assign_cell_type(cells[deactive_cell], Lyt::technology::EMPTY);
                        }

                        temp_params.input_bit = i;

                        critical_temperature(lyt, temp_params, &criticalstats);
                        if (criticalstats.critical_temperature < temp)
                        {
                            temp = criticalstats.critical_temperature;
                        }

                        for (const auto& deactive_cell : deactivating_cell_indices[i])
                        {
                            lyt.assign_cell_type(cells[deactive_cell], Lyt::technology::NORMAL);
                        }
                    }

                    //                    write_sqd_layout(
                    //                        lyt,
                    //                        "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton/layout_found.sqd");
                    //
                    if (temp > 0)
                    {
                        found = true;
                        std::cout << "layout found" << std::endl;
                        std::cout << temp << std::endl;
                        const std::lock_guard lock{mutex};
                        found_gates.push_back(lyt);
                        break;
                    }
                }
            });
    }
    for (auto& thread : threads)
    {
        thread.join();
    }

    return found_gates;
}

}  // namespace fiction

#endif  // FICTION_BESTAGON_GATE_GENERATOR_HPP
