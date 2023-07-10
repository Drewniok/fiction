//
// Created by Jan Drewniok on 19.06.23.
//

#ifndef FICTION_BESTAGON_GATE_GENERATOR_HPP
#define FICTION_BESTAGON_GATE_GENERATOR_HPP

#include "fiction/algorithms/simulation/sidb/add_cells_to_layout.hpp"
#include "fiction/algorithms/simulation/sidb/critical_temperature.hpp"
#include "fiction/algorithms/simulation/sidb/random_layout_generator.hpp"
#include "fiction/io/read_sqd_layout.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/truth_table_utils.hpp"

#include <kitty/dynamic_truth_table.hpp>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <numeric>
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

void makeCombiUtil(std::vector<std::vector<uint64_t>>& ans, std::vector<uint64_t>& tmp, uint64_t n, uint64_t left,
                   uint64_t k)
{
    // Pushing this vector to a vector of vector
    if (k == 0)
    {
        ans.push_back(tmp);
        return;
    }

    // i iterates from left to n. First time
    // left will be 1
    for (uint64_t i = left; i <= n; ++i)
    {
        tmp.push_back(i);
        makeCombiUtil(ans, tmp, n, i + 1, k - 1);

        // Popping out last inserted element
        // from the vector
        tmp.pop_back();
    }
}

std::vector<std::vector<uint64_t>> makeCombi(uint64_t n, uint64_t k)
{
    std::vector<std::vector<uint64_t>> ans;
    std::vector<uint64_t>              tmp;
    makeCombiUtil(ans, tmp, n, 0, k);
    return ans;
}

sidb_cell_clk_lyt_siqad::cell index_to_cell(uint64_t index, cube::coord_t& canvas_size, siqad::coord_t& nordwest)
{
    const auto y     = (index - index % (canvas_size.x)) / (canvas_size.x);
    const auto x     = index % (canvas_size.x);
    auto       siqad = siqad::to_siqad_coord(cube::coord_t{x, y});

    return siqad + nordwest;
}

struct bestagon_gate_generator_params
{
    /**
     * Canvas size.
     */
    cube::coord_t canvas_size{10, 16};
    /**
     * Truth table of the given gate (if layout is simulated in `gate-based` mode).
     */
    tt truth_table{};
    /**
     * Number of generated layouts.
     */
    uint64_t number_of_sidbs = 0;
    /**
     * Number of layouts to find a working gate.
     */

    bestagon_direction gate_direction = bestagon_direction::STRAIGHT;
};

std::vector<sidb_cell_clk_lyt_siqad> exhaustive_gate_generator(bestagon_gate_generator_params& params)
{
    std::vector<sidb_cell_clk_lyt_siqad>   found_gates{};
    sidb_cell_clk_lyt_siqad                layout{};
    typename sidb_cell_clk_lyt_siqad::cell top_left_cell{};
    typename sidb_cell_clk_lyt_siqad::cell bottom_right_cell{};

    const auto canvas_siqad = siqad::to_siqad_coord(params.canvas_size);

    std::vector<std::vector<uint64_t>> deactivating_cell_indices{};

    if (params.truth_table.num_vars() == 3)      // truth table of (double wire, cx, etc.) consists of three variables.
    {
        if (params.truth_table.num_bits() == 8)  // number of bits of truth table.
        {
            layout = read_sqd_layout<sidb_cell_clk_lyt_siqad>(
                "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton/skeleton_hex_inputsdbp_2i2o.sqd");

            top_left_cell     = {19 - (canvas_siqad.x - canvas_siqad.x % 2) / 2, 7, 0};
            bottom_right_cell = {19 + (canvas_siqad.x + canvas_siqad.x % 2) / 2, 7 + canvas_siqad.y, canvas_siqad.z};
            deactivating_cell_indices = {{2, 3}, {1, 2}, {0, 3}, {0, 1}};
        }
    }

    else if (params.truth_table.num_vars() == 1 && params.truth_table.num_bits() == 2)
    {
        if (params.gate_direction == bestagon_direction::STRAIGHT)
        {
            layout =
                read_sqd_layout<sidb_cell_clk_lyt_siqad>("/Users/jandrewniok/CLionProjects/fiction_fork/experiments/"
                                                         "skeleton/skeleton_hex_inputsdbp_1i1o_straight.sqd");
            top_left_cell     = {14 - (canvas_siqad.x - canvas_siqad.x % 2) / 2, 6, 0};
            bottom_right_cell = {14 + (canvas_siqad.x + canvas_siqad.x % 2) / 2, 6 + canvas_siqad.y, canvas_siqad.z};
            deactivating_cell_indices = {{1}, {0}};
        }
        else
        {
            layout =
                read_sqd_layout<sidb_cell_clk_lyt_siqad>("/Users/jandrewniok/CLionProjects/fiction_fork/experiments/"
                                                         "skeleton/skeleton_hex_inputsdbp_1i1o_diagonal.sqd");
            top_left_cell     = {19 - (canvas_siqad.x - canvas_siqad.x % 2) / 2, 7, 0};
            bottom_right_cell = {19 + (canvas_siqad.x + canvas_siqad.x % 2) / 2, 7 + canvas_siqad.y, canvas_siqad.z};
            deactivating_cell_indices = {{1}, {0}};
        }
    }

    else if (params.truth_table.num_vars() == 2)
    {
        if (params.truth_table.num_bits() == 4 && params.truth_table != create_fan_out_tt())  // and, or, nand, etc.
        {
            layout = read_sqd_layout<sidb_cell_clk_lyt_siqad>(
                "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton/skeleton_hex_inputsdbp_2i1o.sqd");
            top_left_cell     = {19 - (canvas_siqad.x - canvas_siqad.x % 2) / 2, 6, 0};
            bottom_right_cell = {19 + (canvas_siqad.x + canvas_siqad.x % 2) / 2, 6 + canvas_siqad.y, canvas_siqad.z};
            deactivating_cell_indices = {{2, 3}, {1, 2}, {0, 3}, {0, 1}};
        }
        else
        {
            layout = read_sqd_layout<sidb_cell_clk_lyt_siqad>(
                "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton/skeleton_hex_inputsdbp_1i2o.sqd");
            top_left_cell     = {19 - (canvas_siqad.x - canvas_siqad.x % 2) / 2, 8, 0};
            bottom_right_cell = {19 + (canvas_siqad.x + canvas_siqad.x % 2) / 2, 8 + canvas_siqad.y, canvas_siqad.z};
            deactivating_cell_indices = {{1}, {0}};
        }
    }

    std::vector<sidb_cell_clk_lyt_siqad::cell> cells{};
    cells.reserve(layout.num_cells());
    layout.foreach_cell([&cells](const auto& cell) { cells.push_back(cell); });
    std::sort(cells.begin(), cells.end());

    const uint64_t canvas_cells_number = static_cast<uint64_t>((params.canvas_size.x) * (params.canvas_size.y));

    std::atomic<bool> found(false);

    const auto ans = makeCombi(canvas_cells_number, params.number_of_sidbs);
    // std::cout << "done" << std::endl;
    // std::cout << ans.size() << std::endl;

    //    for (const auto &conf : ans)
    //    {
    //        for (const auto &idx : conf)
    //        {
    //            std::cout << " | " + std::to_string(idx);
    //        }
    //        std::cout << std::endl;
    //    }

    uint64_t loop_counter = 0;

    uint64_t const           num_threads = 10;
    std::vector<std::thread> threads{};
    threads.reserve(num_threads);
    std::mutex mutex{};  // used to control access to shared resources

    const auto number_per_thread = (ans.size() - (ans.size() % num_threads)) / num_threads;
    const auto number_last       = ans.size() % num_threads;

    for (uint64_t z = 0u; z < num_threads; ++z)
    {
        threads.emplace_back(
            [&, z, ans]
            {
                critical_temperature_params temp_params{
                    simulation_engine::EXACT,
                    critical_temperature_mode::GATE_BASED_SIMULATION,
                    quicksim_params{sidb_simulation_parameters{2, -0.32}, 100, 0.65},
                    0.99,
                    350,
                    params.truth_table};

                for (auto i = z * number_per_thread; i < (z + 1) * number_per_thread; i++)
                {
                    std::vector<sidb_cell_clk_lyt_siqad::cell> placable_cells{};
                    const auto                                 config = ans[i];
                    placable_cells.reserve(config.size());
                    loop_counter += 1;
                    // std::cout << (config).size() << std::endl;
                    for (const auto& idx : config)
                    {
                        // std::cout << " | " << std::to_string(idx);
                        placable_cells.push_back(index_to_cell(idx, params.canvas_size, top_left_cell));
                    }

                    // std::cout << std::endl;

                    bool constraint_violation_positive_sidbs = false;

                    for (const auto& plcell1 : placable_cells)
                    {
                        for (const auto& plcell2 : placable_cells)
                        {
                            if (plcell1 != plcell2)
                            {
                                if (euclidean_distance<sidb_cell_clk_lyt_siqad>(layout, plcell1, plcell2) < 3)
                                {
                                    constraint_violation_positive_sidbs = true;
                                }
                            }
                        }
                    }

                    //                    if (loop_counter % 10000 == 0)
                    //                    {
                    //                        std::cout << loop_counter << std::endl;
                    //                    }

                    auto layout_with_placed = add_cells_to_layout(layout, placable_cells);

                    if (constraint_violation_positive_sidbs)
                    {
                        continue;
                    }

                    double temp = 1000;

                    for (auto i = 0u; i < deactivating_cell_indices.size(); i++)
                    {
                        critical_temperature_stats<sidb_cell_clk_lyt_siqad> criticalstats{};
                        for (const auto& deactive_cell : deactivating_cell_indices[i])
                        {
                            layout_with_placed.assign_cell_type(cells[deactive_cell],
                                                                sidb_cell_clk_lyt_siqad::technology::EMPTY);
                        }

                        temp_params.input_bit = i;

                        critical_temperature(layout_with_placed, temp_params, &criticalstats);
                        if (criticalstats.critical_temperature < temp)
                        {
                            temp = criticalstats.critical_temperature;
                        }

                        for (const auto& deactive_cell : deactivating_cell_indices[i])
                        {
                            layout_with_placed.assign_cell_type(cells[deactive_cell],
                                                                sidb_cell_clk_lyt_siqad::technology::NORMAL);
                        }
                    }
                    //                                        write_sqd_layout(
                    //                                            layout_with_placed,
                    //                                            "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton/layout_found.sqd");
                    if (temp > 0)
                    {
                        //                        for (const auto& idx : config)
                        //                        {
                        //                            std::cout << idx << std::endl;
                        //                        }
                        //                        write_sqd_layout(
                        //                            layout_with_placed,
                        //                            "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton/nand/nand_2/layout_found_temp_"
                        //                            +
                        //                                std::to_string(loop_counter) + " .sqd");
                        found = true;
                        // std::cout << "layout found" << std::endl;
                        // std::cout << temp << std::endl;
                        const std::lock_guard lock{mutex};
                        found_gates.push_back(layout_with_placed);
                        // break;
                    }
                }
            });
    }

    for (auto& thread : threads)
    {
        thread.join();
    }

    //    for (const auto& config : ans)
    //    {
    //        loop_counter += 1;
    //        std::vector<sidb_cell_clk_lyt_siqad::cell> placable_cells{};
    //        placable_cells.reserve(config.size());
    //        for (const auto& idx : config)
    //        {
    //            placable_cells.push_back(index_to_cell(idx, params.canvas_size, top_left_cell));
    //        }
    //
    //        bool constraint_violation_positive_sidbs = false;
    //
    //        for (const auto& plcell1 : placable_cells)
    //        {
    //            for (const auto& plcell2 : placable_cells)
    //            {
    //                if (plcell1 != plcell2)
    //                {
    //                    if (euclidean_distance<sidb_cell_clk_lyt_siqad>(layout, plcell1, plcell2) < 3)
    //                    {
    //                        constraint_violation_positive_sidbs = true;
    //                    }
    //                }
    //            }
    //        }
    //
    //        if (loop_counter % 10000 == 0)
    //        {
    //            std::cout << loop_counter << std::endl;
    //        }
    //
    //        auto layout_with_placed = add_cells_to_layout(layout, placable_cells);
    //
    //        if (constraint_violation_positive_sidbs)
    //        {
    //            continue;
    //        }
    //
    //        double temp = 1000;
    //
    //        for (auto i = 0u; i < deactivating_cell_indices.size(); i++)
    //        {
    //            critical_temperature_stats<sidb_cell_clk_lyt_siqad> criticalstats{};
    //            for (const auto& deactive_cell : deactivating_cell_indices[i])
    //            {
    //                layout_with_placed.assign_cell_type(cells[deactive_cell],
    //                sidb_cell_clk_lyt_siqad::technology::EMPTY);
    //            }
    //
    //            temp_params.input_bit = i;
    //
    //            critical_temperature(layout_with_placed, temp_params, &criticalstats);
    //            if (criticalstats.critical_temperature < temp)
    //            {
    //                temp = criticalstats.critical_temperature;
    //            }
    //
    //            for (const auto& deactive_cell : deactivating_cell_indices[i])
    //            {
    //                layout_with_placed.assign_cell_type(cells[deactive_cell],
    //                sidb_cell_clk_lyt_siqad::technology::NORMAL);
    //            }
    //        }
    //        write_sqd_layout(layout_with_placed,
    //                         "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton/layout_found.sqd");
    //        if (temp > 0)
    //        {
    //            write_sqd_layout(
    //                layout_with_placed,
    //                "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton/layout_found_temp.sqd");
    //            found = true;
    //            std::cout << "layout found" << std::endl;
    //            std::cout << temp << std::endl;
    //            // const std::lock_guard lock{mutex};
    //            found_gates.push_back(layout_with_placed);
    //            break;
    //        }
    //    }

    return found_gates;
}

//    params.random_lyt_params.coordinate_pair = std::make_pair(top_left_cell, bottom_right_cell);
//
//    std::cout << "layouts produced" << std::endl;
//
//    std::vector<sidb_cell_clk_lyt_siqad> found_gates{};
//
//    uint64_t counter = 0;
//
//    uint64_t const           num_threads = 10;
//    std::vector<std::thread> threads{};
//    threads.reserve(num_threads);
//    std::mutex mutex{};  // used to control access to shared resources
//
//    std::atomic<bool> found(false);
//
//    for (uint64_t z = 0u; z < num_threads; z++)
//    {
//        threads.emplace_back(
//            [&]
//            {
//                critical_temperature_params temp_params{
//                    simulation_engine::EXACT,
//                    critical_temperature_mode::GATE_BASED_SIMULATION,
//                    quicksim_params{sidb_simulation_parameters{2, -0.32}, 100, 0.65},
//                    0.99,
//                    350,
//                    params.truth_table};
//
//                while (!found)
//                {
//
//                    auto lyt = generate_random_layout(params.random_lyt_params, layout);
//
//                    double temp = 1000;
//
//                    for (auto i = 0u; i < deactivating_cell_indices.size(); i++)
//                    {
//                        critical_temperature_stats<sidb_cell_clk_lyt_siqad> criticalstats{};
//                        for (const auto& deactive_cell : deactivating_cell_indices[i])
//                        {
//                            lyt.assign_cell_type(cells[deactive_cell],
//                            sidb_cell_clk_lyt_siqad::technology::EMPTY);
//                        }
//
//                        temp_params.input_bit = i;
//
//                        critical_temperature(lyt, temp_params, &criticalstats);
//                        if (criticalstats.critical_temperature < temp)
//                        {
//                            temp = criticalstats.critical_temperature;
//                        }
//
//                        for (const auto& deactive_cell : deactivating_cell_indices[i])
//                        {
//                            lyt.assign_cell_type(cells[deactive_cell],
//                            sidb_cell_clk_lyt_siqad::technology::NORMAL);
//                        }
//                    }
//
//                    //                    write_sqd_layout(
//                    //                        lyt,
//                    // "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton/layout_found.sqd");
//                    //
//                    if (temp > 0)
//                    {
//                        found = true;
//                        std::cout << "layout found" << std::endl;
//                        std::cout << temp << std::endl;
//                        const std::lock_guard lock{mutex};
//                        found_gates.push_back(lyt);
//                        break;
//                    }
//                }
//            });
//    }
//    for (auto& thread : threads)
//    {
//        thread.join();
//    }
//
//    return found_gates;

}  // namespace fiction

#endif  // FICTION_BESTAGON_GATE_GENERATOR_HPP
