//
// Created by Jan Drewniok 01.01.23
//

#include "experiments.hpp"
#include "fiction/algorithms/physical_design/design_sidb_gates.hpp"
#include "fiction/algorithms/physical_design/efficient_gate_design.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/io/read_sqd_layout.hpp"  // reader for SiDB layouts including surface scan data
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/utils/truth_table_utils.hpp"
#include "fiction_experiments.hpp"

#include <fmt/format.h>  // output formatting

#include <array>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

using namespace fiction;

int main()  // NOLINT
{

    experiments::experiment<std::string, uint64_t, double, uint64_t, uint64_t, double, uint64_t, double, uint64_t,
                            double, double, uint64_t, double, double>
        simulation_exp{"benchmark",
                       "gate",
                       " all layouts",
                       "exhaustive runtime",
                       " exhaustive gate number",

                       "L p1",
                       "L p1 / N",

                       "L p2",
                       "L p2 / N",

                       "L p3",
                       "L p3 / N",

                       "runtime for pruning",
                       "number of final gates",
                       "runtime total",
                       "percentual time reduction"};

    const auto truth_tables = std::vector<std::vector<tt>>{
        std::vector<tt>{create_and_tt()}, std::vector<tt>{create_nand_tt()}, std::vector<tt>{create_or_tt()},
        std::vector<tt>{create_nor_tt()}, std::vector<tt>{create_xor_tt()},  std::vector<tt>{create_xnor_tt()},
        std::vector<tt>{create_lt_tt()},  std::vector<tt>{create_gt_tt()},   std::vector<tt>{create_le_tt()},
        std::vector<tt>{create_ge_tt()},  create_crossing_wire_tt(),         create_half_adder_tt(),
        create_double_wire_tt()};

    //    const auto truth_tables = std::vector<std::vector<tt>>{std::vector<tt>{create_nand_tt()}};

    static const std::array<std::string, 13> gate_names = {"and", "nand", "or", "nor", "xor", "xnor",     "lt",
                                                           "gt",  "le",   "ge", "cx",  "ha",  "hourglass"};

    auto lyt = read_sqd_layout<sidb_100_cell_clk_lyt_siqad>(
        "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton_bestagons_with_tags/"
        "skeleton_hex_inputsdbp_2i1o.sqd");

    auto lyt_2_2 = read_sqd_layout<sidb_100_cell_clk_lyt_siqad>(
        "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton_bestagons_with_tags/"
        "skeleton_hex_inputsdbp_2i2o.sqd");

    const design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad> params_2_in_1_out{
        sidb_simulation_parameters{2, -0.32},
        design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::EXHAUSTIVE,
        {{14, 6, 0}, {24, 10, 0}},
        3,
        sidb_simulation_engine::QUICKEXACT};

    const design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad> params_2_in_2_out{
        sidb_simulation_parameters{2, -0.32},
        design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::EXHAUSTIVE,
        {{14, 6, 0}, {24, 14, 0}},
        3,
        sidb_simulation_engine::QUICKEXACT};

    const auto efficent_design_params = efficient_gate_design_params<sidb_100_cell_clk_lyt_siqad>{
        detect_bdl_pairs_params{}, params_2_in_1_out, DESIGN_MODE::PRUNING_AND_OPERATIONAL_CHECK};

    const auto efficent_design_params_2_2 = efficient_gate_design_params<sidb_100_cell_clk_lyt_siqad>{
        detect_bdl_pairs_params{}, params_2_in_2_out, DESIGN_MODE::PRUNING_AND_OPERATIONAL_CHECK};

    double sum_exhaustive_runtime = 0;
    double sum_efficient_runtime  = 0;

    for (auto i = 0u; i < truth_tables.size(); i++)
    {
        const auto& table = truth_tables[i];

        design_sidb_gates_stats stats_exhaustive{};

        std::vector<sidb_100_cell_clk_lyt_siqad> exhaustive_design{};

        if (gate_names[i] == "cx" || gate_names[i] == "ha" || gate_names[i] == "hourglass")
        {
            exhaustive_design = design_sidb_gates(lyt_2_2, table, params_2_in_2_out, &stats_exhaustive);
        }
        else
        {
            exhaustive_design = design_sidb_gates(lyt, table, params_2_in_1_out, &stats_exhaustive);
        }

        std::cout << mockturtle::to_seconds(stats_exhaustive.time_total) << std::endl;

        std::vector<sidb_100_cell_clk_lyt_siqad> efficient_design{};
        efficient_gate_design_stats              efficent_stats{};

        if (gate_names[i] == "cx" || gate_names[i] == "ha" || gate_names[i] == "hourglass")
        {
            // continue;
            efficient_design = efficient_gate_design(lyt_2_2, table, efficent_design_params_2_2, &efficent_stats);
        }
        else
        {
            efficient_design = efficient_gate_design(lyt, table, efficent_design_params, &efficent_stats);
        }

        sum_exhaustive_runtime += mockturtle::to_seconds(stats_exhaustive.time_total);
        sum_efficient_runtime += mockturtle::to_seconds(efficent_stats.time_total);

        const auto all_layouts     = efficent_stats.all_possible_layouts;
        const auto time_exhaustive = mockturtle::to_seconds(stats_exhaustive.time_total);

        const auto layouts_after_p1 = all_layouts - efficent_stats.lp1;
        const auto delta_p1         = static_cast<double>(layouts_after_p1) / all_layouts * 100;

        const auto layouts_after_p2 = layouts_after_p1 - efficent_stats.lp2;
        const auto delta_p2         = static_cast<double>(layouts_after_p2) / all_layouts * 100;

        const auto layouts_after_p3 = layouts_after_p2 - efficent_stats.lp3;
        const auto delta_p3         = static_cast<double>(layouts_after_p3) / all_layouts * 100;

        const auto time_efficient = mockturtle::to_seconds(efficent_stats.time_total);
        const auto time_pruning   = mockturtle::to_seconds(efficent_stats.time_total_pruning);

        const auto time_reduction = time_exhaustive / time_efficient;

        const auto final_number_of_gates = efficient_design.size();

        const auto time_efficient_per_pruning = layouts_after_p3 / all_layouts * 1000;

        simulation_exp(gate_names[i], all_layouts, time_exhaustive, exhaustive_design.size(), layouts_after_p1,
                       delta_p1, layouts_after_p2, delta_p2, layouts_after_p3, delta_p3, time_pruning,
                       final_number_of_gates, time_efficient, time_reduction);

        simulation_exp.save();
        simulation_exp.table();
    }

    const auto total_time_reduction = sum_exhaustive_runtime / sum_efficient_runtime;

    simulation_exp("", 0, sum_exhaustive_runtime, 0, 0, 0, 0, 0, 0, 0, 0, 0, sum_efficient_runtime,
                   total_time_reduction);

    simulation_exp.save();
    simulation_exp.table();
    return EXIT_SUCCESS;
}
