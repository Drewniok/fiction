//
// Created by Jan Drewniok 10.06.24
//

#include "fiction/algorithms/physical_design/design_sidb_gates.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/io/read_sqd_layout.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
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
    experiments::experiment<std::string, double, uint64_t, uint64_t, double, double> simulation_exp{
        "benchmark",
        "gate",
        "runtime (Automatic Exhaustive)",
        "#Gates (Automatic Exhaustive)",
        "#Gates (QuickCell)",
        "runtime (QuickCell)",
        "percentual time reduction"};

    const auto truth_tables = std::vector<std::vector<tt>>{
        std::vector<tt>{create_and_tt()}, std::vector<tt>{create_nand_tt()}, std::vector<tt>{create_or_tt()},
        std::vector<tt>{create_nor_tt()}, std::vector<tt>{create_xor_tt()},  std::vector<tt>{create_xnor_tt()},
        std::vector<tt>{create_lt_tt()},  std::vector<tt>{create_gt_tt()},   std::vector<tt>{create_le_tt()},
        std::vector<tt>{create_ge_tt()},  create_crossing_wire_tt(),         create_half_adder_tt(),
        create_double_wire_tt()};

    static const std::array<std::string, 13> gate_names = {"and", "nand", "or", "nor", "xor", "xnor",     "lt",
                                                           "gt",  "le",   "ge", "cx",  "ha",  "hourglass"};

    static const std::string folder = fmt::format("{}skeleton_bestagons_with_tags/", EXPERIMENTS_PATH);

    auto skeleton_one_input_two_output =
        read_sqd_layout<sidb_100_cell_clk_lyt_siqad>(fmt::format("{}{}", folder, "/skeleton_hex_inputsdbp_2i1o.sqd"));

    auto skeleton_two_input_two_output =
        read_sqd_layout<sidb_100_cell_clk_lyt_siqad>(fmt::format("{}{}", folder, "/skeleton_hex_inputsdbp_2i2o.sqd"));

    design_sidb_gates_params<fiction::cell<sidb_100_cell_clk_lyt_siqad>> params_2_in_1_out{
        is_operational_params{sidb_simulation_parameters{2, -0.32}},
        design_sidb_gates_params<
            fiction::cell<sidb_100_cell_clk_lyt_siqad>>::design_sidb_gates_mode::AUTOMATIC_EXHAUSTIVE_GATE_DESIGNER,
        {{14, 6, 0}, {24, 10, 0}},
        3};

    design_sidb_gates_params<fiction::cell<sidb_100_cell_clk_lyt_siqad>> params_2_in_2_out{
        is_operational_params{sidb_simulation_parameters{2, -0.32}},
        design_sidb_gates_params<
            fiction::cell<sidb_100_cell_clk_lyt_siqad>>::design_sidb_gates_mode::AUTOMATIC_EXHAUSTIVE_GATE_DESIGNER,
        {{14, 6, 0}, {24, 14, 0}},
        3};

    double sum_exhaustive_runtime = 0;
    double sum_quickcell_runtime  = 0;

    for (auto i = 0u; i < truth_tables.size(); i++)
    {
        const auto& table = truth_tables[i];

        design_sidb_gates_stats stats_automatic_exhaustive_design{};

        std::vector<sidb_100_cell_clk_lyt_siqad> automatic_exhaustive_design{};

        if (gate_names[i] == "cx" || gate_names[i] == "ha" || gate_names[i] == "hourglass")
        {
            automatic_exhaustive_design = design_sidb_gates(skeleton_two_input_two_output, table, params_2_in_2_out,
                                                            &stats_automatic_exhaustive_design);
        }
        else
        {
            automatic_exhaustive_design = design_sidb_gates(skeleton_one_input_two_output, table, params_2_in_1_out,
                                                            &stats_automatic_exhaustive_design);
        }

        std::cout << mockturtle::to_seconds(stats_automatic_exhaustive_design.time_total) << '\n';

        std::vector<sidb_100_cell_clk_lyt_siqad> quickcell_design{};
        design_sidb_gates_stats                  stats_quickcell{};

        params_2_in_1_out.design_mode =
            design_sidb_gates_params<fiction::cell<sidb_100_cell_clk_lyt_siqad>>::design_sidb_gates_mode::QUICKCELL;
        params_2_in_2_out.design_mode =
            design_sidb_gates_params<fiction::cell<sidb_100_cell_clk_lyt_siqad>>::design_sidb_gates_mode::QUICKCELL;

        if (gate_names[i] == "cx" || gate_names[i] == "ha" || gate_names[i] == "hourglass")
        {
            quickcell_design =
                design_sidb_gates(skeleton_two_input_two_output, table, params_2_in_2_out, &stats_quickcell);
        }
        else
        {
            quickcell_design =
                design_sidb_gates(skeleton_one_input_two_output, table, params_2_in_1_out, &stats_quickcell);
        }

        const auto runtime_automatic_exhaustive_design =
            mockturtle::to_seconds(stats_automatic_exhaustive_design.time_total);
        const auto runtime_quickcell = mockturtle::to_seconds(stats_quickcell.time_total);

        sum_exhaustive_runtime += runtime_automatic_exhaustive_design;
        sum_quickcell_runtime += runtime_quickcell;

        const auto time_reduction = runtime_automatic_exhaustive_design / runtime_quickcell;

        const auto final_number_of_gates = quickcell_design.size();

        simulation_exp(gate_names[i], runtime_automatic_exhaustive_design, automatic_exhaustive_design.size(),
                       final_number_of_gates, runtime_quickcell, time_reduction);

        simulation_exp.save();
        simulation_exp.table();
    }

    const auto total_time_reduction = sum_exhaustive_runtime / sum_quickcell_runtime;

    simulation_exp("", sum_exhaustive_runtime, 0, 0, sum_quickcell_runtime, total_time_reduction);

    simulation_exp.save();
    simulation_exp.table();

    return EXIT_SUCCESS;
}