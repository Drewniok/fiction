//
// Created by Jan Drewniok 10.06.24
//

#include "fiction_experiments.hpp"

#include <fiction/algorithms/iter/bdl_input_iterator.hpp>
#include <fiction/algorithms/physical_design/design_sidb_gates.hpp>
#include <fiction/algorithms/simulation/sidb/is_operational.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_engine.hpp>
#include <fiction/io/read_sqd_layout.hpp>
#include <fiction/traits.hpp>
#include <fiction/types.hpp>
#include <fiction/utils/truth_table_utils.hpp>

#include <fmt/format.h>
#include <mockturtle/utils/stopwatch.hpp>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

// This script uses the *Automatic Exhaustive Gate Designer* and *QuickCell* to design gate implementations for 2-input
// Boolean functions. It records the number of designed gate implementations and the runtime required for each
// algorithm. The final column displays the runtime reduction factor achieved by *QuickCell*.

using namespace fiction;

int main()  // NOLINT
{

    const auto truth_tables_and_names =
        std::array<std::pair<std::vector<tt>, std::string>, 1>{{{std::vector<tt>{create_and_tt()}, "and"}}};

    static const std::string folder = fmt::format("{}/gate_skeletons/skeleton_bestagons_with_tags", EXPERIMENTS_PATH);

    const auto skeleton_one_input_two_output =
        read_sqd_layout<sidb_100_cell_clk_lyt_siqad>(fmt::format("{}/{}", folder, "skeleton_hex_inputsdbp_2i1o.sqd"));

    design_sidb_gates_params<fiction::cell<sidb_100_cell_clk_lyt_siqad>> params_2_in_1_out{
        is_operational_params{sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT,
                              bdl_input_iterator_params{}, is_operational_params::operational_condition::REJECT_KINKS},
        design_sidb_gates_params<
            fiction::cell<sidb_100_cell_clk_lyt_siqad>>::design_sidb_gates_mode::AUTOMATIC_EXHAUSTIVE_GATE_DESIGNER,
        {{14, 6, 0}, {24, 10, 0}},
        2};

    for (const auto& [truth_table, gate_name] : truth_tables_and_names)
    {
        design_sidb_gates_stats stats_automatic_exhaustive_design{};

        std::vector<sidb_100_cell_clk_lyt_siqad> automatic_exhaustive_design{};

        params_2_in_1_out.design_mode = design_sidb_gates_params<
            fiction::cell<sidb_100_cell_clk_lyt_siqad>>::design_sidb_gates_mode::AUTOMATIC_EXHAUSTIVE_GATE_DESIGNER;
        params_2_in_1_out.operational_params.op_condition = is_operational_params::operational_condition::REJECT_KINKS;

        automatic_exhaustive_design = design_sidb_gates(skeleton_one_input_two_output, truth_table, params_2_in_1_out,
                                                        &stats_automatic_exhaustive_design);
    }

    return EXIT_SUCCESS;
}
