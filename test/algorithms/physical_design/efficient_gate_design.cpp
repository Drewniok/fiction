//
// Created by Jan Drewniok on 08.05.24.
//

#include <catch2/catch_template_test_macros.hpp>

#include "utils/blueprints/layout_blueprints.hpp"

#include <fiction/algorithms/physical_design/design_sidb_gates.hpp>
#include <fiction/algorithms/physical_design/efficient_gate_design.hpp>
#include <fiction/algorithms/simulation/sidb/is_operational.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_engine.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp>
#include <fiction/io/print_layout.hpp>
#include <fiction/io/write_sqd_layout.hpp>
#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/layouts/coordinates.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/technology/sidb_defect_surface.hpp>
#include <fiction/technology/sidb_defects.hpp>
#include <fiction/technology/sidb_lattice.hpp>
#include <fiction/technology/sidb_lattice_orientations.hpp>
#include <fiction/traits.hpp>
#include <fiction/types.hpp>
#include <fiction/utils/layout_utils.hpp>
#include <fiction/utils/truth_table_utils.hpp>

#include <vector>

using namespace fiction;

TEST_CASE("Efficient gate design, AND", "[design-sidb-gates]")
{
    auto lyt = blueprints::two_input_one_output_skeleton_new<sidb_100_cell_clk_lyt_siqad>();
    // const auto chains = detect_wire_bdl_chains(lyt);

    const design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad> params{
        sidb_simulation_parameters{2, -0.32},
        design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::EXHAUSTIVE,
        {{14, 6, 0}, {24, 13, 0}},
        3,
        sidb_simulation_engine::QUICKEXACT};

    efficient_gate_design_params<sidb_100_cell_clk_lyt_siqad> efficient_params{};
    efficient_params.design_params = params;
    //
//    design_sidb_gates_stats st{};
//    const auto              exhaustive_design = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params, &st);
//    std::cout << exhaustive_design.size() << std::endl;
//    std::cout << mockturtle::to_seconds(st.time_total) << std::endl;
//    write_sqd_layout(exhaustive_design[0], "/Users/jandrewniok/Desktop/and_efficient_gate_exh.sqd");

    const auto all_gate_candidates = design_all_efficient_gates(lyt, std::vector<tt>{create_and_tt()}, efficient_params);
    std::cout << all_gate_candidates.size() << std::endl;
    // write_sqd_layout(all_gate_candidates[0], "/Users/jandrewniok/Desktop/and_efficient_gate.sqd");

    uint64_t counter = 0;
    for (const auto& gate : all_gate_candidates)
    {
        if (is_operational(gate, std::vector<tt>{create_and_tt()}, is_operational_params{params.simulation_parameters}).first ==
            operational_status::OPERATIONAL)
        {
            counter++;
            write_sqd_layout(gate, fmt::format("/Users/jandrewniok/Desktop/and_efficient_gate_{}.sqd", counter));
        }
    }
    std::cout << counter << std::endl;
}