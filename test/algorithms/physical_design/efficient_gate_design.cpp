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
#include <fiction/io/read_sqd_layout.hpp>
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

// TEST_CASE("Efficient gate design, AND", "[design-sidb-gates]")
//{
//     auto lyt =
//     read_sqd_layout<sidb_100_cell_clk_lyt_siqad>("/Users/jandrewniok/CLionProjects/fiction_fork/experiments/"
//                                                             "skeleton_reversible/two_input_one_output_long_wire.sqd");
//     //    auto lyt = blueprints::two_input_two_output_skeleton<sidb_100_cell_clk_lyt_siqad>();
//     // const auto chains = detect_bdl_wires(lyt);
//
//     //    const design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad> params{
//     //        sidb_simulation_parameters{2, -0.32},
//     //        design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::EXHAUSTIVE,
//     //        {{14, 6, 0}, {24, 13, 0}},
//     //        3,
//     //        sidb_simulation_engine::QUICKEXACT};
//
//     const design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad> params{
//         sidb_simulation_parameters{2, -0.32},
//         design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::EXHAUSTIVE,
//         {{41, 14, 0}, {52, 18, 0}},
//         3,
//         sidb_simulation_engine::QUICKEXACT};
//
////    const design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad> params{
////        sidb_simulation_parameters{2, -0.32},
////        design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::EXHAUSTIVE,
////        {{15, 8, 0}, {31, 13, 0}},
////        3,
////        sidb_simulation_engine::QUICKEXACT};
//
//    efficient_gate_design_params<sidb_100_cell_clk_lyt_siqad> efficient_params{};
//    efficient_params.design_params = params;
////
////        std::vector<tt> truth_tables = {create_and3_tt(),    create_xor_and_tt(), create_or_and_tt(),
/// create_onehot_tt(), /                                        create_maj_tt(),     create_gamble_tt(),
/// create_dot_tt(),    create_ite_tt(), /                                        create_and_xor_tt(),
/// create_xor3_tt()};
//
//    //
//    //    design_sidb_gates_stats st{};
//    //    const auto              exhaustive_design = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params,
//    //    &st); std::cout << exhaustive_design.size() << std::endl; std::cout << mockturtle::to_seconds(st.time_total)
//    //    << std::endl; write_sqd_layout(exhaustive_design[0],
//    "/Users/jandrewniok/Desktop/and_efficient_gate_exh.sqd");
//
//    efficient_gate_design_stats stats{};
//    const auto all_gate_candidates = design_all_efficient_gates(lyt, std::vector<tt>{create_and_tt()},
//    efficient_params, &stats); CHECK(all_gate_candidates.size() == 1173); std::cout <<
//    mockturtle::to_seconds(stats.time_total) << std::endl;
//    // std::cout << all_gate_candidates.size() << std::endl;
//    // write_sqd_layout(all_gate_candidates[0], "/Users/jandrewniok/Desktop/and_efficient_gate.sqd");
//
//    uint64_t counter = 0;
//        for (const auto& gate : all_gate_candidates)
//        {
//            counter++;
//            if (counter <5)
//            {
//                write_sqd_layout(gate, fmt::format("/Users/jandrewniok/Desktop/large_efficient_gate_{}.sqd",
//                counter));
//            }
//        }
//
//    std::cout << counter << std::endl;
//}

TEST_CASE("Efficient gate design, 3-input", "[design-sidb-gates]")
{
    auto lyt = read_sqd_layout<sidb_100_cell_clk_lyt_siqad>(
        "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/"
        "skeleton_reversible/3_in_1_out_larger_canvas_larger_distance_closer_input.sqd");

    const design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad> params{
        sidb_simulation_parameters{2, -0.30},
        design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::EXHAUSTIVE,
        {{23, 7, 0}, {34, 12, 0}},
        4,
        sidb_simulation_engine::QUICKEXACT};

    //    const design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad> params{
    //        sidb_simulation_parameters{2, -0.32},
    //        design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::EXHAUSTIVE,
    //        {{15, 8, 0}, {31, 13, 0}},
    //        3,
    //        sidb_simulation_engine::QUICKEXACT};

    efficient_gate_design_params<sidb_100_cell_clk_lyt_siqad> efficient_params{};
    efficient_params.design_params = params;
    //
    std::vector<tt> truth_tables = {create_and3_tt(),    create_xor_and_tt(), create_or_and_tt(), create_onehot_tt(),
                                    create_maj_tt(),     create_gamble_tt(),  create_dot_tt(),    create_ite_tt(),
                                    create_and_xor_tt(), create_xor3_tt()};

    std::vector<std::string> names = {"and3",   "xor_and", "or_and", "onehot", "maj",
                                      "gamble", "dot",     "ite",    "and_xor", "xor3"};

    //
    //    design_sidb_gates_stats st{};
    //    const auto              exhaustive_design = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params,
    //    &st); std::cout << exhaustive_design.size() << std::endl; std::cout << mockturtle::to_seconds(st.time_total)
    //    << std::endl; write_sqd_layout(exhaustive_design[0], "/Users/jandrewniok/Desktop/and_efficient_gate_exh.sqd");

    uint64_t truth_counter = 0;
    for (const auto& table : truth_tables)
    {
        efficient_params.design_params.simulation_parameters.mu_minus = -0.30;
        for (auto i = 0u; i < 3; i++)
        {
            efficient_params.design_params.simulation_parameters.mu_minus -= 0.01;
            std::cout << efficient_params.design_params.simulation_parameters.mu_minus << std::endl;
            efficient_gate_design_stats stats{};
            const auto                  all_gate_candidates =
                design_all_efficient_gates(lyt, std::vector<tt>{table}, efficient_params, &stats);
            CHECK(all_gate_candidates.size() == 1173);
            std::cout << mockturtle::to_seconds(stats.time_total) << std::endl;
            // std::cout << all_gate_candidates.size() << std::endl;
            // write_sqd_layout(all_gate_candidates[0], "/Users/jandrewniok/Desktop/and_efficient_gate.sqd");

            uint64_t counter = 0;
            for (const auto& gate : all_gate_candidates)
            {
                counter++;
                if (counter < 5)
                {
                    write_sqd_layout(gate, fmt::format("/Users/jandrewniok/Desktop/{}_mu_{}.sqd", names[truth_counter],
                                                       efficient_params.design_params.simulation_parameters.mu_minus,
                                                       counter));
                }
            }
        }
        truth_counter++;
    }

    // std::cout << counter << std::endl;
}