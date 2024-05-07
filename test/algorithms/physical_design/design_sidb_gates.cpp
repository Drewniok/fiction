//
// Created by Jan Drewniok on 12.09.23.
//

#include <catch2/catch_template_test_macros.hpp>

#include "utils/blueprints/layout_blueprints.hpp"

#include <fiction/algorithms/physical_design/design_sidb_gates.hpp>
#include <fiction/algorithms/simulation/sidb/is_operational.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_engine.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp>
#include <fiction/io/print_layout.hpp>
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
#include <fiction/io/write_sqd_layout.hpp>

#include <vector>

using namespace fiction;

// TEST_CASE("Use SiQAD XNOR skeleton and generate SiQAD XNOR gate, exhaustive", "[design-sidb-gates]")
//{
//     using offset_layout = sidb_100_cell_clk_lyt;
//     using siqad_layout  = sidb_100_cell_clk_lyt_siqad;
//     using cube_layout   = sidb_100_cell_clk_lyt_cube;
//
//     siqad_layout lyt{};
//
//     lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::INPUT);
//     lyt.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::INPUT);
//
//     lyt.assign_cell_type({20, 0, 0}, sidb_technology::cell_type::INPUT);
//     lyt.assign_cell_type({18, 1, 0}, sidb_technology::cell_type::INPUT);
//
//     lyt.assign_cell_type({6, 3, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({14, 3, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({4, 2, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({16, 2, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({10, 6, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({10, 7, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({10, 9, 0}, sidb_technology::cell_type::OUTPUT);
//     lyt.assign_cell_type({10, 10, 0}, sidb_technology::cell_type::OUTPUT);
//
//     lyt.assign_cell_type({10, 12, 1}, sidb_technology::cell_type::NORMAL);
//
//     CHECK(lyt.num_cells() == 13);
//
//     const design_sidb_gates_params<siqad_layout> params{
//         sidb_simulation_parameters{2, -0.32},
//         design_sidb_gates_params<siqad_layout>::design_sidb_gates_mode::EXHAUSTIVE,
//         {{10, 4, 0}, {10, 4, 0}},
//         1,
//         sidb_simulation_engine::QUICKEXACT};
//
//     const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_xnor_tt()}, params);
//
//     REQUIRE(found_gate_layouts.size() == 1);
//     CHECK(found_gate_layouts[0].num_cells() == 14);
//     CHECK(found_gate_layouts[0].get_cell_type({10, 4, 0}) == siqad_layout::technology::NORMAL);
//
//     // using cube coordinates
//     const auto                                  lyt_in_cube_coord = convert_to_fiction_coordinates<cube_layout>(lyt);
//     const design_sidb_gates_params<cube_layout> params_cube{
//         sidb_simulation_parameters{2, -0.32},
//         design_sidb_gates_params<cube_layout>::design_sidb_gates_mode::EXHAUSTIVE,
//         {siqad::to_fiction_coord<cube::coord_t>(siqad::coord_t{10, 4, 0}),
//          siqad::to_fiction_coord<cube::coord_t>(siqad::coord_t{10, 4, 0})},
//         1,
//         sidb_simulation_engine::QUICKEXACT};
//
//     const auto found_gate_layouts_cube =
//         design_sidb_gates(lyt_in_cube_coord, std::vector<tt>{create_xnor_tt()}, params_cube);
//
//     REQUIRE(found_gate_layouts_cube.size() == 1);
//     CHECK(found_gate_layouts_cube[0].num_cells() == 14);
//     CHECK(found_gate_layouts_cube[0].get_cell_type(siqad::to_fiction_coord<cube::coord_t>(siqad::coord_t{10, 4, 0}))
//     ==
//           siqad_layout::technology::NORMAL);
//
//     // using offset coordinates
//     const auto lyt_in_offset_coord = convert_to_fiction_coordinates<offset_layout>(lyt);
//     const design_sidb_gates_params<offset_layout> params_offset{
//         sidb_simulation_parameters{2, -0.32},
//         design_sidb_gates_params<offset_layout>::design_sidb_gates_mode::EXHAUSTIVE,
//         {siqad::to_fiction_coord<offset::ucoord_t>(siqad::coord_t{10, 4, 0}),
//          siqad::to_fiction_coord<offset::ucoord_t>(siqad::coord_t{10, 4, 0})},
//         1,
//         sidb_simulation_engine::QUICKEXACT};
//
//     const auto found_gate_layouts_offset =
//         design_sidb_gates(lyt_in_offset_coord, std::vector<tt>{create_xnor_tt()}, params_offset);
//
//     REQUIRE(found_gate_layouts_offset.size() == 1);
//     CHECK(found_gate_layouts_offset[0].num_cells() == 14);
//     CHECK(found_gate_layouts_offset[0].get_cell_type(
//               siqad::to_fiction_coord<offset::ucoord_t>(siqad::coord_t{10, 4, 0})) ==
//               siqad_layout::technology::NORMAL);
// }
//
// TEST_CASE("Use SiQAD's AND gate skeleton to generate all possible AND gates", "[design-sidb-gates]")
//{
//     sidb_100_cell_clk_lyt_siqad lyt{};
//
//     lyt.assign_cell_type({0, 0, 1}, sidb_technology::cell_type::INPUT);
//     lyt.assign_cell_type({2, 1, 1}, sidb_technology::cell_type::INPUT);
//
//     lyt.assign_cell_type({20, 0, 1}, sidb_technology::cell_type::INPUT);
//     lyt.assign_cell_type({18, 1, 1}, sidb_technology::cell_type::INPUT);
//
//     lyt.assign_cell_type({4, 2, 1}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({6, 3, 1}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({14, 3, 1}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({16, 2, 1}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({10, 6, 0}, sidb_technology::cell_type::OUTPUT);
//     lyt.assign_cell_type({10, 7, 0}, sidb_technology::cell_type::OUTPUT);
//
//     lyt.assign_cell_type({10, 9, 1}, sidb_technology::cell_type::NORMAL);
//
//     design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad> params{
//         sidb_simulation_parameters{2, -0.28},
//         design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::EXHAUSTIVE,
//         {{4, 4, 0}, {14, 5, 1}},
//         1,
//         sidb_simulation_engine::EXGS};
//
//     SECTION("Exhaustive Generation")
//     {
//         const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params);
//         CHECK(!found_gate_layouts.empty());
//     }
//
//     SECTION("Random Generation")
//     {
//         params.design_mode = design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::RANDOM;
//         const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params);
//         CHECK(!found_gate_layouts.empty());
//     }
// }
//
// TEST_CASE("Use FO2 Bestagon gate without SiDB at {17, 11, 0} and generate original one", "[design-sidb-gates]")
//{
//     sidb_100_cell_clk_lyt_siqad lyt{};
//
//     lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::INPUT);
//     lyt.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::INPUT);
//
//     // SiDB, originally part of the Bestagon fo2 gate, is excluded.
//     // lyt.assign_cell_type({17, 11, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({21, 11, 1}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({12, 4, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({18, 13, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({6, 2, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({8, 3, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({19, 7, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({14, 5, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({18, 6, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({24, 15, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({26, 16, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({12, 16, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({14, 15, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({8, 17, 0}, sidb_technology::cell_type::OUTPUT);
//     lyt.assign_cell_type({6, 18, 0}, sidb_technology::cell_type::OUTPUT);
//
//     lyt.assign_cell_type({30, 17, 0}, sidb_technology::cell_type::OUTPUT);
//     lyt.assign_cell_type({32, 18, 0}, sidb_technology::cell_type::OUTPUT);
//
//     lyt.assign_cell_type({36, 19, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({2, 19, 0}, sidb_technology::cell_type::NORMAL);
//
//     SECTION("generate original FO2")
//     {
//         const design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad> params{
//             sidb_simulation_parameters{2, -0.32},
//             design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::EXHAUSTIVE,
//             {{17, 11, 0}, {17, 11, 0}},
//             1,
//             sidb_simulation_engine::QUICKEXACT};
//
//         CHECK(lyt.get_cell_type({17, 11, 0}) == sidb_100_cell_clk_lyt_siqad::technology::EMPTY);
//
//         // generate gate by placing one SiDB
//
//         const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_fan_out_tt()}, params);
//
//         REQUIRE(found_gate_layouts.size() == 1);
//         CHECK(found_gate_layouts[0].num_cells() == 21);
//         CHECK(found_gate_layouts[0].get_cell_type({17, 11, 0}) == sidb_100_cell_clk_lyt_siqad::technology::NORMAL);
//     }
//
//     SECTION("replace the output perturbers by equivalent negatively charged defects")
//     {
//         const design_sidb_gates_params<sidb_defect_surface<sidb_100_cell_clk_lyt_siqad>> params{
//             sidb_simulation_parameters{2, -0.32},
//             design_sidb_gates_params<
//                 sidb_defect_surface<sidb_100_cell_clk_lyt_siqad>>::design_sidb_gates_mode::EXHAUSTIVE,
//             {{17, 11, 0}, {17, 11, 0}},
//             1,
//             sidb_simulation_engine::QUICKEXACT};
//
//         sidb_defect_surface defect_layout{lyt};
//         defect_layout.assign_cell_type({36, 19, 0}, sidb_100_cell_clk_lyt_siqad::cell_type::EMPTY);
//         defect_layout.assign_cell_type({2, 19, 0}, sidb_100_cell_clk_lyt_siqad::cell_type::EMPTY);
//         CHECK(defect_layout.get_cell_type({36, 19, 0}) == sidb_100_cell_clk_lyt_siqad::cell_type::EMPTY);
//         CHECK(defect_layout.get_cell_type({2, 19, 0}) == sidb_100_cell_clk_lyt_siqad::cell_type::EMPTY);
//
//         defect_layout.assign_sidb_defect({36, 19, 0},
//                                          sidb_defect{sidb_defect_type::DB, -1,
//                                          params.simulation_parameters.epsilon_r,
//                                                      params.simulation_parameters.lambda_tf});
//         defect_layout.assign_sidb_defect({2, 19, 0},
//                                          sidb_defect{sidb_defect_type::DB, -1,
//                                          params.simulation_parameters.epsilon_r,
//                                                      params.simulation_parameters.lambda_tf});
//
//         const auto found_gate_layouts = design_sidb_gates(defect_layout, std::vector<tt>{create_fan_out_tt()},
//         params);
//
//         REQUIRE(found_gate_layouts.size() == 1);
//         CHECK(found_gate_layouts[0].num_cells() == 19);
//         CHECK(found_gate_layouts[0].get_cell_type({17, 11, 0}) == sidb_100_cell_clk_lyt_siqad::cell_type::NORMAL);
//     }
// }
//
// TEST_CASE("Design AND Bestagon shaped gate", "[design-sidb-gates]")
//{
//     sidb_100_cell_clk_lyt_siqad lyt{};
//
//     lyt.assign_cell_type({38, 0, 0}, sidb_technology::cell_type::INPUT);
//     lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::INPUT);
//
//     lyt.assign_cell_type({36, 1, 0}, sidb_technology::cell_type::INPUT);
//     lyt.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::INPUT);
//
//     lyt.assign_cell_type({6, 2, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({32, 2, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({30, 3, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({8, 3, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({26, 4, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({12, 4, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({24, 5, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({14, 5, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({24, 15, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({26, 16, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({30, 17, 0}, sidb_technology::cell_type::OUTPUT);
//     lyt.assign_cell_type({32, 18, 0}, sidb_technology::cell_type::OUTPUT);
//     lyt.assign_cell_type({36, 19, 0}, sidb_technology::cell_type::NORMAL);
//
//     SECTION("Random Generation")
//     {
//         const design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad> params{
//             sidb_simulation_parameters{2, -0.32},
//             design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::RANDOM,
//             {{14, 6, 0}, {24, 12, 0}},
//             3,
//             sidb_simulation_engine::QUICKEXACT};
//
//         const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params);
//         REQUIRE(!found_gate_layouts.empty());
//         CHECK(found_gate_layouts.front().num_cells() == lyt.num_cells() + 3);
//     }
//
//     SECTION("Random Generation with defects")
//     {
//         sidb_defect_surface defect_layout{lyt};
//
//         const design_sidb_gates_params<sidb_defect_surface<sidb_100_cell_clk_lyt_siqad>> params{
//             sidb_simulation_parameters{2, -0.32},
//             design_sidb_gates_params<sidb_defect_surface<sidb_100_cell_clk_lyt_siqad>>::design_sidb_gates_mode::RANDOM,
//             {{14, 6, 0}, {24, 12, 0}},
//             3,
//             sidb_simulation_engine::QUICKEXACT};
//
//         defect_layout.assign_sidb_defect({15, 10, 0},
//                                          sidb_defect{sidb_defect_type::DB, -1,
//                                          params.simulation_parameters.epsilon_r,
//                                                      params.simulation_parameters.lambda_tf});
//         defect_layout.assign_sidb_defect({20, 12, 0},
//                                          sidb_defect{sidb_defect_type::DB, -1,
//                                          params.simulation_parameters.epsilon_r,
//                                                      params.simulation_parameters.lambda_tf});
//
//         const auto found_gate_layouts = design_sidb_gates(defect_layout, std::vector<tt>{create_and_tt()}, params);
//         REQUIRE(!found_gate_layouts.empty());
//         CHECK(found_gate_layouts.front().num_defects() == 2);
//         CHECK(found_gate_layouts.front().num_cells() == lyt.num_cells() + 3);
//
//         found_gate_layouts.front().foreach_cell(
//             [](const auto& cell)
//             {
//                 CHECK(cell != siqad::coord_t{15, 10, 0});
//                 CHECK(cell != siqad::coord_t{20, 12, 0});
//             });
//     }
// }
//
// TEST_CASE("Design AND Bestagon shaped gate on H-Si 111", "[design-sidb-gates]")
//{
//     sidb_111_cell_clk_lyt_siqad lyt{};
//
//     lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::INPUT);
//     lyt.assign_cell_type({25, 0, 0}, sidb_technology::cell_type::INPUT);
//
//     lyt.assign_cell_type({23, 1, 1}, sidb_technology::cell_type::INPUT);
//     lyt.assign_cell_type({1, 1, 1}, sidb_technology::cell_type::INPUT);
//
//     lyt.assign_cell_type({4, 4, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({21, 2, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({5, 5, 1}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({19, 5, 1}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({8, 8, 0}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({17, 8, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({9, 9, 1}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({15, 9, 1}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({15, 21, 1}, sidb_technology::cell_type::NORMAL);
//     lyt.assign_cell_type({17, 23, 0}, sidb_technology::cell_type::NORMAL);
//
//     lyt.assign_cell_type({19, 25, 1}, sidb_technology::cell_type::OUTPUT);
//     lyt.assign_cell_type({21, 27, 0}, sidb_technology::cell_type::OUTPUT);
//
//     lyt.assign_cell_type({23, 29, 1}, sidb_technology::cell_type::NORMAL);
//
//     SECTION("Random Generation")
//     {
//         const design_sidb_gates_params<sidb_111_cell_clk_lyt_siqad> params{
//             sidb_simulation_parameters{2, -0.32},
//             design_sidb_gates_params<sidb_111_cell_clk_lyt_siqad>::design_sidb_gates_mode::RANDOM,
//             {{10, 11, 0}, {14, 17, 0}},
//             3,
//             sidb_simulation_engine::QUICKEXACT};
//
//         const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_nor_tt()}, params);
//         REQUIRE(!found_gate_layouts.empty());
//         CHECK(found_gate_layouts.front().num_cells() == lyt.num_cells() + 3);
//     }
//
//     SECTION("Random Generation")
//     {
//         const design_sidb_gates_params<sidb_111_cell_clk_lyt_siqad> params{
//             sidb_simulation_parameters{2, -0.32},
//             design_sidb_gates_params<sidb_111_cell_clk_lyt_siqad>::design_sidb_gates_mode::EXHAUSTIVE,
//             {{10, 11, 0}, {14, 17, 0}},
//             3,
//             sidb_simulation_engine::QUICKEXACT};
//
//         const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_nor_tt()}, params);
//         REQUIRE(found_gate_layouts.size() == 206);
//         CHECK(found_gate_layouts.front().num_cells() == lyt.num_cells() + 3);
//     }
// }

template <typename Lyt>
Lyt determine_canvas_sidbs(const Lyt& lyt_with_canvas_sidbs, const Lyt& lyt_skeleton) noexcept
{
    std::vector<cell<Lyt>> canvas_sidbs{};
    lyt_with_canvas_sidbs.foreach_cell(
        [&](const typename Lyt::cell& c)
        {
            if (lyt_skeleton.get_cell_type(c) == sidb_technology::cell_type::EMPTY)
            {
                canvas_sidbs.push_back(c);
            }
        });
    Lyt lyt{};

    for (const auto& c : canvas_sidbs)
    {
        lyt.assign_cell_type(c, sidb_technology::cell_type::NORMAL);
    }

    return lyt;
}

// TEST_CASE("Efficient gate design, AND", "[design-sidb-gates]")
//{
//     const auto lyt = blueprints::two_input_one_output_skeleton<sidb_100_cell_clk_lyt_siqad>();
//
//     //print_sidb_layout(std::cout, lyt);
//
//     auto lyt_00_skeleton = lyt.clone();
//     lyt_00_skeleton.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::EMPTY);
//     lyt_00_skeleton.assign_cell_type({36, 1, 0}, sidb_technology::cell_type::EMPTY);
//
//     charge_distribution_surface cds_skeleton{lyt_00_skeleton};
//
//     cds_skeleton.assign_charge_state({0, 0, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton.assign_charge_state({38, 0, 0}, sidb_charge_state::NEGATIVE);
//
//     cds_skeleton.assign_charge_state({6, 2, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton.assign_charge_state({32, 2, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton.assign_charge_state({8, 3, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton.assign_charge_state({30, 3, 0}, sidb_charge_state::NEUTRAL);
//
//     cds_skeleton.assign_charge_state({12, 4, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton.assign_charge_state({26, 4, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton.assign_charge_state({14, 5, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton.assign_charge_state({24, 5, 0}, sidb_charge_state::NEUTRAL);
//
//     cds_skeleton.assign_charge_state({24, 15, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton.assign_charge_state({26, 16, 0}, sidb_charge_state::NEUTRAL);
//
//     cds_skeleton.assign_charge_state({30, 17, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton.assign_charge_state({32, 18, 0}, sidb_charge_state::NEUTRAL);
//
//     cds_skeleton.assign_charge_state({36, 19, 0}, sidb_charge_state::NEGATIVE);
//
//     print_sidb_layout(std::cout, cds_skeleton);
//
//     auto lyt_01_skeleton = lyt.clone();
//     lyt_01_skeleton.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::EMPTY);
//     lyt_01_skeleton.assign_cell_type({38, 0, 0}, sidb_technology::cell_type::EMPTY);
//
//     charge_distribution_surface cds_skeleton_01{lyt_01_skeleton};
//
//     cds_skeleton_01.assign_charge_state({0, 0, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_01.assign_charge_state({36, 1, 0}, sidb_charge_state::NEGATIVE);
//
//     cds_skeleton_01.assign_charge_state({6, 2, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_01.assign_charge_state({32, 2, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton_01.assign_charge_state({8, 3, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton_01.assign_charge_state({30, 3, 0}, sidb_charge_state::NEGATIVE);
//
//     cds_skeleton_01.assign_charge_state({12, 4, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_01.assign_charge_state({26, 4, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton_01.assign_charge_state({14, 5, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton_01.assign_charge_state({24, 5, 0}, sidb_charge_state::NEGATIVE);
//
//     cds_skeleton_01.assign_charge_state({24, 15, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_01.assign_charge_state({26, 16, 0}, sidb_charge_state::NEUTRAL);
//
//     cds_skeleton_01.assign_charge_state({30, 17, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_01.assign_charge_state({32, 18, 0}, sidb_charge_state::NEUTRAL);
//
//     cds_skeleton_01.assign_charge_state({36, 19, 0}, sidb_charge_state::NEGATIVE);
//
//     print_sidb_layout(std::cout, cds_skeleton_01);
//
//
//     auto lyt_10_skeleton = lyt.clone();
//     lyt_10_skeleton.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::EMPTY);
//     lyt_10_skeleton.assign_cell_type({36, 1, 0}, sidb_technology::cell_type::EMPTY);
//
//     charge_distribution_surface cds_skeleton_10{lyt_10_skeleton};
//
//     cds_skeleton_10.assign_charge_state({2, 1, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_10.assign_charge_state({38, 0, 0}, sidb_charge_state::NEGATIVE);
//
//     cds_skeleton_10.assign_charge_state({6, 2, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton_10.assign_charge_state({32, 2, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_10.assign_charge_state({8, 3, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_10.assign_charge_state({30, 3, 0}, sidb_charge_state::NEUTRAL);
//
//     cds_skeleton_10.assign_charge_state({12, 4, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton_10.assign_charge_state({26, 4, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_10.assign_charge_state({14, 5, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_10.assign_charge_state({24, 5, 0}, sidb_charge_state::NEUTRAL);
//
//     cds_skeleton_10.assign_charge_state({24, 15, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_10.assign_charge_state({26, 16, 0}, sidb_charge_state::NEUTRAL);
//
//     cds_skeleton_10.assign_charge_state({30, 17, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_10.assign_charge_state({32, 18, 0}, sidb_charge_state::NEUTRAL);
//
//     cds_skeleton_10.assign_charge_state({36, 19, 0}, sidb_charge_state::NEGATIVE);
//
//     print_sidb_layout(std::cout, cds_skeleton_10);
//
//
//     auto lyt_11_skeleton = lyt.clone();
//     lyt_11_skeleton.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::EMPTY);
//     lyt_11_skeleton.assign_cell_type({38, 0, 0}, sidb_technology::cell_type::EMPTY);
//
//     charge_distribution_surface cds_skeleton_11{lyt_11_skeleton};
//
//     cds_skeleton_11.assign_charge_state({2, 1, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_11.assign_charge_state({36, 1, 0}, sidb_charge_state::NEGATIVE);
//
//     cds_skeleton_11.assign_charge_state({6, 2, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton_11.assign_charge_state({32, 2, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton_11.assign_charge_state({8, 3, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_11.assign_charge_state({30, 3, 0}, sidb_charge_state::NEGATIVE);
//
//     cds_skeleton_11.assign_charge_state({12, 4, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton_11.assign_charge_state({26, 4, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton_11.assign_charge_state({14, 5, 0}, sidb_charge_state::NEGATIVE);
//     cds_skeleton_11.assign_charge_state({24, 5, 0}, sidb_charge_state::NEGATIVE);
//
//     cds_skeleton_11.assign_charge_state({24, 15, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton_11.assign_charge_state({26, 16, 0}, sidb_charge_state::NEGATIVE);
//
//     cds_skeleton_11.assign_charge_state({30, 17, 0}, sidb_charge_state::NEUTRAL);
//     cds_skeleton_11.assign_charge_state({32, 18, 0}, sidb_charge_state::NEGATIVE);
//
//     cds_skeleton_11.assign_charge_state({36, 19, 0}, sidb_charge_state::NEGATIVE);
//
//     print_sidb_layout(std::cout, cds_skeleton_11);
//
//     const design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad> params{
//         sidb_simulation_parameters{2, -0.32},
//         design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::EXHAUSTIVE,
//         {{14, 6, 0}, {24, 12, 0}},
//         2,
//         sidb_simulation_engine::QUICKEXACT};
//
//     //print_sidb_layout(std::cout, lyt);
//
//     auto found_gate_layouts = design_sidb_gate_candidates(lyt, std::vector<tt>{create_and_tt()}, params);
//
//     std::vector<sidb_100_cell_clk_lyt_siqad> valid_candidates{};
//
//     for (auto &gate_candidate : found_gate_layouts)
//     {
//         const auto gate_copy = gate_candidate.clone();
////        print_sidb_layout(std::cout, gate_candidate);
//
//        gate_candidate.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::EMPTY);
//        gate_candidate.assign_cell_type({36, 1, 0}, sidb_technology::cell_type::EMPTY);
//
//        charge_distribution_surface candidate_cds{gate_candidate};
//        auto canvas_lyt = determine_canvas_sidbs(gate_candidate, lyt);
//        charge_distribution_surface canvas_cds{canvas_lyt, params.simulation_parameters};
//
//        cds_skeleton.foreach_cell(
//            [&](const auto& c)
//            {
//                candidate_cds.assign_charge_state(c, cds_skeleton.get_charge_state(c));
//            });
//
//        while (canvas_cds.get_charge_index_and_base().first < canvas_cds.get_max_charge_index())
//        {
//            canvas_cds.foreach_cell(
//                [&](const auto& c)
//                {
//                    candidate_cds.assign_charge_state(c, canvas_cds.get_charge_state(c));
//                });
//            candidate_cds.update_after_charge_change();
//            //print_sidb_layout(std::cout, candidate_cds);
//            if (candidate_cds.is_physically_valid())
//            {
//                valid_candidates.push_back(gate_copy);
//            }
//            canvas_cds.increase_charge_index_by_one();
//        }
//    }
//
//    std::vector<sidb_100_cell_clk_lyt_siqad> second_round{};
//
//    //print_sidb_layout(std::cout, valid_candidates[0]);
//
//    for (auto &gate_candidate_second : valid_candidates)
//    {
//        const auto gate_copy_second = gate_candidate_second.clone();
//        //        print_sidb_layout(std::cout, gate_candidate);
//        gate_candidate_second.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::EMPTY);
//        gate_candidate_second.assign_cell_type({38, 0, 0}, sidb_technology::cell_type::EMPTY);
//
//        charge_distribution_surface candidate_cds{gate_candidate_second};
//        auto canvas_lyt = determine_canvas_sidbs(gate_candidate_second, lyt);
//        charge_distribution_surface canvas_cds{canvas_lyt, params.simulation_parameters};
//
//        cds_skeleton_01.foreach_cell(
//            [&](const auto& c)
//            {
//                candidate_cds.assign_charge_state(c, cds_skeleton_01.get_charge_state(c));
//            });
//
//        while (canvas_cds.get_charge_index_and_base().first < canvas_cds.get_max_charge_index())
//        {
//            canvas_cds.foreach_cell(
//                [&](const auto& c)
//                {
//                    candidate_cds.assign_charge_state(c, canvas_cds.get_charge_state(c));
//                });
//            candidate_cds.update_after_charge_change();
//            //print_sidb_layout(std::cout, candidate_cds);
//            if (candidate_cds.is_physically_valid())
//            {
//                second_round.push_back(gate_copy_second);
//            }
//            canvas_cds.increase_charge_index_by_one();
//        }
//    }
//
//
//    std::vector<sidb_100_cell_clk_lyt_siqad> third_round{};
//
//    print_sidb_layout(std::cout, second_round[0]);
//
//    for (auto &gate_candidate_third : second_round)
//    {
//        const auto gate_copy_third = gate_candidate_third.clone();
//        //        print_sidb_layout(std::cout, gate_candidate);
//        gate_candidate_third.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::EMPTY);
//        gate_candidate_third.assign_cell_type({36, 1, 0}, sidb_technology::cell_type::EMPTY);
//
//        charge_distribution_surface candidate_cds{gate_candidate_third};
//        auto canvas_lyt = determine_canvas_sidbs(gate_candidate_third, lyt);
//        charge_distribution_surface canvas_cds{canvas_lyt, params.simulation_parameters};
//
//        cds_skeleton_10.foreach_cell(
//            [&](const auto& c)
//            {
//                candidate_cds.assign_charge_state(c, cds_skeleton_10.get_charge_state(c));
//            });
//
//        while (canvas_cds.get_charge_index_and_base().first < canvas_cds.get_max_charge_index())
//        {
//            canvas_cds.foreach_cell(
//                [&](const auto& c)
//                {
//                    candidate_cds.assign_charge_state(c, canvas_cds.get_charge_state(c));
//                });
//            candidate_cds.update_after_charge_change();
//            //print_sidb_layout(std::cout, candidate_cds);
//            if (candidate_cds.is_physically_valid())
//            {
//                third_round.push_back(gate_copy_third);
//            }
//            canvas_cds.increase_charge_index_by_one();
//        }
//    }
//
//
//
//    std::vector<sidb_100_cell_clk_lyt_siqad> fourth_round{};
//
//    print_sidb_layout(std::cout, third_round[0]);
//
//    for (auto &gate_candidate_forth : third_round)
//    {
//        //        print_sidb_layout(std::cout, gate_candidate);
//        gate_candidate_forth.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::EMPTY);
//        gate_candidate_forth.assign_cell_type({38, 0, 0}, sidb_technology::cell_type::EMPTY);
//
//        charge_distribution_surface candidate_cds{gate_candidate_forth};
//        auto canvas_lyt = determine_canvas_sidbs(gate_candidate_forth, lyt);
//        charge_distribution_surface canvas_cds{canvas_lyt, params.simulation_parameters};
//
//        cds_skeleton_11.foreach_cell(
//            [&](const auto& c)
//            {
//                candidate_cds.assign_charge_state(c, cds_skeleton_11.get_charge_state(c));
//            });
//
//        while (canvas_cds.get_charge_index_and_base().first < canvas_cds.get_max_charge_index())
//        {
//            canvas_cds.foreach_cell(
//                [&](const auto& c)
//                {
//                    candidate_cds.assign_charge_state(c, canvas_cds.get_charge_state(c));
//                });
//            candidate_cds.update_after_charge_change();
//            //print_sidb_layout(std::cout, candidate_cds);
//            if (candidate_cds.is_physically_valid())
//            {
//                fourth_round.push_back(gate_candidate_forth);
//            }
//            canvas_cds.increase_charge_index_by_one();
//        }
//    }
//
//
//    std::cout << found_gate_layouts.size() << std::endl;
//    std::cout << second_round.size() << std::endl;
//    std::cout << third_round.size() << std::endl;
//    std::cout << fourth_round.size() << std::endl;
//    print_sidb_layout(std::cout, fourth_round[0]);
//}

TEST_CASE("Efficient gate design, OR", "[design-sidb-gates]")
{
    const auto lyt = blueprints::two_input_one_output_skeleton_new<sidb_100_cell_clk_lyt_siqad>();



    // print_sidb_layout(std::cout, lyt);

    auto lyt_00_skeleton = lyt.clone();
    lyt_00_skeleton.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::EMPTY);
    lyt_00_skeleton.assign_cell_type({36, 1, 0}, sidb_technology::cell_type::EMPTY);

    charge_distribution_surface cds_skeleton{lyt_00_skeleton};

    cds_skeleton.assign_charge_state({0, 0, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton.assign_charge_state({38, 0, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton.assign_charge_state({6, 2, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton.assign_charge_state({32, 2, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton.assign_charge_state({8, 3, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton.assign_charge_state({30, 3, 0}, sidb_charge_state::NEUTRAL);

    cds_skeleton.assign_charge_state({12, 4, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton.assign_charge_state({26, 4, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton.assign_charge_state({14, 5, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton.assign_charge_state({24, 5, 0}, sidb_charge_state::NEUTRAL);

    cds_skeleton.assign_charge_state({19, 13, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton.assign_charge_state({20, 14, 0}, sidb_charge_state::NEUTRAL);

    cds_skeleton.assign_charge_state({24, 15, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton.assign_charge_state({26, 16, 0}, sidb_charge_state::NEUTRAL);

    cds_skeleton.assign_charge_state({30, 17, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton.assign_charge_state({32, 18, 0}, sidb_charge_state::NEUTRAL);

    cds_skeleton.assign_charge_state({36, 19, 0}, sidb_charge_state::NEGATIVE);

    print_sidb_layout(std::cout, cds_skeleton);

    auto lyt_01_skeleton = lyt.clone();
    lyt_01_skeleton.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::EMPTY);
    lyt_01_skeleton.assign_cell_type({38, 0, 0}, sidb_technology::cell_type::EMPTY);

    charge_distribution_surface cds_skeleton_01{lyt_01_skeleton};

    cds_skeleton_01.assign_charge_state({0, 0, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton_01.assign_charge_state({36, 1, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_01.assign_charge_state({6, 2, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton_01.assign_charge_state({32, 2, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_01.assign_charge_state({8, 3, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_01.assign_charge_state({30, 3, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_01.assign_charge_state({12, 4, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton_01.assign_charge_state({26, 4, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_01.assign_charge_state({14, 5, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_01.assign_charge_state({24, 5, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_01.assign_charge_state({19, 13, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_01.assign_charge_state({20, 14, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_01.assign_charge_state({24, 15, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_01.assign_charge_state({26, 16, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_01.assign_charge_state({30, 17, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_01.assign_charge_state({32, 18, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_01.assign_charge_state({36, 19, 0}, sidb_charge_state::NEGATIVE);

    print_sidb_layout(std::cout, cds_skeleton_01);

    auto lyt_10_skeleton = lyt.clone();
    lyt_10_skeleton.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::EMPTY);
    lyt_10_skeleton.assign_cell_type({36, 1, 0}, sidb_technology::cell_type::EMPTY);

    charge_distribution_surface cds_skeleton_10{lyt_10_skeleton};

    cds_skeleton_10.assign_charge_state({2, 1, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton_10.assign_charge_state({38, 0, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_10.assign_charge_state({6, 2, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_10.assign_charge_state({32, 2, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton_10.assign_charge_state({8, 3, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton_10.assign_charge_state({30, 3, 0}, sidb_charge_state::NEUTRAL);

    cds_skeleton_10.assign_charge_state({12, 4, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_10.assign_charge_state({26, 4, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton_10.assign_charge_state({14, 5, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton_10.assign_charge_state({24, 5, 0}, sidb_charge_state::NEUTRAL);

    cds_skeleton_10.assign_charge_state({19, 13, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_10.assign_charge_state({20, 14, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_10.assign_charge_state({24, 15, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_10.assign_charge_state({26, 16, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_10.assign_charge_state({30, 17, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_10.assign_charge_state({32, 18, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_10.assign_charge_state({36, 19, 0}, sidb_charge_state::NEGATIVE);

    print_sidb_layout(std::cout, cds_skeleton_10);

    auto lyt_11_skeleton = lyt.clone();
    lyt_11_skeleton.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::EMPTY);
    lyt_11_skeleton.assign_cell_type({38, 0, 0}, sidb_technology::cell_type::EMPTY);

    charge_distribution_surface cds_skeleton_11{lyt_11_skeleton};

    cds_skeleton_11.assign_charge_state({2, 1, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton_11.assign_charge_state({36, 1, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_11.assign_charge_state({6, 2, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_11.assign_charge_state({32, 2, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_11.assign_charge_state({8, 3, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton_11.assign_charge_state({30, 3, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_11.assign_charge_state({12, 4, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_11.assign_charge_state({26, 4, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_11.assign_charge_state({14, 5, 0}, sidb_charge_state::NEGATIVE);
    cds_skeleton_11.assign_charge_state({24, 5, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_11.assign_charge_state({19, 13, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_11.assign_charge_state({20, 14, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_11.assign_charge_state({24, 15, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_11.assign_charge_state({26, 16, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_11.assign_charge_state({30, 17, 0}, sidb_charge_state::NEUTRAL);
    cds_skeleton_11.assign_charge_state({32, 18, 0}, sidb_charge_state::NEGATIVE);

    cds_skeleton_11.assign_charge_state({36, 19, 0}, sidb_charge_state::NEGATIVE);

    print_sidb_layout(std::cout, cds_skeleton_11);

    const design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad> params{
        sidb_simulation_parameters{2, -0.32},
        design_sidb_gates_params<sidb_100_cell_clk_lyt_siqad>::design_sidb_gates_mode::EXHAUSTIVE,
        {{14, 6, 0}, {24, 11, 0}},
        3,
        sidb_simulation_engine::QUICKEXACT};

    // print_sidb_layout(std::cout, lyt);
//    auto exh = design_sidb_gates(lyt, std::vector<tt>{create_or_tt()}, params);
//    std::cout << exh.size() << std::endl;

    auto found_gate_layouts = design_sidb_gate_candidates(lyt, std::vector<tt>{create_or_tt()}, params);

    std::vector<sidb_100_cell_clk_lyt_siqad> valid_candidates{};

    for (auto& gate_candidate : found_gate_layouts)
    {
        const auto gate_copy = gate_candidate.clone();
        //        print_sidb_layout(std::cout, gate_candidate);

        gate_candidate.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::EMPTY);
        gate_candidate.assign_cell_type({36, 1, 0}, sidb_technology::cell_type::EMPTY);

        charge_distribution_surface candidate_cds{gate_candidate};
        auto                        canvas_lyt = determine_canvas_sidbs(gate_candidate, lyt);
        charge_distribution_surface canvas_cds{canvas_lyt, params.simulation_parameters};

        cds_skeleton.foreach_cell([&](const auto& c)
                                  { candidate_cds.assign_charge_state(c, cds_skeleton.get_charge_state(c)); });

        while (canvas_cds.get_charge_index_and_base().first < canvas_cds.get_max_charge_index())
        {
            canvas_cds.foreach_cell([&](const auto& c)
                                    { candidate_cds.assign_charge_state(c, canvas_cds.get_charge_state(c)); });
            candidate_cds.update_after_charge_change();
            // print_sidb_layout(std::cout, candidate_cds);
            if (candidate_cds.is_physically_valid())
            {
                valid_candidates.push_back(gate_copy);
                break;
            }
            canvas_cds.increase_charge_index_by_one();
        }
    }

    std::vector<sidb_100_cell_clk_lyt_siqad> second_round{};

    // print_sidb_layout(std::cout, valid_candidates[0]);

    for (auto& gate_candidate_second : valid_candidates)
    {
        const auto gate_copy_second = gate_candidate_second.clone();
        //        print_sidb_layout(std::cout, gate_candidate);
        gate_candidate_second.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::EMPTY);
        gate_candidate_second.assign_cell_type({38, 0, 0}, sidb_technology::cell_type::EMPTY);

        charge_distribution_surface candidate_cds{gate_candidate_second};
        auto                        canvas_lyt = determine_canvas_sidbs(gate_candidate_second, lyt);
        charge_distribution_surface canvas_cds{canvas_lyt, params.simulation_parameters};

        cds_skeleton_01.foreach_cell([&](const auto& c)
                                     { candidate_cds.assign_charge_state(c, cds_skeleton_01.get_charge_state(c)); });

        while (canvas_cds.get_charge_index_and_base().first < canvas_cds.get_max_charge_index())
        {
            canvas_cds.foreach_cell([&](const auto& c)
                                    { candidate_cds.assign_charge_state(c, canvas_cds.get_charge_state(c)); });
            candidate_cds.update_after_charge_change();
            // print_sidb_layout(std::cout, candidate_cds);
            if (candidate_cds.is_physically_valid())
            {
                second_round.push_back(gate_copy_second);
                break;
            }
            canvas_cds.increase_charge_index_by_one();
        }
    }

    std::vector<sidb_100_cell_clk_lyt_siqad> third_round{};

    print_sidb_layout(std::cout, second_round[0]);

    for (auto& gate_candidate_third : second_round)
    {
        const auto gate_copy_third = gate_candidate_third.clone();
        //        print_sidb_layout(std::cout, gate_candidate);
        gate_candidate_third.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::EMPTY);
        gate_candidate_third.assign_cell_type({36, 1, 0}, sidb_technology::cell_type::EMPTY);

        charge_distribution_surface candidate_cds{gate_candidate_third};
        auto                        canvas_lyt = determine_canvas_sidbs(gate_candidate_third, lyt);
        charge_distribution_surface canvas_cds{canvas_lyt, params.simulation_parameters};

        cds_skeleton_10.foreach_cell([&](const auto& c)
                                     { candidate_cds.assign_charge_state(c, cds_skeleton_10.get_charge_state(c)); });

        while (canvas_cds.get_charge_index_and_base().first < canvas_cds.get_max_charge_index())
        {
            canvas_cds.foreach_cell([&](const auto& c)
                                    { candidate_cds.assign_charge_state(c, canvas_cds.get_charge_state(c)); });
            candidate_cds.update_after_charge_change();
            // print_sidb_layout(std::cout, candidate_cds);
            if (candidate_cds.is_physically_valid())
            {
                third_round.push_back(gate_copy_third);
                write_sqd_layout(gate_copy_third, "/Users/jandrewniok/Desktop/and_efficient_gate.sqd");
                break;
            }
            canvas_cds.increase_charge_index_by_one();
        }
    }

    std::vector<sidb_100_cell_clk_lyt_siqad> fourth_round{};

    print_sidb_layout(std::cout, third_round[0]);

    for (auto& gate_candidate_forth : third_round)
    {
        const auto gate_copy_forth = gate_candidate_forth.clone();
        //        print_sidb_layout(std::cout, gate_candidate);
        gate_candidate_forth.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::EMPTY);
        gate_candidate_forth.assign_cell_type({38, 0, 0}, sidb_technology::cell_type::EMPTY);

        charge_distribution_surface candidate_cds{gate_candidate_forth};
        auto                        canvas_lyt = determine_canvas_sidbs(gate_candidate_forth, lyt);
        charge_distribution_surface canvas_cds{canvas_lyt, params.simulation_parameters};

        cds_skeleton_11.foreach_cell([&](const auto& c)
                                     { candidate_cds.assign_charge_state(c, cds_skeleton_11.get_charge_state(c)); });

        while (canvas_cds.get_charge_index_and_base().first < canvas_cds.get_max_charge_index())
        {
            canvas_cds.foreach_cell([&](const auto& c)
                                    { candidate_cds.assign_charge_state(c, canvas_cds.get_charge_state(c)); });
            candidate_cds.update_after_charge_change();
            // print_sidb_layout(std::cout, candidate_cds);
            if (candidate_cds.is_physically_valid())
            {
                fourth_round.push_back(gate_copy_forth);
                break;
            }
            canvas_cds.increase_charge_index_by_one();
        }
    }

    std::cout << found_gate_layouts.size() << std::endl;
    std::cout << second_round.size() << std::endl;
    std::cout << third_round.size() << std::endl;
    std::cout << fourth_round.size() << std::endl;

    print_sidb_layout(std::cout, fourth_round[0]);

    uint64_t counter = 0;

    for (const auto& gate : fourth_round)
    {
        if (is_operational(gate, std::vector<tt>{create_or_tt()}, is_operational_params{params.simulation_parameters}).first ==
            operational_status::OPERATIONAL)
        {
            counter++;
            if (counter % 10 == 0)
            {
                std::cout << counter << std::endl;
                write_sqd_layout(gate, "/Users/jandrewniok/Desktop/and_efficient_gate.sqd");
            }
        };
    }
    std::cout << counter << std::endl;
    print_sidb_layout(std::cout, fourth_round[0]);
}
