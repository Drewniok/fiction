//
// Created by Jan Drewniok on 12.09.23.
//

#include <catch2/catch_template_test_macros.hpp>

#include <fiction/algorithms/physical_design/design_sidb_gates.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_engine.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/traits.hpp>
#include <fiction/types.hpp>
#include <fiction/utils/truth_table_utils.hpp>

using namespace fiction;

TEST_CASE("Use SiQAD XNOR skeleton and generate SiQAD XNOR gate, exhaustive", "[design-sidb-gates]")
{
    using layout = sidb_cell_clk_lyt_siqad;

    layout lyt{};

    lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({20, 0, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({18, 1, 0}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({6, 3, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({14, 3, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({4, 2, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({16, 2, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({10, 6, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({10, 7, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({10, 9, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({10, 10, 0}, sidb_technology::cell_type::OUTPUT);

    lyt.assign_cell_type({10, 12, 1}, sidb_technology::cell_type::NORMAL);

    CHECK(lyt.num_cells() == 13);

    const design_sidb_gates_params params{sidb_simulation_parameters{2, -0.32},
                                          design_sidb_gates_params::design_sidb_gates_mode::EXHAUSTIVE,
                                          {{10, 4, 0}, {10, 4, 0}},
                                          1,
                                          sidb_simulation_engine::QUICKEXACT};

    const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_xnor_tt()}, params);

    REQUIRE(found_gate_layouts.size() == 1);
    CHECK(found_gate_layouts[0].num_cells() == 14);
    CHECK(found_gate_layouts[0].get_cell_type({10, 4, 0}) == layout::technology::NORMAL);
}

TEST_CASE("Use SiQAD's AND gate skeleton to generate all possible AND gates", "[design-sidb-gates]")
{
    using layout = sidb_cell_clk_lyt_siqad;

    layout lyt{};

    lyt.assign_cell_type({0, 0, 1}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({2, 1, 1}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({20, 0, 1}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({18, 1, 1}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({4, 2, 1}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({6, 3, 1}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({14, 3, 1}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({16, 2, 1}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({10, 6, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({10, 7, 0}, sidb_technology::cell_type::OUTPUT);

    lyt.assign_cell_type({10, 9, 1}, sidb_technology::cell_type::NORMAL);

    design_sidb_gates_params params{sidb_simulation_parameters{2, -0.28},
                                    design_sidb_gates_params::design_sidb_gates_mode::EXHAUSTIVE,
                                    {{4, 4, 0}, {14, 5, 1}},
                                    1,
                                    sidb_simulation_engine::QUICKEXACT};

    SECTION("Exhaustive Generation")
    {
        const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params);
        CHECK(!found_gate_layouts.empty());
    }

    SECTION("Random Generation")
    {
        params.design_mode            = design_sidb_gates_params::design_sidb_gates_mode::RANDOM;
        const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params);
        CHECK(!found_gate_layouts.empty());
    }
}

TEST_CASE("Use FO2 Bestagon gate without SiDB at {17, 11, 0} and generate original one", "[design-sidb-gates]")
{
    using layout = sidb_cell_clk_lyt_siqad;

    layout lyt{};

    lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::INPUT);

    // SiDB, originally part of the Bestagon fo2 gate, is excluded.
    // lyt.assign_cell_type({17, 11, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({21, 11, 1}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({12, 4, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({18, 13, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({6, 2, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({8, 3, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({19, 7, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({14, 5, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({18, 6, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({24, 15, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({26, 16, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({12, 16, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({14, 15, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({8, 17, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({6, 18, 0}, sidb_technology::cell_type::OUTPUT);

    lyt.assign_cell_type({30, 17, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({32, 18, 0}, sidb_technology::cell_type::OUTPUT);

    lyt.assign_cell_type({36, 19, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({2, 19, 0}, sidb_technology::cell_type::NORMAL);

    const design_sidb_gates_params params{sidb_simulation_parameters{2, -0.32},
                                          design_sidb_gates_params::design_sidb_gates_mode::EXHAUSTIVE,
                                          {{17, 11, 0}, {17, 11, 0}},
                                          1,
                                          sidb_simulation_engine::QUICKEXACT};

    SECTION("generate original FO2")
    {
        CHECK(lyt.get_cell_type({17, 11, 0}) == layout::technology::EMPTY);

        // generate gate by placing one SiDB

        const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_fan_out_tt()}, params);

        REQUIRE(found_gate_layouts.size() == 1);
        CHECK(found_gate_layouts[0].num_cells() == 21);
        CHECK(found_gate_layouts[0].get_cell_type({17, 11, 0}) == layout::technology::NORMAL);
    }

    SECTION("replace the output perturbers by equivalent negatively charged defects")
    {
        sidb_defect_cell_clk_lyt_siqad defect_layout{lyt};
        defect_layout.assign_cell_type({36, 19, 0}, technology<sidb_defect_cell_clk_lyt_siqad>::cell_type::EMPTY);
        defect_layout.assign_cell_type({2, 19, 0}, technology<sidb_defect_cell_clk_lyt_siqad>::cell_type::EMPTY);
        CHECK(defect_layout.get_cell_type({36, 19, 0}) == technology<sidb_defect_cell_clk_lyt_siqad>::cell_type::EMPTY);
        CHECK(defect_layout.get_cell_type({2, 19, 0}) == technology<sidb_defect_cell_clk_lyt_siqad>::cell_type::EMPTY);

        defect_layout.assign_sidb_defect(
            {36, 19, 0},
            sidb_defect{sidb_defect_type::DB, -1, params.phys_params.epsilon_r, params.phys_params.lambda_tf});
        defect_layout.assign_sidb_defect({2, 19, 0}, sidb_defect{sidb_defect_type::DB, -1, params.phys_params.epsilon_r,
                                                                 params.phys_params.lambda_tf});

        const auto found_gate_layouts = design_sidb_gates(defect_layout, std::vector<tt>{create_fan_out_tt()}, params);

        REQUIRE(found_gate_layouts.size() == 1);
        CHECK(found_gate_layouts[0].num_cells() == 19);
        CHECK(found_gate_layouts[0].get_cell_type({17, 11, 0}) ==
              technology<sidb_defect_cell_clk_lyt_siqad>::cell_type::NORMAL);
    }
}

TEST_CASE("Design AND Bestagon shaped gate", "[design-sidb-gates]")
{
    using layout = sidb_cell_clk_lyt_siqad;

    layout lyt{};

    lyt.assign_cell_type({38, 0, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({36, 1, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({6, 2, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({32, 2, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({30, 3, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({8, 3, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({26, 4, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({12, 4, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({24, 5, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({14, 5, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({24, 15, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({26, 16, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({30, 17, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({32, 18, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({36, 19, 0}, sidb_technology::cell_type::NORMAL);

    // generate gate by placing one SiDB
    const design_sidb_gates_params params{sidb_simulation_parameters{2, -0.32},
                                          design_sidb_gates_params::design_sidb_gates_mode::RANDOM,
                                          {{14, 6, 0}, {24, 12, 0}},
                                          3,
                                          sidb_simulation_engine::QUICKEXACT};

    SECTION("Random Generation")
    {
        const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params);
        REQUIRE(!found_gate_layouts.empty());
        CHECK(found_gate_layouts.front().num_cells() == lyt.num_cells() + 3);
    }

    SECTION("Random Generation with defects")
    {
        sidb_defect_cell_clk_lyt_siqad defect_layout{lyt};

        defect_layout.assign_sidb_defect(
            {15, 10, 0},
            sidb_defect{sidb_defect_type::DB, -1, params.phys_params.epsilon_r, params.phys_params.lambda_tf});
        defect_layout.assign_sidb_defect(
            {20, 12, 0},
            sidb_defect{sidb_defect_type::DB, -1, params.phys_params.epsilon_r, params.phys_params.lambda_tf});
        defect_layout.assign_sidb_defect({23, 12, 0}, sidb_defect{sidb_defect_type::GUNK});

        const auto found_gate_layouts = design_sidb_gates(defect_layout, std::vector<tt>{create_and_tt()}, params);
        REQUIRE(!found_gate_layouts.empty());
        CHECK(found_gate_layouts.front().num_defects() == 3);
        CHECK(found_gate_layouts.front().num_cells() == lyt.num_cells() + 3);

        found_gate_layouts.front().foreach_cell(
            [](const auto& cell)
            {
                CHECK(cell != siqad::coord_t{15, 10, 0});
                CHECK(cell != siqad::coord_t{20, 12, 0});
            });
    }
}

TEST_CASE("Design AND Bestagon shaped gate with several metrics", "[design-sidb-gates]")
{
    using layout = sidb_cell_clk_lyt_siqad;

    layout lyt{};

    lyt.assign_cell_type({38, 0, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({36, 1, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({6, 2, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({32, 2, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({30, 3, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({8, 3, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({26, 4, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({12, 4, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({24, 5, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({14, 5, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({24, 15, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({26, 16, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({30, 17, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({32, 18, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({36, 19, 0}, sidb_technology::cell_type::NORMAL);

    SECTION("Temperature as metric")
    {
        // generate gate by placing one SiDB
        design_sidb_gates_params params{sidb_simulation_parameters{2, -0.32},
                                        design_sidb_gates_params::design_sidb_gates_mode::RANDOM,
                                        {{14, 6, 0}, {24, 12, 0}},
                                        3,
                                        sidb_simulation_engine::QUICKEXACT};

        SECTION("2.0 K")
        {
            const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params,
                                                              {gate_metric::CRITICAL_TEMPERATURE, 2.0});
            REQUIRE(!found_gate_layouts.empty());
            CHECK(found_gate_layouts.front().num_cells() == lyt.num_cells() + 3);

            // check if gained layouts fulfill set conditions
            for (const auto& designed : found_gate_layouts)
            {
                const auto critical_temperature = critical_temperature_gate_based(
                    designed, std::vector<tt>{create_and_tt()}, params.metric_params.critical_temp_params);
                CHECK(critical_temperature > 2.0);
            }
        }

        SECTION("5.0 K")
        {
            const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params,
                                                              {gate_metric::CRITICAL_TEMPERATURE, 5.0});
            REQUIRE(!found_gate_layouts.empty());
            CHECK(found_gate_layouts.front().num_cells() == lyt.num_cells() + 3);

            // check if gained layouts fulfill set conditions
            for (const auto& designed : found_gate_layouts)
            {
                const auto critical_temperature = critical_temperature_gate_based(
                    designed, std::vector<tt>{create_and_tt()}, params.metric_params.critical_temp_params);
                CHECK(critical_temperature > 5.0);
            }
        }
    }

    SECTION("Population stability as metric")
    {
        // generate gate by placing one SiDB
        design_sidb_gates_params params{sidb_simulation_parameters{2, -0.32},
                                        design_sidb_gates_params::design_sidb_gates_mode::RANDOM,
                                        {{14, 6, 0}, {24, 12, 0}},
                                        3,
                                        sidb_simulation_engine::QUICKEXACT};

        SECTION("5 nm")
        {
            const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params,
                                                              {gate_metric::POPULATION_STABILITY, 5.0});
            REQUIRE(!found_gate_layouts.empty());
            CHECK(found_gate_layouts.front().num_cells() == lyt.num_cells() + 3);

            // check if gained layouts fulfill set conditions
            for (const auto& designed : found_gate_layouts)
            {
                const auto result =
                    assess_physical_population_stability(designed, params.metric_params.population_stability_params);
                REQUIRE(!result.empty());
                CHECK(result[0].distance_corresponding_to_potential > 5.0);
            }
        }

        SECTION("10 nm")
        {
            const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params,
                                                              {gate_metric::POPULATION_STABILITY, 10.0});
            REQUIRE(!found_gate_layouts.empty());
            CHECK(found_gate_layouts.front().num_cells() == lyt.num_cells() + 3);

            // check if gained layouts fulfill set conditions
            for (const auto& designed : found_gate_layouts)
            {
                const auto result =
                    assess_physical_population_stability(designed, params.metric_params.population_stability_params);
                REQUIRE(!result.empty());
                CHECK(result[0].distance_corresponding_to_potential > 10.0);
            }
        }
    }
}

TEST_CASE("Use SiQAD's AND gate skeleton", "[design-sidb-gates]")
{
    using layout = sidb_cell_clk_lyt_siqad;

    layout lyt{};

    lyt.assign_cell_type({0, 0, 1}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({2, 1, 1}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({20, 0, 1}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({18, 1, 1}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({4, 2, 1}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({6, 3, 1}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({14, 3, 1}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({16, 2, 1}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({10, 6, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({10, 7, 0}, sidb_technology::cell_type::OUTPUT);

    lyt.assign_cell_type({10, 9, 1}, sidb_technology::cell_type::NORMAL);

    design_sidb_gates_params params{sidb_simulation_parameters{2, -0.28},
                                    design_sidb_gates_params::design_sidb_gates_mode::EXHAUSTIVE,
                                    {{4, 4, 0}, {14, 5, 1}},
                                    1,
                                    sidb_simulation_engine::QUICKEXACT};

    SECTION("Operational domain as metric")
    {
        SECTION("1 % operational domain area")
        {
            const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params,
                                                              {gate_metric::OPERATIONAL_DOMAIN, 0.01});
            REQUIRE(!found_gate_layouts.empty());
            CHECK(found_gate_layouts.front().num_cells() == lyt.num_cells() + 1);

            // check if gained layouts fulfill set conditions
            for (const auto& designed : found_gate_layouts)
            {
                operational_domain_stats op_domain_stats{};
                const auto               op_domain = operational_domain_flood_fill(
                    designed, std::vector<tt>{create_and_tt()}, 1, params.metric_params.op_domain_params,
                    operational_domain::parameter_point(params.phys_params.epsilon_r, params.phys_params.lambda_tf),
                    &op_domain_stats);
                CHECK(op_domain_stats.percentual_operational_area > 0.01);
            }
        }

        SECTION("10 % operational domain area")
        {
            const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params,
                                                              {gate_metric::OPERATIONAL_DOMAIN, 0.1});
            REQUIRE(!found_gate_layouts.empty());
            CHECK(found_gate_layouts.front().num_cells() == lyt.num_cells() + 1);

            // check if gained layouts fulfill set conditions
            for (const auto& designed : found_gate_layouts)
            {
                operational_domain_stats op_domain_stats{};
                const auto               op_domain = operational_domain_flood_fill(
                    designed, std::vector<tt>{create_and_tt()}, 1, params.metric_params.op_domain_params,
                    operational_domain::parameter_point(params.phys_params.epsilon_r, params.phys_params.lambda_tf),
                    &op_domain_stats);
                CHECK(op_domain_stats.percentual_operational_area > 0.1);
            }
        }
    }

    SECTION("Avoidance distance as metric")
    {
        SECTION("5 nm")
        {
            const auto found_gate_layouts = design_sidb_gates(lyt, std::vector<tt>{create_and_tt()}, params,
                                                              {gate_metric::MAXIMUM_DEFECT_INFLUENCE_DISTANCE, 5.0});
            REQUIRE(!found_gate_layouts.empty());
            CHECK(found_gate_layouts.front().num_cells() == lyt.num_cells() + 1);

            // check if gained layouts fulfill set conditions
            const auto [defect_pos, distance] = maximum_defect_influence_position_and_distance(
                found_gate_layouts[0], params.metric_params.influence_defect_distance_param);
            CHECK(distance < 5.0);
        }
    }
}
