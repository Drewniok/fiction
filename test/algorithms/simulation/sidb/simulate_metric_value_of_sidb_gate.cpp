//
// Created by Jan Drewniok on 20.11.23.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fiction/algorithms/simulation/sidb/simulate_metric_value_of_sidb_gate.hpp>
#include <fiction/technology/physical_constants.hpp>
#include <fiction/types.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace fiction;

TEST_CASE("BDL wire operational domain computation", "[simulate-etric-value-of-sidb-gate]")
{
    using layout = sidb_cell_clk_lyt_siqad;

    layout lyt{{24, 0}, "BDL wire"};

    lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({3, 0, 0}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({6, 0, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({8, 0, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({12, 0, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({14, 0, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({18, 0, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({20, 0, 0}, sidb_technology::cell_type::OUTPUT);

    // output perturber
    lyt.assign_cell_type({24, 0, 0}, sidb_technology::cell_type::NORMAL);

    SECTION("Simulate critical temperature")
    {
        auto temperature_params                               = critical_temperature_params{};
        temperature_params.simulation_params.phys_params.base = 2;
        const auto critical_temperature =
            critical_temperature_gate_based<layout>(lyt, std::vector<tt>{create_id_tt()}, temperature_params);

        auto metric_params                 = simulate_metric_value_of_sidb_gate_params{};
        metric_params.critical_temp_params = temperature_params;

        const auto critical_temperature_via_metric_information = simulate_metric_value_of_sidb_gate(
            lyt, std::vector<tt>{create_id_tt()}, gate_metric::CRITICAL_TEMPERATURE, metric_params);

        CHECK_THAT(std::abs(critical_temperature_via_metric_information - critical_temperature),
                   Catch::Matchers::WithinAbs(0.0, fiction::physical_constants::POP_STABILITY_ERR));
    }

    SECTION("Simulate population stability")
    {
        const auto pop_stability_params =
            assess_physical_population_stability_params{sidb_simulation_parameters{2, -0.32}, 2};

        // input 0
        lyt.assign_cell_type({3, 0, 0}, sidb_technology::cell_type::EMPTY);
        const auto result_input_0 = assess_physical_population_stability(lyt, pop_stability_params);
        REQUIRE(!result_input_0.empty());

        // input 1
        lyt.assign_cell_type({3, 0, 0}, sidb_technology::cell_type::NORMAL);
        lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::EMPTY);
        const auto result_input_1 = assess_physical_population_stability(lyt, pop_stability_params);
        REQUIRE(!result_input_1.empty());

        const auto& total_population_stability_distance =
            std::max(result_input_0.front().distance_corresponding_to_potential,
                     result_input_1.front().distance_corresponding_to_potential);

        auto metric_params                        = simulate_metric_value_of_sidb_gate_params{};
        metric_params.population_stability_params = pop_stability_params;

        const auto population_stability_via_metric_information = simulate_metric_value_of_sidb_gate(
            lyt, std::vector<tt>{create_id_tt()}, gate_metric::POPULATION_STABILITY, metric_params);

        CHECK_THAT(std::abs(population_stability_via_metric_information - total_population_stability_distance),
                   Catch::Matchers::WithinAbs(0.0, fiction::physical_constants::POP_STABILITY_ERR));
    }

    SECTION("Simulate maximum defect influence distance")
    {
        const sidb_defect defect{sidb_defect_type::UNKNOWN, -1, sidb_simulation_parameters{}.epsilon_r,
                                 sidb_simulation_parameters{}.lambda_tf};
        const auto        defect_influence_params =
            maximum_defect_influence_distance_params{defect, sidb_simulation_parameters{}, {50, 6}};

        // input 0
        lyt.assign_cell_type({3, 0, 0}, sidb_technology::cell_type::EMPTY);
        const auto [defect_pos_input_0, distance_input_0] =
            maximum_defect_influence_position_and_distance(lyt, defect_influence_params);

        // input 1
        lyt.assign_cell_type({3, 0, 0}, sidb_technology::cell_type::NORMAL);
        lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::EMPTY);
        const auto [defect_pos_input_1, distance_input_1] =
            maximum_defect_influence_position_and_distance(lyt, defect_influence_params);

        const auto& maximum_defect_infleunce_distance_for_gate = std::max(distance_input_0, distance_input_1);

        auto metric_params = simulate_metric_value_of_sidb_gate_params{
            sidb_simulation_parameters{}, sidb_simulation_engine::QUICKEXACT, critical_temperature_params{},
            assess_physical_population_stability_params{sidb_simulation_parameters{}, 2}, defect_influence_params};

        const auto population_stability_via_metric_information = simulate_metric_value_of_sidb_gate(
            lyt, std::vector<tt>{create_id_tt()}, gate_metric::MAXIMUM_DEFECT_INFLUENCE_DISTANCE, metric_params);

        CHECK_THAT(std::abs(population_stability_via_metric_information - maximum_defect_infleunce_distance_for_gate),
                   Catch::Matchers::WithinAbs(0.0, fiction::physical_constants::POP_STABILITY_ERR));
    }

    SECTION("Simulate operational domain")
    {
        operational_domain_stats op_domain_stats{};
        auto                     op_domain_params = operational_domain_params{};
        op_domain_params.sim_params.base          = 2;

        const auto op_domain =
            operational_domain_flood_fill(lyt, std::vector<tt>{create_id_tt()}, 0, op_domain_params,
                                          operational_domain::parameter_point{op_domain_params.sim_params.epsilon_r,
                                                                              op_domain_params.sim_params.lambda_tf},
                                          &op_domain_stats);

        auto metric_params             = simulate_metric_value_of_sidb_gate_params{};
        metric_params.op_domain_params = op_domain_params;

        const auto operational_domain_via_metric_information = simulate_metric_value_of_sidb_gate(
            lyt, std::vector<tt>{create_id_tt()}, gate_metric::OPERATIONAL_DOMAIN, metric_params);

        CHECK_THAT(std::abs(operational_domain_via_metric_information - op_domain_stats.percentual_operational_area),
                   Catch::Matchers::WithinAbs(0.0, fiction::physical_constants::POP_STABILITY_ERR));
    }
}
