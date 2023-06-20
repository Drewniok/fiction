//
// Created by Jan Drewniok on 19.06.23.
//

#include <catch2/catch_template_test_macros.hpp>

#include <fiction/algorithms/simulation/sidb/bestagon_gate_generator.hpp>

using namespace fiction;

TEST_CASE("Test critical_temperature function", "[bestagon_gate_generator]")
{
    bestagon_gate_generator_params<sidb_cell_clk_lyt_siqad> params{
        std::pair(15, 7), random_layout_params<sidb_cell_clk_lyt_siqad>{{}, 3, true}, create_crossing_wire_tt(), 5000};
    auto gates = bestagon_gate_generator(params);
    // write_sqd_layout(gates[0], "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/skeleton/test.sqd");
    CHECK(gates.size() == 1);
}