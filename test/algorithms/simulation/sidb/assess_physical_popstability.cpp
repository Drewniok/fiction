//
// Created by Jan Drewniok on 02.11.23.
//

#include <catch2/catch_template_test_macros.hpp>

#include <fiction/algorithms/simulation/sidb/assess_physical_popstability.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/types.hpp>

using namespace fiction;

using layout = sidb_cell_clk_lyt_siqad;

TEST_CASE("Single SiDB", "[assess-physical-popstability]")
{
    layout lyt{};

    lyt.assign_cell_type({1, 1, 0}, sidb_technology::cell_type::NORMAL);

    const auto result = assess_physical_popstability(lyt, sidb_simulation_parameters{2, -0.32});
    REQUIRE(result.size() == 1);
    const auto [cell, status, potential_difference] = result.at(0);
    CHECK(cell == siqad::coord_t{1, 1, 0});
    CHECK(status == transition_type::NEGATIVE_TO_NEUTRAL);
    CHECK(potential_difference == 0.32);
}

TEST_CASE("Three SiDBs with positive charge states", "[assess-physical-popstability]")
{
    layout lyt{};

    lyt.assign_cell_type({1, 1, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({1, 1, 1}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::NORMAL);

    const auto result = assess_physical_popstability(lyt, sidb_simulation_parameters{2, -0.32});
    REQUIRE(result.size() == 3);

    SECTION("Ground state")
    {
        const auto [cell, status, potential_difference] = result.at(0);
        CHECK(cell == siqad::coord_t{2, 1, 0});
        CHECK(status == transition_type::NEGATIVE_TO_NEUTRAL);
        CHECK(potential_difference < 0.43);
    }
    SECTION("1st excited state")
    {
        const auto [cell, status, potential_difference] = result.at(1);
        CHECK(cell == siqad::coord_t{2, 1, 0});
        CHECK(status == transition_type::NEGATIVE_TO_NEUTRAL);
        CHECK(potential_difference < 0.23);
    }

    SECTION("2nd excited state")
    {
        const auto [cell, status, potential_difference] = result.at(2);
        CHECK(cell == siqad::coord_t{1, 1, 1});
        CHECK(status == transition_type::NEUTRAL_TO_NEGATIVE);
        CHECK(potential_difference < 0.21);
    }
}

TEST_CASE("Bestagon AND gate", "[assess-physical-popstability]")
{
    layout lyt{};

    lyt.assign_cell_type({36, 1, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({38, 0, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({23, 9, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({18, 11, 1}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({18, 9, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({19, 8, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({20, 14, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({19, 13, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({26, 16, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({24, 15, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({32, 2, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({30, 3, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({26, 4, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({24, 5, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({12, 4, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({14, 5, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({6, 2, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({8, 3, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({32, 18, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({30, 17, 0}, sidb_technology::cell_type::OUTPUT);

    lyt.assign_cell_type({36, 19, 0}, sidb_technology::cell_type::NORMAL);

    const auto result = assess_physical_popstability(lyt, sidb_simulation_parameters{2, -0.32});
    REQUIRE(result.size() == 8);
    const auto [cell, status, potential_difference] = result.at(0);
    CHECK(cell == siqad::coord_t{2, 1, 0});
    CHECK(status == transition_type::NEUTRAL_TO_NEGATIVE);
    CHECK(potential_difference < 0.021);
}

TEST_CASE("Bestagon CROSSING gate input 11", "[assess-physical-popstability]")
{
    layout lyt{};

    lyt.assign_cell_type({36, 1, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({6, 2, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({20, 12, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({8, 3, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({14, 5, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({14, 11, 1}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({12, 4, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({14, 15, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({26, 4, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({14, 9, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({24, 15, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({12, 16, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({18, 9, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({26, 16, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({24, 13, 1}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({24, 5, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({30, 3, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({16, 13, 1}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({32, 2, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({20, 8, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({30, 17, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({6, 18, 0}, sidb_technology::cell_type::OUTPUT);

    lyt.assign_cell_type({32, 18, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({8, 17, 0}, sidb_technology::cell_type::OUTPUT);

    lyt.assign_cell_type({2, 19, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({36, 19, 0}, sidb_technology::cell_type::NORMAL);

    CHECK(lyt.num_cells() == 27);

    const auto result = assess_physical_popstability(lyt, sidb_simulation_parameters{2, -0.32});
    REQUIRE(result.size() == 20);
    const auto [cell, status, potential_difference] = result.at(0);
    CHECK(cell == siqad::coord_t{14, 9, 0});
    CHECK(status == transition_type::NEUTRAL_TO_NEGATIVE);
    CHECK(potential_difference < 0.01);
}