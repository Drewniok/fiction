//
// Created by Jan Drewniok on 13.12.22.
//
#include <catch2/catch_template_test_macros.hpp>
#include <fiction/algorithms/simulation_sidb/distance_matrix.hpp>
#include <fiction/algorithms/simulation_sidb/local_potential.hpp>
#include <fiction/algorithms/simulation_sidb/potential_matrix.hpp>
#include <fiction/algorithms/simulation_sidb/validity_check.hpp>
#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/layouts/hexagonal_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>

using namespace fiction;

TEMPLATE_TEST_CASE(
    "Physical validity check, far distance", "[system-energy]",
    (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_column_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_column_hex>>>))
{
    TestType                    lyt{{5, 5}};
    charge_distribution_surface charge_layout{lyt};

    charge_layout.assign_cell_type({0, 0, 0}, TestType::cell_type::NORMAL);
    charge_layout.assign_cell_type({0, 2, 0}, TestType::cell_type::NORMAL);
    charge_layout.assign_cell_type({4, 1, 1}, TestType::cell_type::NORMAL);
    charge_layout.assign_charge_state({0, 0, 0}, sidb_charge_state::NEGATIVE);
    charge_layout.assign_charge_state({0, 2, 0}, sidb_charge_state::NEGATIVE);
    charge_layout.assign_charge_state({4, 1, 1}, sidb_charge_state::NEGATIVE);

    auto distance  = distance_SiDBs(charge_layout);
    auto potential = potential_SiDBs<charge_distribution_surface<TestType>>(distance);
    auto local_pot = local_potential<charge_distribution_surface<TestType>>(charge_layout, potential);
    auto valid = validity_check(charge_layout,local_pot,potential);
    CHECK(valid == 1);
}

TEMPLATE_TEST_CASE(
    "Physical validity check, close distance", "[system-energy]",
    (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_column_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_column_hex>>>))
{
    TestType                    lyt{{5, 5}};
    charge_distribution_surface charge_layout{lyt};

    charge_layout.assign_cell_type({0, 0, 0}, TestType::cell_type::NORMAL);
    charge_layout.assign_cell_type({0, 2, 0}, TestType::cell_type::NORMAL);
    charge_layout.assign_cell_type({0, 3, 0}, TestType::cell_type::NORMAL);
    charge_layout.assign_charge_state({0, 0, 0}, sidb_charge_state::NEGATIVE);
    charge_layout.assign_charge_state({0, 2, 0}, sidb_charge_state::NEGATIVE);
    charge_layout.assign_charge_state({0, 3, 1}, sidb_charge_state::NEGATIVE);

    auto distance  = distance_SiDBs(charge_layout);
    auto potential = potential_SiDBs<charge_distribution_surface<TestType>>(distance);
    auto local_pot = local_potential<charge_distribution_surface<TestType>>(charge_layout, potential);
    auto valid = validity_check(charge_layout,local_pot,potential);
    CHECK(valid == 0);

}



