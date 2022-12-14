//
// Created by Jan Drewniok on 13.12.22.
//
#include <catch2/catch_template_test_macros.hpp>
#include <fiction/algorithms/simulation_sidb/distance_matrix.hpp>
#include <fiction/algorithms/simulation_sidb/local_potential.hpp>
#include <fiction/algorithms/simulation_sidb/potential_matrix.hpp>
#include <fiction/algorithms/simulation_sidb/system_energy.hpp>
#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/layouts/hexagonal_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>

using namespace fiction;

TEMPLATE_TEST_CASE(
    "system energy calculation", "[system-energy]",
    (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_column_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_column_hex>>>))
{
    TestType                    lyt{{5, 5}};
    charge_distribution_surface charge_layout{lyt};

    charge_layout.assign_cell_type({0, 0, 0}, TestType::cell_type::NORMAL);
    charge_layout.assign_cell_type({1, 0, 0}, TestType::cell_type::NORMAL);
    charge_layout.assign_cell_type({1, 1, 1}, TestType::cell_type::NORMAL);
    charge_layout.assign_charge_state({0, 0, 0}, sidb_charge_state::POSITIVE);
    charge_layout.assign_charge_state({1, 0, 0}, sidb_charge_state::POSITIVE);
    charge_layout.assign_charge_state({1, 1, 1}, sidb_charge_state::POSITIVE);

    // system energy is zero when all SiDBs are positively charged.
    auto distance  = distance_sidbs(charge_layout);
    auto potential = potential_sidbs<charge_distribution_surface<TestType>>(distance);
    auto local_pot = local_potential<charge_distribution_surface<TestType>>(charge_layout, potential);
    auto energy = system_energy<TestType>(charge_layout, local_pot);
    CHECK(energy>0);

    // system energy is zero when all SiDBs are neutrally charged.
    charge_layout.foreach_charge_state([&charge_layout](const auto& cd)
                                       { charge_layout.assign_charge_state(cd.first, sidb_charge_state::NEUTRAL); });
    auto local_pot_neutral = local_potential<charge_distribution_surface<TestType>>(charge_layout, potential);
    auto energy_neutral = system_energy<TestType>(charge_layout, local_pot);
    CHECK(energy_neutral==0);

    // system energy is zero when all SiDBs are positively charged.
    charge_layout.foreach_charge_state([&charge_layout](const auto& cd)
                                       { charge_layout.assign_charge_state(cd.first, sidb_charge_state::POSITIVE); });
    auto local_pot_positive = local_potential<charge_distribution_surface<TestType>>(charge_layout, potential);
    auto energy_positive = system_energy<TestType>(charge_layout, local_pot);
    CHECK(energy_positive>0);

}