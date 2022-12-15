//
// Created by Jan Drewniok on 12.12.22.
//

#include <catch2/catch_template_test_macros.hpp>

#include <fiction/algorithms/simulation_sidb/distance_matrix.hpp>
#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/layouts/hexagonal_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/algorithms/simulation_sidb/chargeconf_to_index.hpp>


using namespace fiction;

TEMPLATE_TEST_CASE(
    "Distance Matrix calculation", "[distance-matrix]",
    (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_column_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_column_hex>>>))
{
    TestType lyt{{10,10}};

    charge_distribution_surface charge_layout{lyt};

    charge_layout.assign_cell_type({0, 0, 1}, TestType::cell_type::NORMAL);
    charge_layout.assign_cell_type({1, 3, 0}, TestType::cell_type::NORMAL);
    charge_layout.assign_cell_type({10, 5, 1}, TestType::cell_type::NORMAL);
    charge_layout.assign_charge_state({0, 0, 1}, sidb_charge_state::NEGATIVE);
    charge_layout.assign_charge_state({1, 3, 0}, sidb_charge_state::NEGATIVE);
    charge_layout.assign_charge_state({10, 5, 1}, sidb_charge_state::NEGATIVE);

    CHECK(chargeconf_to_index(charge_layout, 3)== 0);
    CHECK(chargeconf_to_index(charge_layout, 2)== 0);

    charge_layout.foreach_charge_state([&charge_layout](const auto& cd)
                                       { charge_layout.assign_charge_state(cd.first, sidb_charge_state::POSITIVE); });

    CHECK(chargeconf_to_index(charge_layout, 3)== 26);

    charge_layout.foreach_charge_state([&charge_layout](const auto& cd)
                                       { charge_layout.assign_charge_state(cd.first, sidb_charge_state::NEUTRAL); });

    CHECK(chargeconf_to_index(charge_layout, 3)== 13);
    CHECK(chargeconf_to_index(charge_layout, 2)== 7);



} // namespace fiction
