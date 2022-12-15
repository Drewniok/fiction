//
// Created by Jan Drewniok on 13.12.22.
//

#include <catch2/catch_template_test_macros.hpp>

#include <fiction/algorithms/simulation_sidb/distance_matrix.hpp>
#include <fiction/algorithms/simulation_sidb/potential_matrix.hpp>
#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/layouts/hexagonal_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/technology/charge_distribution_surface.hpp>

using namespace fiction;

TEMPLATE_TEST_CASE(
    "Potential Matrix calculation", "[potential-matrix]",
    (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_column_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_column_hex>>>))
{
    TestType lyt{{3, 3}};

    charge_distribution_surface charge_layout{lyt};

    charge_layout.assign_cell_type({0, 0, 0}, TestType::cell_type::NORMAL);
    charge_layout.assign_cell_type({1, 0, 0}, TestType::cell_type::NORMAL);
    charge_layout.assign_cell_type({1, 1, 1}, TestType::cell_type::NORMAL);

    auto distance  = initialize_sidb_distance_matrix(charge_layout);
    auto potential = potential_sidbs<charge_distribution_surface<TestType>>(distance);

    CHECK(potential.size() == 9);
    CHECK(potential.at({{0, 0, 0}, {0, 0, 0}}) == 0.0);
    CHECK(potential.at({{1, 0, 0}, {1, 0, 0}}) == 0.0);
    CHECK(potential.at({{1, 1, 1}, {1, 1, 1}}) == 0.0);
    CHECK((potential.at({{1, 1, 1}, {0, 0, 0}}) - 0.19518121) < 0.00000001);
    CHECK(potential.at({{0, 0, 0}, {1, 1, 1}}) == potential.at({{1, 1, 1}, {0, 0, 0}}));
    CHECK(potential.at({{0, 0, 0}, {1, 1, 1}}) < potential.at({{1, 0, 0}, {0, 0, 0}}));
}
