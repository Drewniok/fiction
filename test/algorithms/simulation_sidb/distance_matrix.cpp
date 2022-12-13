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

using namespace fiction;

TEMPLATE_TEST_CASE(
    "Distance Matrix calculation", "[distance-matrix]",
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

    auto distance = distance_SiDBs(charge_layout);
    CHECK(distance.at({{0, 0, 0}, {0, 0, 0}}) == 0.0);
    CHECK(distance.at({{0, 0, 0}, {1, 0, 0}}) == simulation_parameter{}.lat_a);
    CHECK(distance.at({{1, 0, 0}, {0, 0, 0}}) == simulation_parameter{}.lat_a);
    CHECK(distance.at({{1, 0, 0}, {1, 0, 0}}) == 0.0);
    CHECK(distance.at({{0, 0, 0}, {1, 1, 1}}) ==
          std::hypot(simulation_parameter{}.lat_a, simulation_parameter{}.lat_b + simulation_parameter{}.lat_c));
    CHECK(distance.at({{1, 1, 1}, {1, 1, 1}}) == 0);

} // namespace fiction
