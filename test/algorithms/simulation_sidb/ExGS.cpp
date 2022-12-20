//
// Created by Jan Drewniok on 18.12.22.
//

#include <catch2/catch_template_test_macros.hpp>
#include <fiction/algorithms/simulation_sidb/distance_matrix.hpp>
#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/layouts/hexagonal_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include "fiction/algorithms/simulation_sidb/ExGS.hpp"

using namespace fiction;
using namespace detail;

TEMPLATE_TEST_CASE(
    "Local potential", "[local-potential]",
    (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_column_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_column_hex>>>))
{
    TestType                    lyt{{100, 22}};

    lyt.assign_cell_type({5, 0, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({1, 3, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({1, 5, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({1, 1, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({1, 6, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({1, 7, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({50, 0, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({10, 3, 0}, TestType::cell_type::NORMAL);


    charge_distribution_surface charge_layout{lyt};

    CHECK(charge_layout.get_charge_state({5, 0,1}) == sidb_charge_state::NEGATIVE);

   std::vector<std::pair<charge_distribution_surface<TestType>, double>> output = GS<TestType>(charge_layout);



}
