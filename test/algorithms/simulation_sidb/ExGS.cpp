//
// Created by Jan Drewniok on 18.12.22.
//

#include "fiction/algorithms/simulation_sidb/ExGS.hpp"

#include <catch2/catch_template_test_macros.hpp>

#include "fiction/algorithms/simulation_sidb/new_approach.hpp"

#include <fiction/algorithms/simulation_sidb/TTS.hpp>
#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/layouts/hexagonal_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>

using namespace fiction;
using namespace detail;

TEMPLATE_TEST_CASE(
    "exhaustive ground state search", "[ExGS]",
    (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_column_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_column_hex>>>))
{
    TestType lyt{{20, 10}};

    lyt.assign_cell_type({1, 3, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({3, 3, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({4, 3, 0}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({6, 3, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({7, 3, 0}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({6, 10, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({7, 10, 0}, TestType::cell_type::NORMAL);

    charge_distribution_surface                                       charge_layout{lyt};
    std::pair<double, std::vector<charge_distribution_surface<TestType>>> output = exgs<TestType>(charge_layout);
    CHECK(!output.second.empty());

    const physical_params     params{2, 5.6, 5.0 * 1E-9, -0.32, 3.84 * 1E-10, 7.68 * 1E-10, 2.25 * 1E-10};

    charge_distribution_surface charge_layout_new{lyt, params};
    std::pair<double, std::vector<charge_distribution_surface<TestType>>> output_new = exgs<TestType>(charge_layout_new);
    CHECK(!output_new.second.empty());

    for (const auto& it : output_new.second)
    {
        it.foreach_cell([&it](const auto& c)
                                       { CHECK(it.get_charge_state_cell(c) != sidb_charge_state::POSITIVE); });
    }

}
