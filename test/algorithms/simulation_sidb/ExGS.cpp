//
// Created by Jan Drewniok on 18.12.22.
//

#include <catch2/catch_template_test_macros.hpp>
#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/layouts/hexagonal_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include "fiction/algorithms/simulation_sidb/ExGS.hpp"
#include "fiction/algorithms/simulation_sidb/new_approach.hpp"

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
    lyt.assign_cell_type({5, 0, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({1, 3, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({2, 3, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({10, 3, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({10, 5, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({5, 10, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({1, 9, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({2, 4, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({7, 10, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({5, 9, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({6, 4, 0}, TestType::cell_type::NORMAL);


    charge_distribution_surface                        charge_layout{lyt};
    std::vector<charge_distribution_surface<TestType>> output = metastable_layouts<TestType>(charge_layout);
    CHECK(output.size() > 0);

    const simulation_params     params{5.6, 5.0 * 1E-9, -0.32, 3.84 * 1E-10, 7.68 * 1E-10, 2.25 * 1E-10, 2};
    charge_distribution_surface charge_layout_new{lyt, params};
    std::vector<charge_distribution_surface<TestType>> output_new = metastable_layouts<TestType>(charge_layout_new);
    CHECK(output_new.size() > 0);
    for (auto& it : output_new)
    {
        it.foreach_charge_state([&it](const auto& c)
                                                { CHECK(it.get_charge_state(c.first) != sidb_charge_state::POSITIVE); });
    }

    //std::vector<charge_distribution_surface<TestType>> output_sim = Sim<TestType>(charge_layout_new);
    std::unordered_map<double, charge_distribution_surface<TestType>> output_sim = Sim<TestType>(charge_layout_new);
}
