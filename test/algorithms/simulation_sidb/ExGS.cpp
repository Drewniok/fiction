//
// Created by Jan Drewniok on 18.12.22.
//

#include <catch2/catch_template_test_macros.hpp>
#include <fiction/algorithms/simulation_sidb/distance_matrix.hpp>
#include <fiction/algorithms/simulation_sidb/local_potential.hpp>
#include <fiction/algorithms/simulation_sidb/potential_matrix.hpp>
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
    TestType                    lyt{{10, 22}};
    charge_distribution_surface charge_layout{lyt};


    charge_layout.assign_cell_type({5, 0, 1}, TestType::cell_type::NORMAL);
    charge_layout.assign_cell_type({1, 3, 0}, TestType::cell_type::NORMAL);
    charge_layout.assign_cell_type({10, 5, 1}, TestType::cell_type::NORMAL);





//    charge_layout.assign_cell_type({12, 3, 0}, TestType::cell_type::NORMAL);

//    charge_layout.assign_charge_state({0, 0, 1}, sidb_charge_state::NEGATIVE);
//    charge_layout.assign_charge_state({1, 3, 0}, sidb_charge_state::NEUTRAL);
//    charge_layout.assign_charge_state({10, 5, 1}, sidb_charge_state::NEGATIVE);

    std::vector<std::pair<charge_distribution_surface<TestType>, double>> output = GS<TestType>(charge_layout);

//    for (auto &it: output)
//    {
//        CHECK(it.second == 1);
//     }
//    simulation.initialize_sidb_potential_matrix();
//    simulation.local_potential();
//    simulation.system_energy();
//    CHECK(simulation.get_system_energy() == -1);
//
//    simulation.chargeconf_to_index();
//    CHECK(simulation.get_charge_index().first == 0);



//    simulation.run();


}
