//
// Created by Jan Drewniok on 05.12.22.
//

#include <catch2/catch_template_test_macros.hpp>

#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/layouts/hexagonal_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/technology/sidb_surface.hpp>

using namespace fiction;

TEMPLATE_TEST_CASE(
    "charge distribution surface traits and construction (layer on cell_level_layout)", "[charge-distribution-surface]",
    (sidb_surface<cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>>),
    (sidb_surface<cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_row_hex>>>>),
    (sidb_surface<cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_row_hex>>>>),
    (sidb_surface<
        cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_column_hex>>>>),
    (sidb_surface<
        cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_column_hex>>>>))
{
    REQUIRE(is_cell_level_layout_v<TestType>);
    CHECK(!has_assign_charge_state_v<TestType>);
    CHECK(!has_get_charge_state_v<TestType>);
    CHECK(!has_foreach_charge_state_v<TestType>);

    TestType lyt{};

    using charge_layout = charge_distribution_surface<TestType>;
    CHECK(is_cell_level_layout_v<charge_layout>);

    const charge_layout defect_lyt{};
    const charge_layout charge_lyt_from_lyt{lyt};

    using charge_charge_layout = sidb_surface<charge_layout>;
    CHECK(is_cell_level_layout_v<charge_charge_layout>);

}

TEMPLATE_TEST_CASE(
    "charge distribution surface traits (layer on cell_level_layout)", "[charge-distribution-surface]",
    (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_column_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_column_hex>>>))
{
    REQUIRE(is_cell_level_layout_v<TestType>);
    CHECK(!has_assign_charge_state_v<TestType>);
    CHECK(!has_get_charge_state_v<TestType>);
    CHECK(!has_foreach_charge_state_v<TestType>);
    CHECK(!has_assign_sidb_defect_v<TestType>);
    CHECK(!has_get_sidb_defect_v<TestType>);
    CHECK(!has_foreach_sidb_defect_v<TestType>);

    TestType lyt{};
    using charge_layout = charge_distribution_surface<TestType>;
    CHECK(is_cell_level_layout_v<charge_layout>);
    CHECK(!has_assign_sidb_defect_v<TestType>);
    CHECK(!has_get_sidb_defect_v<TestType>);
    CHECK(!has_foreach_sidb_defect_v<TestType>);
}

TEMPLATE_TEST_CASE(
    "assign and delete charge states", "[charge-distribution-surface]",
    (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_column_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_column_hex>>>),

    (sidb_surface<cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>>),
    (sidb_surface<cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_row_hex>>>>),
    (sidb_surface<cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_row_hex>>>>),
    (sidb_surface<
        cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, odd_column_hex>>>>),
    (sidb_surface<
        cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<siqad::coord_t, even_column_hex>>>>))

{
    TestType                    lyt{{11, 11}};
    charge_distribution_surface charge_layout{lyt, physical_params{}};

    SECTION("assign and read out charge states")
    {
        // assign SiDBs and charge states to three different cells and read the charge state
        charge_layout.assign_cell_type({5, 4}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({5, 5}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({5, 6}, TestType::cell_type::NORMAL);
        charge_layout.initialize();
        charge_layout.assign_charge_state_cell({5,4}, sidb_charge_state::POSITIVE);
        charge_layout.assign_charge_state_cell({5, 5}, sidb_charge_state::NEUTRAL);
        charge_layout.assign_charge_state_cell({5, 6}, sidb_charge_state::NEGATIVE);
        CHECK(charge_layout.get_charge_state_cell({5, 4}) == sidb_charge_state::POSITIVE);
        CHECK(charge_layout.get_charge_state_cell({5, 5}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_layout.get_charge_state_cell({5, 6}) == sidb_charge_state::NEGATIVE);

        CHECK(charge_layout.get_charge_state_cell({7, 6}) == sidb_charge_state::NONE);

        charge_layout.assign_charge_state_index(0, sidb_charge_state::NEUTRAL);
        charge_layout.assign_charge_state_index(1, sidb_charge_state::POSITIVE);
        charge_layout.assign_charge_state_index(2, sidb_charge_state::POSITIVE);
        CHECK(charge_layout.get_charge_state_index(0) == sidb_charge_state::NEUTRAL);
        CHECK(charge_layout.get_charge_state_index(1) == sidb_charge_state::POSITIVE);
        CHECK(charge_layout.get_charge_state_index(2) == sidb_charge_state::POSITIVE);

        // check if charge state for empty cell is 'NONE'
       CHECK(charge_layout.get_charge_state_index(3) == sidb_charge_state::NONE);
    }

    SECTION("assign different charge states to the same cell")
    {
        // check if the charge state can be overwritten
        charge_layout.assign_cell_type({5, 1}, TestType::cell_type::NORMAL);
        charge_layout.initialize();
        charge_layout.assign_charge_state_cell({5, 1}, sidb_charge_state::NEUTRAL);
        charge_layout.assign_charge_state_cell({5, 1}, sidb_charge_state::NONE);
        charge_layout.assign_charge_state_cell({5, 1}, sidb_charge_state::POSITIVE);
        charge_layout.assign_charge_state_cell({5, 1}, sidb_charge_state::NEGATIVE);
        CHECK(charge_layout.get_charge_state_cell({5, 1}) == sidb_charge_state::NEGATIVE);
    }


    SECTION("overwrite an assigned charge state")
    {
        // assign SiDBs and charge states to three different cells
        charge_layout.assign_cell_type({5, 4}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({5, 5}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({5, 6}, TestType::cell_type::NORMAL);
        charge_layout.initialize();
        charge_layout.assign_charge_state_cell({5, 4}, sidb_charge_state::POSITIVE);
        charge_layout.assign_charge_state_cell({5, 5}, sidb_charge_state::NEUTRAL);
        charge_layout.assign_charge_state_cell({5, 6}, sidb_charge_state::NEGATIVE);
//
        // all SiDBs' charge states are set to positive
        charge_layout.set_charge_states(sidb_charge_state::POSITIVE);
//
//        // read SiDBs' charge states
        CHECK(charge_layout.get_charge_state_cell({5, 4}) == sidb_charge_state::POSITIVE);
        CHECK(charge_layout.get_charge_state_cell({5, 5}) == sidb_charge_state::POSITIVE);
        CHECK(charge_layout.get_charge_state_cell({5, 6}) == sidb_charge_state::POSITIVE);
        CHECK(charge_layout.get_charge_state_cell({5, 1}) == sidb_charge_state::NONE);
//
        charge_layout.set_charge_states(sidb_charge_state::POSITIVE);

        // all SiDBs' charge states are set to neutral
        charge_layout.set_charge_states(sidb_charge_state::NEUTRAL);
//
//        // read SiDBs' charge states
        CHECK(charge_layout.get_charge_state_cell({5, 4}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_layout.get_charge_state_cell({5, 5}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_layout.get_charge_state_cell({5, 6}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_layout.get_charge_state_cell({5, 1}) == sidb_charge_state::NONE);
//
        charge_layout.set_charge_states(sidb_charge_state::NEUTRAL);
//
//        // all SiDBs' charge states are set to negative
        charge_layout.set_charge_states(sidb_charge_state::NEGATIVE);
//
//        // read SiDBs' charge states
        CHECK(charge_layout.get_charge_state_cell({5, 4}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_layout.get_charge_state_cell({5, 5}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_layout.get_charge_state_cell({5, 6}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_layout.get_charge_state_cell({5, 1}) == sidb_charge_state::NONE);
//
        charge_layout.set_charge_states(sidb_charge_state::NEGATIVE);
     }

    SECTION("overwrite the charge state")
    {
        // assign SiDBs and charge states to three different cells
        charge_layout.assign_cell_type({5, 4}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({5, 5}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({5, 6}, TestType::cell_type::NORMAL);
        charge_layout.initialize();
        charge_layout.assign_charge_state_cell({5, 4}, sidb_charge_state::POSITIVE);
        charge_layout.assign_charge_state_cell({5, 5}, sidb_charge_state::NEUTRAL);
        charge_layout.assign_charge_state_cell({5, 6}, sidb_charge_state::NEGATIVE);

        // remove previously assigned charge states
        charge_layout.assign_charge_state_cell({5, 4}, sidb_charge_state::NONE);
        charge_layout.assign_charge_state_cell({5, 5}, sidb_charge_state::NONE);
        charge_layout.assign_charge_state_cell({5, 6}, sidb_charge_state::POSITIVE);
        CHECK(charge_layout.get_charge_state_cell({5, 4}) == sidb_charge_state::NONE);
        CHECK(charge_layout.get_charge_state_cell({5, 5}) == sidb_charge_state::NONE);
        CHECK(charge_layout.get_charge_state_cell({5, 6}) == sidb_charge_state::POSITIVE);
    }

    SECTION("Distance matrix")
    {

        charge_layout.assign_cell_type({0, 0, 0}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({1, 0, 0}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({1, 1, 1}, TestType::cell_type::NORMAL);
        charge_layout.initialize();

        CHECK(charge_layout.get_distance_cell({0, 0, 0}, {0, 0, 0}) == 0.0);
        CHECK(charge_layout.get_distance_cell({0, 0, 0}, {1, 0, 0}) == physical_params{}.lat_a);
        CHECK(charge_layout.get_distance_cell({1, 0, 0}, {0, 0, 0}) == physical_params{}.lat_a);
        CHECK(charge_layout.get_distance_cell({1, 0, 0}, {1, 0, 0}) == 0.0);
        CHECK(charge_layout.get_distance_cell({0, 0, 0}, {1, 1, 1}) ==
              std::hypot(physical_params{}.lat_a, physical_params{}.lat_b + physical_params{}.lat_c));
        CHECK(charge_layout.get_distance_cell({1, 1, 1}, {1, 1, 1}) == 0);
    }
//
    SECTION("Potential matrix")
    {
        charge_layout.assign_cell_type({0, 0, 0}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({1, 8, 0}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({1, 10, 1}, TestType::cell_type::NORMAL);
        charge_layout.initialize();

        CHECK(charge_layout.get_potential_cell({0, 0, 0}, {0, 0, 0}) == 0.0);
        CHECK(charge_layout.get_potential_cell({1, 8, 0}, {1, 8, 0}) == 0.0);
        CHECK(charge_layout.get_potential_cell({1, 10, 1}, {1, 10, 1}) == 0.0);
        CHECK((charge_layout.get_potential_cell({1, 8, 1}, {0, 0, 0}) - 0.0121934043) < 0.00000001);
        CHECK(charge_layout.get_potential_cell({0, 0, 0}, {1, 10, 1}) == charge_layout.get_potential_cell({1, 10, 1}, {0, 0, 0}));
        CHECK(charge_layout.get_potential_cell({0, 0, 0}, {1, 8, 0}) > charge_layout.get_potential_cell({1, 10, 1}, {0, 0, 0}));
    }
//
    SECTION("Local Potential")
    {
        charge_layout.assign_cell_type({0, 0, 0}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({1, 8, 0}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({1, 10, 1}, TestType::cell_type::NORMAL);
        charge_layout.initialize();
        charge_layout.assign_charge_state_cell({0, 0, 0}, sidb_charge_state::POSITIVE);
        charge_layout.assign_charge_state_cell({1, 8, 0}, sidb_charge_state::POSITIVE);
        charge_layout.assign_charge_state_cell({1, 10, 1}, sidb_charge_state::POSITIVE);

        charge_layout.local_potential();

        charge_layout.foreach_cell([&charge_layout](const auto& c)
                                           { CHECK(*charge_layout.get_loc_pot_cell(c) > 0.0); });

        charge_layout.assign_charge_state_cell({0, 0, 0}, sidb_charge_state::NEGATIVE);
        charge_layout.assign_charge_state_cell({1, 8, 0}, sidb_charge_state::NEGATIVE);
        charge_layout.assign_charge_state_cell({1, 10, 1}, sidb_charge_state::NEGATIVE);

        charge_layout.local_potential();

        charge_layout.foreach_cell([&charge_layout](const auto& c)
                                           { CHECK(*charge_layout.get_loc_pot_cell(c) < 0.0); });

        charge_layout.set_charge_states(sidb_charge_state::NEUTRAL);

        charge_layout.local_potential();

        charge_layout.foreach_cell([&charge_layout](const auto& c)
                                           { CHECK(*charge_layout.get_loc_pot_cell(c) == 0); });
    }
//
    SECTION("Electrostatic potential energy of the charge configuration")
    {
        charge_layout.assign_cell_type({0, 0, 0}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({1, 0, 0}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({1, 1, 1}, TestType::cell_type::NORMAL);
        charge_layout.initialize();
        charge_layout.assign_charge_state_cell({0, 0, 0}, sidb_charge_state::NEGATIVE);
        charge_layout.assign_charge_state_cell({1, 0, 0}, sidb_charge_state::NEGATIVE);
        charge_layout.assign_charge_state_cell({1, 1, 1}, sidb_charge_state::NEGATIVE);

        // system energy is zero when all SiDBs are positively charged.
        charge_layout.local_potential();
        charge_layout.system_energy();
        CHECK(charge_layout.get_system_energy() > 0);

        // system energy is zero when all SiDBs are neutrally charged.
        charge_layout.set_charge_states(sidb_charge_state::NEUTRAL);
        charge_layout.local_potential();
        charge_layout.system_energy();
        CHECK(charge_layout.get_system_energy() == 0);

        // system energy is zero when all SiDBs are positively charged.
        charge_layout.set_charge_states(sidb_charge_state::POSITIVE);
        charge_layout.local_potential();
        charge_layout.system_energy();
        CHECK(charge_layout.get_system_energy() > 0);
    }
//
    SECTION("Physical validity check, far distance of SIDBs, all NEGATIVE")
    {

        TestType                    layout{{11, 11}};
        layout.assign_cell_type({0, 0, 0}, TestType::cell_type::NORMAL);
        layout.assign_cell_type({0, 2, 0}, TestType::cell_type::NORMAL);
        layout.assign_cell_type({4, 1, 1}, TestType::cell_type::NORMAL);

        charge_distribution_surface charge_layout_five{layout, physical_params{}};
        CHECK(charge_layout_five.get_charge_state_cell({0, 0,0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_layout_five.get_charge_state_cell({0, 2,0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_layout_five.get_charge_state_cell({4, 1,1}) == sidb_charge_state::NEGATIVE);


        charge_layout_five.local_potential();
        charge_layout_five.system_energy();
        charge_layout_five.validity_check();
        CHECK(charge_layout_five.get_charge_index().first == 0);

        charge_layout_five.assign_charge_state_cell({4, 1, 1}, sidb_charge_state::POSITIVE);
        CHECK(charge_layout_five.get_charge_state_cell({4, 1,1}) == sidb_charge_state::POSITIVE);
        charge_layout_five.chargeconf_to_index();
        CHECK(charge_layout_five.get_charge_index().first == 2);

        charge_layout_five.increase_charge_index();
        CHECK(charge_layout_five.get_charge_index().first == 3);
        charge_layout_five.increase_charge_index();
        CHECK(charge_layout_five.get_charge_index().first == 4);
        charge_layout_five.increase_charge_index();
        CHECK(charge_layout_five.get_charge_index().first == 5);
        charge_layout_five.increase_charge_index();
        CHECK(charge_layout_five.get_charge_index().first == 6);
        charge_layout_five.increase_charge_index();
        CHECK(charge_layout_five.get_charge_index().first == 7);
        charge_layout_five.increase_charge_index();
        CHECK(charge_layout_five.get_charge_index().first == 8);
        charge_layout_five.increase_charge_index();
        CHECK(charge_layout_five.get_charge_index().first == 9);
        charge_layout_five.increase_charge_index();
        CHECK(charge_layout_five.get_charge_index().first == 10);

        CHECK(charge_layout_five.get_charge_state_cell({0, 0,0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_layout_five.get_charge_state_cell({0, 2,0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_layout_five.get_charge_state_cell({4, 1,1}) == sidb_charge_state::NEUTRAL);
    }
//
    SECTION("Physical validity check, small distance, not all can be negatively charged anymore")
    {
        charge_layout.assign_cell_type({1, 0, 0}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({0, 2, 0}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({0, 2, 1}, TestType::cell_type::NORMAL);
        charge_layout.assign_charge_state_cell({1, 0, 0}, sidb_charge_state::NEGATIVE);
        charge_layout.assign_charge_state_cell({0, 2, 0}, sidb_charge_state::NEGATIVE);
        charge_layout.assign_charge_state_cell({0, 2, 1}, sidb_charge_state::NEGATIVE);
        charge_layout.initialize();

        // closely arranged SiDBs cannot be all negatively charged
        charge_layout.local_potential();
        charge_layout.system_energy();
        charge_layout.validity_check();
        CHECK(charge_layout.get_validity() == 0);

        charge_layout.assign_cell_type({0, 0, 1}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({1, 3, 0}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({10, 5, 1}, TestType::cell_type::NORMAL);
        charge_layout.assign_charge_state_cell({0, 0, 1}, sidb_charge_state::NEGATIVE);
        charge_layout.assign_charge_state_cell({1, 3, 0}, sidb_charge_state::NEGATIVE);
        charge_layout.assign_charge_state_cell({10, 5, 1}, sidb_charge_state::NEGATIVE);

        charge_layout.chargeconf_to_index();
        CHECK(charge_layout.get_charge_index().first == 0);
    }

    SECTION("increase charge index")
    {
        TestType                lyt_new{{11, 11}};
        const physical_params params{3, 5.0 * 1E-9, -0.32, 3.84 * 1E-10, 7.68 * 1E-10, 2.25 * 1E-10};

        lyt_new.assign_cell_type({0, 0, 1}, TestType::cell_type::NORMAL);
        lyt_new.assign_cell_type({1, 3, 0}, TestType::cell_type::NORMAL);
        lyt_new.assign_cell_type({10, 5, 1}, TestType::cell_type::NORMAL);

        charge_distribution_surface charge_layout_new{lyt_new, params};
        CHECK(charge_layout_new.get_charge_index().first == 0);

        charge_layout_new.set_charge_states(sidb_charge_state::POSITIVE);
        CHECK(charge_layout_new.get_charge_state_cell({0,0,1}) == sidb_charge_state::POSITIVE);
        CHECK(charge_layout_new.get_charge_state_cell({1,3,0}) == sidb_charge_state::POSITIVE);
        CHECK(charge_layout_new.get_charge_state_cell({10,5,1}) == sidb_charge_state::POSITIVE);

        charge_layout_new.chargeconf_to_index();
        CHECK(charge_layout_new.get_charge_index().first == 26);

        charge_layout_new.set_charge_states(sidb_charge_state::NEUTRAL);

        charge_layout_new.chargeconf_to_index();
        CHECK(charge_layout_new.get_charge_index().first == 13);

        charge_layout_new.increase_charge_index();
        charge_layout_new.chargeconf_to_index();
        CHECK(charge_layout_new.get_charge_index().first == 14);
        charge_layout_new.increase_charge_index();
        CHECK(charge_layout_new.get_charge_index().first == 15);
    }
}
