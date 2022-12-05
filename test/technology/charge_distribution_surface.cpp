//
// Created by Jan Drewniok on 05.12.22.
//
//
// Created by marcel on 07.03.22.
//

#include "catch.hpp"

#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/layouts/hexagonal_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/technology/sidb_surface.hpp>
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/traits.hpp>

#include <set>
#include <type_traits>

using namespace fiction;


TEMPLATE_TEST_CASE(
    "Charged and neutral defect extent", "[sidb-surface]",
    (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<offset::ucoord_t>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<offset::ucoord_t, odd_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<offset::ucoord_t, even_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<offset::ucoord_t, odd_column_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<offset::ucoord_t, even_column_hex>>>),

    (sidb_surface<cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<offset::ucoord_t>>>>),
    (sidb_surface<cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<offset::ucoord_t, odd_row_hex>>>>),
    (sidb_surface<cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<offset::ucoord_t, even_row_hex>>>>),
    (sidb_surface<cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<offset::ucoord_t, odd_column_hex>>>>),
    (sidb_surface<cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<offset::ucoord_t, even_column_hex>>>>))

{
    TestType lyt{{11, 9}};
    charge_distribution_surface                 charge_layout{lyt};

    SECTION("charged defects")
    {
        // assign defects
        charge_layout.assign_cell_type({5, 4}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({5, 5}, TestType::cell_type::NORMAL);
        charge_layout.assign_cell_type({5, 6}, TestType::cell_type::NORMAL);
        charge_layout.assign_charge_state({5, 4}, sidb_charge{sidb_charge_states::POSITIVE});
        charge_layout.assign_charge_state({5, 5}, sidb_charge{sidb_charge_states::NEUTRAL});
        charge_layout.assign_charge_state({5, 6}, sidb_charge{sidb_charge_states::NEGATIVE});
        CHECK(charge_layout.get_chargestate({5,4}).charge_state == sidb_charge_states::POSITIVE);
        CHECK(charge_layout.get_chargestate({5,5}).charge_state == sidb_charge_states::NEUTRAL);
        CHECK(charge_layout.get_chargestate({5,6}).charge_state == sidb_charge_states::NEGATIVE);

        CHECK(charge_layout.get_chargestate({5,7}).charge_state == sidb_charge_states::NONE);
    }
}

TEMPLATE_TEST_CASE(
    "Charged and neutral defect extent at layout edges", "[sidb-surface]",
    (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<offset::ucoord_t>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<offset::ucoord_t, odd_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<offset::ucoord_t, even_row_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<offset::ucoord_t, odd_column_hex>>>),
    (cell_level_layout<sidb_technology, clocked_layout<hexagonal_layout<offset::ucoord_t, even_column_hex>>>))
{
    TestType lyt{aspect_ratio<TestType>{11, 9}};

    sidb_surface<TestType> defect_layout{lyt};



}




