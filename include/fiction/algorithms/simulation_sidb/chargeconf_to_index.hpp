//
// Created by Jan Drewniok on 15.12.22.
//

#ifndef FICTION_CHARGECONF_TO_INDEX_HPP
#define FICTION_CHARGECONF_TO_INDEX_HPP

#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/traits.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/technology/sidb_charge_state.hpp>

namespace fiction
{

template <typename Lyt>
uint64_t chargeconf_to_index(charge_distribution_surface<Lyt>& lyt, const uint8_t & base)
{
    uint64_t chargeindex = 0;
    uint64_t counter_start = 0;
    uint64_t counter = 0;
    lyt.foreach_charge_state([&chargeindex, &base, &lyt, &counter, &counter_start](const auto& cs) {
                                 chargeindex += static_cast<unsigned int>((transform_to_sign(cs.second) + 1)
                                                 * std::pow(base, lyt.assigned_charges() - static_cast<int>(counter-counter_start) - 1));
                             counter +=1;}
    );
    return chargeindex;
    }

}  // namespace fiction

#endif  // FICTION_CHARGECONF_TO_INDEX_HPP

