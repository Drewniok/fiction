//
// Created by Jan Drewniok on 15.12.22.
//

#ifndef FICTION_CHARGECONF_TO_INDEX_HPP
#define FICTION_CHARGECONF_TO_INDEX_HPP

#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/technology/sidb_charge_state.hpp>
#include <fiction/traits.hpp>

namespace fiction
{

template <typename Lyt>
std::pair<uint64_t, uint8_t> chargeconf_to_index(charge_distribution_surface<Lyt>& lyt, const uint8_t& base)
{
    uint64_t chargeindex   = 0;
    uint64_t counter_start = 0;
    uint64_t counter       = 0;
    lyt.foreach_charge_state(
        [&chargeindex, &base, &lyt, &counter, &counter_start](const auto& cs)
        {
            chargeindex += static_cast<unsigned int>(
                (transform_to_sign(cs.second) + 1) *
                std::pow(base, lyt.num_charges() - static_cast<int>(counter - counter_start) - 1));
            counter += 1;
        });
    return std::make_pair(chargeindex, base);
}

template <typename Lyt>
void index_to_chargeconf(charge_distribution_surface<Lyt>& lyt, const std::pair<uint64_t, uint8_t>& cp)
{
    uint64_t chargeindex = 0;

    auto charge_quot = cp.first;

    while (charge_quot > 0)
    {
        auto  num_charges = lyt.num_charges() - 1;
        div_t d;  // Structure to represent the result value of an integral division performed by function div.
        d                = div(static_cast<int>(charge_quot), cp.second);
        charge_quot      = static_cast<uint64_t>(d.quot);
        uint64_t counter = 0;
        lyt.foreach_charge_state(
            [&chargeindex, &lyt, &counter, &d, &num_charges](const auto& cs)
            {
                if (counter == num_charges)
                {
                    lyt.assign_charge_state(cs.first,sign_to_label(d.rem - 1));
                }
                counter += 1;
            });
        num_charges -= 1;
    }
}

}  // namespace fiction

#endif  // FICTION_CHARGECONF_TO_INDEX_HPP
