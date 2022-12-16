//
// Created by Jan Drewniok on 15.12.22.
//

#ifndef FICTION_CHARGECONF_TO_INDEX_HPP
#define FICTION_CHARGECONF_TO_INDEX_HPP

#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/technology/sidb_charge_state.hpp>
#include <fiction/traits.hpp>
#include <cassert>

namespace fiction
{
/**
 * The charge distribution of the charge distribution surface is converted to a unique index. It is used to map every possible charge distribution in a SiDB layout.
 *
 * @tparam Lyt cell-level layout.
 * @tparam base number of charge states per SiDB. Base = 2 when only neutrally and negatively charged SiDBs are taken into account, otherwise base = 3.
 * @return a pair with the calculated charge distribution's corresponding charge index and the chosen base as the first and second element, respectively.
 */
template <typename Lyt>
std::pair<uint64_t, uint8_t> chargeconf_to_index(charge_distribution_surface<Lyt>& lyt, const uint8_t& base)
{
    assert(base == 2 || base == 3 && "base must be 2 or 3");
    // TODO check if in uint64
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is cell-level layout");
    uint64_t chargeindex   = 0;
    lyt.foreach_charge_state(
        [&chargeindex, &base, &lyt, i = 0u](const auto& cs) mutable
        {
            chargeindex += static_cast<unsigned int>(
                (transform_to_sign(cs.second) + 1) *
                std::pow(base, lyt.num_charges() - i - 1));
            i++;
        });
    return std::make_pair(chargeindex, base);
}
/**
 *  The unique index is converted to the charge distribution of the charge distribution surface.
 *
 * @tparam Lyt cell-level layout.
 * @param cp a pair of the charge index and the chosen base number.
 * @return a pair with the calculated charge distribution's corresponding charge index and the chosen base as the first and second element, respectively.
 */
template <typename Lyt>
void index_to_chargeconf(charge_distribution_surface<Lyt>& lyt, const std::pair<uint64_t, uint8_t>& cp)
{
    uint64_t chargeindex = 0;

    auto charge_quot = cp.first;

    while (charge_quot > 0)
    {
        auto  num_charges = lyt.num_charges() - 1;
        div_t d;
        d                = div(static_cast<int>(charge_quot), cp.second);
        charge_quot      = static_cast<uint64_t>(d.quot);
        lyt.foreach_charge_state(
            [&chargeindex, &lyt, i = 0u, &d, &num_charges](const auto& cs) mutable
            {
                if (i == num_charges)
                {
                    lyt.assign_charge_state(cs.first,sign_to_label(d.rem - 1));
                }
                i++;
            });
        num_charges -= 1;
    }
}

}  // namespace fiction

#endif  // FICTION_CHARGECONF_TO_INDEX_HPP
