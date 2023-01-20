//
// Created by Jan Drewniok on 20.01.23.
//

#ifndef FICTION_SIGN_TO_LABEL_HPP
#define FICTION_SIGN_TO_LABEL_HPP

#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/traits.hpp"
#include <iostream>

namespace fiction
{

/**
 * Converts the charge state index (-1,0,1) into enum.
 *
 * @param int sidb charge state as index (-1,0,1).
 * @return enum representing the SiDB's charge state.
 */
[[nodiscard]] sidb_charge_state sign_to_label(const int& sg) noexcept
{
    if (sg == 1)
    {
        return sidb_charge_state::POSITIVE;
    }
    if (sg == -1)
    {
        return sidb_charge_state::NEGATIVE;
    }

    if (sg == 0)
    {
        return sidb_charge_state::NEUTRAL;
    }

    return sidb_charge_state::NONE;
}

} // namespace fiction

#endif  // FICTION_SIGN_TO_LABEL_HPP
