//
// Created by Jan Drewniok on 20.01.23.
//

#ifndef FICTION_TRANSFORM_TO_SIGN_HPP
#define FICTION_TRANSFORM_TO_SIGN_HPP

#include "fiction/algorithms/path_finding/distance.hpp"
#include "fiction/algorithms/simulation_sidb/simulation_parameters.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/technology/sign_to_label.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"

namespace fiction
{

/**
 * Converts the charge state into int (-1,0,1)
 *
 * @param cs sidb charge state.
 * @return int representing the SiDB's charge state.
 */
[[nodiscard]] int transform_to_sign(const sidb_charge_state& cs) noexcept
{
    if (cs == sidb_charge_state::POSITIVE)
    {
        return +1;
    }
    if (cs == sidb_charge_state::NEGATIVE)
    {
        return -1;
    }

    if (cs == sidb_charge_state::NEUTRAL)
    {
        return 0;
    }
    return 0;
}

} // namespace fiction

#endif  // FICTION_TRANSFORM_TO_SIGN_HPP
