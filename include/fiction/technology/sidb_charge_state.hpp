//
// Created by Jan Drewniok on 08.12.22.
//

#ifndef FICTION_SIDB_CHARGE_STATE_HPP
#define FICTION_SIDB_CHARGE_STATE_HPP

namespace fiction
{
/**
 * Possible SiDB charges.
 */
enum class sidb_charge_state
{
    NONE,  // assigned when layout cell is empty
    POSITIVE,
    NEUTRAL,
    NEGATIVE
};

/**
 * Converts the charge state into int (-1,0,1).
 *
 * @param cs sidb charge state.
 * @return int representing the SiDB's charge state.
 */
[[nodiscard]] static constexpr int transform_to_sign(const sidb_charge_state& cs) noexcept
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

/**
 * Converts the charge state (-1,0,1) into enum.
 *
 * @param sq charge state as integer (-1,0,1).
 * @return sidb_charge_state.
 */
[[nodiscard]] static constexpr sidb_charge_state sign_to_label(const int& sg) noexcept
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

#endif  // FICTION_SIDB_CHARGE_STATE_HPP
