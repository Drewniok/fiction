//
// Created by Jan Drewniok on 08.12.22.
//

#ifndef FICTION_SIDB_CHARGE_STATE_HPP
#define FICTION_SIDB_CHARGE_STATE_HPP

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
 * Convert the charge state into int (-1,0,1)
 *
 * @param cs sidb charge state
 * @return int representing the SiDB's charge state
 */
[[nodiscard]] int transform_to_sign(const sidb_charge_state& cs) noexcept
{
    if (cs == sidb_charge_state::POSITIVE)
    {
        return +1;
    }
    else if (cs == sidb_charge_state::NEGATIVE)
    {
        return -1;
    }

    else if (cs == sidb_charge_state::NEUTRAL)
    {
        return 0;
    }

    else
    {
        return 0;
    }

}


[[nodiscard]] sidb_charge_state sign_to_label(const int& sg) noexcept
{
    if (sg == 1)
    {
        return sidb_charge_state::POSITIVE;
    }
    else if (sg == -1)
    {
        return sidb_charge_state::NEGATIVE;
    }

    else if (sg == 0)
    {
        return sidb_charge_state::NEUTRAL;
    }

    else
    {
        return sidb_charge_state::NONE;
    }

}


#endif  // FICTION_SIDB_CHARGE_STATE_HPP