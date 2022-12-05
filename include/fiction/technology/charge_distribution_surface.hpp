//
// Created by Jan Drewniok on 23.11.22.
//

#ifndef FICTION_EXGS_HPP
#define FICTION_EXGS_HPP


#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/sidb_defects.hpp"
#include "fiction/traits.hpp"

#include <unordered_map>

namespace fiction
{

enum class sidb_charge_states
{
    NONE,
    POSITIVE,
    NEUTRAL,
    NEGATIVE
};

struct sidb_charge
{
    /**
     * Standard constructor.
     */
    constexpr explicit sidb_charge(const sidb_charge_states defect_type = sidb_charge_states::NONE) noexcept :
            charge_state{defect_type}
    {}
    /**
     * Type of SiDB charge.
     */
    const sidb_charge_states charge_state;



};

[[nodiscard]] static constexpr bool is_positive(const sidb_charge &si_ch) noexcept
{
    return si_ch.charge_state == sidb_charge_states::POSITIVE;
}

[[nodiscard]] static constexpr bool is_negative(const sidb_charge &si_ch) noexcept
{
    return si_ch.charge_state == sidb_charge_states::NEGATIVE;
}

[[nodiscard]] static constexpr bool is_neutral(const sidb_charge &si_ch) noexcept
{
    return si_ch.charge_state == sidb_charge_states::NEUTRAL;
}

[[nodiscard]] static constexpr bool is_none(const sidb_charge &si_ch) noexcept
{
    return si_ch.charge_state == sidb_charge_states::NONE;
}


template <typename Lyt>
class charge_distribution_surface: public Lyt
{
  public:

    struct charge_distribution_storage
    {
        //explicit charge_distribution_storage() : charge_coordinates{} {}
        std::unordered_map<coordinate<Lyt>, sidb_charge> charge_coordinates{};
    };

    using storage = std::shared_ptr<charge_distribution_storage>;

    /**
     * Standard constructor for empty layouts.
     */
    explicit charge_distribution_surface() : Lyt(), strg{std::make_shared<charge_distribution_storage>()}
    {

        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(std::is_same_v<technology<Lyt>, sidb_technology>, "Lyt is not an SiDB layout");

    }

    explicit charge_distribution_surface(const Lyt& lyt) : Lyt(lyt), strg{std::make_shared<charge_distribution_storage>()}
    {
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(std::is_same_v<technology<Lyt>, sidb_technology>, "Lyt is not an SiDB layout");
    };
    /**
     * Assigns a given charge state to the given coordinate.
     *
     * @param c Coordinate to assign charge state cs.
     * @param d Defect to assign to coordinate c.
     */
    void assign_charge_state(const coordinate<Lyt>& c, const sidb_charge& cs) noexcept
    {
        if (!Lyt::is_empty_cell(c))
        {
            strg->charge_coordinates.insert({c, cs});
        }
    }
    /**
     * Returns the given coordinate's assigned charge state.
     *
     * @param c Coordinate to check.
     * @return charge state previously assigned to c or NONE if no defect was yet assigned.
     */
    [[nodiscard]] sidb_charge get_chargestate(const coordinate<Lyt>& c) const noexcept
    {
        if (const auto it = strg->charge_coordinates.find(c); it != strg->charge_coordinates.cend())
        {
            return it->second;
        }
        return sidb_charge{sidb_charge_states::NONE};
    }



  private:
    storage strg;
};

template <class T>
charge_distribution_surface(const T&) -> charge_distribution_surface<T>;

}  // namespace fiction

#endif  // FICTION_EXGS_HPP
