//
// Created by Jan Drewniok on 23.11.22.
//

#ifndef FICTION_CHARGE_DISTRIBUTION_SURFACE_HPP
#define FICTION_CHARGE_DISTRIBUTION_SURFACE_HPP

#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/sidb_defects.hpp"
#include "fiction/traits.hpp"

#include <unordered_map>

namespace fiction
{

/**
 * Possible SiDB charges.
 */
enum class sidb_charge_states
{
    NONE,  // assigned when layout cell is empty
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

    const sidb_charge_states charge_state;
};

template <typename Lyt>
class charge_distribution_surface : public Lyt
{
  public:
    struct charge_distribution_storage
    {
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
    /**
     * Standard constructor for existing layouts.
     */
    explicit charge_distribution_surface(const Lyt& lyt) :
            Lyt(lyt),
            strg{std::make_shared<charge_distribution_storage>()}
    {
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(std::is_same_v<technology<Lyt>, sidb_technology>, "Lyt is not an SiDB layout");
    };
    /**
     * Assigns a given charge state to the given coordinate.
     *
     * @param c Coordinate to assign charge state cs.
     * @param cs Charge state to assign to coordinate c.
     */
    void assign_charge_state(const coordinate<Lyt>& c, const sidb_charge& cs) noexcept
    {
        if (!Lyt::is_empty_cell(c))
        {
            if (const auto it = strg->charge_coordinates.find(c); it != strg->charge_coordinates.cend())
            {
                strg->charge_coordinates.erase(c);
                strg->charge_coordinates.insert({c, cs});
            }
            strg->charge_coordinates.insert({c, cs});
        }
    }
    /**
     * Returns the given coordinate's assigned charge state.
     *
     * @param c Coordinate to check.
     * @return Charge state previously assigned to c or NONE if cell owns emtpy cell_type.
     */
    [[nodiscard]] sidb_charge get_charge_state(const coordinate<Lyt>& c) const noexcept
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

#endif  // CHARGE_DISTRIBUTION_SURFACE
