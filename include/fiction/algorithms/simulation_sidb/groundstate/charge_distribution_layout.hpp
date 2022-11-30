//
// Created by Jan Drewniok on 23.11.22.
//

#ifndef FICTION_EXGS_HPP
#define FICTION_EXGS_HPP

#include "fiction/algorithms/simulation_sidb/groundstate/constants.hpp"
#include "fiction/algorithms/simulation_sidb/groundstate/simulation_parameter.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/sidb_defects.hpp"
#include "fiction/traits.hpp"

#include <unordered_map>

namespace fiction
{

struct exgs_params
{
    constants            constants{};
    simulation_parameter physical_params{};
};

enum class sidb_chargestate
{
    positive = 1,
    neutral  = 0,
    negative = -1
};

template <typename Lyt>
class charge_distribution_surface<Lyt, false> : public Lyt
{
  public:
    struct charge_distribution_storage
    {
        explicit charge_distribution_storage();
        std::unordered_map<coordinate<Lyt>, sidb_chargestate> charge_coordinates{};
    };

    using storage = std::shared_ptr<charge_distribution_storage>;

    explicit charge_distribution_surface(const Lyt& lyt) : Lyt(lyt)
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
    void assign_charge_state(const coordinate<Lyt>& c, const sidb_chargestate& cs) noexcept
    {
        strg->defective_coordinates.insert({c, cs});
    }
    /**
     * Returns the given coordinate's assigned charge state.
     *
     * @param c Coordinate to check.
     * @return charge state previously assigned to c or NONE if no defect was yet assigned.
     */
    [[nodiscard]] sidb_chargestate get_sidb_chargestate(const coordinate<Lyt>& c) const noexcept
    {
        if (const auto  it = strg->charge_coordinates.find(c); it != strg->defective_coordinates.cend())
        {
            return it->second;
        }

    }

  private:
    storage strg;
};

template <class T>
charge_distribution_surface(const T&) -> charge_distribution_surface<T>;

namespace detail
{}
}  // namespace fiction

#endif  // FICTION_EXGS_HPP
