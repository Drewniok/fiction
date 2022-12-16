//
// Created by Jan Drewniok on 23.11.22.
//

#ifndef FICTION_CHARGE_DISTRIBUTION_SURFACE_HPP
#define FICTION_CHARGE_DISTRIBUTION_SURFACE_HPP

#include "fiction/algorithms/path_finding/distance.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/electrostatic_potential.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/traits.hpp"

#include <cassert>

namespace fiction
{

/**
 * A layout type to layer on top of any SiDB cell-level layout. It implements an interface to store and access
 * SiDBs' charge states.
 *
 * @tparam Lyt SiDB cell-level layout based in SiQAD-coordinates.
 * @tparam has_sidb_charge_distribution Automatically determines whether a charge distribution interface is already
 * present.
 */
template <typename Lyt, bool has_charge_distribution_interface =
                            std::conjunction_v<has_assign_charge_state<Lyt>, has_get_charge_state<Lyt>>>
class charge_distribution_surface : public Lyt
{};

template <typename Lyt>
class charge_distribution_surface<Lyt, true> : public Lyt
{
  public:
    explicit charge_distribution_surface(const Lyt& lyt) : Lyt(lyt) {}
};

template <typename Lyt>
class charge_distribution_surface<Lyt, false> : public Lyt
{
  public:
    struct charge_distribution_storage
    {
      private:
        /**
         * The distance matrix is an unordered map with pairs of cells as key and the corresponding euclidean distance
         * as value.
         */
        using distance_matrix = std::unordered_map<std::pair<const cell<Lyt>, const cell<Lyt>>, double>;
        /**
         * The potential matrix is an unordered map with pairs of cells as key and the corresponding electrostatic
         * potential as value.
         */
        using potential_matrix = std::unordered_map<std::pair<const cell<Lyt>, const cell<Lyt>>, double>;
        /**
         * The local electrostatic potential matrix is an unordered map with cells as key and the corresponding local
         * electrostatic potential as value.
         */
        using local_potential = std::unordered_map<cell<Lyt>, double>;

      public:
        std::unordered_map<typename Lyt::coordinate, sidb_charge_state> charge_coordinates{};
        distance_matrix                                                 dist_mat{};
        potential_matrix                                                pot_mat{};
        local_potential                                                 loc_pot{};
        simulation_params                                               sim_params{};
        double                                                          system_energy;
        bool validity;
    };

    using storage = std::shared_ptr<charge_distribution_storage>;

    /**
     * Standard constructor for empty layouts.
     */
    explicit charge_distribution_surface() : Lyt(), strg{std::make_shared<charge_distribution_storage>()}
    {
        static_assert(std::is_same_v<cell<Lyt>, siqad::coord_t>, "Lyt is not based on SiQAD coordinates");
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
    }
    /**
     * Standard constructor for existing layouts and simulation parameter as input.
     */
    explicit charge_distribution_surface(const Lyt& lyt, const simulation_params& sim_params) :
            Lyt(lyt),
            strg{std::make_shared<charge_distribution_storage>(sim_params)}
    {
        static_assert(std::is_same_v<cell<Lyt>, siqad::coord_t>, "Lyt is not based on SiQAD coordinates");
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
    };
    /**
     * Standard constructor for existing layouts.
     */
    explicit charge_distribution_surface(const Lyt& lyt) :
            Lyt(lyt),
            strg{std::make_shared<charge_distribution_storage>()}
    {
        static_assert(std::is_same_v<cell<Lyt>, siqad::coord_t>, "Lyt is not based on SiQAD coordinates");
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
    };
    /**
     * Assigns a given charge state to the given coordinate.
     *
     * @param c Coordinate to assign charge state cs.
     * @param cs Charge state to assign to coordinate c.
     */
    void assign_charge_state(const coordinate<Lyt>& c, const sidb_charge_state& cs) noexcept
    {
        if (!Lyt::is_empty_cell(c))
        {
            if (cs != sidb_charge_state::NONE)
            {
                strg->charge_coordinates.insert_or_assign(c, cs);
            }
            else
            {
                strg->charge_coordinates.erase(c);
            }
        }
    }
    /**
     * Returns the given coordinate's assigned charge state.
     *
     * @param c Coordinate to check.
     * @return Charge state previously assigned to c or NONE if cell owns emtpy cell_type.
     */
    [[nodiscard]] sidb_charge_state get_charge_state(const coordinate<Lyt>& c) const noexcept
    {
        if (const auto it = strg->charge_coordinates.find(c); it != strg->charge_coordinates.cend())
        {
            return it->second;
        }
        return sidb_charge_state::NONE;
    }
    /**
     * Applies a function to all SiDBs' charge states on the surface. Since the charge states are fetched directly from
     * the storage map, the given function has to receive a pair of a coordinate and a charge state as its parameter.
     *
     * @tparam Fn Functor type that has to comply with the restrictions imposed by mockturtle::foreach_element.
     * @param fn Functor to apply to each charge coordinate.
     */
    template <typename Fn>
    void foreach_charge_state(Fn&& fn) const
    {
        mockturtle::detail::foreach_element(strg->charge_coordinates.cbegin(), strg->charge_coordinates.cend(),
                                            std::forward<Fn>(fn));
    }

    /**
     * The Euclidean distance for every cell (needs to exhibit an assigned cell type) pair in the layout is calculated
     * and returned as a matrix (unordered map).
     *
     * @tparam Lyt Cell-level layout type (SiQAD coordinates are required).
     * @tparam Dist Floating-point type for the distance.
     * @param lyt Layout.
     * @return Distance matrix
     */
    void initialize_sidb_distance_matrix()
    {
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(std::is_same_v<cell<Lyt>, siqad::coord_t>, "Lyt is not based on SiQAD coordinates");
        this->foreach_cell(
            [this](const auto& c1)
            {
                this->foreach_cell(
                    [c1, this](const auto& c2)
                    { strg->dist_mat.insert(std::make_pair(std::make_pair(c1, c2), distance_sidb_pair(c1, c2))); });
            });
    };

    /**
     * SiQAD coordinates are converted to a nm location on the Si-substrate by taking Silicon's lattice constants into
     * account (see simulation_parameters.hpp).
     *
     * @tparam Lyt cell layout type (SiQAD coordinates are required).
     * @tparam Dist Data type for the distance.
     * @param c cell of the layout Lyt.
     * @return nm position
     */
    std::pair<double, double> nm_position(const cell<Lyt>& c)
    {
        const auto x = static_cast<double>(c.x * strg->sim_params.lat_a);
        const auto y = static_cast<double>(c.y * strg->sim_params.lat_b + c.z * strg->sim_params.lat_c);
        return std::make_pair(x, y);
    }

    /**
     * All electrostatic inter-potentials between two cells are calculated. The Euclidean distance is provided by the
     * distance matrix. Electrostatic potential between identical cells is set to 0.
     *
     * @tparam Lyt Cell-level layout type (SiQAD coordinates are required).
     * @tparam Dist Floating-point type for the distance.
     * @tparam Potential Floating-point type for the electrostatic potential.
     * @param lyt Layout.
     * @return Potential matrix
     */
    void initialize_sidb_potential_matrix()
    {
        for (const auto& it : strg->dist_mat)
        {
            strg->pot_mat.insert(std::make_pair(it.first, potential_sidb_pair(it.first.first, it.first.second)));
        }
    };

    /**
     * The Euclidean distance in nm between two SiDBs on the H-Si surface is calculated (SiQAD coordinates are
     * required). In the first step, SiQAD coordinates are converted to a nm position on the Si-substrate by taking
     * Silicon's lattice constants into account (see simulation_parameters.hpp). Afterward, the Euclidean distance is
     * calculated.
     *
     * @tparam Lyt cell-level layout type.
     * @tparam Dist Floating-point type for the distance.
     * @param lyt Layout.
     * @param c1 cell coordinate.
     * @param c2 cell coordinate.
     * @return Euclidean distance between c1 and c2.
     */

    [[nodiscard]] constexpr double distance_sidb_pair(const cell<Lyt>& c1, const cell<Lyt>& c2)
    {
        const auto pos_c1 = nm_position(c1);
        const auto pos_c2 = nm_position(c2);
        const auto x      = static_cast<double>(pos_c1.first) - static_cast<double>(pos_c2.first);
        const auto y      = static_cast<double>(pos_c1.second) - static_cast<double>(pos_c2.second);
        return std::hypot(x, y);
    }

    [[nodiscard]] constexpr double dist(const cell<Lyt>& c1, const cell<Lyt>& c2)
    {
        for (auto& it : strg->dist_mat)
        {
            if (it.first.first == c1 && it.first.second == c2)
            {
                strg->dist_mat.at({c1, c2});
                return it.second;
            }
        }
    }

    /**
     * Returns the potential between two cells in the layout.
     *
     * @return number (uint64_t)
     */
    [[nodiscard]] constexpr std::optional<double> pot(const cell<Lyt>& c1, const cell<Lyt>& c2)
    {
        if (auto it = strg->pot_mat.find(std::make_pair(c1, c2)); it != strg->pot_mat.end())
        {
            return it->second;
        }
        return std::nullopt;
    }

    /**
     * Returns the number of SiDBs assigned to a charge state
     *
     * @return number (uint64_t)
     */
    [[nodiscard]] uint64_t num_charges() const noexcept
    {
        return strg->charge_coordinates.size();
    }

    /**
     * The electrostatic potential for a given Euclidean distance is calculated.
     *
     * @tparam Potential Data type for the electrostatic potential.
     * @tparam Dist Euclidean distance.
     * @param dist Euclidean distance value.
     * @param k Coulomb constant (default value).
     * @param lambda_tf Thomas-Fermi screening distance (default value).
     * @return Electrostatic potential value.
     */

    [[nodiscard]] double potential_sidb_pair(const cell<Lyt>& c1, const cell<Lyt>& c2)
    {
        if (auto it = strg->dist_mat.find(std::make_pair(c1, c2)); it != strg->dist_mat.end())
        {
            auto dis = it->second;
            if (dis == 0.0)
            {
                return 0.0;
            }
            return (strg->sim_params.k / dis * std::exp(-dis / strg->sim_params.lambda_tf) *
                    physical_sim_constants::ELECTRIC_CHARGE);
        }
    }

    void local_potential()
    {
        this->foreach_charge_state(
            [this](const auto& cs)
            {
                double collect = 0;
                for (auto& it : strg->pot_mat)
                {

                    if (it.first.first == cs.first)
                    {
                        collect += it.second * transform_to_sign(get_charge_state(cs.first));
                    }
                }
                strg->loc_pot.insert_or_assign(cs.first, collect);
            });
    }

    std::optional<double> get_loc_pot(const cell<Lyt>& c1)
    {
        if (auto it = strg->loc_pot.find(c1); it != strg->loc_pot.end())
        {
            return strg->loc_pot[c1];
        }
        return std::nullopt;
    };

    /**
     * Calculate the system's total electrostatic potential energy.
     *
     * @param lyt charge distribution layout
     * @param loc_pot local electrostatic potential
     */
    void system_energy()
    {
        double total_energy = 0;
        for (auto& it : strg->loc_pot)
        {
            total_energy += 0.5 * it.second * transform_to_sign(get_charge_state(it.first));
        }
        strg->system_energy = total_energy;
    }

    [[nodiscard]] double get_system_energy()
    {
        return strg->system_energy;
    }


    void validity_check()
    {
        bool valid = false;
        for (auto& it : strg->loc_pot)
        {
            valid = (((this->get_charge_state(it.first) == sidb_charge_state::NEGATIVE) &&
                      ((-it.second + strg->sim_params.mu) < physical_sim_constants::POP_STABILITY_ERR)) ||
                     ((this->get_charge_state(it.first) == sidb_charge_state::POSITIVE) &&
                      ((-it.second + strg->sim_params.mu_p) > physical_sim_constants::POP_STABILITY_ERR)) ||
                     ((this->get_charge_state(it.first) == sidb_charge_state::NEUTRAL) &&
                      ((-it.second + strg->sim_params.mu) > physical_sim_constants::POP_STABILITY_ERR) &&
                      (-it.second + strg->sim_params.mu_p) < physical_sim_constants::POP_STABILITY_ERR));

            if (!valid)
            {
                strg->validity = false;
                break;
            }

            else
            {
                strg->validity = true;
            }
        }

        auto hopDel = [this](const cell<Lyt>& c1, const cell<Lyt>& c2)
        {
            int dn_i = (this->get_charge_state(c1) == sidb_charge_state::NEGATIVE) ? 1 : -1;
            int dn_j = -dn_i;

            return strg->loc_pot.at(c1) * dn_i + strg->loc_pot.at(c1) * dn_j - strg->pot_mat.at(std::make_pair(c1, c2)) * 1;
        };

        for (auto& it : strg->loc_pot)
        {
            if (this->get_charge_state(it.first) == sidb_charge_state::POSITIVE)
            {
                continue;
            }

            for (auto& it_second : strg->loc_pot)
            {
                auto E_del = hopDel(it.first, it_second.first);
                if ((transform_to_sign(get_charge_state(it.first)) >
                     transform_to_sign(get_charge_state(it.first))) &&
                    (E_del < -physical_sim_constants::POP_STABILITY_ERR))
                {
                    strg->validity = false;
                }
            }
        }
        //strg->validity = true;
    }

    [[nodiscard]] bool get_validity()
    {
        return strg->validity;
    }

  private:
    storage strg;
};

template <class T>
charge_distribution_surface(const T&) -> charge_distribution_surface<T>;

}  // namespace fiction

#endif  // CHARGE_DISTRIBUTION_SURFACE
