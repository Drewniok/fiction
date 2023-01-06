//
// Created by Jan Drewniok on 23.11.22.
//

#ifndef FICTION_CHARGE_DISTRIBUTION_SURFACE_HPP
#define FICTION_CHARGE_DISTRIBUTION_SURFACE_HPP

#include "fiction/algorithms/path_finding/distance.hpp"
#include "fiction/algorithms/simulation_sidb/simulation_parameters.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"

#include <cassert>
#include <limits>

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
        using distance_matrix = std::vector<std::vector<double>>;
        /**
         * The potential matrix is an unordered map with pairs of cells as key and the corresponding electrostatic
         * potential as value.
         */
        using potential_matrix = std::vector<std::vector<double>>;
        /**
         * The local electrostatic potential matrix is an unordered map with cells as key and the corresponding local
         * electrostatic potential as value.
         */
        using local_potential = std::vector<double>;

      public:
        explicit charge_distribution_storage(const simulation_params& sim_param_default = simulation_params{}) :
                sim_params{sim_param_default} {};

        simulation_params              sim_params{};
        std::vector<sidb_charge_state> charge_coordinates{};
        std::vector<cell<Lyt>>         coord_vec;
        distance_matrix                dist_mat{};
        potential_matrix               pot_mat{};
        local_potential                loc_pot{};
        double                         system_energy{};
        bool                           validity = false;
        std::pair<uint64_t, uint8_t>   charge_index{};
        uint64_t                       max_charge_index{};
    };

    using storage = std::shared_ptr<charge_distribution_storage>;

    /**
     * Standard constructor for empty layouts.
     */
    explicit charge_distribution_surface(const sidb_charge_state& cs                = sidb_charge_state::NEGATIVE,
                                         const simulation_params& sim_param_default = simulation_params{}) :
            Lyt(),
            strg{std::make_shared<charge_distribution_storage>(sim_param_default)}
    {
        initialize_new(cs);
        static_assert(std::is_same_v<cell<Lyt>, siqad::coord_t>, "Lyt is not based on SiQAD coordinates");
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
    }
    /**
     * Standard constructor for existing layouts and simulation parameter as input.
     */
    explicit charge_distribution_surface(const Lyt&               lyt,
                                         const simulation_params& sim_param_default = simulation_params{},
                                         const sidb_charge_state& cs                = sidb_charge_state::NEGATIVE) :
            Lyt(lyt),
            strg{std::make_shared<charge_distribution_storage>(sim_param_default)}
    {
        initialize_new(cs);
        static_assert(std::is_same_v<cell<Lyt>, siqad::coord_t>, "Lyt is not based on SiQAD coordinates");
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
    };
    /**
     * Copy constructor.
     */
    explicit charge_distribution_surface(const charge_distribution_surface<Lyt>& lyt) :
            strg{std::make_shared<charge_distribution_storage>(*lyt.strg)}
    {}
    /**
     * initialisation function used for the construction of the charge distribution layout.
     *
     * @param cs charge state assigned to all SiDBs.
     */
    void initialize(const sidb_charge_state& cs)
    {
        strg->coord_vec.reserve(this->num_cells());
        this->foreach_cell([this](const auto& c) { strg->coord_vec.push_back(c); });
        this->foreach_cell([this, &cs](const auto& c1) { strg->charge_coordinates.insert(std::make_pair(c1, cs)); });
        assert((std::pow(strg->sim_params.base, this->num_cells()) - 1) < std::numeric_limits<uint64_t>::max() &&
               "number of SiDBs is too large");
        strg->max_charge_index = std::pow(strg->sim_params.base, this->num_cells()) - 1;
        this->chargeconf_to_index();
        this->initialize_sidb_distance_matrix();
        this->initialize_sidb_potential_matrix();
        this->local_potential();
        this->system_energy();
        this->validity_check();
    };


    void initialize_new(const sidb_charge_state& cs)
    {
        strg->coord_vec.reserve(this->num_cells());
        strg->charge_coordinates.reserve(this->num_cells());
        this->foreach_cell([this](const auto& c) { strg->coord_vec.push_back(c); });
        this->foreach_cell([this, &cs](const auto& c1) { strg->charge_coordinates.push_back(cs); });
        assert((std::pow(strg->sim_params.base, this->num_cells()) - 1) < std::numeric_limits<uint64_t>::max() &&
               "number of SiDBs is too large");
        strg->max_charge_index = std::pow(strg->sim_params.base, this->num_cells()) - 1;
        this->chargeconf_to_index();
        this->initialize_sidb_distance_matrix_new();
        this->initialize_sidb_potential_matrix_new();
        this->local_potential();
        this->system_energy_new();
        this->validity_check_new();
    };

    /**
     * Assigns a given charge state to the given coordinate.
     *
     * @param c Coordinate to assign charge state cs.
     * @param cs Charge state to assign to coordinate c.
     */
    void assign_charge_state(const coordinate<Lyt>& c, const sidb_charge_state& cs) const
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

    void assign_charge_state_new(const coordinate<Lyt>& c, const sidb_charge_state& cs) const
    {
        if (!Lyt::is_empty_cell(c))
        {
            if (cs != sidb_charge_state::NONE)
            {
                if (auto it = std::find(strg->coord_vec.begin(), strg->coord_vec.end(), c); it != strg->coord_vec.end())
                {
                    int index = std::distance(strg->coord_vec.begin(), it);
                    strg->charge_coordinates[index] = cs;
                }
                else
                {
                    std::cout << "Element not found" << std::endl;
                }

            }
        }
    }

    /**
     * Returns the given coordinate's assigned charge state.
     *
     * @param c Coordinate to check.
     * @return Charge state previously assigned to c or NONE if cell owns emtpy cell_type.
     */

//    [[nodiscard]] sidb_charge_state get_charge_state(const coordinate<Lyt>& c) const noexcept
//    {
//        if (const auto it = strg->charge_coordinates.find(c); it != strg->charge_coordinates.cend())
//        {
//            return it->second;
//        }
//        return sidb_charge_state::NONE;
//    }
//
//    [[nodiscard]] sidb_charge_state get_charge_state_new(const coordinate<Lyt>& c) const noexcept
//    {
//        if (const auto it = strg->charge_coordinates.find(c); it != strg->charge_coordinates.cend())
//        {
//            return it->second;
//        }
//        return sidb_charge_state::NONE;
//    }
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
//    void initialize_sidb_distance_matrix() const
//    {
//        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
//        static_assert(std::is_same_v<cell<Lyt>, siqad::coord_t>, "Lyt is not based on SiQAD coordinates");
//        this->foreach_cell(
//            [this](const auto& c1)
//            {
//                this->foreach_cell(
//                    [c1, this](const auto& c2)
//                    { strg->dist_mat.insert(std::make_pair(std::make_pair(c1, c2), distance_sidb_pair(c1, c2))); });
//            });
//    };

    void initialize_sidb_distance_matrix_new() const
    {

        for (unsigned int i = 0; i < strg->coord_vec.size(); i++)
        {

            std::vector<double> vec{};
            for (unsigned int j = 0; j < strg->coord_vec.size(); j++)
            {
                if (i==j)
                {
                    vec.push_back(0);
                }
                else
                {
                    vec.push_back(distance_sidb_pair(strg->coord_vec[i], strg->coord_vec[j]));
                }
            }
            strg->dist_mat.push_back(vec);
        }
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
    std::pair<double, double> nm_position(const cell<Lyt>& c) const
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
    void initialize_sidb_potential_matrix() const
    {
        for (const auto& it : strg->dist_mat)
        {
            strg->pot_mat.insert(std::make_pair(it.first, potential_sidb_pair(it.first.first, it.first.second)));
        }
    };

    void initialize_sidb_potential_matrix_new() const
    {

        for (unsigned int i = 0; i < strg->coord_vec.size(); i++)
        {
            std::vector<double> vec{};
            for (unsigned int j = 0; j < strg->coord_vec.size(); j++)
            {
                if (i==j)
                {
                    vec.push_back(0);
                }
                else
                {
                    vec.push_back(potential_sidb_pair_new(i, j));
                }
            }
            strg->pot_mat.push_back(vec);
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
    [[nodiscard]] constexpr double distance_sidb_pair(const cell<Lyt>& c1, const cell<Lyt>& c2) const
    {
        const auto pos_c1 = nm_position(c1);
        const auto pos_c2 = nm_position(c2);
        const auto x      = static_cast<double>(pos_c1.first) - static_cast<double>(pos_c2.first);
        const auto y      = static_cast<double>(pos_c1.second) - static_cast<double>(pos_c2.second);
        return std::hypot(x, y);
    }

    /**
     * The euclidean distance of two cells (in this case two SiDBs) is returned.
     *
     * @tparam Lyt cell-level layout type.
     * @param c1 cell coordinate.
     * @param c2 cell coordinate.
     * @return Euclidean distance between c1 and c2.
     */
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
     * Returns the number of SiDBs to which a charge state is assigned.
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
    [[nodiscard]] double potential_sidb_pair(const cell<Lyt>& c1, const cell<Lyt>& c2) const
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

    [[nodiscard]] double potential_sidb_pair_new(const int& c1, const int& c2) const
    {
        if (strg->dist_mat[c1][c2] == 0)
        {
            return 0.0;
        }
    else
    {
        return (strg->sim_params.k / strg->dist_mat[c1][c2] * std::exp(-strg->dist_mat[c1][c2] / strg->sim_params.lambda_tf) *
                physical_sim_constants::ELECTRIC_CHARGE);
    }
}



    /**
     * The local electrostatic potential is calculated for each SiDB position.
     * It depends on the charge configuration, hence, the local electrostatic potential has to be updated when the
     * charge distribution is updated.
     */
//    void local_potential()
//    {
//        this->foreach_charge_state(
//            [this](const auto& cs)
//            {
//                double collect = 0;
//                for (auto& it : strg->pot_mat)
//                {
//                    if (it.first.second == cs.first)
//                    {
//                        collect += it.second * transform_to_sign(get_charge_state(it.first.first));
//                    }
//                }
//                strg->loc_pot[cs.first] = collect;
//            });
//    }

    void local_potential()
    {
        strg->loc_pot.resize(this->num_cells());
        std::fill(strg->loc_pot.begin(), strg->loc_pot.end(), 0);
        for (int i = 0; i <strg->coord_vec.size(); i++)
        {
                double collect = 0;
                for (int j= 0; j <  strg->pot_mat[i].size(); j++)
                {
                    collect += strg->pot_mat[i][j] * transform_to_sign(strg->charge_coordinates[j]);
                }
                strg->loc_pot[i] = collect;
            }
    }
    /**
     * The local electrostatic potential is calculated for each SiDB position.
     * It depends on the charge configuration, hence, the local electrostatic potential has to be updated when the
     * charge distribution is updated.
     * @param c1 cell for which the value of the local electrostatic potential is wanted.
     * @return local potential is returned if an SiDB is at position c1. Otherwise, std::nullopt is returned.
     */
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
//    void system_energy()
//    {
//        double total_energy = 0;
//        for (auto& it : strg->loc_pot)
//        {
//            total_energy += 0.5 * it.second * transform_to_sign(get_charge_state(it.first));
//        }
//        strg->system_energy = total_energy;
//    }

    void system_energy_new()
    {
        double total_energy = 0;
        for (int i = 0; i < strg->loc_pot.size(); i++)
        {
            //total_energy += 0.5 * strg->loc_pot[i] * transform_to_sign(strg->charge_coordinates[i]);
            total_energy += 0.5 * strg->loc_pot[i];
        }
        strg->system_energy = total_energy;
    }
    /**
     * Return the system's total electrostatic potential energy.
     *
     * @return system's total electrostatic potential energy.
     */
    [[nodiscard]] double get_system_energy() const
    {
        return strg->system_energy;
    }
    /**
     * The physically validity of the current charge distribution is evaluated and stored in the storage struct.
     * A charge distribution is valid if the *Population Stability* and the *Configuration Stability* is given.
     *
     * @return bool indicating the validity.
     */
    void validity_check()
    {
        int  flag_fd = 0;
        for (auto& it : strg->loc_pot)
        {
            bool valid = (((this->get_charge_state(it.first) == sidb_charge_state::NEGATIVE) &&
                      ((-it.second + strg->sim_params.mu) < physical_sim_constants::POP_STABILITY_ERR)) ||
                     ((this->get_charge_state(it.first) == sidb_charge_state::POSITIVE) &&
                      ((-it.second + strg->sim_params.mu_p) > -physical_sim_constants::POP_STABILITY_ERR)) ||
                     ((this->get_charge_state(it.first) == sidb_charge_state::NEUTRAL) &&
                      ((-it.second + strg->sim_params.mu) > -physical_sim_constants::POP_STABILITY_ERR) &&
                      (-it.second + strg->sim_params.mu_p) < physical_sim_constants::POP_STABILITY_ERR));

            if (!valid)
            {
                strg->validity = false;
                flag_fd += 1;
                break;
            }
        }

        if (flag_fd == 0)
        {

            auto hopDel = [this](const cell<Lyt>& c1, const cell<Lyt>& c2)
            {
                int dn_i = (this->get_charge_state(c1) == sidb_charge_state::NEGATIVE) ? 1 : -1;
                int dn_j = -dn_i;

                return strg->loc_pot.at(c1) * dn_i + strg->loc_pot.at(c2) * dn_j -
                       strg->pot_mat.at(std::make_pair(c1, c2)) * 1;
            };

            int flag_value = 0;
            for (auto& it : strg->loc_pot)
            {
                if (this->get_charge_state(it.first) == sidb_charge_state::POSITIVE)
                {
                    continue;
                }

                for (auto& it_second : strg->loc_pot)
                {
                    if (flag_value == 1)
                    {
                        break;
                    }

                    auto E_del = hopDel(it.first, it_second.first);
                    if ((transform_to_sign(get_charge_state(it_second.first)) >
                         transform_to_sign(get_charge_state(it.first))) &&
                        (E_del < -physical_sim_constants::POP_STABILITY_ERR))
                    {
                        flag_value = 1;
                        break;
                    }
                }
            }

            if (flag_value == 0)
            {
                strg->validity = true;
            }
            else
            {
                strg->validity = false;
            }
        }
    }

    void validity_check_new()
    {
        int  flag_fd = 0;
        int counter_loop = 0;
        for (auto& it : strg->loc_pot)
        {
            bool valid = (((strg->charge_coordinates[counter_loop] == sidb_charge_state::NEGATIVE) &&
                           ((-it+ strg->sim_params.mu) < physical_sim_constants::POP_STABILITY_ERR)) ||
                          ((strg->charge_coordinates[counter_loop] == sidb_charge_state::POSITIVE) &&
                           ((-it + strg->sim_params.mu_p) > -physical_sim_constants::POP_STABILITY_ERR)) ||
                          ((strg->charge_coordinates[counter_loop] == sidb_charge_state::NEUTRAL) &&
                           ((-it + strg->sim_params.mu) > -physical_sim_constants::POP_STABILITY_ERR) &&
                           (-it + strg->sim_params.mu_p) < physical_sim_constants::POP_STABILITY_ERR));
            counter_loop +=1;
            if (!valid)
            {
                strg->validity = false;
                flag_fd += 1;
                break;
            }
        }

        if (flag_fd == 0)
        {

            auto hopDel = [this](const int& c1, const int& c2)
            {
                int dn_i = (strg->charge_coordinates[c1] == sidb_charge_state::NEGATIVE) ? 1 : -1;
                int dn_j = -dn_i;

                return strg->loc_pot[c1] * dn_i + strg->loc_pot[c2] * dn_j -
                       strg->pot_mat[c1][c2] * 1;
            };

            int flag_value = 0;
            int counter_new = 0;
            for (int i =0; i < strg->loc_pot.size(); i++)
            {
                if (strg->charge_coordinates[counter_new] == sidb_charge_state::POSITIVE)
                {
                    continue;
                }

                for (int j =0; j < strg->loc_pot.size(); j++)
                {
                    if (flag_value == 1)
                    {
                        break;
                    }

                    auto E_del = hopDel(i, j);
                    if ((transform_to_sign(strg->charge_coordinates[j]) >
                         transform_to_sign(strg->charge_coordinates[i])) &&
                        (E_del < -physical_sim_constants::POP_STABILITY_ERR))
                    {
                        flag_value = 1;
                        break;
                    }
                }
            }

            if (flag_value == 0)
            {
                strg->validity = true;
            }
            else
            {
                strg->validity = false;
            }
        }
    }

    [[nodiscard]] bool get_validity()
    {
        return strg->validity;
    }

    /**
     * The charge distribution of the charge distribution surface is converted to a unique index. It is used to map
     * every possible charge distribution in a SiDB layout.
     *
     * @tparam Lyt cell-level layout.
     * @tparam base number of charge states per SiDB. Base = 2 when only neutrally and negatively charged SiDBs are
     * taken into account, otherwise base = 3.
     * @return a pair with the calculated charge distribution's corresponding charge index and the chosen base as the
     * first and second element, respectively.
     */
    void chargeconf_to_index() const
    {
        uint8_t base = strg->sim_params.base;
        assert(base == 2 || base == 3 && "base must be 2 or 3");
        uint64_t chargeindex = 0;
        int      counter     = 0;
        for (const auto& c : strg->charge_coordinates)
        {
            chargeindex += static_cast<unsigned int>((transform_to_sign(c) + 1) *
                                                     std::pow(base, this->num_charges() - counter - 1));
            counter += 1;
        }
        strg->charge_index = {chargeindex, base};
    }

    /**
     * The charge index of the current charge distribution is returned.
     *
     * @return a pair with the charge index and the used base is returned.
     */
    [[nodiscard]] std::pair<uint64_t, uint8_t> get_charge_index()
    {
        return strg->charge_index;
    }

    std::unordered_map<cell<Lyt>, sidb_charge_state> get_charge_coordinates()
    {
        return strg->charge_coordinates;
    }
    /**
     *  The unique index is converted to the charge distribution of the charge distribution surface.
     *
     * @tparam Lyt cell-level layout.
     * @param cp a pair of the charge index and the chosen base number.
     * @return a pair with the calculated charge distribution's corresponding charge index and the chosen base as the
     * first and second element, respectively.
     */
    void index_to_chargeconf() const
    {

        int charge_quot = strg->charge_index.first;
        int base        = strg->charge_index.second;
        int num_charges = this->num_charges() - 1;

        while (charge_quot > 0)
        {
            div_t d;
            d           = div(charge_quot, base);
            charge_quot = d.quot;
            int counter = 0;
            for (auto& it : strg->coord_vec)
            {

                if (counter == num_charges)
                {
                    this->assign_charge_state_new(it, sign_to_label(d.rem - 1));
                }
                counter++;
            }
            num_charges -= 1;
        }
    }

    /**
     *  The unique index is increased by one, but only if it is less than or equal to the maximum charge index after it has been increased. If that's the case, it is increased and
     *  afterward, the charge configuration is updated by invoking the index_to_chargeconf() function.
     *
     */
    void increase_charge_index()
    {
        if (strg->charge_index.first == strg->max_charge_index)
        {
            std::cout << "maximal index is exceeded";
        }
        else
        {
            strg->charge_index.first += 1;
            this->index_to_chargeconf();
        }
    }

    uint64_t get_max_charge_index()
    {
        return strg->max_charge_index;
    }

    void assign_charge_index(const uint64_t &index)
    {
        strg->charge_index.first = index;
        this->index_to_chargeconf();
    }

    void next_N(const double &alpha)
    {
        auto                                 count    = 0;
        float                                dist_max = 0;
        cell<Lyt>                            coord{};
        std::unordered_map<cell<Lyt>, float> max_dist_unocc_occ{};

        for (auto& unocc : strg->charge_coordinates)
        {
            if (unocc.second != sidb_charge_state::NEUTRAL)
            {
                continue;
            }
            count += 1;

            auto                                 dist_min = MAXFLOAT;
            for (auto& occ : strg->charge_coordinates)
            {
                if (((occ.second == sidb_charge_state::NEGATIVE) || (occ.second == sidb_charge_state::POSITIVE)) &&
                    (distance_sidb_pair(unocc.first, occ.first) < dist_min))
                {
                    dist_min = distance_sidb_pair(unocc.first, occ.first);
                    coord    = unocc.first;
                }
            }

            if (count == 1)
            {
                max_dist_unocc_occ.insert(std::pair(coord, dist_min));
                dist_max = dist_min;
            }
            else if (count > 1)
            {
                max_dist_unocc_occ.insert(std::pair(coord, dist_min));

                if (dist_min > dist_max)
                {
                    dist_max = dist_min;
                }

                else if ((dist_min == dist_max) && (sum_distance(coord, unocc.first)))
                {
                    dist_max = dist_min;
                }
            }
        }

        auto it = max_dist_unocc_occ.begin();
        while (it != max_dist_unocc_occ.end())
        {
            if (it->second < alpha * dist_max)
            {
                it = max_dist_unocc_occ.erase(it);
            }
            else
            {
                it++;
            }
        }

        auto candidate = max_dist_unocc_occ.begin();
        auto end       = max_dist_unocc_occ.end();

        // Generate a random number of steps to advance from the beginning iterator
        std::random_device              rd;
        std::mt19937                    gen(rd());
        std::uniform_int_distribution<> dis(0, std::distance(candidate, end) - 1);
        int                             steps = dis(gen);

        // Advance the iterator by the random number of steps
        std::advance(candidate, steps);
        // candidate += steps;

        this->assign_charge_state(candidate->first, sidb_charge_state::NEGATIVE);

        strg->system_energy += -this->get_loc_pot(candidate->first).value();

        for (auto& loop : strg->pot_mat)
        {
            if (loop.first.second== candidate->first)
            {
//                strg->loc_pot.insert_or_assign(
//                    std::pair(loop.first.second, -this->pot(loop.first.first, candidate->first).value()));
            strg->loc_pot[loop.first.first] += -this->pot(loop.first.first, candidate->first).value();
            }
        }




    }

    bool sum_distance(const cell<Lyt>& c1, const cell<Lyt>& c2)
    {
        float sum_old = 0;
        float sum_now = 0;
        for (auto& it : strg->charge_coordinates)
        {
            if (it.second == sidb_charge_state::NEUTRAL)
            {
                continue;
            }
            sum_old += distance_sidb_pair(it.first, c1);
            sum_now += distance_sidb_pair(it.first, c2);
        }
        return sum_now > sum_old;
    }

  private:
    storage strg;
};

template <class T>
charge_distribution_surface(const T&) -> charge_distribution_surface<T>;

template <class T>
charge_distribution_surface(const T&, const simulation_params&) -> charge_distribution_surface<T>;

}  // namespace fiction

#endif  // CHARGE_DISTRIBUTION_SURFACE
