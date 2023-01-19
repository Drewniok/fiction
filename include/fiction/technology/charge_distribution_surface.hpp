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
#include "fiction/utils/hash.hpp"

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
         * The distance matrix is a vector of vectors storing the euclidean distance.
         */
        using distance_matrix = std::vector<std::vector<double>>;

        /**
         * The distance matrix is a vector of vectors storing the euclidean distance.
         */
        using potential_matrix = std::vector<std::vector<double>>;

        /**
         * It is a vector that stores the electrostatic potential
         * electrostatic potential as value.
         */
        using local_potential = std::vector<double>;

      public:
        explicit charge_distribution_storage(const physical_params& sim_param_default = physical_params{}) :
                phys_params{sim_param_default} {};

        physical_params        phys_params{};  // all physical parameters used for the simulation are stored in a struct
        std::vector<cell<Lyt>> sidb_order{};   // all cells that are occupied by an SiDB are stored in this vector
        std::vector<sidb_charge_state>
            cell_charge{};  // the SiDBs' charge states are stored. Corresponding cells are stored in "sidb_order"
        distance_matrix dist_mat{};  // distance between SiDBs are stored as matrix
        potential_matrix
            pot_mat{};  // electrostatic potential between SiDBs are stored as matrix (here, still charge independent)
        local_potential loc_pot{};  // electrostatic potential at each SiDB postion (has to be updated when charge
                                    // distribution is changed)
        double system_energy{0.0};  // electrostatic energy of a given charge distribution
        bool   validity = false;    // labels if given charge distribution is valid (compare Paper XXX)
        std::pair<uint64_t, uint8_t>
            charge_index{};  // each charge distribution is assigned a unique index (first entry of pair), second one
                             // stores the base number (2 or 3 state simulation)
        uint64_t max_charge_index{};  // depending on the number of SiDBs and the base number, a maximal number of
                                      // possible charge distributions exists
    };

    using storage = std::shared_ptr<charge_distribution_storage>;

    /**
     * Standard constructor for empty layouts.
     *
     * @param sidb_charge_state charge state used for the initialization of all SiDBs, default is negative charge
     * @param physical_params physical parameters used for the simulation (µ_minus, base number,...).
     */
    explicit charge_distribution_surface(const physical_params&   sim_param_default = physical_params{},
                                         const sidb_charge_state& cs                = sidb_charge_state::NEGATIVE) :
            Lyt(),
            strg{std::make_shared<charge_distribution_storage>(sim_param_default)}
    {
        static_assert(std::is_same_v<cell<Lyt>, siqad::coord_t>, "Lyt is not based on SiQAD coordinates");
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
        initialize(cs);
    }
    /**
     * Standard constructor for existing layouts.
     *
     * @param sidb_charge_state charge state used for the initialization of all SiDBs, default is negative charge
     * @param physical_params physical parameters used for the simulation (µ_minus, base number,...).
     */
    explicit charge_distribution_surface(const Lyt& lyt, const physical_params& sim_param_default = physical_params{},
                                         const sidb_charge_state& cs = sidb_charge_state::NEGATIVE) :
            Lyt(lyt),
            strg{std::make_shared<charge_distribution_storage>(sim_param_default)}
    {
        static_assert(std::is_same_v<cell<Lyt>, siqad::coord_t>, "Lyt is not based on SiQAD coordinates");
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
        initialize(cs);
    };

    /**
     * Copy constructor.
     *
     * @param charge distribution layout that wants to be copied
     */
    explicit charge_distribution_surface(const charge_distribution_surface<Lyt>& lyt) :
            strg{std::make_shared<charge_distribution_storage>(*lyt.strg)}
    {}

    /**
     * Initialisation function used for the construction of the charge distribution layout.
     *
     * @param cs charge state assigned to all SiDBs.
     */
    void initialize(const sidb_charge_state& cs = sidb_charge_state::NEGATIVE)
    {
        strg->sidb_order.reserve(this->num_cells());
        strg->cell_charge.reserve(this->num_cells());
        this->foreach_cell([this](const auto& c1) { strg->sidb_order.push_back(c1); });
        this->foreach_cell([this, &cs](const auto&) { strg->cell_charge.push_back(cs); });
        if (strg->phys_params.base == 3)
        {
            assert(this->num_cells() < 41 && "number of SiDBs is too large");
        }
        else
        {
            assert(this->num_cells() < 64 && "number of SiDBs is too large");
        }
        this->chargeconf_to_index();
        this->initialize_distance_matrix();
        this->initialize_potential_matrix();
        strg->max_charge_index = static_cast<uint64_t>(std::pow(static_cast<double>(strg->phys_params.base), this->num_cells()) - 1);
        this->local_potential();
        this->system_energy();
        this->validity_check();
    };

    /**
     * Function to set the physical parameters for the simulation
     *
     * @param sim_param Physical parameters to be set
     */
    void set_physical_parameters(const physical_params& sim_param)
    {
        strg->phys_params         = sim_param;
        strg->charge_index.second = sim_param.base;
        strg->max_charge_index    = static_cast<uint64_t>(std::pow(strg->phys_params.base, this->num_cells())) - 1;
        this->local_potential();
        this->system_energy();
        this->validity_check();
    }

    [[nodiscard]] bool charge_exists(const sidb_charge_state& cs) const
    {
        return std::accumulate(strg->cell_charge.begin(), strg->cell_charge.end(), true,
                               [&cs](bool acc, const sidb_charge_state& c) { return acc && c == cs; });
    }

    /**
     * Retrieves the physical parameters of the simulation.
     *
     * @return physical_params struct containing the physical parameters of the simulation.
     */
    [[nodiscard]] physical_params get_phys_params() const
    {
        return strg->phys_params;
    }

    /**
     * This function assigns the given charge state to the cell of the layout at the specified index. It updates the
     `cell_charge` member of `strg` object with the new charge state of the specified cell.
     *
     * @param i The index of the cell.
     * @param cs The charge state to be assigned to the cell.
     *

     */
    void state_index(const uint64_t& i, const sidb_charge_state& cs) const
    {
        strg->cell_charge[i] = cs;
    }

    /**
     * This function assigns the given charge state to the cell of the layout.
     *
     * @param cell The cell whose charge state is to be assigned.
     * @param cs The charge state to be assigned to the cell.
     */
    void assign_charge_state_cell(const cell<Lyt>& cell, const sidb_charge_state& cs) const
    {
        auto index = cell_to_index(cell);
        if (index != -1)
        {
            strg->cell_charge[index] = cs;
        }
        else
        {
            std::cout << "not part of the charge distribution surface" << std::endl;
        }
    }

    /**
     * This function assigns the given charge state to the cell (accessed by index) of the layout.
     *
     * @param index The cell's index whose charge state is to be assigned.
     * @param cs The charge state to be assigned to the cell.
     */
    void assign_charge_state_index(const uint64_t& index, const sidb_charge_state& cs) const
    {
        strg->cell_charge[index] = cs;
    }

    /**
     * Sets the charge state of all NOT "EMPTY" cells in the layout to a given state.
     *
     * @param cs The charge state to be assigned to all the cells.
     */
    void set_charge_states(const sidb_charge_state& cs) const
    {
        for (uint64_t i = 0u; i < strg->cell_charge.size(); i++)
        {
            strg->cell_charge[i] = cs;
        }
    }

    /**
     * Returns the charge state of a cell of the layout at a given index.
     *
     * @param index The index of the cell.
     * @return The charge state of the cell at the given index.
     */
    [[nodiscard]] sidb_charge_state get_charge_state_index(const uint64_t& index) const noexcept
    {
        if (index < (strg->cell_charge.size()))
        {
            return strg->cell_charge[index];
        }
        return sidb_charge_state::NONE;
    }

    /**
     * Returns the charge state of a cell of the layout at a given coordinate.
     *
     * @param c The coordinate of the cell.
     * @return The charge state of the cell at the given coordinate.
     */
    [[nodiscard]] sidb_charge_state get_charge_state_cell(const coordinate<Lyt>& c) const noexcept
    {
        if (auto index = cell_to_index(c); index != -1)
        {
            return strg->cell_charge[index];
        }
        return sidb_charge_state::NONE;
    }

    /**
     * @brief Initializes the distance matrix between all the cells of the layout.
     */
    void initialize_distance_matrix() const
    {
        strg->dist_mat = std::vector<std::vector<double>>(this->num_cells(), std::vector<double>(this->num_cells(), 0));

        for (uint64_t i = 0u; i < strg->sidb_order.size(); i++)
        {
            for (uint64_t j = 0u; j < strg->sidb_order.size(); j++)
            {
                strg->dist_mat[i][j] = distance_sidb_pair(strg->sidb_order[i], strg->sidb_order[j]);
            }
        }
    };

    /**
     * Computes the position of a cell in nanometers.
     *
     * @param c The cell to compute the position for.
     * @return A pair of double values representing the x,y position of the cell in nanometers.
     */
    std::pair<double, double> nm_position(const cell<Lyt>& c) const
    {
        const auto x = static_cast<double>(c.x * strg->phys_params.lat_a);
        const auto y = static_cast<double>(c.y * strg->phys_params.lat_b + c.z * strg->phys_params.lat_c);
        return std::make_pair(x, y);
    }

    /**
     * Initializes the potential matrix between all the cells of the layout.
     */
    void initialize_potential_matrix() const
    {
        strg->pot_mat = std::vector<std::vector<double>>(this->num_cells(), std::vector<double>(this->num_cells(), 0));
        for (uint64_t i = 0u; i < strg->sidb_order.size(); i++)
        {
            for (uint64_t j = 0u; j < strg->sidb_order.size(); j++)
            {
                strg->pot_mat[i][j] = potential_sidb_pair_index(i, j);
            }
        }
    };

    /**
     * Computes the distance between two cells in nanometers.
     *
     * @param c1 The first cell.
     * @param c2 The second cell.
     * @return The distance between the two cells in nanometers.
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
     * @brief Finds the index of a cell in the layout.
     *
     * @param cell The cell to find the index of.
     * @return The index of the cell in the layout. Returns -1 if the cell is not in the layout.
     */
    [[nodiscard]] constexpr int64_t cell_to_index(const cell<Lyt>& cell) const
    {
        if (auto it = std::find(strg->sidb_order.begin(), strg->sidb_order.end(), cell); it != strg->sidb_order.end())
        {
            return static_cast<int64_t>(std::distance(strg->sidb_order.begin(), it));
        }
        return -1;
    }

    /**
     *  Calculates the distance between two cells
     *  @param cell1 the first cell to compare
     *  @param cell2 the second cell to compare
     *  @return a constexpr double representing the distance between the two cells
     */
    [[nodiscard]] constexpr double get_distance_cell(const cell<Lyt>& cell1, const cell<Lyt>& cell2)
    {
        return strg->dist_mat[cell_to_index(cell1)][cell_to_index(cell2)];
    }

    /**
     * Calculates and returns the distance index between two input values.
     *
     * @param input1 The first input value
     * @param input2 The second input value
     * @return The distance index between input1 and input2
     */
    [[nodiscard]] constexpr double get_distance_index(const uint64_t& input1, const uint64_t& input2)
    {
        return strg->dist_mat[input1][input2];
    }

    /**
     * Calculates and returns the potential of two cells.
     *
     * @param input1 The first cell
     * @param input2 The second cell
     * @return The potential between input1 and input2
     */
    [[nodiscard]] constexpr double get_potential_cell(const cell<Lyt>& cell1, const cell<Lyt>& cell2)
    {
        return strg->pot_mat[cell_to_index(cell1)][cell_to_index(cell2)];
    }

    /**
     * Calculates and returns the potential of two indices.
     *
     * @param input1 The first index
     * @param input2 The second index
     * @return The potential between input1 and input2
     */
    [[nodiscard]] constexpr double get_potential_index(const uint64_t& input1, const uint64_t& input2)
    {
        return strg->pot_mat[input1][input2];
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
    [[nodiscard]] double potential_sidb_pair_index(const uint64_t& c1, const uint64_t& c2) const
    {
        if (strg->dist_mat[c1][c2] == 0)
        {
            return 0.0;
        }

        return (strg->phys_params.k / strg->dist_mat[c1][c2] *
                std::exp(-strg->dist_mat[c1][c2] / strg->phys_params.lambda_tf) *
                physical_sim_constants::ELECTRIC_CHARGE);
    }

    /**
     * Calculates and returns the potential of a pair of cells based on their distance and simulation parameters.
     *
     * @param c1 The first cell
     * @param c2 The second cell
     * @return The potential between c1 and c2
     */
    [[nodiscard]] double potential_sidb_pair_cell(const cell<Lyt>& c1, const cell<Lyt>& c2) const
    {
        auto index1 = cell_to_index(c1);
        auto index2 = cell_to_index(c2);

        if (strg->dist_mat[index1][index2] == 0)
        {
            return 0.0;
        }
        return (strg->sim_params.k / strg->dist_mat[index1][index2] *
                std::exp(-strg->dist_mat[index1][index2] / strg->sim_params.lambda_tf) *
                physical_sim_constants::ELECTRIC_CHARGE);
    }

    /**
     * The function calculates the electrostatic potential for each SiDB position (local).
     */
    void local_potential()
    {
        strg->loc_pot.resize(this->num_cells(), 0);
        for (uint64_t i = 0u; i < strg->sidb_order.size(); i++)
        {
            double collect = 0;
            for (uint64_t j = 0u; j < strg->sidb_order.size(); j++)
            {
                collect += strg->pot_mat[i][j] * transform_to_sign(strg->cell_charge[j]);
            }
            strg->loc_pot[i] = collect;
        }
    }

    std::optional<double> get_loc_pot_cell(const cell<Lyt>& c1)
    {
        if (auto index = cell_to_index(c1); index != -1)
        {
            return strg->loc_pot[index];
        }
        return std::nullopt;
    };

    std::optional<double> get_loc_pot_index(const uint64_t& c1)
    {
        if (c1 < strg->sidb_order.size())
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
        for (uint64_t i = 0; i < strg->loc_pot.size(); i++)
        {
            total_energy += 0.5 * strg->loc_pot[i] * transform_to_sign(strg->cell_charge[i]);
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
     */
    void validity_check()
    {
        uint64_t flag_fd      = 0;
        uint64_t counter_loop = 0;
        for (auto& it : strg->loc_pot)
        {
            bool valid = (((strg->cell_charge[counter_loop] == sidb_charge_state::NEGATIVE) &&
                           ((-it + strg->phys_params.mu) < physical_sim_constants::POP_STABILITY_ERR)) ||
                          ((strg->cell_charge[counter_loop] == sidb_charge_state::POSITIVE) &&
                           ((-it + strg->phys_params.mu_p) > -physical_sim_constants::POP_STABILITY_ERR)) ||
                          ((strg->cell_charge[counter_loop] == sidb_charge_state::NEUTRAL) &&
                           ((-it + strg->phys_params.mu) > -physical_sim_constants::POP_STABILITY_ERR) &&
                           (-it + strg->phys_params.mu_p) < physical_sim_constants::POP_STABILITY_ERR));
            counter_loop += 1;
            if (!valid)
            {
                strg->validity = false;
                flag_fd += 1;
                break;
            }
        }

        if (flag_fd == 0)
        {

            auto hopDel = [this](const uint64_t& c1, const uint64_t& c2)
            {
                int dn_i = (strg->cell_charge[c1] == sidb_charge_state::NEGATIVE) ? 1 : -1;
                int dn_j = -dn_i;

                return strg->loc_pot[c1] * dn_i + strg->loc_pot[c2] * dn_j - strg->pot_mat[c1][c2] * 1;
            };

            uint64_t flag_value  = 0;
            uint64_t counter_new = 0;
            for (uint64_t i = 0u; i < strg->loc_pot.size(); i++)
            {
                if (strg->cell_charge[counter_new] == sidb_charge_state::POSITIVE)
                {
                    continue;
                }

                for (uint64_t j = 0u; j < strg->loc_pot.size(); j++)
                {
                    if (flag_value == 1)
                    {
                        break;
                    }

                    auto E_del = hopDel(i, j);
                    if ((transform_to_sign(strg->cell_charge[j]) > transform_to_sign(strg->cell_charge[i])) &&
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

    /**
     * Returns the validity of the present charge distribution layout.
     *
     * @param strg A pointer to a storage object.
     * @returns bool The validity charge distribution.
     */
    [[nodiscard]] bool get_validity()
    {
        return strg->validity;
    }

    /**
     * The charge distribution of the charge distribution surface is converted to a unique index. It is used to map
     * every possible charge distribution of a SiDB layout to an index.
     *
     * @tparam Lyt cell-level layout.
     * @tparam base number of charge states per SiDB. Base = 2 when only neutrally and negatively charged SiDBs are
     * taken into account, otherwise base = 3.
     * @return a pair with the calculated charge distribution's corresponding charge index and the chosen base as the
     * first and second element, respectively.
     */
    void chargeconf_to_index() const
    {
        uint8_t base = strg->phys_params.base;
        assert(base == 2 || base == 3 && "base must be 2 or 3");
        uint64_t chargeindex = 0;
        uint64_t counter     = 0;
        for (const auto& c : strg->cell_charge)
        {
            chargeindex +=
                static_cast<uint64_t>((transform_to_sign(c) + 1) * std::pow(base, this->num_cells() - counter - 1));
            counter += 1;
        }
        strg->charge_index = {chargeindex, base};
    }

    /**
     * The charge index of the current charge distribution is returned.
     *
     * @return a pair with the charge index and the used base is returned.
     */
    [[nodiscard]] std::pair<uint64_t, uint8_t> get_charge_index() const
    {
        return strg->charge_index;
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
        auto charge_quot = strg->charge_index.first;
        auto base        = strg->charge_index.second;
        auto num_charges = this->num_cells() - 1;
        auto counter     = num_charges;

        while (charge_quot > 0)
        {
            div_t d;
            d           = div(charge_quot, base);
            charge_quot = static_cast<uint64_t>(d.quot);

            this->assign_charge_state_index(counter, sign_to_label(d.rem - 1));
            counter -= 1;
        }
    }

    /**
     * The unique index is increased by one, but only if it is less than or equal to the maximum charge index after it
     * has been increased. If that's the case, it is increased and afterward, the charge configuration is updated by
     * invoking the index_to_chargeconf() function.
     *
     */
    void increase_charge_index()
    {
        if (strg->charge_index.first == strg->max_charge_index)
        {
            std::cout << "maximal index is exceeded \n";
        }
        else
        {
            strg->charge_index.first += 1;
            this->index_to_chargeconf();
        }
    }

    /**
     * Returns the maximum index of the cell-level layout.
     *
     * @returns uint64_t The maximal index.
     */
    uint64_t get_max_charge_index()
    {
        return strg->max_charge_index;
    }

    /**
     * The funcion assigns a certain charge state to the layout and charge distribution is updated correspondingly.
     *
     * @returns uint64_t The maximal index.
     */
    void assign_charge_index(const uint64_t& index)
    {
        strg->charge_index.first = index;
        this->index_to_chargeconf();
    }

    /**
     * This method is used for the quicksim algorithm (see quicksim.hpp).
     * It gets a vector with indices representing negatively charged SiDBs as input. Afterward, a distant and neutrally
     * charged SiDB is localized using a min-max diversity algorithm. This selected SiDB is set to negativ and the index
     * is added to the input vector such that the next iteration works correctly.
     *
     * @param alpha a parameter for the algorithm (default: 0.7).
     * @param index_db vector of cell indices that are already negatively charged.
     */
    void adjacent_search(const double alpha, std::vector<uint64_t>& index_db)
    {
        uint64_t              count    = 0;
        double                dist_max = 0;
        int                   coord    = -1;
        std::vector<uint64_t> index_vector{};
        std::vector<double>   distance{};
        auto                  reservesize = this->num_cells() - index_db.size();
        index_vector.reserve(reservesize);
        distance.reserve(reservesize);

        for (uint64_t unocc = 0u; unocc < strg->cell_charge.size(); unocc++)
        {
            if (strg->cell_charge[unocc] != sidb_charge_state::NEUTRAL)
            {
                continue;
            }
            count += 1;

            auto dist_min = std::accumulate(index_db.begin(), index_db.end(), std::numeric_limits<double>::max(),
                                            [&](double acc, uint64_t occ)
                                            { return std::min(acc, this->get_distance_index(unocc, occ)); });

            index_vector.push_back(unocc);
            distance.push_back(dist_min);

            if (dist_min > dist_max)
            {
                dist_max = dist_min;
                coord    = static_cast<int>(unocc);
            }

            else if ((dist_min == dist_max) && sum_distance(index_db, coord, unocc))
            {
                dist_max = dist_min;
                coord    = static_cast<int>(unocc);
            }
        }

        std::vector<uint64_t> candidates{};
        candidates.reserve(reservesize);

        for (uint64_t i = 0u; i < distance.size(); i++)
        {
            if (distance[i] >= (alpha * dist_max))
            {
                candidates.push_back(i);
            }
        }

        if (!candidates.empty())
        {
            auto random_index = static_cast<uint64_t>(rand() % candidates.size());

            auto random_element               = index_vector[candidates[random_index]];
            strg->cell_charge[random_element] = sidb_charge_state::NEGATIVE;
            index_db.push_back(random_element);
            strg->system_energy += -this->get_loc_pot_index(random_element).value();

            for (uint64_t i = 0u; i < strg->pot_mat.size(); i++)
            {
                strg->loc_pot[i] += -this->get_potential_index(i, random_element);
            }
        }
    }

    /**
     * Compares the sum of distances between elements of a given input vector and two indices (old and now) in a
     * distance matrix.
     * @param input A vector of integers representing the elements to compare the distances with.
     * @param old An integer representing the old index to compare the distances with.
     * @param now An integer representing the now index to compare the distances with.
     * @returns bool A boolean value indicating if the sum of distances with "now index" is greater than the sum of
     * distances with "old index".
     */
    [[nodiscard]] bool sum_distance(const std::vector<uint64_t>& input, const uint64_t& old, const uint64_t& now) const
    {
        double sum_old = 0;
        double sum_now = 0;

        for (uint64_t j = 0; j < input.size(); j++)
        {
            sum_old += strg->dist_mat[j][old];
            sum_now += strg->dist_mat[j][now];
        }
        return sum_now > sum_old;
    };

  private:
    storage strg;
};

template <class T>
charge_distribution_surface(const T&) -> charge_distribution_surface<T>;

template <class T>
charge_distribution_surface(const T&, const physical_params&, const sidb_charge_state& cs)
    -> charge_distribution_surface<T>;

template <class T>
charge_distribution_surface(const T&, const physical_params&) -> charge_distribution_surface<T>;

}  // namespace fiction

#endif  // CHARGE_DISTRIBUTION_SURFACE
