//
// Created by Jan Drewniok on 23.11.22.
//

#ifndef FICTION_CHARGE_DISTRIBUTION_SURFACE_HPP
#define FICTION_CHARGE_DISTRIBUTION_SURFACE_HPP

#include "fiction/algorithms/simulation_sidb/simulation_parameters.hpp"
#include "fiction/layouts/cell_level_layout.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/types.hpp"

#include <cassert>
#include <cstdint>
#include <limits>
#include <random>

namespace fiction
{

/**
 * A layout type to layer on top of any SiDB cell-level layout. It implements an interface to store and access
 * SiDBs' charge states.
 *
 * @tparam Lyt cell-level layout based in SiQAD-coordinates.
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
         * The potential matrix is a vector of vectors storing the euclidean distance.
         */
        using potential_matrix = std::vector<std::vector<double>>;

        /**
         * It is a vector that stores the electrostatic potential as value.
         */
        using local_potential = std::vector<double>;

      public:
        explicit charge_distribution_storage(const physical_params& sim_param_default = physical_params{}) :
                phys_params{sim_param_default} {};

        physical_params phys_params{};        // all physical parameters used for the simulation are stored in a struct.
        std::vector<typename Lyt::cell> sidb_order{};  // all cells that are occupied by an SiDB are stored in this vector.
        std::vector<sidb_charge_state>
            cell_charge{};  // the SiDBs' charge states are stored. Corresponding cells are stored in "sidb_order".
        distance_matrix dist_mat{};  // distance between SiDBs are stored as matrix.
        potential_matrix
            pot_mat{};  // electrostatic potential between SiDBs are stored as matrix (here, still charge-independent).
        local_potential loc_pot{};  // electrostatic potential at each SiDB postion (has to be updated when charge
                                    // distribution is changed).
        double system_energy{0.0};  // electrostatic energy of a given charge distribution.
        bool   validity = false;    // labels if given charge distribution is physically valid (compare Paper XXX).
        std::pair<uint64_t, uint8_t>
            charge_index{};  // each charge distribution is assigned a unique index (first entry of pair), second one
                             // stores the base number (2 or 3 state simulation).
        uint64_t max_charge_index{};  // depending on the number of SiDBs and the base number, a maximal number of
                                      // possible charge distributions exists.
    };

    using storage = std::shared_ptr<charge_distribution_storage>;

    /**
     * Standard constructor for empty layouts.
     *
     * @param sidb_charge_state charge state used for the initialization of all SiDBs, default is negative charge for
     * all SiDBs.
     * @param physical_params physical parameters used for the simulation (µ_minus, base number, ...).
     */
    explicit charge_distribution_surface(const physical_params&   sim_param_default = physical_params{},
                                         const sidb_charge_state& cs                = sidb_charge_state::NEGATIVE) :
            Lyt(),
            strg{std::make_shared<charge_distribution_storage>(sim_param_default)}
    {
        static_assert(std::is_same_v<typename Lyt::cell, siqad::coord_t>, "Lyt is not based on SiQAD coordinates");
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
        initialize(cs);
    }
    /**
     * Standard constructor for existing layouts.
     *
     * @param sidb_charge_state charge state used for the initialization of all SiDBs, default is negative charge.
     * @param physical_params physical parameters used for the simulation (µ_minus, base number, ...).
     */
    explicit charge_distribution_surface(const Lyt& lyt, const physical_params& sim_param_default = physical_params{},
                                         const sidb_charge_state& cs = sidb_charge_state::NEGATIVE) :
            Lyt(lyt),
            strg{std::make_shared<charge_distribution_storage>(sim_param_default)}
    {
        static_assert(std::is_same_v<typename Lyt::cell, siqad::coord_t>, "Lyt is not based on SiQAD coordinates");
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
        initialize(cs);
    };

    /**
     * Copy constructor.
     *
     * @param lyt charge_distribution_surface that wants to be copied.
     */
    explicit charge_distribution_surface(const charge_distribution_surface<Lyt>& lyt) :
            strg{std::make_shared<charge_distribution_storage>(*lyt.strg)}
    {}

    /**
     * Initialisation function used for the construction of the charge distribution surface.
     *
     * @param cs charge state assigned to all SiDBs.
     */
    void initialize(const sidb_charge_state& cs = sidb_charge_state::NEGATIVE)
    {
        strg->sidb_order.reserve(this->num_cells());
        strg->cell_charge.reserve(this->num_cells());
        this->foreach_cell([this](const auto& c1) { strg->sidb_order.push_back(c1); });
        this->foreach_cell([this, &cs](const auto&) { strg->cell_charge.push_back(cs); });
        assert(((this->num_cells() < 41) && (strg->phys_params.base == 3)) && "number of SiDBs is too large");
        assert(((strg->phys_params.base == 3) && (this->num_cells() < 64)) && "number of SiDBs is too large");
        this->chargeconf_to_index();
        this->initialize_distance_matrix();
        this->initialize_potential_matrix();
        strg->max_charge_index =
            static_cast<uint64_t>(std::pow(static_cast<double>(strg->phys_params.base), this->num_cells()) - 1);
        this->local_potential();
        this->system_energy();
        this->validity_check();
    };

    /**
     * Set the physical parameters for the simulation.
     *
     * @param sim_param Physical parameters to be set.
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

    /**
     * Check if any SiDB is eqaul the given charge state.
     *
     * @param cs charge state for which it is checked whether any SiDB has it.
     */
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
     */
    void state_index(const uint64_t& i, const sidb_charge_state& cs) const
    {
        strg->cell_charge[i] = cs;
    }

    /**
     * This function assigns the given charge state to the given cell of the layout.
     *
     * @param cell The cell to which a state of charge is to be assigned.
     * @param cs The charge state to be assigned to the cell.
     */
    void assign_charge_state_cell(const typename Lyt::cell& cell, const sidb_charge_state& cs) const
    {
        auto index = cell_to_index(cell);
        if (index != -1)
        {
            strg->cell_charge[static_cast<uint64_t>(index)] = cs;
        }
        else
        {
            std::cout << "not part of the charge distribution surface" << std::endl;
        }
    }

    /**
     * This function assigns the given charge state to the cell (accessed by the input index) of the layout.
     *
     * @param index The index of the cell to which a state of charge is to be assigned.
     * @param cs The charge state to be assigned to the cell.
     */
    void assign_charge_state_index(const uint64_t& index, const sidb_charge_state& cs) const
    {
        strg->cell_charge[index] = cs;
    }

    /**
     * Sets the charge state of all SiDBs in the layout to a given charge state.
     *
     * @param cs The charge state to be assigned to all the SiDBs.
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
     * Returns the charge state of a given cell.
     *
     * @param c cell.
     * @return The charge state of the given cell.
     */
    [[nodiscard]] sidb_charge_state get_charge_state_cell(const typename Lyt::cell& c) const noexcept
    {
        if (auto index = cell_to_index(c); index != -1)
        {
            return strg->cell_charge[static_cast<uint64_t>(index)];
        }
        return sidb_charge_state::NONE;
    }

    /**
     * Initializes the distance matrix between all the cells of the layout.
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
    std::pair<double, double> nm_position(const typename Lyt::cell& c) const
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
    [[nodiscard]] constexpr double distance_sidb_pair(const typename Lyt::cell& c1, const typename Lyt::cell& c2) const
    {
        // @Marcel: I checked if I can use the euclidean_distance (see distance.hpp) function, but as it is now, I guess I cannot use it.
        const auto pos_c1 = nm_position(c1);
        const auto pos_c2 = nm_position(c2);
        const auto x      = static_cast<double>(pos_c1.first) - static_cast<double>(pos_c2.first);
        const auto y      = static_cast<double>(pos_c1.second) - static_cast<double>(pos_c2.second);
        return std::hypot(x, y);
    }

    /**
     * Finds the index of an SiDB.
     *
     * @param cell The cell to find the index of.
     * @return The index of the cell in the layout. Returns -1 if the cell is not part of the layout.
     */
    [[nodiscard]] constexpr int64_t cell_to_index(const typename Lyt::cell& cell) const
    {
        if (auto it = std::find(strg->sidb_order.begin(), strg->sidb_order.end(), cell); it != strg->sidb_order.end())
        {
            return static_cast<int64_t>(std::distance(strg->sidb_order.begin(), it));
        }
        return -1;
    }

    /**
     *  Returns the distance between two cells
     *
     *  @param cell1 the first cell to compare
     *  @param cell2 the second cell to compare
     *  @return a constexpr double representing the distance between the two cells
     */
    [[nodiscard]] constexpr double get_distance_cell(const typename Lyt::cell& cell1, const typename Lyt::cell& cell2)
    {
        auto index1 = cell_to_index(cell1);
        if (auto index2 = cell_to_index(cell2); (index1 != -1) && (index2 != -1))
        {
            return strg->dist_mat[static_cast<uint64_t>(index1)][static_cast<uint64_t>(index2)];
        }
        return 0;
    }

    /**
     * Calculates and returns the distance between two input indices.
     *
     * @param input1 The first input index.
     * @param input2 The second input index.
     * @return The distance index between index1 and index2 (indices correspond to unique SiDBs).
     */
    [[nodiscard]] constexpr double get_distance_index(const uint64_t& input1, const uint64_t& input2)
    {
        return strg->dist_mat[input1][input2];
    }

    /**
     * Calculates and returns the electrostatic potential between two cells.
     *
     * @param input1 The first cell
     * @param input2 The second cell
     * @return The potential between input1 and input2
     */
    [[nodiscard]] constexpr double get_potential_cell(const typename Lyt::cell& cell1, const typename Lyt::cell& cell2)
    {
        auto index1 = cell_to_index(cell1);
        if (auto index2 = cell_to_index(cell2); (index1 != -1) && (index2 != -1))
        {
            return strg->pot_mat[static_cast<uint64_t>(index1)][static_cast<uint64_t>(index2)];
        }
        return 0;
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
     * The electrostatic potential between two cells (SiDBs) is calculated.
     *
     * @param input1 The first index.
     * @param input2 The second index.
     * @return The potential between input1 and input2.
     */
    [[nodiscard]] double potential_sidb_pair_index(const uint64_t& index1, const uint64_t& index2) const
    {
        if (strg->dist_mat[index1][index2] == 0)
        {
            return 0.0;
        }

        return (strg->phys_params.k / strg->dist_mat[index1][index2] *
                std::exp(-strg->dist_mat[index1][index2] / strg->phys_params.lambda_tf) *
                physical_sim_constants::ELECTRIC_CHARGE);
    }

    /**
     * Calculates and returns the potential of a pair of cells based on their distance and simulation parameters.
     *
     * @param c1 The first cell.
     * @param c2 The second cell.
     * @return The potential between c1 and c2.
     */
    [[nodiscard]] double potential_sidb_pair_cell(const typename Lyt::cell& c1, const typename Lyt::cell& c2) const
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
                collect += strg->pot_mat[i][j] * static_cast<double>(transform_to_sign(strg->cell_charge[j]));
            }
            strg->loc_pot[i] = collect;
        }
    }

    /**
     * The function returns the local electrostatic potential at a given SiDB position.
     *
     * @param cell.
     * @return local potential at given cell position. If there is no SiDB at the given cell, the null pointer is
     * returned.
     */
    std::optional<double> get_loc_pot_cell(const typename Lyt::cell& c1)
    {
        if (auto index = cell_to_index(c1); index != -1)
        {
            return strg->loc_pot[static_cast<uint64_t>(index)];
        }
        return std::nullopt;
    };

    /**
     * The function returns the local electrostatic potential at a given index position.
     *
     * @param cell.
     * @return local potential at given index position. If there is no SiDB at the given index (which corresponds to a
     * unique cell), the null pointer is returned.
     */
    std::optional<double> get_loc_pot_index(const uint64_t& c1)
    {
        if (c1 < strg->sidb_order.size())
        {
            return strg->loc_pot[c1];
        }
        return std::nullopt;
    };

    /**
     * Calculates the system's total electrostatic potential energy.
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
     * The physically validity of the current charge distribution is evaluated and stored in the storage struct. A
     * charge distribution is valid if the *Population Stability* and the *Configuration Stability* is fulfilled.
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

            auto hop_del = [this](const uint64_t& c1, const uint64_t& c2)
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

                    auto e_del = hop_del(i, j);
                    if ((transform_to_sign(strg->cell_charge[j]) > transform_to_sign(strg->cell_charge[i])) &&
                        (e_del < -physical_sim_constants::POP_STABILITY_ERR))
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
     * @returns The validity of the present charge distribution.
     */
    [[nodiscard]] bool get_validity()
    {
        return strg->validity;
    }

    /**
     * The charge distribution of the charge distribution surface is converted to a unique index. It is used to map
     * every possible charge distribution of an SiDB layout to a unique index.
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
            d           = div(static_cast<int>(charge_quot), static_cast<int>(base));
            charge_quot = static_cast<uint64_t>(d.quot);

            this->assign_charge_state_index(counter, sign_to_label(d.rem - 1));
            counter -= 1;
        }
    }

    /**
     * The unique index is increased by one, but only if it is less than or equal to the maximum charge index after it
     has been increased. If that's the case, it is increased and afterward, the charge configuration is updated by
     invoking the index_to_chargeconf() function.
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
     * @returns The maximal possible charge distribution index.
     */
    uint64_t get_max_charge_index()
    {
        return strg->max_charge_index;
    }

    /**
     * Assigns a certain charge state to a given index (which corresponds to a certain SiDB) and the charge distribution
     * is updated correspondingly.
     */
    void assign_charge_index(const uint64_t& index)
    {
        strg->charge_index.first = index;
        this->index_to_chargeconf();
    }

    /**
     * This function is used for the quicksim algorithm (see quicksim.hpp).
     It gets a vector with indices representing negatively charged SiDBs as input. Afterward, a distant and a neutrally
     charged SiDB is localized using a min-max diversity algorithm. This selected SiDB is set to "negativ" and the index
     is added to the input vector such that the next iteration works correctly.
     *
     * @param alpha a parameter for the algorithm (default: 0.7).
     * @param index_db vector of SiDBs indices that are already negatively charged (double occupied).
     */
    void adjacent_search(const double alpha, std::vector<uint64_t>& index_db)
    {
        double                dist_max = 0;
        int                   coord    = -1;
        std::vector<uint64_t> index_vector{};
        std::vector<double>   distance{};
        auto                  reserve_size = this->num_cells() - index_db.size();
        index_vector.reserve(reserve_size);
        distance.reserve(reserve_size);

        for (uint64_t unocc = 0u; unocc < strg->cell_charge.size(); unocc++)
        {
            if (strg->cell_charge[unocc] != sidb_charge_state::NEUTRAL)
            {
                continue;
            }

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
        }

        std::vector<uint64_t> candidates{};
        candidates.reserve(reserve_size);

        for (uint64_t i = 0u; i < distance.size(); i++)
        {
            if (distance[i] >= (alpha * dist_max))
            {
                candidates.push_back(i);
            }
        }

        if (!candidates.empty())
        {

            auto random_index =
                static_cast<uint64_t>(rand()) %
                candidates.size();  // Yes, it is wright what clang-tidy says, but the other way is much slower.

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
