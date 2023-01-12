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

        physical_params              phys_params{};
        std::vector<cell<Lyt>>         cells{};
        std::vector<sidb_charge_state> cell_charge{};
        distance_matrix                dist_mat{};
        potential_matrix               pot_mat{};
        local_potential                loc_pot{};
        double                         system_energy{0.0};
        bool                           validity = false;
        std::pair<uint64_t, uint8_t>   charge_index{};
        uint64_t                       max_charge_index{};
    };

    using storage = std::shared_ptr<charge_distribution_storage>;

    /**
     * Standard constructor for empty layouts.
     */
    explicit charge_distribution_surface(const sidb_charge_state& cs                = sidb_charge_state::NEGATIVE,
                                         const physical_params& sim_param_default = physical_params{}) :
            Lyt(),
            strg{std::make_shared<charge_distribution_storage>(sim_param_default)}
    {
        initialize(cs);
        static_assert(std::is_same_v<cell<Lyt>, siqad::coord_t>, "Lyt is not based on SiQAD coordinates");
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
    }
    /**
     * Standard constructor for existing layouts and simulation parameter as input.
     */
    explicit charge_distribution_surface(const Lyt&               lyt,
                                         const physical_params& sim_param_default = physical_params{},
                                         const sidb_charge_state& cs                = sidb_charge_state::NEGATIVE) :
            Lyt(lyt),
            strg{std::make_shared<charge_distribution_storage>(sim_param_default)}
    {
        initialize(cs);
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

    void initialize(const sidb_charge_state& cs = sidb_charge_state::NEGATIVE)
    {
        strg->cells.reserve(this->num_cells());
        strg->cell_charge.reserve(this->num_cells());
        this->foreach_cell([this](const auto& c1) { strg->cells.push_back(c1); });
        this->foreach_cell([this, &cs](const auto& c1) { strg->cell_charge.push_back(cs); });
        assert((std::pow(strg->phys_params.base, this->num_cells()) - 1) < std::numeric_limits<uint64_t>::max() &&
               "number of SiDBs is too large");
        strg->max_charge_index = std::pow(strg->phys_params.base, this->num_cells()) - 1;
        this->chargeconf_to_index();
        this->initialize_distance_matrix();
        this->initialize_potential_matrix();
        this->local_potential();
        this->system_energy();
        this->validity_check();
    };

    bool operator==(const charge_distribution_surface<Lyt>& other) const {
        return this->get_charge_index().first == other->get_charge_index().first;
    }

    /**
     *
     * @brief Assigns a certain charge state to a particular cell (assigned via an index) of the layout.
     *
     * @param i The index of the cell.
     * @param cs The charge state to be assigned to the cell.
     *
     * This function assigns the given charge state to the cell of the layout at the specified index.
     * It updates the `cell_charge` member of `strg` object with the new charge state of the specified cell.
     */

    void
    state_index(const uint64_t& i, const sidb_charge_state& cs) const
    {
        strg->cell_charge[i] = cs;
    }

    /**
     * @brief Assigns the charge state of a particular cell of the layout.
     *
     * @param cell The cell whose charge state is to be assigned.
     * @param cs The charge state to be assigned to the cell.
     *
     * This function assigns the given charge state to the cell of the layout.
     * It first uses the `cell_to_index` function to convert the cell to its corresponding index
     * and then updates the `cell_charge` member of `strg` object with the new charge state of the specified cell.
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

    void assign_charge_state_index(const uint64_t& index, const sidb_charge_state& cs) const
    {
        strg->cell_charge[index] = cs;
    }
    /**
     * @brief Sets the charge state of all NOT "EMPTY" cells in the layout to a given state.
     *
     * @param cs The charge state to be assigned to all the cells.
     *
     * This function sets the charge state of all the cells in the layout to the given state.
     * It iterates over the `cell_charge` member of `strg` object, which is a vector
     * of charge states for each cell, and assigns the given charge state to each element.
     */

    void set_charge_states(const sidb_charge_state& cs) const
    {
        for (int i = 0u; i < strg->cell_charge.size(); i++)
        {
            strg->cell_charge[i] = cs;
        }
    }

    /**
     * @brief Returns the charge state of a cell of the layout at a given index.
     *
     * @param index The index of the cell.
     * @return The charge state of the cell at the given index.
     *
     * This function returns the charge state of a cell of the layout at a given index.
     * It checks if the index is valid by comparing it with the size of the `cell_charge` vector.
     * If the index is valid, it accesses the `cell_charge` member of `strg` object and
     * returns the charge state at the given index.
     * If the index is invalid, it returns NONE.
     * This function is marked with [[nodiscard]] and noexcept, indicating that it is safe to use and
     * does not throw exceptions.
     */

    [[nodiscard]] sidb_charge_state get_charge_state_index(const uint64_t &index) const noexcept
    {
        if (index < (strg->cell_charge.size()))
        {
            return strg->cell_charge[index];
        }
        else
        {
            return sidb_charge_state::NONE;
        }
    }

    /**
     * @brief Returns the charge state of a cell of the layout at a given coordinate.
     *
     * @param c The coordinate of the cell.
     * @return The charge state of the cell at the given coordinate.
     *
     * This function returns the charge state of a cell of the layout at a given coordinate.
     * It uses the `cell_to_index` function to convert the coordinate of the cell to its corresponding index.
     * Then it checks if the index is valid and if it is, it accesses the `cell_charge` member of `strg` object and
     * returns the charge state at the given index.
     * If the index is invalid, it returns NONE.
     * This function is marked with [[nodiscard]] and noexcept, indicating that it is safe to use and
     * does not throw exceptions.
     */

    [[nodiscard]] sidb_charge_state get_charge_state_cell(const coordinate<Lyt>& c) const noexcept
    {
        auto index = cell_to_index(c);
        if (index != -1)
        {
            return strg->cell_charge[index];
        }
        else
        {
            return sidb_charge_state::NONE;
        }
    }

    /**
     * @brief Applies a given function to each charge state in the layout.
     *
     * @param fn A function to be applied to each charge state.
     *
     * This function applies a given function to each charge state in the layout.
     * It uses mockturtle::detail::foreach_element function to iterate over the elements of the
     * `cell_charge` member of the `strg` object and applies the given function `fn` to each element.
     * This allows to perform an operation on all the elements, it can be useful if you want to apply
     * a certain operation to all charge states.
     *
     * @tparam Fn Type of the function object, should be callable with the signature
     * void(const sidb_charge_state&)
     */
    template <typename Fn>
    void foreach_charge_state(Fn&& fn) const
    {
        mockturtle::detail::foreach_element(strg->cell_charge.cbegin(), strg->cell_charge.cend(),
                                            std::forward<Fn>(fn));
    }

    /**
     * @brief Initializes the distance matrix between all the cells of the layout.
     *
     * This function initializes the distance matrix between all the cells of the layout.
     * First, it creates an empty 2D vector (strg->dist_mat) of the size of the number of cells with all elements set to
     * 0. Then, it iterates over the `cells` member of `strg` object, which is a vector of cell coordinates, for each
     * pair of cells, it computes the distance between them using the `distance_sidb_pair` function and stores the
     * distance in the corresponding position of the distance matrix.
     *
     */

    void initialize_distance_matrix() const
    {
        strg->dist_mat = std::vector<std::vector<double>>(this->num_cells(), std::vector<double>(this->num_cells(), 0));

        for (int i = 0u; i < strg->cells.size(); i++)
        {
            for (int j = 0u; j < strg->cells.size(); j++)
            {
                strg->dist_mat[i][j] = distance_sidb_pair(strg->cells[i], strg->cells[j]);
            }
        }
    };

    /**
     * @brief Computes the position of a cell in nanometers.
     *
     * @param c The cell to compute the position for.
     * @return A pair of double values representing the x,y position of the cell in nanometers.
     *
     * This function computes the position of a cell in nanometers.
     * It uses the x,y,z coordinates of the input cell `c` and the lattice constants stored in
     * the `sim_params` member of `strg` object to compute the position of the cell in nanometers.
     * It returns a pair of double values representing the x,y position of the cell in nanometers.
     */
    std::pair<double, double> nm_position(const cell<Lyt>& c) const
    {
        const auto x = static_cast<double>(c.x * strg->phys_params.lat_a);
        const auto y = static_cast<double>(c.y * strg->phys_params.lat_b + c.z * strg->phys_params.lat_c);
        return std::make_pair(x, y);
    }

    /**
     * @brief Initializes the potential matrix between all the cells of the layout.
     *
     * This function initializes the potential matrix between all the cells of the layout.
     * First, it creates an empty 2D vector (strg->pot_mat) of the size of the number of cells with all elements set to
     * 0. Then, it iterates over the `cells` member of `strg` object, which is a vector of cell coordinates, for each
     * pair of cells, it computes the potential between them using the `potential_sidb_pair_index` function and stores
     * the potential in the corresponding position of the potential matrix.
     *
     */

    void initialize_potential_matrix() const
    {
        strg->pot_mat = std::vector<std::vector<double>>(this->num_cells(), std::vector<double>(this->num_cells(), 0));
        for (int i = 0u; i < strg->cells.size(); i++)
        {
            for (int j = 0u; j < strg->cells.size(); j++)
            {
                strg->pot_mat[i][j] = potential_sidb_pair_index(i, j);
            }
        }
    };

    /**
     * @brief Computes the distance between two cells in nanometers.
     *
     * @param c1 The first cell.
     * @param c2 The second cell.
     * @return The distance between the two cells in nanometers.
     *
     * This function computes the distance between two cells in nanometers.
     * It uses the `nm_position` function to compute the positions of the cells in nanometers
     * and then it computes the Euclidean distance between the two positions using the hypot function.
     * The distance is returned as a double.
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
     *
     * This function finds the index of a cell in the layout.
     * It uses the std::find function to search for the input cell in the 'cells' member of the 'strg' object
     * which is a vector of cell coordinates. If the cell is found in the vector, the function returns the index
     * of the cell in the vector using the std::distance function. Otherwise, it returns -1 indicating that the
     * cell is not part of the layout.
     */

    [[nodiscard]] constexpr int64_t cell_to_index(const cell<Lyt>& cell) const
    {
        auto it = std::find(strg->cells.begin(), strg->cells.end(), cell);
        if (it != strg->cells.end())
        {
            std::cout << std::distance(strg->cells.begin(), it) << std::endl;
            return static_cast<int64_t>(std::distance(strg->cells.begin(), it));
        }
        else
        {
            return -1;
        }
    }

    /**
     *  @brief calculates the distance between two cells
     *  @param cell1 the first cell to compare
     *  @param cell2 the second cell to compare
     *  @return a constexpr double representing the distance between the two cells
     *
     *  This function uses the pre-calculated distance matrix stored in the strg pointer
     *  to look up the distance between the two cells passed in as arguments. The cell_to_index
     *  function is used to convert the cell objects into indices that can be used to look up
     *  the distance in the matrix.
     *  The [[nodiscard]] attribute indicate that the returned value should not be discarded by the caller.
     */

    [[nodiscard]] constexpr double get_distance_cell(const cell<Lyt>& cell1, const cell<Lyt>& cell2)
    {
        return strg->dist_mat[cell_to_index(cell1)][cell_to_index(cell2)];
    }

    /**
     * @brief Calculates and returns the distance index between two input values.
     *
     * This function calculates and returns the distance index between two input values, represented by input1 and
     * input2. The distance index is retrieved from a 2D distance matrix stored in a struct pointed to by the strg
     * pointer.
     *
     * @param input1 The first input value
     * @param input2 The second input value
     *
     * @return The distance index between input1 and input2
     */
    [[nodiscard]] constexpr double get_distance_index(const uint64_t& input1, const uint64_t& input2)
    {
        return strg->dist_mat[input1][input2];
    }

    /**
     * @brief Calculates and returns the potential of two cells.
     *
     * This function calculates and returns the potential of two cells, represented by input1 and input2.
     * The potential is retrieved from a 2D potential matrix stored in a struct pointed to by the strg pointer.
     * The cell_to_index function is used to convert the cells to corresponding indices in the potential matrix.
     *
     * @param input1 The first cell
     * @param input2 The second cell
     *
     * @return The potential between input1 and input2
     */

    [[nodiscard]] constexpr double get_potential_cell(const cell<Lyt>& input1, const cell<Lyt>& input2)
    {
        return strg->pot_mat[cell_to_index(input1)][cell_to_index(input2)];
    }

    /**
     * @brief Calculates and returns the potential of two indices.
     *
     * This function calculates and returns the potential of two indices, represented by input1 and input2.
     * The potential is retrieved from a 2D potential matrix stored in a struct pointed to by the strg pointer.
     *
     * @attribution The constexpr keyword indicates that the function's return value is a constant expression and can be
     * evaluated at compile-time if its inputs are known at compile-time.
     *
     * @param input1 The first index
     * @param input2 The second index
     *
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
    [[nodiscard]] double potential_sidb_pair_index(const int& c1, const int& c2) const
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
     * @brief Calculates and returns the potential of a pair of cells based on their distance and simulation parameters.
     *
     * This function calculates and returns the potential of a pair of cells, represented by c1 and c2.
     * The potential is calculated based on the distance between the cells and simulation parameters stored in a struct
     * pointed to by the strg pointer. The cell_to_index function is used to convert the cells to corresponding indices
     * in the distance matrix. If the distance between the cells is zero, the function returns 0.0.
     *
     * @param c1 The first cell
     * @param c2 The second cell
     *
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
    * The function calculates the electrostatic potential for each SiDB position (therefore local).
    */
    void local_potential()
    {
        strg->loc_pot.resize(this->num_cells(), 0);
        for (int i = 0u; i <  strg->cells.size(); i++)
        {
                double collect = 0;
                for (int j = 0u; j <  strg->cells.size(); j++)
                {
                    collect += strg->pot_mat[i][j] * transform_to_sign(strg->cell_charge[j]);
                }
                strg->loc_pot[i] = collect;
            }
    }

    std::optional<double> get_loc_pot_cell(const cell<Lyt>& c1)
    {
        auto index = cell_to_index(c1);
        if (index != -1)
        {
            return strg->loc_pot[index];
        }

        return std::nullopt;
    };

    std::optional<double> get_loc_pot_index(const int& c1)
    {
        if (c1 < strg->cells.size())
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
        for (int i = 0; i < strg->loc_pot.size(); i++)
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
     *
     */

    void validity_check()
    {
        int  flag_fd = 0;
        int counter_loop = 0;
        for (auto& it : strg->loc_pot)
        {
            bool valid = (((strg->cell_charge[counter_loop] == sidb_charge_state::NEGATIVE) &&
                           ((-it+ strg->phys_params.mu) < physical_sim_constants::POP_STABILITY_ERR)) ||
                          ((strg->cell_charge[counter_loop] == sidb_charge_state::POSITIVE) &&
                           ((-it + strg->phys_params.mu_p) > -physical_sim_constants::POP_STABILITY_ERR)) ||
                          ((strg->cell_charge[counter_loop] == sidb_charge_state::NEUTRAL) &&
                           ((-it + strg->phys_params.mu) > -physical_sim_constants::POP_STABILITY_ERR) &&
                           (-it + strg->phys_params.mu_p) < physical_sim_constants::POP_STABILITY_ERR));
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
                int dn_i = (strg->cell_charge[c1] == sidb_charge_state::NEGATIVE) ? 1 : -1;
                int dn_j = -dn_i;

                return strg->loc_pot[c1] * dn_i + strg->loc_pot[c2] * dn_j -
                       strg->pot_mat[c1][c2] * 1;
            };

            int flag_value = 0;
            int counter_new = 0;
            for (int i =0; i < strg->loc_pot.size(); i++)
            {
                if (strg->cell_charge[counter_new] == sidb_charge_state::POSITIVE)
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
                    if ((transform_to_sign(strg->cell_charge[j]) >
                         transform_to_sign(strg->cell_charge[i])) &&
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
    * @brief Returns the validity of the present charge distribution layout.
    * The function accesses the 'validity' member variable of the 'strg' pointer, which is a boolean value
    indicating whether the storage object is valid or not. The function returns the value of 'validity' without
    modifying it.
    * @param strg A pointer to a storage object.
    * @returns bool The validity of
    the storage object.
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
        int      counter     = 0;
        for (const auto& c : strg->cell_charge)
        {
            chargeindex += static_cast<unsigned int>((transform_to_sign(c) + 1) *
                                                     std::pow(base, this->num_cells() - counter - 1));
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
        int charge_quot = strg->charge_index.first;
        int base        = strg->charge_index.second;
        int num_charges = this->num_cells() - 1;
        int counter = num_charges;

        while (charge_quot > 0)
        {
            div_t d;
            d           = div(charge_quot, base);
            charge_quot = d.quot;

                    this->assign_charge_state_index(counter, sign_to_label(d.rem - 1));
                    counter-= 1;
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

    /**
    * @brief Returns the maximum index of the cell-level layout.
    * The function accesses the 'max_charge_index' member variable of the 'strg' pointer, which is an unsigned 64-bit
    integer. The function returns the value of 'max_charge_index' without modifying it.
    * @returns uint64_t The maximal index.
     */
    uint64_t get_max_charge_index()
    {
        return strg->max_charge_index;
    }

    /**
* The funcion assigns a certain charge state to the layout and charge distribution is updated correspondingly.
* @returns uint64_t The maximal index.
 */
    void assign_charge_index(const uint64_t& index)
    {
        strg->charge_index.first = index;
        this->index_to_chargeconf();
    }

    /**
     * This method is used for the faccusim algorithm (see new_approach.hpp).
     * It gets a vector with indices representing negatively charged SiDBs as input. Afterward, a distant and neutrlly charged SiDB is localized using a new method.
     * This selected SiDB is set to negativ and the index is added to the input vector such that the next iteration works correctly.
     * @param alpha
     * @param index_db
     */

    void next_N(const double alpha, std::vector<int>& index_db)
    {
        auto               count    = 0;
        float              dist_max = 0;
        auto               coord    = -1;
        std::vector<int>   index_vector{};
        std::vector<float> distance{};
        auto               reservesize = this->num_cells() - index_db.size();
        index_vector.reserve(reservesize);
        distance.reserve(reservesize);

        for (int unocc = 0u; unocc < strg->cell_charge.size(); unocc++)
        {
            if (strg->cell_charge[unocc] != sidb_charge_state::NEUTRAL)
            {
                continue;
            }
            count += 1;

            auto dist_min = MAXFLOAT;
            for (int occ : index_db)
            {
                if (this->get_distance_index(unocc, occ) < dist_min)
                {
                    dist_min = this->get_distance_index(unocc, occ);
                }
            }

            index_vector.push_back(unocc);
            distance.push_back(dist_min);

            if (dist_min > dist_max)
            {
                dist_max = dist_min;
                coord    = unocc;
            }

            else if (dist_min == dist_max)
            {
                if (sum_distance(index_db, coord, unocc))
                {
                    dist_max = dist_min;
                    coord    = unocc;
                }
            }
        }

        std::vector<int> candidates{};
        candidates.reserve(reservesize);

        for (int i = 0u; i < distance.size(); i++)
        {
            if (distance[i] >= (alpha * dist_max))
            {
                candidates.push_back(i);
            }
        }

        int random_index = rand() % candidates.size();

        int random_element                = index_vector[candidates[random_index]];
        strg->cell_charge[random_element] = sidb_charge_state::NEGATIVE;
        index_db.push_back(random_element);
        strg->system_energy += -this->get_loc_pot_index(random_element).value();

        for (int i = 0u; i < strg->pot_mat.size(); i++)
        {
            strg->loc_pot[i] += -this->get_potential_index(i, random_element);
        }
    }

    /**
    * Compares the sum of distances between elements of a given input vector and two indices (old and now) in a distance matrix.
    * @param input A vector of integers representing the elements to compare the distances with.
    * @param old An integer representing the old index to compare the distances with.
    * @param now An integer representing the now index to compare the distances with.
    * @returns bool A boolean value indicating if the sum of distances with now index is greater than the sum of distances with old index.
    */

    bool sum_distance(std::vector<int>& input, int& old, int& now)
    {
        float sum_old = 0;
        float sum_now = 0;
        for (int j = 0; j < input.size(); j++)
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
charge_distribution_surface(const T&, const physical_params&) -> charge_distribution_surface<T>;

}  // namespace fiction

#endif  // CHARGE_DISTRIBUTION_SURFACE
