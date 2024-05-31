//
// Created by Jan Drewniok on 11.09.23.
//

#ifndef FICTION_IS_OPERATIONAL_HPP
#define FICTION_IS_OPERATIONAL_HPP

#include "fiction/algorithms/iter/bdl_input_iterator.hpp"
#include "fiction/algorithms/simulation/sidb/can_positive_charges_occur.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_pairs.hpp"
#include "fiction/algorithms/simulation/sidb/energy_distribution.hpp"
#include "fiction/algorithms/simulation/sidb/exhaustive_ground_state_simulation.hpp"
#include "fiction/algorithms/simulation/sidb/quickexact.hpp"
#include "fiction/algorithms/simulation/sidb/quicksim.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_engine.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_result.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/traits.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_wires.hpp"

#include <kitty/bit_operations.hpp>
#include <kitty/traits.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

namespace fiction
{
/**
 * Possible operational status of a layout.
 */
enum class operational_status
{
    /**
     * The layout is operational.
     */
    OPERATIONAL,
    /**
     * The layout is non-operational.
     */
    NON_OPERATIONAL
};

enum class operation_ctiterium
{
    /**
     * The layout is operational.
     */
    ALLOW_KINKS,
    /**
     * The layout is non-operational.
     */
    FORBID_KINKS
};
/**
 * Parameters for the `is_operational` algorithm.
 */
struct is_operational_params
{
    /**
     * The simulation parameters for the physical simulation of the ground state.
     */
    sidb_simulation_parameters simulation_parameters{};
    /**
     * The simulation engine to be used for the operational domain computation.
     */
    sidb_simulation_engine sim_engine{sidb_simulation_engine::QUICKEXACT};
    /**
     * Parameters for the BDL pair detection algorithms.
     */
    detect_bdl_pairs_params bdl_params{};

    operation_ctiterium condition = operation_ctiterium::FORBID_KINKS;
};

namespace detail
{
/**
 * Implementation of the `is_operational` algorithm for a given gate layout.
 *
 * This class provides an implementation of the `is_operational` algorithm for
 * a specified gate layout and parameters. It checks whether the gate layout is operational
 * by simulating its behavior for different input combinations and comparing the results
 * to expected outputs from a truth table.
 *
 * @tparam Lyt SiDB cell-level layout type.
 * @tparam TT The type of the truth table specifying the gate behavior.
 */
template <typename Lyt, typename TT>
class is_operational_impl
{
  public:
    /**
     * Constructor to initialize the algorithm with a layout and parameters.
     *
     * @param lyt The SiDB cell-level layout to be checked.
     * @param spec Expected Boolean function of the layout given as a multi-output truth table.
     * @param params Parameters for the `is_operational` algorithm.
     */
    is_operational_impl(const Lyt& lyt, const Lyt& skeleton_lyt, const std::vector<TT>& tt,
                        const is_operational_params& params) :
            layout{lyt.clone()},
            skeleton{skeleton_lyt.clone()},
            truth_table{tt},
            parameters{params},
            output_bdl_pairs(detect_bdl_pairs(skeleton, sidb_technology::cell_type::OUTPUT, parameters.bdl_params)),
            bii(bdl_input_iterator<Lyt>{layout, skeleton, parameters.bdl_params}),
            input_wires{detect_bdl_wires(skeleton, WIRE::INPUT_WO_INPUT_BDL)},
            output_wires{detect_bdl_wires(skeleton, WIRE::OUTPUT)},
            all_wires{detect_bdl_wires(skeleton)},
            wire_directions{determine_wire_direction(detect_bdl_pairs(skeleton, sidb_technology::cell_type::INPUT),
                                                     detect_bdl_wires(skeleton, WIRE::INPUT_W_INPUT_BDL))}
    {}

    /**
     * Run the `is_operational` algorithm.
     *
     * This function executes the operational status checking algorithm for the gate layout
     * and parameters provided during initialization.
     *
     * @return The operational status of the gate layout (either `OPERATIONAL` or `NON_OPERATIONAL`).
     */
    [[nodiscard]] operational_status run() noexcept
    {
        assert(!output_bdl_pairs.empty() && "No output cell provided.");
        assert((truth_table.size() == output_bdl_pairs.size()) &&
               "Number of truth tables and output BDL pairs does not match");

        // number of different input combinations
        for (auto i = 0u; i < truth_table.front().num_bits(); ++i, ++bii)
        {
            ++simulator_invocations;

            // if positively charged SiDBs can occur, the SiDB layout is considered as non-operational
            if (can_positive_charges_occur(*bii, parameters.simulation_parameters))
            {
                return operational_status::NON_OPERATIONAL;
            }

            // performs physical simulation of a given SiDB layout at a given input combination
            const auto simulation_results = physical_simulation_of_layout(bii);

            // if no physically valid charge distributions were found, the layout is non-operational
            if (simulation_results.charge_distributions.empty())
            {
                return operational_status::NON_OPERATIONAL;
            }

            // find the ground state, which is the charge distribution with the lowest energy
            const auto ground_state = std::min_element(
                simulation_results.charge_distributions.cbegin(), simulation_results.charge_distributions.cend(),
                [](const auto& lhs, const auto& rhs) { return lhs.get_system_energy() < rhs.get_system_energy(); });

            // ground state is degenerate
            if ((energy_distribution(simulation_results.charge_distributions).begin()->second) > 1)
            {
                return operational_status::NON_OPERATIONAL;
            }

            // fetch the charge states of the output BDL pair
            for (auto output = 0u; output < output_bdl_pairs.size(); output++)
            {
                const auto charge_state_output_upper = ground_state->get_charge_state(output_bdl_pairs[output].upper);
                const auto charge_state_output_lower = ground_state->get_charge_state(output_bdl_pairs[output].lower);

                // if the output charge states are equal, the layout is not operational
                if (charge_state_output_lower == charge_state_output_upper)
                {
                    return operational_status::NON_OPERATIONAL;
                }

                // if the expected output is 1, the expected charge states are (upper, lower) = (0, -1)
                if (kitty::get_bit(truth_table[output], i))
                {
                    if (charge_state_output_upper != sidb_charge_state::NEUTRAL ||
                        charge_state_output_lower != sidb_charge_state::NEGATIVE)
                    {
                        return operational_status::NON_OPERATIONAL;
                    }
                }
                // if the expected output is 0, the expected charge states are (upper, lower) = (-1, 0)
                else
                {
                    if (charge_state_output_upper != sidb_charge_state::NEGATIVE ||
                        charge_state_output_lower != sidb_charge_state::NEUTRAL)
                    {
                        return operational_status::NON_OPERATIONAL;
                    }
                }
            }
        }

        // if we made it here, the layout is operational
        return operational_status::OPERATIONAL;
    }

    void set_charge_distribution_based_on_logic(charge_distribution_surface<Lyt>& layout,
                                                const uint64_t                    current_input_index)
    {
        layout.assign_all_charge_states(sidb_charge_state::NEGATIVE, false);

        for (auto i = 0u; i < input_wires.size(); i++)
        {
            if (wire_directions[input_wires.size() - 1 - i] == bdl_wire_direction::TOP_DOWN)
            {
                if ((current_input_index & (uint64_t{1ull} << i)) != 0ull)
                {
                    for (const auto& bdl : input_wires[input_wires.size() - 1 - i])
                    {
                        layout.assign_charge_state(bdl.upper, sidb_charge_state::NEUTRAL, false);
                        layout.assign_charge_state(bdl.lower, sidb_charge_state::NEGATIVE, false);
                    }
                }
                else
                {
                    for (const auto& bdl : input_wires[input_wires.size() - 1 - i])
                    {
                        layout.assign_charge_state(bdl.upper, sidb_charge_state::NEGATIVE, false);
                        layout.assign_charge_state(bdl.lower, sidb_charge_state::NEUTRAL, false);
                    }
                }
            }
            else if (wire_directions[input_wires.size() - 1 - i] == bdl_wire_direction::DOWN_TOP)
            {
                if ((current_input_index & (uint64_t{1ull} << i)) != 0ull)
                {
                    for (const auto& bdl : input_wires[input_wires.size() - 1 - i])
                    {
                        layout.assign_charge_state(bdl.upper, sidb_charge_state::NEGATIVE, false);
                        layout.assign_charge_state(bdl.lower, sidb_charge_state::NEUTRAL, false);
                    }
                }
                else
                {
                    for (const auto& bdl : input_wires[input_wires.size() - 1 - i])
                    {
                        layout.assign_charge_state(bdl.upper, sidb_charge_state::NEUTRAL, false);
                        layout.assign_charge_state(bdl.lower, sidb_charge_state::NEGATIVE, false);
                    }
                }
            }
        }

        for (auto i = 0u; i < output_wires.size(); i++)
        {
            for (const auto& bdl : output_wires[i])
            {
                if (kitty::get_bit(truth_table[i], current_input_index))
                {
                    layout.assign_charge_state(bdl.upper, sidb_charge_state::NEUTRAL, false);
                    layout.assign_charge_state(bdl.lower, sidb_charge_state::NEGATIVE, false);
                }
                else
                {
                    layout.assign_charge_state(bdl.upper, sidb_charge_state::NEGATIVE, false);
                    layout.assign_charge_state(bdl.lower, sidb_charge_state::NEUTRAL, false);
                }
            }
        }
    }

    [[nodiscard]] operational_status run_kinks_are_forbidden() noexcept
    {
        assert(!output_bdl_pairs.empty() && "No output cell provided.");
        assert((truth_table.size() == output_bdl_pairs.size()) &&
               "Number of truth tables and output BDL pairs does not match");

        // number of different input combinations
        for (auto i = 0u; i < truth_table.front().num_bits(); ++i, ++bii)
        {
            ++simulator_invocations;

            // if positively charged SiDBs can occur, the SiDB layout is considered as non-operational
            if (can_positive_charges_occur(*bii, parameters.simulation_parameters))
            {
                return operational_status::NON_OPERATIONAL;
            }

            // performs physical simulation of a given SiDB layout at a given input combination
            const auto simulation_results = physical_simulation_of_layout(bii);

            // if no physically valid charge distributions were found, the layout is non-operational
            if (simulation_results.charge_distributions.empty())
            {
                return operational_status::NON_OPERATIONAL;
            }

            // find the ground state, which is the charge distribution with the lowest energy
            const auto ground_state = std::min_element(
                simulation_results.charge_distributions.cbegin(), simulation_results.charge_distributions.cend(),
                [](const auto& lhs, const auto& rhs) { return lhs.get_system_energy() < rhs.get_system_energy(); });

            // ground state is degenerate
            if ((energy_distribution(simulation_results.charge_distributions).begin()->second) > 1)
            {
                return operational_status::NON_OPERATIONAL;
            }

            auto charge_layout = charge_distribution_surface<Lyt>{skeleton.clone()};
            set_charge_distribution_based_on_logic(charge_layout, i);

            //print_sidb_layout(std::cout, *ground_state);
            //print_sidb_layout(std::cout, charge_layout);

            bool non_operational = false;
            charge_layout.foreach_cell(
                [this, &non_operational, &charge_layout, &ground_state](const auto& c)
                {
                    if ((*bii).get_cell_type(c) != Lyt::technology::cell_type::INPUT &&
                        (*bii).get_cell_type(c) != Lyt::technology::cell_type::EMPTY)
                    {
                        if (ground_state->get_charge_state(c) != charge_layout.get_charge_state(c))
                        {
                            non_operational = true;
                            return;
                        }
                    }
                });
            if (non_operational)
            {
                //std::cout << "non operational" << std::endl;
                return operational_status::NON_OPERATIONAL;
            }
        }

        // if we made it here, the layout is operational
        return operational_status::OPERATIONAL;
    }
    /**
     * Determines the input combinations yielding the correct output.
     *
     * @return All inputs (e.g. 2-input Boolean function: 00 ^= 0; 10 ^= 2) for which the correct output is computed.
     */
    [[nodiscard]] std::set<uint64_t> determine_operational_input_patterns() noexcept
    {
        assert(!output_bdl_pairs.empty() && "No output cell provided.");
        assert((truth_table.size() == output_bdl_pairs.size()) &&
               "Number of truth tables and output BDL pairs does not match");

        std::set<uint64_t> operational_inputs{};

        // number of different input combinations
        for (auto i = 0u; i < truth_table.front().num_bits(); ++i, ++bii)
        {
            ++simulator_invocations;

            // if positively charged SiDBs can occur, the SiDB layout is considered as non-operational
            if (can_positive_charges_occur(*bii, parameters.simulation_parameters))
            {
                continue;
            }

            // performs physical simulation of a given SiDB layout at a given input combination
            const auto simulation_results = physical_simulation_of_layout(bii);

            // if no physically valid charge distributions were found, the layout is non-operational
            if (simulation_results.charge_distributions.empty())
            {
                continue;
            }

            // find the ground state, which is the charge distribution with the lowest energy
            const auto ground_state = std::min_element(
                simulation_results.charge_distributions.cbegin(), simulation_results.charge_distributions.cend(),
                [](const auto& lhs, const auto& rhs) { return lhs.get_system_energy() < rhs.get_system_energy(); });

            // ground state is degenerate
            if ((energy_distribution(simulation_results.charge_distributions).begin()->second) > 1)
            {
                continue;
            }

            bool correct_output = true;
            // fetch the charge states of the output BDL pair
            for (auto output = 0u; output < output_bdl_pairs.size(); output++)
            {
                auto charge_state_output_upper = ground_state->get_charge_state(output_bdl_pairs[output].upper);
                auto charge_state_output_lower = ground_state->get_charge_state(output_bdl_pairs[output].lower);

                // if the output charge states are equal, the layout is not operational
                if (charge_state_output_lower == charge_state_output_upper)
                {
                    correct_output = false;
                    break;
                }

                // if the expected output is 1, the expected charge states are (upper, lower) = (0, -1)
                if (kitty::get_bit(truth_table[output], i))
                {
                    if (charge_state_output_lower != sidb_charge_state::NEGATIVE ||
                        charge_state_output_upper != sidb_charge_state::NEUTRAL)
                    {
                        correct_output = false;
                    }
                }
                // if the expected output is 0, the expected charge states are (upper, lower) = (-1, 0)
                else
                {
                    if (charge_state_output_lower != sidb_charge_state::NEUTRAL ||
                        charge_state_output_upper != sidb_charge_state::NEGATIVE)
                    {
                        correct_output = false;
                    }
                }
            }
            if (correct_output)
            {
                operational_inputs.insert(i);
            }
        }

        // if we made it here, the layout is operational
        return operational_inputs;
    }
    /**
     * Returns the total number of simulator invocations.
     *
     * @return The number of simulator invocations.
     */
    [[nodiscard]] std::size_t get_number_of_simulator_invocations() const noexcept
    {
        return simulator_invocations;
    }

  private:
    /**
     * SiDB cell-level layout.
     */
    const Lyt layout;

    const Lyt skeleton;
    /**
     * The specification of the layout.
     */
    const std::vector<TT>&
        truth_table;  // TODO implement the matching of multi-input truth table inputs and BDL pair ordering
    /**
     * Parameters for the `is_operational` algorithm.
     */
    is_operational_params parameters;

    std::vector<bdl_wire<Lyt>>            input_wires;
    std::vector<bdl_wire<Lyt>>            output_wires;
    std::vector<bdl_wire<Lyt>>            all_wires;
    const std::vector<bdl_wire_direction> wire_directions;
    /**
     * Output BDL pairs.
     */
    std::vector<bdl_pair<Lyt>> output_bdl_pairs;
    /**
     * Iterator that iterates over all possible input states.
     */
    bdl_input_iterator<Lyt> bii;
    /**
     * Number of simulator invocations.
     */
    std::size_t simulator_invocations{0};
    /**
     * This function conducts physical simulation of the given layout (gate layout with certain input combination). The
     * simulation results are stored in the `sim_result` variable.
     *
     * @param bdl_iterator A reference to a BDL input iterator representing the gate layout at a given input
     * combination. The simulation is performed based on the configuration represented by the iterator.
     * @return Simulation results.
     */
    [[nodiscard]] sidb_simulation_result<Lyt>
    physical_simulation_of_layout(const bdl_input_iterator<Lyt>& bdl_iterator) noexcept
    {
        assert(parameters.simulation_parameters.base == 2 && "base number is set to 3");
        if (parameters.sim_engine == sidb_simulation_engine::EXGS)
        {
            // perform an exhaustive ground state simulation
            return exhaustive_ground_state_simulation(*bdl_iterator, parameters.simulation_parameters);
        }
        if (parameters.sim_engine == sidb_simulation_engine::QUICKSIM)
        {
            // perform a heuristic simulation
            const quicksim_params qs_params{parameters.simulation_parameters, 500, 0.6};
            return quicksim(*bdl_iterator, qs_params);
        }
        if (parameters.sim_engine == sidb_simulation_engine::QUICKEXACT)
        {
            // perform exact simulation
            const quickexact_params<cell<Lyt>> quickexact_params{
                parameters.simulation_parameters,
                fiction::quickexact_params<cell<Lyt>>::automatic_base_number_detection::OFF};
            return quickexact(*bdl_iterator, quickexact_params);
        }

        assert(false && "unsupported simulation engine");

        return sidb_simulation_result<Lyt>{};
    }
};

}  // namespace detail

/**
 * Determine the operational status of an SiDB layout.
 *
 * This function checks the operational status of a given gate layout using the `is_operational` algorithm. It
 * determines whether the gate layout is operational and returns the correct result for all \f$2^n\f$ input
 * combinations.
 *
 * @tparam Lyt SiDB cell-level layout type.
 * @tparam TT The type of the truth table specifying the layout behavior.
 * @param lyt The SiDB cell-level layout to be checked.
 * @param spec Expected Boolean function of the layout given as a multi-output truth table.
 * @param params Parameters for the `is_operational` algorithm.
 * @return A pair containing the operational status of the gate layout (either `OPERATIONAL` or `NON_OPERATIONAL`) and
 * the number of input combinations tested.
 */
template <typename Lyt, typename TT>
[[nodiscard]] std::pair<operational_status, std::size_t>
is_operational(const Lyt& lyt, const Lyt& skeleton, const std::vector<TT>& spec,
               const is_operational_params& params = {}) noexcept
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
    static_assert(kitty::is_truth_table<TT>::value, "TT is not a truth table");

    assert(lyt.num_pis() > 0 && "lyt needs input cells");
    assert(lyt.num_pos() > 0 && "lyt needs output cells");

    assert(!spec.empty());
    // all elements in tts must have the same number of variables
    assert(std::adjacent_find(spec.cbegin(), spec.cend(),
                              [](const auto& a, const auto& b)
                              { return a.num_vars() != b.num_vars(); }) == spec.cend());

    detail::is_operational_impl<Lyt, TT> p{lyt, skeleton, spec, params};

    if (params.condition == operation_ctiterium::ALLOW_KINKS)
    {
        return {p.run(), p.get_number_of_simulator_invocations()};
    }

    return {p.run_kinks_are_forbidden(), p.get_number_of_simulator_invocations()};
}
/**
 * This function determines the input combinations for which the SiDB-based logic, represented by the
 * provided layout (`lyt`) and truth table specifications (`spec`), produces the correct output.
 *
 * @tparam Lyt Type of the cell-level layout.
 * @tparam TT Type of the truth table.
 * @param lyt The SiDB layout.
 * @param spec Vector of truth table specifications.
 * @param params Parameters to simualte if a input combination is operational.
 * @return The count of operational input combinations.
 */
template <typename Lyt, typename TT>
[[nodiscard]] std::set<uint64_t> operational_input_patterns(const Lyt& lyt, const std::vector<TT>& spec,
                                                            const is_operational_params& params = {}) noexcept
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
    static_assert(kitty::is_truth_table<TT>::value, "TT is not a truth table");

    assert(lyt.num_pis() > 0 && "skeleton needs input cells");
    assert(lyt.num_pos() > 0 && "skeleton needs output cells");

    assert(!spec.empty());
    // all elements in tts must have the same number of variables
    assert(std::adjacent_find(spec.begin(), spec.end(),
                              [](const auto& a, const auto& b) { return a.num_vars() != b.num_vars(); }) == spec.end());

    detail::is_operational_impl<Lyt, TT> p{lyt, spec, params};

    return p.determine_operational_input_patterns();
}

}  // namespace fiction

#endif  // FICTION_IS_OPERATIONAL_HPP
