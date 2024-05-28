//
// Created by Jan Drewniok on 08.05.24.
//

#ifndef FICTION_EFFICIENT_GATE_DESIGN_HPP
#define FICTION_EFFICIENT_GATE_DESIGN_HPP

#include "fiction/algorithms/iter/bdl_input_iterator.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_pairs.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_wires.hpp"
#include "fiction/io/print_layout.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/technology/physical_constants.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/traits.hpp"

#include <cstdint>
#include <future>
#include <mutex>
#include <optional>
#include <set>
#include <vector>

namespace fiction
{

/**
 * Statistics for the design of SiDB gates.
 */
struct efficient_gate_design_stats
{
    /**
     * The total runtime of SiDB gate design process.
     */
    mockturtle::stopwatch<>::duration time_total{0};

    mockturtle::stopwatch<>::duration time_total_pruning{0};

    mockturtle::stopwatch<>::duration time_total_operational_check{0};

    uint64_t all_possible_layouts{};

    uint64_t lp1{};

    uint64_t lp2{};

    uint64_t lp3{};
};

enum class DESIGN_MODE
{
    PRUNING,
    PRUNING_AND_OPERATIONAL_CHECK,
};

template <typename Lyt>
struct efficient_gate_design_params
{
    detect_bdl_pairs_params       bdl_params{};
    design_sidb_gates_params<Lyt> design_params{};
    DESIGN_MODE                   mode = DESIGN_MODE::PRUNING_AND_OPERATIONAL_CHECK;
};

template <typename Lyt>
using bdl_wire = std::vector<bdl_pair<Lyt>>;

namespace detail
{

template <typename Lyt, typename TT>
class efficient_gate_design_impl
{
  public:
    /**
     * Standard constructor. Initializes the layout, the truth table, the parameters and the statistics. Also
     * detects the output BDL pair, which is necessary for the operational domain computation. The layout must
     * have exactly one output BDL pair.
     *
     * @param lyt SiDB cell-level layout to be evaluated.
     * @param spec Expected Boolean function of the layout given as a multi-output truth table.
     * @param ps Parameters for the operational domain computation.
     * @param st Statistics of the process.
     */
    efficient_gate_design_impl(Lyt& lyt, const std::vector<TT>& tt, const efficient_gate_design_params<Lyt>& ps,
                               efficient_gate_design_stats& st) noexcept :
            skeleton{lyt},
            truth_table{tt},
            params{ps},
            input_wires{detect_bdl_wires(skeleton, WIRE::INPUT_WO_INPUT_BDL)},
            output_wires{detect_bdl_wires(skeleton, WIRE::OUTPUT)},
            all_wires{detect_bdl_wires(lyt)},
            all_canvas_layouts{design_sidb_gate_candidates(skeleton, truth_table, params.design_params)},
            stats{st},
            wire_directions{determine_wire_direction(detect_bdl_pairs(skeleton, sidb_technology::cell_type::INPUT),
                                                     detect_bdl_wires(skeleton, WIRE::INPUT_W_INPUT_BDL))}
    {
        std::cout << all_canvas_layouts.size() << std::endl;
    }

    void set_charge_distribution(charge_distribution_surface<Lyt>& layout, const uint64_t current_input_index)
    {
        layout.assign_all_charge_states(sidb_charge_state::NEGATIVE, false);

        for (auto i = 0u; i < all_wires.size(); i++)
        {
            if ((current_input_index & (uint64_t{1ull} << i)) != 0ull)
            {
                for (const auto& bdl : all_wires[i])
                {
                    layout.assign_charge_state(bdl.upper, sidb_charge_state::NEUTRAL, false);
                    layout.assign_charge_state(bdl.lower, sidb_charge_state::NEGATIVE, false);
                }
            }
            else
            {
                for (const auto& bdl : all_wires[i])
                {
                    layout.assign_charge_state(bdl.upper, sidb_charge_state::NEGATIVE, false);
                    layout.assign_charge_state(bdl.lower, sidb_charge_state::NEUTRAL, false);
                }
            }
        }
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

    [[nodiscard]] bool is_physical_validity_feasible(const Lyt& canvas_lyt)
    {
        std::mutex mutex_to_protect_stats;
        auto       current_layout = skeleton.clone();
        cell<Lyt>  dependent_cell{};
        canvas_lyt.foreach_cell(
            [&](const auto& c)
            {
                current_layout.assign_cell_type(c, Lyt::technology::cell_type::NORMAL);
                dependent_cell = c;
            });

        charge_distribution_surface cds_canvas{canvas_lyt, params.design_params.simulation_parameters,
                                               sidb_charge_state::NEGATIVE, false};
        cds_canvas.assign_dependent_cell(dependent_cell);

        // print_sidb_layout(std::cout, current_layout);

        auto bii = bdl_input_iterator<Lyt>{current_layout, all_wires, wire_directions, params.bdl_params};

        for (auto i = 0u; i < truth_table.front().num_bits(); ++i, ++bii)
        {
            charge_distribution_surface cds_layout{*bii, params.design_params.simulation_parameters};

            if (can_positive_charges_occur(cds_layout, params.design_params.simulation_parameters))
            {
                const std::lock_guard lock(mutex_to_protect_stats);
                stats.lp1++;
                return false;
            }

            cds_layout.assign_dependent_cell(dependent_cell);

            // print_sidb_layout(std::cout, cds_layout);

            set_charge_distribution_based_on_logic(cds_layout, i);

            bool physical_valid = false;

            cds_canvas.assign_charge_index(0);
            uint64_t counter_cds_first = 0;
            while (cds_canvas.get_charge_index_and_base().first < cds_canvas.get_max_charge_index())
            {
                cds_canvas.foreach_cell([&](const auto& c)
                                        { cds_layout.assign_charge_state(c, cds_canvas.get_charge_state(c), false); });
                cds_layout.update_after_charge_change(dependent_cell_mode::VARIABLE,
                                                      energy_calculation::KEEP_OLD_ENERGY_VALUE);

                // print_sidb_layout(std::cout, cds_layout);

                if (cds_layout.is_physically_valid())
                {
                    // print_sidb_layout(std::cout, cds_layout);

                    cds_layout.recompute_system_energy();
                    const auto energy = cds_layout.get_system_energy();
                    for (auto kink_states = 0u; kink_states < std::pow(2, (input_wires.size() + output_wires.size()));
                         ++kink_states)
                    {
                        set_charge_distribution(cds_layout, kink_states);
                        cds_canvas.assign_charge_index(0);

                        // vprint_sidb_layout(std::cout, cds_layout);

                        uint64_t charge_index_counter = 0;
                        while (cds_canvas.get_charge_index_and_base().first < cds_canvas.get_max_charge_index())
                        {
                            cds_canvas.foreach_cell(
                                [&](const auto& c)
                                { cds_layout.assign_charge_state(c, cds_canvas.get_charge_state(c), false); });

                            cds_layout.update_after_charge_change(dependent_cell_mode::VARIABLE,
                                                                  energy_calculation::KEEP_OLD_ENERGY_VALUE,
                                                                  charge_distribution_history::NEGLECT);
                            if (cds_layout.is_physically_valid())
                            {
                                cds_layout.recompute_system_energy();
                                if (cds_layout.get_system_energy() + physical_constants::POP_STABILITY_ERR < energy)
                                {
                                    const std::lock_guard lock(mutex_to_protect_stats);
                                    stats.lp3++;
                                    return false;
                                }
                            }
                            charge_index_counter++;
                            cds_canvas.assign_charge_index(charge_index_counter);
                            //                            cds_canvas.increase_charge_index_by_one(dependent_cell_mode::FIXED,
                            //                                                                    energy_calculation::KEEP_OLD_ENERGY_VALUE,
                            //                                                                    charge_distribution_history::NEGLECT);
                        }
                    }
                    physical_valid = true;
                    break;
                }
                counter_cds_first++;
                //                cds_canvas.increase_charge_index_by_one(dependent_cell_mode::FIXED,
                //                                                        energy_calculation::KEEP_OLD_ENERGY_VALUE,
                //                                                        charge_distribution_history::NEGLECT);
                cds_canvas.assign_charge_index(counter_cds_first);
            }
            if (!physical_valid)
            {
                const std::lock_guard lock(mutex_to_protect_stats);
                stats.lp2++;
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] std::vector<Lyt> design() noexcept
    {
        mockturtle::stopwatch stop{stats.time_total_pruning};
        stats.all_possible_layouts = all_canvas_layouts.size();
        std::vector<Lyt> all_designs{};
        all_designs.reserve(all_canvas_layouts.size());

        std::vector<std::future<void>> futures{};
        futures.reserve(all_canvas_layouts.size());

        std::mutex mutex_to_protect_designer_gate_layouts;  // Mutex for protecting shared resources

        // Function to check validity and add layout to all_designs
        auto add_valid_layout = [&](const Lyt& canvas_lyt)
        {
            if (is_physical_validity_feasible(canvas_lyt))
            {
                Lyt modified_skeleton = skeleton.clone();  // Make a copy of skeleton_lyt to avoid data race
                canvas_lyt.foreach_cell([&](const auto& c)
                                        { modified_skeleton.assign_cell_type(c, Lyt::technology::cell_type::NORMAL); });

                // Lock mutex before modifying shared resource
                const std::lock_guard lock(mutex_to_protect_designer_gate_layouts);
                all_designs.push_back(modified_skeleton);
                // std::cout << "FOUND" << std::endl;
            }
        };

        //        Launch async tasks
        for (const auto& combination : all_canvas_layouts)
        {
            futures.emplace_back(std::async(std::launch::async, add_valid_layout, combination));
        }

        //                        for (const auto& combination : all_canvas_layouts)
        //                        {
        //                            add_valid_layout(combination);
        //                        }

        // Wait for all tasks to finish
        for (auto& future : futures)
        {
            future.wait();
        }

        return all_designs;
    }

    [[nodiscard]] std::vector<Lyt> select_all_operational_gates(const std::vector<Lyt>& all_left_over_layouts) noexcept
    {
        mockturtle::stopwatch stop{stats.time_total_operational_check};
        std::vector<Lyt>      all_design_after_operational_check{};
        all_design_after_operational_check.reserve(all_left_over_layouts.size());

        std::mutex mutex_to_protect_designer_gate_layouts;  // Mutex for protecting shared resources

        const auto op_params = is_operational_params{params.design_params.simulation_parameters};

        std::vector<std::future<void>> futures{};
        futures.reserve(all_left_over_layouts.size());

        auto add_valid_layout = [&](const Lyt& canvas_lyt)
        {
            if (is_operational(canvas_lyt, skeleton, truth_table, op_params).first == operational_status::OPERATIONAL)
            {
                // Lock mutex before modifying shared resource
                const std::lock_guard lock(mutex_to_protect_designer_gate_layouts);
                all_design_after_operational_check.push_back(canvas_lyt);
                // std::cout << "FOUND" << std::endl;
            }
        };

        //        Launch async tasks
        for (const auto& combination : all_canvas_layouts)
        {
            futures.emplace_back(std::async(std::launch::async, add_valid_layout, combination));
        }

        // Wait for all tasks to finish
        for (auto& future : futures)
        {
            future.wait();
        }

        return all_design_after_operational_check;
    }

  private:
    Lyt&                                    skeleton;
    const std::vector<TT>                   truth_table;
    const efficient_gate_design_params<Lyt> params;
    std::vector<bdl_wire<Lyt>>              input_wires;
    std::vector<bdl_wire<Lyt>>              output_wires;
    std::vector<bdl_wire<Lyt>>              all_wires;
    std::vector<Lyt>                        all_canvas_layouts{};
    efficient_gate_design_stats&            stats;
    const std::vector<bdl_wire_direction>   wire_directions;
};

}  // namespace detail

template <typename Lyt, typename TT>
[[nodiscard]] std::vector<Lyt> efficient_gate_design(Lyt& skeleton, const std::vector<TT>& truth_table,
                                                     const efficient_gate_design_params<Lyt>& params,
                                                     efficient_gate_design_stats*             stats = nullptr)
{
    efficient_gate_design_stats st{};
    std::vector<Lyt>            result{};
    {
        mockturtle::stopwatch                       stop{st.time_total};
        detail::efficient_gate_design_impl<Lyt, TT> gate_design{skeleton, truth_table, params, st};
        result = gate_design.design();
        if (params.mode == DESIGN_MODE::PRUNING_AND_OPERATIONAL_CHECK)
        {
            result = gate_design.select_all_operational_gates(result);
        }
    }

    if (stats)
    {
        *stats = st;
    }
    return result;
};

}  // namespace fiction

#endif  // FICTION_EFFICIENT_GATE_DESIGN_HPP
