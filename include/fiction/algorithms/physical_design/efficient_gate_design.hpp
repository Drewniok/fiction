//
// Created by Jan Drewniok on 08.05.24.
//

#ifndef FICTION_EFFICIENT_GATE_DESIGN_HPP
#define FICTION_EFFICIENT_GATE_DESIGN_HPP

#include "fiction/algorithms/iter/bdl_input_iterator.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_pairs.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/traits.hpp"

#include <atomic>
#include <cstdint>
#include <future>
#include <mutex>
#include <optional>
#include <set>
#include <vector>

namespace fiction
{

template <typename Lyt>
std::optional<bdl_pair<Lyt>> bdl_parter_upper(const bdl_pair<Lyt>& given_bdl, const std::set<bdl_pair<Lyt>>& bdl_pairs)
{
    for (const auto& bdl : bdl_pairs)
    {
        if (sidb_nm_distance<Lyt>(Lyt{}, given_bdl.upper, bdl.lower) < 4 && given_bdl != bdl &&
            given_bdl.upper.y >= bdl.lower.y)
        {
            return std::optional<bdl_pair<Lyt>>(bdl);
        }
    }
    return std::nullopt;
};

template <typename Lyt>
std::optional<bdl_pair<Lyt>> bdl_parter_lower(const bdl_pair<Lyt>& given_bdl, const std::set<bdl_pair<Lyt>>& bdl_pairs)
{
    for (const auto& bdl : bdl_pairs)
    {
        if (sidb_nm_distance<Lyt>(Lyt{}, given_bdl.lower, bdl.upper) < 4 && given_bdl != bdl &&
            given_bdl.lower.y <= bdl.upper.y)
        {
            return std::optional<bdl_pair<Lyt>>(bdl);
        }
    }
    return std::nullopt;
};

enum class WIRE
{
    ALL,
    INPUT,
    OUTPUT
};

template <typename Lyt>
struct efficient_gate_design_params
{
    detect_bdl_pairs_params       bdl_params{};
    design_sidb_gates_params<Lyt> design_params{};
};

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
    efficient_gate_design_impl(Lyt& lyt, const std::vector<TT>& tt,
                               const efficient_gate_design_params<Lyt>& ps) noexcept :
            skeleton{lyt},
            truth_table{tt},
            params{ps},
            chains_input{detect_wire_bdl_chains(skeleton, WIRE::INPUT)},
            chains_output{detect_wire_bdl_chains(skeleton, WIRE::OUTPUT)},
            all_canvas_layouts{design_sidb_gate_candidates(skeleton, truth_table, params.design_params)},
            all_chains{detect_wire_bdl_chains(lyt)}
    {
        std::cout << all_canvas_layouts.size() << std::endl;
    }

    std::vector<std::vector<bdl_pair<Lyt>>> detect_wire_bdl_chains(const Lyt& layout,
                                                                   const WIRE wire_selection = WIRE::ALL)
    {
        std::set<bdl_pair<Lyt>> bdl_pairs{};
        const auto              normal_bdls = detect_bdl_pairs(layout, sidb_technology::cell_type::NORMAL);
        auto                    input_bdls  = detect_bdl_pairs(layout, sidb_technology::cell_type::INPUT);
        std::sort(input_bdls.begin(), input_bdls.end());
        auto output_bdls = detect_bdl_pairs(layout, sidb_technology::cell_type::OUTPUT);
        std::sort(output_bdls.begin(), output_bdls.end());
        for (const auto& bdl : normal_bdls)
        {
            bdl_pairs.insert(bdl);
        }

        for (auto& bdl : output_bdls)
        {
            bdl.type = sidb_technology::cell_type::NORMAL;
            bdl_pairs.insert(bdl);
        }

        for (auto& bdl : input_bdls)
        {
            bdl.type = sidb_technology::cell_type::NORMAL;
            bdl_pairs.insert(bdl);
        }

        std::vector<bdl_pair<Lyt>> bdl_pairs_vector(bdl_pairs.begin(), bdl_pairs.end());

        // Sort the vector
        std::sort(bdl_pairs_vector.begin(), bdl_pairs_vector.end());

        // Construct a new set from the sorted vector
        std::set<bdl_pair<Lyt>> sorted_bdl_pairs(bdl_pairs_vector.begin(), bdl_pairs_vector.end());

        // std::sort(bdl_pairs.begin(), bdl_pairs.end());

        std::vector<std::vector<bdl_pair<Lyt>>> chains{};

        while (!sorted_bdl_pairs.empty())
        {
            std::vector<bdl_pair<Lyt>> chain{};

            bool       neighbor_bdl_found = true;
            auto       front_bdl_pair     = *sorted_bdl_pairs.begin();
            const auto start_pair         = front_bdl_pair;

            chain.push_back(front_bdl_pair);
            sorted_bdl_pairs.erase(front_bdl_pair);

            while (neighbor_bdl_found)
            {
                const auto lower_partner = bdl_parter_lower<Lyt>(front_bdl_pair, sorted_bdl_pairs);
                if (lower_partner.has_value())
                {
                    chain.push_back(lower_partner.value());
                    sorted_bdl_pairs.erase(lower_partner.value());
                    front_bdl_pair = lower_partner.value();
                }
                else
                {
                    front_bdl_pair = start_pair;
                    if (bdl_parter_upper<Lyt>(front_bdl_pair, sorted_bdl_pairs).has_value())
                    {
                        chain.push_back(lower_partner.value());
                        sorted_bdl_pairs.erase(lower_partner.value());
                    }
                    else
                    {
                        neighbor_bdl_found = false;
                        chains.push_back(chain);
                    }
                }
            }
        }

        if (wire_selection == WIRE::INPUT)
        {
            std::vector<std::vector<bdl_pair<Lyt>>> result{};
            for (auto& bdl : input_bdls)
            {
                bdl.type = sidb_technology::cell_type::NORMAL;
                for (auto& chain : chains)
                {
                    auto it = std::find(chain.begin(), chain.end(), bdl);
                    if (it != chain.end())
                    {
                        chain.erase(it);
                        result.push_back(chain);
                        break;
                    }
                }
            }
            return result;
        }
        if (wire_selection == WIRE::OUTPUT)
        {
            std::vector<std::vector<bdl_pair<Lyt>>> result{};
            for (auto& bdl : output_bdls)
            {
                bdl.type = sidb_technology::cell_type::NORMAL;
                for (auto& chain : chains)
                {
                    auto it = std::find(chain.begin(), chain.end(), bdl);
                    if (it != chain.end())
                    {
                        result.push_back(chain);
                        break;
                    }
                }
            }
            return result;
        }

        std::vector<std::vector<bdl_pair<Lyt>>> result{};

        for (auto& chain : chains)
        {
            for (auto& bdl : input_bdls)
            {
                bdl.type = sidb_technology::cell_type::NORMAL;
                auto it  = std::find(chain.begin(), chain.end(), bdl);
                if (it != chain.end())
                {
                    chain.erase(it);
                    result.push_back(chain);
                    break;
                }
                result.push_back(chain);
                break;
            }
        }

        return result;
    }

    void set_charge_distribution(charge_distribution_surface<Lyt>& layout, uint64_t current_input_index)
    {
        layout.assign_all_charge_states(sidb_charge_state::NEGATIVE);

        for (auto i = 0u; i < all_chains.size(); i++)
        {
            if ((current_input_index & (uint64_t{1ull} << i)) != 0ull)
            {
                for (const auto& bdl : all_chains[i])
                {
                    layout.assign_charge_state(bdl.upper, sidb_charge_state::NEUTRAL);
                    layout.assign_charge_state(bdl.lower, sidb_charge_state::NEGATIVE);
                }
            }
            else
            {
                for (const auto& bdl : all_chains[i])
                {
                    layout.assign_charge_state(bdl.upper, sidb_charge_state::NEGATIVE);
                    layout.assign_charge_state(bdl.lower, sidb_charge_state::NEUTRAL);
                }
            }
        }
    }

    void set_charge_distribution_based_on_logic(charge_distribution_surface<Lyt>& layout,
                                                const uint64_t                    current_input_index)
    {
        layout.assign_all_charge_states(sidb_charge_state::NEGATIVE);

        for (auto i = 0u; i < chains_input.size(); i++)
        {
            if ((current_input_index & (uint64_t{1ull} << i)) != 0ull)
            {
                for (const auto& bdl : chains_input[chains_input.size() - 1 - i])
                {
                    layout.assign_charge_state(bdl.upper, sidb_charge_state::NEUTRAL);
                    layout.assign_charge_state(bdl.lower, sidb_charge_state::NEGATIVE);
                }
            }
            else
            {
                for (const auto& bdl : chains_input[chains_input.size() - 1 - i])
                {
                    layout.assign_charge_state(bdl.upper, sidb_charge_state::NEGATIVE);
                    layout.assign_charge_state(bdl.lower, sidb_charge_state::NEUTRAL);
                }
            }
        }

        for (auto i = 0u; i < chains_output.size(); i++)
        {
            for (const auto& bdl : chains_output[i])
            {
                if (kitty::get_bit(truth_table[i], current_input_index))
                {
                    layout.assign_charge_state(bdl.upper, sidb_charge_state::NEUTRAL);
                    layout.assign_charge_state(bdl.lower, sidb_charge_state::NEGATIVE);
                }
                else
                {
                    layout.assign_charge_state(bdl.upper, sidb_charge_state::NEGATIVE);
                    layout.assign_charge_state(bdl.lower, sidb_charge_state::NEUTRAL);
                }
            }
        }
    }

    bool physically_validity_is_fulfilled(const Lyt& canvas_lyt)
    {
        auto current_layout = skeleton.clone();
        canvas_lyt.foreach_cell([&](const auto& c)
                                { current_layout.assign_cell_type(c, Lyt::technology::cell_type::NORMAL); });

        charge_distribution_surface cds_canvas{canvas_lyt, params.design_params.simulation_parameters};

        auto bii = bdl_input_iterator<Lyt>{current_layout, params.bdl_params};

        for (auto i = 0u; i < truth_table.front().num_bits(); ++i, ++bii)
        {
            charge_distribution_surface cds_layout{*bii, params.design_params.simulation_parameters};

            set_charge_distribution_based_on_logic(cds_layout, i);

            bool physical_valid = false;

            cds_canvas.assign_charge_index(0);
            while (cds_canvas.get_charge_index_and_base().first < cds_canvas.get_max_charge_index())
            {
                cds_canvas.foreach_cell([&](const auto& c)
                                        { cds_layout.assign_charge_state(c, cds_canvas.get_charge_state(c)); });
                cds_layout.update_after_charge_change(dependent_cell_mode::FIXED,
                                                      energy_calculation::KEEP_OLD_ENERGY_VALUE);
                if (cds_layout.is_physically_valid())
                {
                    cds_layout.recompute_system_energy();
                    const auto energy = cds_layout.get_system_energy();
                    for (auto kink_states = 0u; kink_states < std::pow(2, (chains_input.size() + chains_output.size()));
                         ++kink_states)
                    {
                        set_charge_distribution(cds_layout, kink_states);
                        cds_canvas.assign_charge_index(0);

                        //print_sidb_layout(std::cout, cds_layout);

                        while (cds_canvas.get_charge_index_and_base().first < cds_canvas.get_max_charge_index())
                        {
                            cds_canvas.foreach_cell(
                                [&](const auto& c)
                                { cds_layout.assign_charge_state(c, cds_canvas.get_charge_state(c)); });
                            // print_sidb_layout(std::cout, cds_layout);
                            cds_layout.update_after_charge_change(dependent_cell_mode::FIXED,
                                                                  energy_calculation::KEEP_OLD_ENERGY_VALUE);
                            if (cds_layout.is_physically_valid())
                            {
                                cds_layout.recompute_system_energy();
                                if (cds_layout.get_system_energy() + physical_constants::POP_STABILITY_ERR < energy)
                                {
                                    return false;
                                }
                            }
                            cds_canvas.increase_charge_index_by_one(dependent_cell_mode::FIXED,
                                                                    energy_calculation::KEEP_OLD_ENERGY_VALUE,
                                                                    charge_distribution_history::NEGLECT);
                        }
                    }
                    physical_valid = true;
                    break;
                }
                cds_canvas.increase_charge_index_by_one(dependent_cell_mode::FIXED,
                                                        energy_calculation::KEEP_OLD_ENERGY_VALUE,
                                                        charge_distribution_history::NEGLECT);
            }
            if (!physical_valid)
            {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] std::vector<Lyt> design() noexcept
    {
        std::vector<Lyt>               all_designs{};
        std::vector<std::future<void>> futures{};
        futures.reserve(all_canvas_layouts.size());

        std::mutex mutex_to_protect_designer_gate_layouts;  // Mutex for protecting shared resources

        // Function to check validity and add layout to all_designs
        auto add_valid_layout = [&](const Lyt& canvas_lyt)
        {
            if (physically_validity_is_fulfilled(canvas_lyt))
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

        //        // Launch async tasks
                for (const auto& combination : all_canvas_layouts)
                {
                    futures.emplace_back(std::async(std::launch::async, add_valid_layout, combination));
                }

        // Launch async tasks
//        for (const auto& combination : all_canvas_layouts)
//        {
//            add_valid_layout(combination);
//        }

        // Wait for all tasks to finish
        for (auto& future : futures)
        {
            future.wait();
        }

        return all_designs;
    }

  private:
    Lyt&                                     skeleton;
    const std::vector<TT>&                   truth_table;
    const efficient_gate_design_params<Lyt>& params;
    std::vector<std::vector<bdl_pair<Lyt>>>  chains_input;
    std::vector<std::vector<bdl_pair<Lyt>>>  chains_output;
    std::vector<std::vector<bdl_pair<Lyt>>>  all_chains;
    std::vector<Lyt>                         all_canvas_layouts{};
};

}  // namespace detail

template <typename Lyt, typename TT>
std::vector<Lyt> design_all_efficient_gates(Lyt& skeleton, const std::vector<TT>& truth_table,
                                            const efficient_gate_design_params<Lyt>& params)
{
    detail::efficient_gate_design_impl<Lyt, TT> gate_design{skeleton, truth_table, params};
    return gate_design.design();
};

}  // namespace fiction

#endif  // FICTION_EFFICIENT_GATE_DESIGN_HPP
