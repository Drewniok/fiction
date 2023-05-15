//
// Created by Jan Drewniok on 18.12.22.
//

#ifndef FICTION_EXHAUSTIVE_GROUND_STATE_SIMULATION_HPP
#define FICTION_EXHAUSTIVE_GROUND_STATE_SIMULATION_HPP

#include "fiction/algorithms/simulation/sidb/energy_distribution.hpp"
#include "fiction/algorithms/simulation/sidb/minimum_energy.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"

#include <fmt/format.h>
#include <mockturtle/utils/stopwatch.hpp>

#include <iostream>
#include <vector>

namespace fiction
{

template <typename Lyt>
struct exgs_stats
{
    mockturtle::stopwatch<>::duration             time_total{0};
    std::vector<charge_distribution_surface<Lyt>> valid_lyts{};

    void report(std::ostream& out = std::cout) const
    {
        out << fmt::format("total time  = {:.2f} secs\n", mockturtle::to_seconds(time_total));
        if (!valid_lyts.empty())
        {
            for (const auto& [energy, count] : energy_distribution<Lyt>(valid_lyts))
            {
                out << fmt::format("energy: {} | occurance: {} \n", energy, count);
            }
            out << fmt::format("the ground state energy is  = {:.4f} \n", minimum_energy(valid_lyts));
        }
        else
        {
            std::cout << "no state found | if two state simulation is used, continue with three state" << std::endl;
        }
        out << fmt::format("{} physically valid charge states were found \n", valid_lyts.size());
        std::cout << "_____________________________________________________ \n";
    }
};

/**
 *  All metastable and physically valid charge distribution layouts are computed, stored in a vector and returned.
 *
 * @tparam Lyt Cell-level layout type.
 * @param lyt The layout to simulate.
 * @param params Simulation parameters.
 * @param ps Simulation statistics.
 */
template <typename Lyt>
void exhaustive_ground_state_simulation(
    Lyt& lyt, const sidb_simulation_parameters& params = sidb_simulation_parameters{}, exgs_stats<Lyt>* ps = nullptr,
    std::vector<std::unordered_map<typename Lyt::cell, const sidb_defect>> defects                  = {},
    const std::unordered_map<typename Lyt::cell, double>&                  local_external_potential = {},
    const double&                                                          global_potential         = 0) noexcept
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");

    exgs_stats<Lyt> st{};

    {
        mockturtle::stopwatch stop{st.time_total};

        charge_distribution_surface charge_lyt{lyt, params, sidb_charge_state::NEGATIVE};
        charge_lyt.set_physical_parameters(params);
        charge_lyt.set_all_charge_states(sidb_charge_state::NEGATIVE);
        charge_lyt.update_after_charge_change();
        uint64_t counter_conf = 0;
        for (const auto& defect_conf : defects)
        {
            counter_conf += 1;
            if (counter_conf % 3 == 0)
            {
                std::cout << counter_conf << std::endl;
            }
            for (const auto& [cell, defect] : defect_conf)
            {
                if (defect.epsilon_r == 0 && defect.lambda_tf == 0)
                {
                    charge_lyt.assign_defect(cell, sidb_defect{defect.type, defect.charge,
                                                               charge_lyt.get_phys_params().epsilon_r,
                                                               charge_lyt.get_phys_params().lambda_tf});
                }
                else if (defect.epsilon_r == 0 && defect.lambda_tf != 0)
                {
                    charge_lyt.assign_defect(
                        cell, sidb_defect{defect.type, defect.charge, charge_lyt.get_phys_params().epsilon_r});
                }
                else if (defect.epsilon_r != 0 && defect.lambda_tf == 0)
                {
                    charge_lyt.assign_defect(cell, sidb_defect{defect.type, defect.charge, defect.epsilon_r,
                                                               charge_lyt.get_phys_params().lambda_tf});
                }
                else
                {
                    charge_lyt.assign_defect(cell, defect);
                }
            }
            auto real_placed_defects = charge_lyt.get_defects();
            charge_lyt.set_local_external_potential(local_external_potential);
            charge_lyt.set_global_external_potential(global_potential, false);
            const auto sidbs_charge_lyt = charge_lyt.get_sidb_order();

            auto all_sidbs_in_lyt_without_detected_ones = sidbs_charge_lyt;
            // determine all SiDBs that have to be negatively charged to fulfill the population stability. This is an
            // efficient way to prune the search space by 2^k with k being the number of detected negatively charged
            // SiDBs.
            const auto detected_negative_sidb_indices = charge_lyt.negative_sidb_detection();
            // determine if three state simulation (i.e. positively charged SiDBs can occur) is required.
            // const auto three_state_simulation_required = charge_lyt.three_state_sim_required();
            const auto                      three_state_simulation_required = false;
            std::vector<typename Lyt::cell> detected_negative_sidbs{};
            detected_negative_sidbs.reserve(detected_negative_sidb_indices.size());

            // if layout has at least two SiDBs, the code inside this if-statement is executed.
            if (sidbs_charge_lyt.size() > 1)
            {
                detected_negative_sidbs.reserve(detected_negative_sidb_indices.size());
                for (const auto& index : detected_negative_sidb_indices)
                {
                    const auto cell = charge_lyt.index_to_cell(index);
                    detected_negative_sidbs.push_back(cell);
                    lyt.assign_cell_type(cell, Lyt::cell_type::EMPTY);
                }

                // all detected negatively charged SiDBs are erased from the all_sidbs_in_lyt_without_detected_ones
                // vector.
                all_sidbs_in_lyt_without_detected_ones.erase(
                    std::remove_if(all_sidbs_in_lyt_without_detected_ones.begin(),
                                   all_sidbs_in_lyt_without_detected_ones.end(),
                                   [&detected_negative_sidbs](const typename Lyt::cell& n)
                                   {
                                       return std::find(detected_negative_sidbs.begin(), detected_negative_sidbs.end(),
                                                        n) != detected_negative_sidbs.end();
                                   }),
                    all_sidbs_in_lyt_without_detected_ones.end());
            }
            if (!all_sidbs_in_lyt_without_detected_ones.empty() && sidbs_charge_lyt.size() > 1)
            {
                // the first cell from all_sidbs_in_lyt_without_detected_ones is chosen as the dependent cell to
                // initialize the layout (detected negatively charged SiDBs were erased in the step before).
                charge_distribution_surface charge_lyt_new{lyt, params, sidb_charge_state::NEUTRAL,
                                                           all_sidbs_in_lyt_without_detected_ones[0]};
                charge_lyt_new.set_local_external_potential(local_external_potential);
                charge_lyt_new.set_global_external_potential(global_potential);
                // IMPORTANT: The detected negatively charged SiDBs (they have to be negatively charged to fulfill the
                // population stability) are considered as negatively charged defects in the layout. Hence, there are no
                // "real" defects assigned but in order to set some SiDBs with a fixed negative charge, this way of
                // implementation is chosen.
                for (const auto& cell : detected_negative_sidbs)
                {
                    charge_lyt_new.assign_defect(cell, sidb_defect{sidb_defect_type::UNKNOWN, -1,
                                                                   charge_lyt_new.get_phys_params().epsilon_r,
                                                                   charge_lyt_new.get_phys_params().lambda_tf});
                }
                for (const auto& [cell, defect] : real_placed_defects)
                {
                    charge_lyt_new.assign_defect(cell, defect);
                }
                // update all local potentials, system energy and physically validity. Flag is set to "false" to allow
                // dependent cell to change its charge state based on the N-1 SiDBs to fulfill the local population
                // stability at its position.

                // false declares that the dependent cell is updated based on the local potential at the position.
                charge_lyt_new.update_after_charge_change(false);

                // if no positively charged DB can occur in the layout, this scope is executed.
                if (!three_state_simulation_required)
                {
                    charge_lyt_new.set_base_num(2);
                    uint64_t current_charge_index = 0;
                    uint64_t old_charge_index     = 0;
                    for (uint64_t i = 0; i <= charge_lyt_new.get_max_charge_index(); i++)
                    {
                        current_charge_index = (i ^ (i >> 1));
                        charge_lyt_new.set_charge_index_by_gray(
                            current_charge_index, old_charge_index, false, false,
                            true);  // "false" allows that the charge state of the dependent cell is
                                    // automatically changed based on the new charge distribution.
                        old_charge_index = current_charge_index;

                        if (charge_lyt_new.is_physically_valid())
                        {
                            charge_distribution_surface<Lyt> charge_lyt_copy{charge_lyt_new};
                            charge_lyt_copy.recompute_system_energy();
                            // the previously detected negatively charged SiDBs are added to the final layout.
                            for (const auto& cell : detected_negative_sidbs)
                            {
                                charge_lyt_copy.adding_sidb_to_layout(cell, -1);
                            }
                            charge_lyt_copy.charge_distr_to_index();
                            st.valid_lyts.push_back(charge_lyt_copy);
                        }
                    }

                    // the cells of the previously detected negatively charged SiDBs are added to the cell level layout.
                    for (const auto& cell : detected_negative_sidbs)
                    {
                        lyt.assign_cell_type(cell, Lyt::cell_type::NORMAL);
                    }
                }

                // if positively charged DBs can occur in the layout, 3-state simulation is conducted.
                else
                {
                    charge_lyt_new.set_all_charge_states(sidb_charge_state::NEGATIVE);
                    charge_lyt_new.update_after_charge_change();
                    const auto three_state_simulation_required_true = charge_lyt_new.three_state_sim_required();
                    charge_lyt_new.base_to_three();
                    charge_lyt_new.update_after_charge_change(false);
                    // std::cout << charge_lyt_new.get_charge_index_sub_layout().first << std::endl;
                    while (charge_lyt_new.get_charge_index().first < charge_lyt_new.get_max_charge_index())
                    {
                        if (charge_lyt_new.is_physically_valid())
                        {
                            charge_distribution_surface<Lyt> charge_lyt_copy{charge_lyt_new};
                            charge_lyt_copy.recompute_system_energy();
                            // the previously detected negatively charged SiDBs are added to the final layout.
                            for (const auto& cell : detected_negative_sidbs)
                            {
                                charge_lyt_copy.adding_sidb_to_layout(cell, -1);
                            }
                            st.valid_lyts.push_back(charge_lyt_copy);
                        }

                        while (charge_lyt_new.get_charge_index_sub_layout().first <
                               charge_lyt_new.get_max_charge_index_sub_layout())
                        {
                            if (charge_lyt_new.is_physically_valid())
                            {
                                charge_distribution_surface<Lyt> charge_lyt_copy{charge_lyt_new};
                                charge_lyt_copy.recompute_system_energy();
                                // the previously detected negatively charged SiDBs are added to the final layout.
                                for (const auto& cell : detected_negative_sidbs)
                                {
                                    charge_lyt_copy.adding_sidb_to_layout(cell, -1);
                                }
                                st.valid_lyts.push_back(charge_lyt_copy);
                            }
                            charge_lyt_new.increase_charge_index_by_sub_layout(
                                false, false, true,
                                true);  // "false" allows that the charge state of the dependent cell is
                                        // automatically changed based on the new charge distribution.
                        }

                        if (charge_lyt_new.is_physically_valid())
                        {
                            charge_distribution_surface<Lyt> charge_lyt_copy{charge_lyt_new};
                            for (const auto& cell : detected_negative_sidbs)
                            {
                                charge_lyt_copy.adding_sidb_to_layout(cell, -1);
                            }
                            st.valid_lyts.push_back(charge_lyt_copy);
                        }

                        charge_lyt_new.reset_charge_index_sub_layout(false, false, true, true);

                        charge_lyt_new.increase_charge_index_by_one(
                            false, false, true, true);  // "false" allows that the charge state of the dependent cell is
                                                        // automatically changed based on the new charge distribution.
                    }

                    while (charge_lyt_new.get_charge_index_sub_layout().first <
                           charge_lyt_new.get_max_charge_index_sub_layout())
                    {
                        if (charge_lyt_new.is_physically_valid())
                        {
                            charge_distribution_surface<Lyt> charge_lyt_copy{charge_lyt_new};
                            charge_lyt_copy.recompute_system_energy();
                            // the previously detected negatively charged SiDBs are added to the final layout.
                            for (const auto& cell : detected_negative_sidbs)
                            {
                                charge_lyt_copy.adding_sidb_to_layout(cell, -1);
                            }
                            st.valid_lyts.push_back(charge_lyt_copy);
                        }
                        charge_lyt_new.increase_charge_index_by_sub_layout(
                            false, false, true,
                            true);  // "false" allows that the charge state of the dependent cell is
                                    // automatically changed based on the new charge distribution.
                    }

                    if (charge_lyt_new.is_physically_valid())
                    {
                        charge_distribution_surface<Lyt> charge_lyt_copy{charge_lyt_new};
                        for (const auto& cell : detected_negative_sidbs)
                        {
                            charge_lyt_copy.adding_sidb_to_layout(cell, -1);
                        }
                        st.valid_lyts.push_back(charge_lyt_copy);
                    }

                    for (const auto& cell : detected_negative_sidbs)
                    {
                        lyt.assign_cell_type(cell, Lyt::cell_type::NORMAL);
                    }
                }
            }

            // in the case with only one SiDB in the layout (due to external potentials or defects, this single SiDB can
            // neutrally or even positively charged.)
            else if (sidbs_charge_lyt.size() == 1)
            {
                if (three_state_simulation_required)
                {
                    charge_lyt.set_base_num(3);
                }
                else
                {
                    charge_lyt.set_base_num(2);
                }

                while (charge_lyt.get_charge_index().first < charge_lyt.get_max_charge_index())
                {

                    if (charge_lyt.is_physically_valid())
                    {
                        charge_distribution_surface<Lyt> charge_lyt_copy{charge_lyt};
                        st.valid_lyts.push_back(charge_lyt_copy);
                    }
                    charge_lyt.increase_charge_index_by_one(
                        false);  // "false" allows that the charge state of the dependent cell is automatically changed
                                 // based on the new charge distribution.
                }
                if (charge_lyt.is_physically_valid())
                {
                    charge_distribution_surface<Lyt> charge_lyt_copy{charge_lyt};
                    charge_lyt_copy.charge_distr_to_index();
                    st.valid_lyts.push_back(charge_lyt_copy);
                }
            }
            // if the layout consists of only detected negatively charged SiDBs, this scope is
            // executed.
            else if (all_sidbs_in_lyt_without_detected_ones.empty() && sidbs_charge_lyt.size() > 1)
            {
                charge_distribution_surface<Lyt> charge_lyt_copy{charge_lyt};
                for (const auto& cell : detected_negative_sidbs)
                {
                    charge_lyt.adding_sidb_to_layout(cell, -1);
                }
                charge_lyt_copy.charge_distr_to_index();
                st.valid_lyts.push_back(charge_lyt_copy);
            }
        }
    }

    if (ps)
    {
        *ps = st;
    }
}

}  // namespace fiction

#endif  // FICTION_EXHAUSTIVE_GROUND_STATE_SIMULATION_HPP
