//
// Created by Jan Drewniok on 11.01.23.
//

#ifndef FICTION_QUICKSIM_HPP
#define FICTION_QUICKSIM_HPP

#include "fiction/algorithms/simulation/sidb/energy_distribution.hpp"
#include "fiction/algorithms/simulation/sidb/minimum_energy.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/traits.hpp"

#include <fmt/format.h>
#include <mockturtle/utils/stopwatch.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <mutex>
#include <thread>
#include <vector>

namespace fiction
{

/**
 * This struct stores the parameters for the *QuickSim* algorithm.
 */
struct quicksim_params
{
    /**
     * General parameters for the simulation of the physical SiDB system.
     */
    sidb_simulation_parameters phys_params{};
    /**
     * Number of iterations to run the simulation for.
     */
    uint64_t interation_steps{80};
    /**
     * `alpha` parameter for the *QuickSim* algorithm (should be reduced if no result is found).
     */
    double alpha{0.7};
    /**
     * Number of threads to spawn. By default the number of threads is set to the number of available hardware threads.
     */
    uint64_t number_threads{std::thread::hardware_concurrency()};
};

/**
 * This struct stores the simulation runtime and all physically valid charge layouts gained by the *QuickSim* algorithm.
 *
 * @paramt Lyt Cell-level layout type.
 */
template <typename Lyt>
struct quicksim_stats
{
    /**
     * Total simulation runtime.
     */
    mockturtle::stopwatch<>::duration time_total{0};
    /**
     * Vector of all physically valid charge layouts.
     */
    std::vector<charge_distribution_surface<Lyt>> valid_lyts{};

    std::vector<std::pair<charge_distribution_surface<Lyt>, uint64_t>> defect_iter_num_valid_lyts{};
    /**
     * Report the simulation statistics in a human-readable fashion.
     *
     * @param out Output stream to write to.
     */
    void report(std::ostream& out = std::cout)
    {
        out << fmt::format("[i] total runtime: {:.2f} secs\n", mockturtle::to_seconds(time_total));

        if (!energy_distribution<Lyt>(valid_lyts).empty())
        {
            for (auto [energy, count] : energy_distribution<Lyt>(valid_lyts))
            {
                out << fmt::format("[i] lowest energy state: {:.4f} meV \n", minimum_energy(valid_lyts));
                out << fmt::format("[i] energy: {} | occurrence: {} \n", energy, count);
            }
        }
        else
        {
            std::cout << "no state found" << std::endl;
        }

        std::cout << "_____________________________________________________ \n";
    }
};

/**
 * The *QuickSim* algorithm is an electrostatic ground state simulation algorithm for SiDB layouts. It determines
 * physically valid charge configurations (with minimal energy) of a given (already initialized) charge distribution
 * layout. Depending on the simulation parameters, the ground state is found with a certain probability after one run.
 *
 * @tparam Lyt Cell-level layout type.
 * @param lyt The layout to simulate.
 * @param ps Physical parameters. They are material-specific and may vary from experiment to experiment.
 * @param pst Statistics. They store the simulation results (simulation runtime as well as all physically valid charge
 * distribution layouts).
 */
template <typename Lyt>
void quicksim(Lyt& lyt, const quicksim_params& ps = quicksim_params{}, quicksim_stats<Lyt>* pst = nullptr,
              std::vector<std::unordered_map<typename Lyt::cell, const sidb_defect>> defects = {})
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<Lyt>, "Lyt must be an SiDB layout");

    quicksim_stats<Lyt> st{};
    st.valid_lyts.reserve(ps.interation_steps);

    // measure run time (artificial scope)
    {
        mockturtle::stopwatch stop{st.time_total};
        uint64_t              counter_conf = 0;
        for (const auto& defect_conf : defects)
        {
            for (const auto& defect_cells : defects[0])
            {
                lyt.assign_cell_type(defect_cells.first, sidb_technology::cell_type::NORMAL);
            }
            charge_distribution_surface charge_lyt{lyt};

            // set the given physical parameters
            charge_lyt.set_physical_parameters(ps.phys_params);
            charge_lyt.set_base_num(2);
            charge_lyt.set_all_charge_states(sidb_charge_state::NEGATIVE);
            for (const auto& defect : defect_conf)
            {
                if (defect.second.charge == 0)
                {
                    charge_lyt.assign_charge_state(defect.first, sidb_charge_state::NEUTRAL);
                }
            }
            charge_lyt.update_after_charge_change(false);
            auto negative_sidb_indices = charge_lyt.negative_sidb_detection();

            for (const auto& defect : defect_conf)
            {
                if (defect.second.charge == -1)
                {
                    auto index = static_cast<uint64_t>(charge_lyt.cell_to_index(defect.first));
                    if (std::find(negative_sidb_indices.begin(), negative_sidb_indices.end(), index) ==
                        negative_sidb_indices.end())
                    {
                        negative_sidb_indices.push_back(charge_lyt.cell_to_index(defect.first));
                    }
                }
            }

            if (charge_lyt.is_physically_valid())
            {
                st.valid_lyts.push_back(charge_distribution_surface<Lyt>{charge_lyt});
            }

            charge_lyt.set_all_charge_states(sidb_charge_state::NEUTRAL);
            charge_lyt.update_after_charge_change();

            if (!negative_sidb_indices.empty())
            {
                if (charge_lyt.is_physically_valid())
                {
                    st.valid_lyts.push_back(charge_distribution_surface<Lyt>{charge_lyt});
                }
            }

            // If the number of threads is initially set to zero, the simulation is run with one thread.
            const uint64_t num_threads = std::max(ps.number_threads, uint64_t{1});

            // split the iterations among threads
            const auto iter_per_thread =
                std::max(ps.interation_steps / num_threads,
                         uint64_t{1});  // If the number of set threads is greater than the number of iterations, the
                                        // number of threads defines how many times QuickSim is repeated

            std::vector<std::thread> threads{};
            threads.reserve(num_threads);
            std::mutex                                    mutex{};  // used to control access to shared resources
            std::vector<charge_distribution_surface<Lyt>> charge_lyt_coll{};
            std::vector<std::pair<charge_distribution_surface<Lyt>, uint64_t>> defect_coll{};
            for (uint64_t z = 0ul; z < num_threads; z++)
            {
                threads.emplace_back(
                    [&]
                    {
                        charge_distribution_surface<Lyt> charge_lyt_copy{charge_lyt};

                        for (uint64_t l = 0ul; l < iter_per_thread; ++l)
                        {
                            for (uint64_t i = 0ul; i < charge_lyt_copy.num_cells(); ++i)
                            {
                                {
                                    const std::lock_guard lock{mutex};
                                    if (std::find(negative_sidb_indices.cbegin(), negative_sidb_indices.cend(), i) !=
                                        negative_sidb_indices.cend())
                                    {
                                        continue;
                                    }
                                }

                                std::vector<uint64_t> index_start{i};

                                charge_lyt_copy.set_all_charge_states(sidb_charge_state::NEUTRAL);

                                for (const auto& index : negative_sidb_indices)
                                {
                                    charge_lyt_copy.assign_charge_state_by_cell_index(static_cast<uint64_t>(index),
                                                                                      sidb_charge_state::NEGATIVE);
                                    index_start.push_back(static_cast<uint64_t>(index));
                                }

                                charge_lyt_copy.assign_charge_state_by_cell_index(i, sidb_charge_state::NEGATIVE);
                                charge_lyt_copy.update_after_charge_change();

                                if (charge_lyt_copy.is_physically_valid())
                                {
                                    const std::lock_guard lock{mutex};
                                    for (const auto& defect_cells : defects[0])
                                    {
                                        lyt.assign_cell_type(defect_cells.first, sidb_technology::cell_type::EMPTY);
                                    }
                                    charge_distribution_surface charge_lyt_copy_final{lyt, ps.phys_params};
                                    charge_lyt_copy_final.foreach_cell(
                                        [&charge_lyt_copy_final, &charge_lyt_copy](const auto& c1) {
                                            charge_lyt_copy_final.assign_charge_state(
                                                c1, charge_lyt_copy.get_charge_state(c1), false);
                                        });
                                    charge_lyt_copy_final.charge_distr_to_index();
                                    charge_lyt_copy_final.update_after_charge_change();
                                    charge_lyt_coll.push_back(charge_lyt_copy_final);
                                    defect_coll.push_back(std::make_pair(charge_lyt_copy_final, counter_conf));
                                    //                                    st.valid_lyts.push_back(charge_distribution_surface<Lyt>{charge_lyt_copy_final});
                                    for (const auto& defect_cells : defects[0])
                                    {
                                        lyt.assign_cell_type(defect_cells.first, sidb_technology::cell_type::NORMAL);
                                    }
                                }

                                const auto upper_limit = std::min(
                                    static_cast<uint64_t>(static_cast<double>(charge_lyt_copy.num_cells()) / 1.5),
                                    charge_lyt_copy.num_cells() - negative_sidb_indices.size());

                                for (uint64_t num = 0ul; num < upper_limit; num++)
                                {
                                    charge_lyt_copy.adjacent_search(ps.alpha, index_start);
                                    charge_lyt_copy.validity_check();

                                    if (charge_lyt_copy.is_physically_valid())
                                    {
                                        const std::lock_guard lock{mutex};
                                        for (const auto& defect_cells : defects[0])
                                        {
                                            lyt.assign_cell_type(defect_cells.first, sidb_technology::cell_type::EMPTY);
                                        }
                                        charge_distribution_surface charge_lyt_copy_final{lyt, ps.phys_params};
                                        charge_lyt_copy_final.foreach_cell(
                                            [&charge_lyt_copy_final, &charge_lyt_copy](const auto& c1) {
                                                charge_lyt_copy_final.assign_charge_state(
                                                    c1, charge_lyt_copy.get_charge_state(c1), false);
                                            });
                                        charge_lyt_copy_final.charge_distr_to_index();
                                        charge_lyt_copy_final.update_after_charge_change();
                                        charge_lyt_coll.push_back(charge_lyt_copy_final);
                                        defect_coll.push_back(std::make_pair(charge_lyt_copy_final, counter_conf));
                                        //                                        st.valid_lyts.push_back(
                                        //                                            charge_distribution_surface<Lyt>{charge_lyt_copy_final});
                                        for (const auto& defect_cells : defects[0])
                                        {
                                            lyt.assign_cell_type(defect_cells.first,
                                                                 sidb_technology::cell_type::NORMAL);
                                        }
                                    }
                                }
                            }
                        }
                    });
            }

            for (auto& thread : threads)
            {
                thread.join();
            }

            if (charge_lyt_coll.size() == 0)
            {
                std::cout << "no layout found!" << std::endl;
            }

            auto min_energy = minimum_energy(charge_lyt_coll);
            for (const auto& lyts : charge_lyt_coll)
            {
                if (lyts.get_system_energy() == min_energy)
                {
                    st.valid_lyts.push_back(lyts);
                    break;
                }
            }
            for (const auto& coll : defect_coll)
            {
                if (coll.first.get_system_energy() == min_energy)
                {
                    st.defect_iter_num_valid_lyts.push_back(std::make_pair(coll.first, coll.second));
                    break;
                }
            }
            counter_conf += 1;
            for (const auto& defect_cells : defects[0])
            {
                lyt.assign_cell_type(defect_cells.first, sidb_technology::cell_type::EMPTY);
            }
        }

        if (pst)
        {
            *pst = st;
        }
    }
}

}  // namespace fiction

#endif  // FICTION_QUICKSIM_HPP
