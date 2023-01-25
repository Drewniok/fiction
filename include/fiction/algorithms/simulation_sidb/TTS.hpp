//
// Created by Jan Drewniok on 23.12.22.
//

#ifndef FICTION_TTS_HPP
#define FICTION_TTS_HPP

#include "fiction/algorithms//simulation_sidb/check_groundstate.hpp"
#include "fiction/algorithms//simulation_sidb/minimum_energy.hpp"
#include "fiction/algorithms/simulation_sidb/exhaustive_ground_state_simulation.hpp"
#include "fiction/algorithms/simulation_sidb/quicksim.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"

#include <chrono>


namespace fiction
{
/**
 * This struct stores the time-to-solution, the simulation accuracy and the average single simulation runtime of quicksim (see quicksim.hpp).
 *
 */
struct tts_stats
{
    double time_to_solution{};
    double acc{};
    double mean_single_runtime{};

    void report(std::ostream& out = std::cout)
    {
        out << fmt::format("time_to_solution: {} | acc: {} | t_(s): {} \n", time_to_solution, acc, mean_single_runtime);
    }
};

/**
 * This function determines the time-to-solution (time_to_solution) and the accuracy (acc) of the quicksim-algorithm.
 *
 * @paramt Lyt cell-level layout.
 * @param charge_distribution_surface<Lyt> charge layout that is used for the simulation.
 * @param tts_stats struct where the results (time_to_solution, acc, single runtime) are stored.
 * @param pp number of repetitions to determine the simulation accuracy (pp = 100 ==> accuracy is precise to 1 %).
 * @param iteration_steps simulation parameter (see quicksim.hpp).
 * @param alpha simulation parameter (see quicksim.hpp).
 * @param convlevel the time-to-solution also depends one the given confidence level which can be set here.
 */
template <typename Lyt>
void sim_acc_tts(charge_distribution_surface<Lyt>& lyt, tts_stats& ts, exgs_stats<Lyt>& result_exact,
                 const int& pp = 100, const int iteration_steps = 100, const double alpha = 0.7,
                 const double convlevel = 0.997)
{
    int                 count = 0;
    std::vector<double> time;
    time.reserve(pp);

    for (uint64_t i = 0; i < pp; i++)
    {
        quicksim_stats<Lyt> stats_quick{};
        const auto          t_start = std::chrono::high_resolution_clock::now();
        quicksim<Lyt>(lyt, stats_quick, lyt.get_phys_params(), iteration_steps, alpha);
        const auto t_end      = std::chrono::high_resolution_clock::now();
        const auto elapsed    = t_end - t_start;
        auto       diff_first = std::chrono::duration<double>(elapsed).count();
        time.push_back(diff_first);

        if (check_groundstate(stats_quick, result_exact))
        {
            count += 1;
        }
    }

    auto single_runtime = std::accumulate(time.begin(), time.end(), 0.0) / pp;
    std::cout << single_runtime << std::endl;
    auto acc = static_cast<double>(count) / static_cast<double>(pp);

    double tts = single_runtime;

    if (acc == 1)
    {
        tts = single_runtime;
    }
    else
    {
        tts = (single_runtime * log(1.0 - convlevel) / log(1.0 - acc));
    }

    ts.time_to_solution    = tts;
    ts.acc                 = acc * 100;
    ts.mean_single_runtime = single_runtime;
}
}  // namespace fiction

#endif  // FICTION_TTS_HPP
