//
// Created by Jan Drewniok on 23.12.22.
//

#ifndef FICTION_TTS_HPP
#define FICTION_TTS_HPP

#include <fiction/algorithms/simulation_sidb/quicksim.hpp>
#include <fiction/technology/charge_distribution_surface.hpp>

#include <chrono>
#include <limits>

namespace fiction
{

template <typename Lyt>
double minimum_energy(const std::vector<charge_distribution_surface<Lyt>>& result)
{
    auto min_energy = std::numeric_limits<double>::max();
    for (const auto &lyt : result)
    {
        if (lyt.get_system_energy() < min_energy)
        {
            min_energy = lyt.get_system_energy();
        }
    }
    return min_energy;
}

template <typename Lyt>
bool found_groundstate(const std::vector<charge_distribution_surface<Lyt>>& result_new_ap, const std::vector<charge_distribution_surface<Lyt>>& result_exact)
{
    auto min_energy_exact  = minimum_energy(result_exact);
    auto min_energy_new_ap = minimum_energy(result_new_ap);

    return std::abs(min_energy_exact - min_energy_new_ap) / min_energy_exact < 0.00001;
}

template <typename Lyt>
[[nodiscard]] std::pair<double, double> sim_acc_tts(charge_distribution_surface<Lyt>&                                   lyt,
                                  const std::vector<charge_distribution_surface<Lyt>>& result_exact,
                               const int& pp = 100, const int iteration_steps = 100, const double alpha = 0.7, const double convlevel = 0.997)
{
    int                 count = 0;
    std::vector<double> time;
    time.reserve(pp);

    for (int i = 0; i < pp; i++)
    {
        const auto t_start    = std::chrono::high_resolution_clock::now();
        auto       output_ap  = detail::quicksim<Lyt>(lyt, iteration_steps, alpha);
        const auto t_end      = std::chrono::high_resolution_clock::now();
        const auto elapsed    = t_end - t_start;
        auto       diff_first = std::chrono::duration<double>(elapsed).count() * 1000;
        time.push_back(diff_first);

        if (found_groundstate(output_ap, result_exact))
        {
            count += 1;
        }
    }

    auto single_runtime     = std::accumulate(time.begin(), time.end(), 0.0) / pp;
    std::cout << single_runtime << std::endl;
    auto acc            = static_cast<double>(count) / static_cast<double>(pp);

    double tts = single_runtime;

    if (acc == 1)
    {
        tts = single_runtime;
    }

    else
    {
        tts = (single_runtime * log(1.0 - convlevel) / log(1.0 - acc));
    }

    return std::make_pair(acc*100.0,tts);
}

} // namespace fiction
#endif  // FICTION_TTS_HPP
