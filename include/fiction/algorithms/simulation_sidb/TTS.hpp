//
// Created by Jan Drewniok on 23.12.22.
//

#ifndef FICTION_TTS_HPP
#define FICTION_TTS_HPP

#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/algorithms/simulation_sidb/new_approach.hpp>
#include <limits>
#include <chrono>


namespace fiction
{
template <typename Lyt>
double minimum_energy(const std::unordered_map<double, charge_distribution_surface<Lyt>>& result)
{
    auto min_energy = std::numeric_limits<double>::max();
    for (auto &it : result)
    {
        if (it.second.get_system_energy() < min_energy)
        {
            min_energy = it.second.get_system_energy();
        }
    }
    return min_energy;
}

template <typename Lyt>
bool found_groundstate(const std::unordered_map<double, charge_distribution_surface<Lyt>>& result_new_ap, const std::unordered_map<double, charge_distribution_surface<Lyt>>& result_exact)
{
    auto min_energy_exact = std::numeric_limits<double>::max();
    for (auto &it : result_exact)
    {
        if (it.second.get_system_energy() < min_energy_exact)
        {
            min_energy_exact = it.second.get_system_energy();
        }
    }

    auto min_energy_new_ap = std::numeric_limits<double>::max();
    for (auto &it : result_new_ap)
    {
        if (it.second.get_system_energy() < min_energy_new_ap)
        {
            min_energy_new_ap = it.second.get_system_energy();
        }
    }

    return std::abs(min_energy_exact - min_energy_new_ap) / min_energy_exact < 0.00001;
}

template <typename Lyt>
[[nodiscard]] std::pair<float,uint64_t> sim_acc_tts(charge_distribution_surface<Lyt>&                                   lyt,
                               const std::unordered_map<double, charge_distribution_surface<Lyt>>& result_exact,
                               const int& pp = 1000, const double& convlevel = 0.997, const int iteration_steps = 10, const double alpha = 0.7)
{
    int                                                          count = 0;
    std::unordered_map<double, charge_distribution_surface<Lyt>> output_ap{};
    auto                                                         t_start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < pp; i++)
    {
        std::unordered_map<double, charge_distribution_surface<Lyt>> output_ap = detail::Sim<Lyt>(lyt, iteration_steps, alpha);

        if (found_groundstate(output_ap, result_exact))
        {
            count += 1;
        }
    }
    auto t_end          = std::chrono::high_resolution_clock::now();
    auto elapsed        = t_end - t_start;
    auto diff_first     = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
    auto single_runtime = static_cast<double>(diff_first) / static_cast<double>(pp);
    auto acc            = static_cast<float>(count / pp);

    auto tts = std::numeric_limits<uint64_t>::max();

    if (acc == 1)
    {
        tts = static_cast<uint64_t>(single_runtime);
    }

    else
    {
        tts = static_cast<uint64_t>(single_runtime * log(1.0 - convlevel) / log(1.0 - static_cast<double>(acc)));
    }

    return std::make_pair(acc*100,tts);
}

} // namespace fiction
#endif  // FICTION_TTS_HPP
