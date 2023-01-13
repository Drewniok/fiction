//
// Created by Jan Drewniok on 11.01.23.
//

#ifndef FICTION_TEMPERATURE_HPP
#define FICTION_TEMPERATURE_HPP

#include "fiction/algorithms/simulation_sidb/TTS.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/algorithms/simulation_sidb/occupation_function.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/hash.hpp"
#include <cmath>

namespace fiction
{

template <typename Lyt>
double critical_temp(const std::vector<charge_distribution_surface<Lyt>>& valid_lyts, const double convlevel = 0.997, const uint64_t temp_limit = 400)
{
    std::vector<double> temp_values{};
    for (uint64_t i = 0; i <= temp_limit * 5; i++)
    {
        temp_values.push_back(static_cast<double>(i) / 5.0);
    }

    std::vector<double>        energies{};
    std::map<double, uint64_t> energy_deg{};

    energies.reserve(valid_lyts.size());

    for (auto& it : valid_lyts)
    {
        energies.push_back(std::round(it.get_system_energy() * 10000) / 10000);
    }

    auto min_energy = minimum_energy(valid_lyts);

    for (double & energy : energies)
    {

        energy      = energy - min_energy;

        uint64_t counter = 0;
        for (auto& it : valid_lyts)
        {
            if (std::abs(energy - (it.get_system_energy() - min_energy)) < 0.0001)
            {
                counter += 1;
            }
            else
            {
                continue;
            }
        }
        energy_deg.insert(std::make_pair(energy, counter));
    }


    for (const auto &temp: temp_values)
    {

        if (occu_prop(energy_deg,temp) < convlevel)
        {
            return temp;
        }
        else if (std::abs(temp - static_cast<double>(temp_limit)) < 0.001)
        {
            return static_cast<double>(temp_limit);
        }

        else
        {
            continue;
        }
    }
}

}  // namespace fiction

#endif  // FICTION_TEMPERATURE_HPP
