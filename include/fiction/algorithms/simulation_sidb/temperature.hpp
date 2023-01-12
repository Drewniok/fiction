//
// Created by Jan Drewniok on 11.01.23.
//

#ifndef FICTION_TEMPERATURE_HPP
#define FICTION_TEMPERATURE_HPP

#include "fiction/algorithms/simulation_sidb/TTS.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/hash.hpp"
#include <cmath>

namespace fiction
{

double occu_prop(std::map<double, uint64_t>& energy_deg, const double &temp)
{
    double part_func = 0;

    for (const auto& it : energy_deg)
    {
        std::cout << static_cast<double>(it.second) << " | " << it.first << std::endl;
        part_func += static_cast<double>(it.second) * std::exp(-it.first * 12000 / temp);
    }

    auto it = energy_deg.begin();

    std::cout << "T: " << temp << " | " << part_func << " | " << static_cast<double>(it->second) * std::exp(-it->first * 12000 / temp) / part_func << std::endl;

    return static_cast<double>(it->second) * std::exp(-it->first * 12000 / temp) / part_func;

};

template <typename Lyt>
std::vector<charge_distribution_surface<Lyt>> removeDuplicates(std::vector<charge_distribution_surface<Lyt>>& lyts)
{

    std::vector<charge_distribution_surface<Lyt>> collect{};

    std::set<uint64_t> index_set{};
    for (auto &it: lyts)
    {
        it.chargeconf_to_index();
        index_set.insert(it.get_charge_index().first);
    }

    for (const auto it_set : index_set)
    {
        for (auto &it: lyts)
        {
            if (it_set == it.get_charge_index().first)
            {
                collect.push_back(it);
                break;
            }
            else
            {
                continue;
            }
        }
    }

    return collect;

}

template <typename Lyt>
double critical_temp(const std::vector<charge_distribution_surface<Lyt>>& valid_lyts, const double convlevel = 0.997, const double temp_limit = 400.0)
{
    std::vector<double> temp_values{};
    double start = 0.1;
    double end = temp_limit;
    for (int i = start; i <= end*5; i ++) {
        temp_values.push_back(static_cast<double>(i)/5.0);
    }

    auto valid_lyts_unique = valid_lyts;
    std::vector<double>        energies{};
    std::map<double, uint64_t> energy_deg{};

    energies.reserve(valid_lyts_unique.size());

    for (auto& it : valid_lyts_unique)
    {
       energies.push_back(std::round(it.get_system_energy()*10000)/10000);
    }

    auto min_energy = minimum_energy(valid_lyts_unique);

    for (int i = 0u; i < energies.size(); i++)
    {

        energies[i]      = energies[i] - min_energy;

        uint64_t counter = 0;
        for (auto& it : valid_lyts_unique)
        {
            if (std::abs(energies[i] - (it.get_system_energy() - min_energy))  < 0.0001)
            {
                counter += 1;
            }
            else
            {
                continue;
            }
        }
        energy_deg.insert(std::make_pair(energies[i], counter));
    }


    for (const auto &temp: temp_values)
    {

        if (occu_prop(energy_deg,temp) < convlevel)
        {
            return temp;
        }
        else if (std::abs(temp - temp_limit) < 0.001)
        {
            return temp_limit;
        }

        else
        {
            continue;
        }
    }
}

}  // namespace fiction

#endif  // FICTION_TEMPERATURE_HPP
