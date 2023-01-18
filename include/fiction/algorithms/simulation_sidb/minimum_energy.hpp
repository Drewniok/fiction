//
// Created by Jan Drewniok on 18.01.23.
//

#ifndef FICTION_MINIMUM_ENERGY_HPP
#define FICTION_MINIMUM_ENERGY_HPP

#include <fiction/technology/charge_distribution_surface.hpp>

namespace fiction
{

/**
* @brief Computes the minimum energy of a vector of charge_distribution_surface objects
* @tparam Lyt template parameter for the charge_distribution_surface class
* @param result vector of charge_distribution_surface objects
* @return double value of the minimum energy found in the input vector
*/
template <typename Lyt>
double minimum_energy(const std::vector<charge_distribution_surface<Lyt>>& result)
{
    auto min_energy = std::numeric_limits<double>::max();
    for (const auto& lyt : result)
    {
        if (lyt.get_system_energy() < min_energy)
        {
            min_energy = lyt.get_system_energy();
        }
    }
    return min_energy;
}
}  // namespace fiction

#endif  // FICTION_MINIMUM_ENERGY_HPP
