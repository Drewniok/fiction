//
// Created by Jan Drewniok on 18.01.23.
//

#ifndef FICTION_MINIMUM_ENERGY_HPP
#define FICTION_MINIMUM_ENERGY_HPP

#include <fiction/technology/charge_distribution_surface.hpp>

namespace fiction
{

/**
 * Computes the minimum energy of a vector of charge_distribution_surface objects.
 *
 * @tparam Lyt template parameter for the charge_distribution_surface class.
 * @param result vector of charge_distribution_surface objects.
 * @return double value of the minimum energy found in the input vector.
 */
template <typename Lyt>
double minimum_energy(const std::vector<charge_distribution_surface<Lyt>>& result)
{
    return std::accumulate(result.begin(), result.end(), std::numeric_limits<double>::max(),
                           [](double a, const auto& lyt) { return std::min(a, lyt.get_system_energy()); });
}

}  // namespace fiction

#endif  // FICTION_MINIMUM_ENERGY_HPP
