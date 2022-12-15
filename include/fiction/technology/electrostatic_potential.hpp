#ifndef FICTION_ELECTROSTATIC_POTENTIAL_HPP
#define FICTION_ELECTROSTATIC_POTENTIAL_HPP

//
// Created by Jan Drewniok on 12.12.22.
//

#include "fiction/algorithms/simulation_sidb/simulation_parameters.hpp"
#include "fiction/traits.hpp"

namespace fiction
{
/**
 * The electrostatic potential for a given Euclidean distance is calculated.
 *
 * @tparam Potential Data type for the electrostatic potential.
 * @tparam Dist Euclidean distance.
 * @param dist Euclidean distance value.
 * @param k Coulomb constant (default value).
 * @param lambda_tf Thomas-Fermi screening distance (default value).
 * @return Electrostatic potential value.
 */
template <typename Potential = double, typename Dist = double>
Potential potential_sidb_pair(const Dist& dist, const simulation_params &sim_params = simulation_params{})
{
    static_assert(std::is_floating_point_v<Dist>, "Dist is not a floating-point type");
    static_assert(std::is_floating_point_v<Potential>, "Potential is not a floating-point type");
    if (dist == static_cast<Dist>(0.0))
    {
        return static_cast<Potential>(0.0);
    }

    return static_cast<Potential>(sim_params.k / dist * std::exp(-dist / sim_params.lambda_tf) * physical_sim_constants::ELECTRIC_CHARGE);
}

}  // namespace fiction

#endif  // FICTION_ELECTROSTATIC_POTENTIAL_HPP
