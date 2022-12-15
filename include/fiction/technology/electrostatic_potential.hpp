#ifndef FICTION_ELECTROSTATIC_POTENTIAL_HPP
#define FICTION_ELECTROSTATIC_POTENTIAL_HPP

//
// Created by Jan Drewniok on 12.12.22.
//

#include "fiction/algorithms/simulation_sidb/simulation_parameter.hpp"
#include "fiction/traits.hpp"

namespace fiction
{
/**
 * The electrostatic potential for a given Euclidean distance is calculated.
 *
 * @tparam Dist Euclidean distance, floating-point type for the distance.
 * @param Potential Floating-point type for the electrostatic potential.
 * @return electrostatic potential value.
 */
template <typename Potential = double, typename Dist = double>
Potential potential_sidb_pair(const Dist& dist, decltype(simulation_params{}.k) k = simulation_params{}.k,
                              decltype(simulation_params{}.lambda_tf) lambda_tf = simulation_params{}.lambda_tf)
{
    if (dist == static_cast<Dist>(0.0))
    {
        return static_cast<Potential>(0.0);
    }

    return static_cast<Potential>(-k / dist * std::exp(-dist / lambda_tf) * physical_sim_constants::ELECTRIC_CHARGE);
}

}  // namespace fiction

#endif  // FICTION_ELECTROSTATIC_POTENTIAL_HPP
