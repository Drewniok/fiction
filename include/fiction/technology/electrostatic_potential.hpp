#ifndef FICTION_ELECTROSTATIC_POTENTIAL_HPP
#define FICTION_ELECTROSTATIC_POTENTIAL_HPP

//
// Created by Jan Drewniok on 12.12.22.
//

#include "fiction/algorithms/simulation_sidb/simulation_parameter.hpp"
#include "fiction/traits.hpp"

namespace fiction
{

template <typename Potential = double, typename Dist = double>
Potential potential_sidb_pair(const Dist& dist, decltype(simulation_parameter{}.k) k = simulation_parameter{}.k,
                              decltype(simulation_parameter{}.lambda_tf) lambda_tf = simulation_parameter{}.lambda_tf)
{
    if (dist == static_cast<Dist>(0.0))
    {
        return static_cast<Potential>(0.0);
    }

    return static_cast<Potential>(k / dist * std::exp(-dist / lambda_tf));
}

}  // namespace fiction

#endif  // FICTION_ELECTROSTATIC_POTENTIAL_HPP
