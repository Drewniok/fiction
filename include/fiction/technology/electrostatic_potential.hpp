#ifndef FICTION_ELECTROSTATIC_POTENTIAL_HPP
#define FICTION_ELECTROSTATIC_POTENTIAL_HPP

//
// Created by Jan Drewniok on 12.12.22.
//

#include "fiction/algorithms/simulation_sidb/simulation_parameter.hpp"
#include "fiction/traits.hpp"

namespace fiction
{

template <typename potential = double, typename Dist = double>
potential potential_SiDB_pair(const Dist &dist)
{
    if (dist == static_cast<Dist>(0.0))
    {
        return static_cast<Dist>(0.0);
    }
    else
    {
        const simulation_parameter simulation_parameter {};
    return simulation_parameter.k / dist * std::exp(-dist / simulation_parameter.lambda_tf);
};
}


} // namespace fiction
#endif  // FICTION_ELECTROSTATIC_POTENTIAL_HPP
