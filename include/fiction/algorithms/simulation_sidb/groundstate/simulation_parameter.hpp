//
// Created by Jan Drewniok on 24.11.22.
//

#ifndef FICTION_SIMULATION_PARAMETER_HPP
#define FICTION_SIMULATION_PARAMETER_HPP

#include "constants.hpp"

namespace fiction
{
struct SimulationParameter
{
    static constexpr float epsilon_screen = 5.6f;
    static constexpr float tf             = 5.0f * 1E-9f;
    static constexpr float k              = 1.0f / (4.0f * 3.141592653f * fiction::Constants::epsilon * epsilon_screen);
    static constexpr float mu             = -0.32f;
    static constexpr float mu_p           = mu - 0.59f;
};
}  // namespace fiction

#endif  // FICTION_SIMULATION_PARAMETER_HPP
