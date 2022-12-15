//
// Created by Jan Drewniok on 24.11.22.
//

#ifndef FICTION_CONSTANTS_HPP
#define FICTION_CONSTANTS_HPP

namespace fiction
{
namespace physical_sim_constants
{
/**
 * epsilon is the vacuum permittivity [F/m].
 */
constexpr double EPSILON = 8.854 * 1E-12;
/**
 * electric charge is the elementary charge [C].
 */
constexpr double ELECTRIC_CHARGE = 1.602 * 1E-19;
/**
 * The pop stability err is used for the physical simulation to avoid floating-point errors.
 */
constexpr double POP_STABILITY_ERR = 1E-6;

constexpr double PI = 3.14159265359;

}  // namespace physical_sim_constants

}  // namespace fiction

#endif  // FICTION_CONSTANTS_HPP
