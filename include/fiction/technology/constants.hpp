//
// Created by Jan Drewniok on 24.11.22.
//

#ifndef FICTION_TECHNOLOGY_PARAMETER_HPP
#define FICTION_TECHNOLOGY_PARAMETER_HPP

namespace fiction
{
namespace physical_sim_constants
{
    /**
     * vacuum permittivity
     */
    const double EPSILON           = 8.854 * 1E-12;
    /**
     * electric charge
     */
    const double ELECTRIC_CHARGE                 = 1.602 * 1E-19;
    /**
     * stability threashold
     */
    const double POP_STABILITY_ERR = 1E-6;

}; // namespace physical_sim_constants

}  // namespace fiction

#endif  // FICTION_TECHNOLOGY_PARAMETER_HPP
