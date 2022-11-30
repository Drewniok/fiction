//
// Created by Jan Drewniok on 24.11.22.
//

#ifndef FICTION_SIMULATION_PARAMETER_HPP
#define FICTION_SIMULATION_PARAMETER_HPP

#include "constants.hpp"

namespace fiction
{
struct simulation_parameter
{
    constexpr explicit simulation_parameter(const double relative_permittivity = 5.6, const double screening_distance = 5.0 * 1E-9, const double mu_ = -0.32) noexcept:
            epsilon_r{relative_permittivity},
            k{1.0 / (4.0 * 3.141592653 * fiction::constants().epsilon * epsilon_r)},
            lambda_tf{screening_distance},
            mu{mu_},
            mu_p{mu - 0.59}
    {}
    /**
     * Electric permittivity.
     */
    const double epsilon_r;
    /**
     * Coulomb constant.
     */
    const double k;
    /**
     * Thomas-Fermi screening distance in nm.
     */
    const double lambda_tf;
    /**
     * µ-
     */
    const double mu;
    /**
     * µ+
     */
    const double mu_p;
};
}  // namespace fiction

#endif  // FICTION_SIMULATION_PARAMETER_HPP
