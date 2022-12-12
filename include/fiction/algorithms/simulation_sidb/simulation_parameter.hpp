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
     explicit simulation_parameter(const double relative_permittivity = 5.6, const double screening_distance = 5.0 * 1E-9, const double mu_ = -0.32, const double a = 3.84 * 1E-10, const double b = 7.68 * 1E-10,
                                            const double c = 2.25 * 1E-10) noexcept:
            lat_a{a},
            lat_b{b},
            lat_c{c},
            epsilon_r{relative_permittivity},
            k{1.0 / (4.0 * 3.141592653 * fiction::constants::epsilon * epsilon_r)},
            lambda_tf{screening_distance},
            mu{mu_},
            mu_p{mu - 0.59}

    {}


    /**
     * lattice vector in x, angstroms (intra dimer row)
     */
    const double lat_a;
    /**
     * lattice vector in y, angstroms (inter dimer row)
     */
    const double lat_b;
    /**
     * dimer pair separation, angstroms
     */
    const double lat_c;
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
