//
// Created by Jan Drewniok on 24.11.22.
//

#ifndef FICTION_SIMULATION_PARAMETERS_HPP
#define FICTION_SIMULATION_PARAMETERS_HPP

#include "fiction/technology/constants.hpp"

namespace fiction
{
struct physical_params
{
    explicit physical_params(const uint8_t base_number = 3, const double relative_permittivity = 5.6, const double screening_distance = 5.0 * 1E-9,
                               const double mu_minus = -0.32, const double a = 3.84 * 1E-10,
                               const double b = 7.68 * 1E-10, const double c = 2.25 * 1E-10) noexcept :
            lat_a{a},
            lat_b{b},
            lat_c{c},
            epsilon_r{relative_permittivity},
            k{1.0 / (4.0 * fiction::physical_sim_constants::PI * fiction::physical_sim_constants::EPSILON * epsilon_r)},
            lambda_tf{screening_distance},
            mu{mu_minus},
            mu_p{mu - 0.59},
            base{base_number}

    {}

    /**
     * lat_a is the lattice vector in x-direction.
     */
    const double lat_a;
    /**
     * lat_b is the lattice vector in y-direction.
     */
    const double lat_b;
    /**
     * lat_c is the dimer pair separation.
     */
    const double lat_c;
    /**
     * epsilon_r is the electric permittivity. It is a material specific number.
     */
    const double epsilon_r;
    /**
     * k is the Coulomb constant and is inversely proportinal to the electric permittivity.
     */
    const double k;
    /**
     * lambda_tf is the Thomas-Fermi screening distance.
     */
    const double lambda_tf;
    /**
     * µ- is the energy transition level (0/-)
     */
    const double mu;
    /**
     * µ+ is the energy transition level (+/0)
     */
    const double mu_p;
    /**
     * base can be either 2 or 3 and describes the assumed number of charge states of one SiDB.
     * It often makes sense to assume only negatively and neutrally charged SiDBs.
     */
    const uint8_t base;
};
}  // namespace fiction

#endif  // FICTION_SIMULATION_PARAMETERS_HPP
