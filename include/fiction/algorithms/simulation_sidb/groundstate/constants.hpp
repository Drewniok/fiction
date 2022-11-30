//
// Created by Jan Drewniok on 24.11.22.
//

#ifndef FICTION_TECHNOLOGY_PARAMETER_HPP
#define FICTION_TECHNOLOGY_PARAMETER_HPP

namespace fiction
{
struct constants
{
    constexpr explicit constants(const double a = 3.84 * 1E-10, const double b = 7.68 * 1E-10,
                                 const double c = 2.25 * 1E-10) noexcept :

            lat_a{a},
            lat_b{b},
            lat_c{c}
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
     * vacuum permittivity
     */
    const double epsilon           = 8.854 * 1E-12;
    /**
     * electric charge
     */
    const double e                 = 1.602 * 1E-19;
    /**
     * stability threashold
     */
    const double POP_STABILITY_ERR = 1E-6;
};

}  // namespace fiction

#endif  // FICTION_TECHNOLOGY_PARAMETER_HPP
