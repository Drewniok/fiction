//
// Created by Jan Drewniok on 24.11.22.
//

#ifndef FICTION_TECHNOLOGY_PARAMETER_HPP
#define FICTION_TECHNOLOGY_PARAMETER_HPP

namespace fiction
{
struct Constants
{
    static constexpr float epsilon           = 8.854f * 1E-12f;
    static constexpr float e                 = 1.602f * 1E-19f;
    static constexpr float POP_STABILITY_ERR = 1E-6f;
    static constexpr float lat_a = 3.84f * 1E-10f;  // lattice vector in x, angstroms (intra dimer row)
    static constexpr float lat_b = 7.68f * 1E-10f;  // lattice vector in y, angstroms (inter dimer row)
    static constexpr float lat_c = 2.25f * 1E-10f;  // dimer pair separation, angstroms
};

}

#endif  // FICTION_TECHNOLOGY_PARAMETER_HPP
