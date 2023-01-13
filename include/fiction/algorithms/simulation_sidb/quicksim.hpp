//
// Created by Jan Drewniok on 22.12.22.
//

#ifndef FICTION_QUICKSIM_HPP
#define FICTION_QUICKSIM_HPP

#include "fiction/technology/charge_distribution_surface.hpp"

namespace fiction::detail

{

template <typename Lyt>
std::vector<charge_distribution_surface<Lyt>> quicksim(charge_distribution_surface<Lyt>& lyt,
                                                       const int iteration_steps = 70, const double alpha = 0.7)

{
    std::vector<charge_distribution_surface<Lyt>> collect{};

    lyt.set_charge_states(sidb_charge_state::NEUTRAL);
    lyt.local_potential();
    lyt.system_energy();
    lyt.validity_check();

    if (lyt.get_validity())
    {
        charge_distribution_surface<Lyt> lyt_new{lyt};
        collect.push_back(lyt_new);
    }

    float best_energy = MAXFLOAT;
    auto  bound       = static_cast<int>(0.6 * lyt.num_cells());
    for (int z = 0; z < iteration_steps; z++)
    {
        for (int i = 0u; i < bound; i++)
        {
            std::vector<int> index_start = {i};
            lyt.set_charge_states(sidb_charge_state::NEUTRAL);
            lyt.assign_charge_state_index(i, sidb_charge_state::NEGATIVE);
            lyt.local_potential();
            lyt.system_energy();

            for (int num = 0; num < lyt.num_cells() / 2 + 4; num++)
            {
                lyt.next_N(alpha, index_start);
                lyt.validity_check();
                lyt.chargeconf_to_index();

                if (lyt.get_validity() && (lyt.get_system_energy() <= best_energy))
                {
                    charge_distribution_surface<Lyt> lyt_new{lyt};
                    collect.push_back(lyt_new);
                }
            }
        }
    }

    return collect;
}

}  // namespace fiction::detail

#endif  // FICTION_QUICKSIM_HPP
