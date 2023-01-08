//
// Created by Jan Drewniok on 22.12.22.
//

#ifndef FICTION_NEW_APPROACH_HPP
#define FICTION_NEW_APPROACH_HPP

#include "fiction/technology/charge_distribution_surface.hpp"

namespace fiction::detail

{

template <typename Lyt>

std::unordered_map<double, charge_distribution_surface<Lyt>> Sim(charge_distribution_surface<Lyt>& lyt, const int iteration_steps = 80, const double alpha = 0.7)
{
    std::unordered_map<double, charge_distribution_surface<Lyt>> collect{};
    lyt.initialize_sidb_distance_matrix_new();
    lyt.initialize_sidb_potential_matrix_new();

    lyt.set_all_neutral();
    lyt.local_potential();
    lyt.system_energy();
    lyt.validity_check();

    if (lyt.get_validity())
    {
        charge_distribution_surface<Lyt> lyt_new{lyt};
        collect.insert(std::pair(lyt_new.get_charge_index().first, lyt_new));
    }

    float best_energy = MAXFLOAT;
    for (int z = 0; z < iteration_steps; z++)
    {
        for (int i = 0u; i < lyt.num_cells()-15; i++)
        {
            std::vector<int> index_start = {i};
            lyt.set_all_neutral();
            lyt.assign_charge_state(i, sidb_charge_state::NEGATIVE);
            lyt.local_potential();
            lyt.system_energy();

            for (int num = 0; num < lyt.num_cells()/2+4; num++)
            {
                lyt.next_N_new(alpha, index_start);
                //lyt.chargeconf_to_index();
                lyt.local_potential();
                lyt.system_energy();
                lyt.validity_check();


                if (lyt.get_validity() && (lyt.get_system_energy() < best_energy))
                {
                    charge_distribution_surface<Lyt> lyt_new{lyt};
                    collect.insert(std::pair(z, lyt_new));
                    break;
                }
            }
        }
    }
    return collect;
}

} //namespace fiction::detail



#endif  // FICTION_NEW_APPROACH_HPP
