//
// Created by Jan Drewniok on 22.12.22.
//

#ifndef FICTION_NEW_APPROACH_HPP
#define FICTION_NEW_APPROACH_HPP

#include "fiction/technology/charge_distribution_surface.hpp"

namespace fiction::detail

{

template <typename Lyt>

std::unordered_map<double, charge_distribution_surface<Lyt>> Sim(charge_distribution_surface<Lyt>& lyt, const int iteration_steps = 3, const double alpha = 0.7)
{
    std::unordered_map<double, charge_distribution_surface<Lyt>> collect{};
    lyt.initialize_sidb_distance_matrix();
    lyt.initialize_sidb_potential_matrix();
    lyt.foreach_charge_state([&lyt](const auto& cd) { lyt.assign_charge_state(cd.first, sidb_charge_state::NEUTRAL); });
    lyt.local_potential();
    lyt.system_energy();
    lyt.validity_check();

    if (lyt.get_validity())
    {
        charge_distribution_surface<Lyt> lyt_new{lyt};
        collect.insert(std::pair(lyt_new.get_charge_index().first, lyt_new));
    }

    for (int z = 0; z < iteration_steps; z++)
    {
        for (auto& it : lyt.get_charge_coordinates())
        {
            lyt.foreach_charge_state([&lyt](const auto& cd)
                                     { lyt.assign_charge_state(cd.first, sidb_charge_state::NEUTRAL); });
            lyt.assign_charge_state(it.first, sidb_charge_state::NEGATIVE);
            lyt.local_potential();
            lyt.system_energy();

            for (int num = 0; num < lyt.num_cells()-1; num++)
            {
                lyt.next_N(alpha);
                lyt.chargeconf_to_index();
                lyt.local_potential();
                lyt.system_energy();
                lyt.validity_check();


                if (lyt.get_validity())
                {
                    charge_distribution_surface<Lyt> lyt_new{lyt};
                    collect.insert(std::pair(lyt_new.get_charge_index().first, lyt_new));
                }
            }
        }
    }
    return collect;
}

} //namespace fiction::detail



#endif  // FICTION_NEW_APPROACH_HPP