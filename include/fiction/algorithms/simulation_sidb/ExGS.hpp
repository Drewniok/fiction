//
// Created by Jan Drewniok on 18.12.22.
//

#ifndef FICTION_EXGS_HPP
#define FICTION_EXGS_HPP

#include "fiction/technology/charge_distribution_surface.hpp"

namespace fiction::detail

{

template <typename Lyt>
std::vector<std::pair<charge_distribution_surface<Lyt>, double>> GS(charge_distribution_surface<Lyt>& lyt)

{
    uint64_t                                                         max_charge_index{};
    std::vector<std::pair<charge_distribution_surface<Lyt>, double>> collect{};

    lyt.initialize_sidb_distance_matrix();
    lyt.initialize_sidb_potential_matrix();

    max_charge_index = static_cast<uint64_t>(std::pow(3, lyt.num_cells()) - 1);

    while (lyt.get_charge_index().first <= max_charge_index)
    {
        lyt.local_potential();
        lyt.system_energy();
        lyt.validity_check();

        if (lyt.get_validity() == 1)
        {
            charge_distribution_surface new_lyt{lyt};
            collect.push_back(std::pair(new_lyt, lyt.get_system_energy()));
        }
        lyt.increase_charge_index();
    }
    return collect;
};

};  // namespace fiction::detail

#endif  // FICTION_EXGS_HPP
