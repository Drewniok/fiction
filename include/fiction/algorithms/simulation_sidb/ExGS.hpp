//
// Created by Jan Drewniok on 18.12.22.
//

#ifndef FICTION_EXGS_HPP
#define FICTION_EXGS_HPP

#include "fiction/technology/charge_distribution_surface.hpp"

namespace fiction::detail

{
/**
 *  All metastable and physically valid charge distribution layouts are computed, stored in a vector and returned.
 *
 * @tparam Lyt cell-level layout.
 * @param lyt charge distribution layout.
 * @return a vector of different charge distribution layouts, all of which satisfy the validity test.
 */
template <typename Lyt>
std::unordered_map<double, charge_distribution_surface<Lyt>> metastable_layouts(charge_distribution_surface<Lyt>& lyt)

{
    std::unordered_map<double, charge_distribution_surface<Lyt>> collect{};

    lyt.initialize_sidb_distance_matrix();
    lyt.initialize_sidb_potential_matrix();

    lyt.local_potential();
    lyt.system_energy();
    lyt.validity_check();

    if (lyt.get_validity() == 1)
    {
        charge_distribution_surface<Lyt> lyt_new{lyt};
        collect.insert(std::pair(lyt_new.get_charge_index().first, lyt_new));
    }

    while (lyt.get_charge_index().first <= lyt.get_max_charge_index() - 1)
    {
        lyt.increase_charge_index();
        lyt.local_potential();
        lyt.system_energy();
        lyt.validity_check();

        if (lyt.get_charge_index().first % 1000 == 0)
        {
            std::cout << lyt.get_charge_index().first << std::endl;
        }
        if (lyt.get_validity() == 1)
        {
            charge_distribution_surface<Lyt> lyt_new{lyt};
            collect.insert(std::pair(lyt_new.get_charge_index().first, lyt_new));
        }
    }
    return collect;
};

};  // namespace fiction::detail

#endif  // FICTION_EXGS_HPP
