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
std::vector<charge_distribution_surface<Lyt>> metastable_layouts(charge_distribution_surface<Lyt>& lyt)

{
    std::vector<charge_distribution_surface<Lyt>> collect{};

    lyt.initialize_sidb_distance_matrix();
    lyt.initialize_sidb_potential_matrix();

    lyt.increase_charge_index();
    lyt.local_potential();
    lyt.system_energy();
    lyt.validity_check();

    if (lyt.get_validity() == 1)
    {
        collect.emplace_back(lyt);
    }

    auto max_charge_index = static_cast<uint64_t>(std::pow(lyt.get_charge_index().second, lyt.num_cells()) - 1);

    while (lyt.get_charge_index().first <= max_charge_index - 1)
    {
        lyt.increase_charge_index();
        lyt.local_potential();
        lyt.system_energy();
        lyt.validity_check();

        if (lyt.get_validity() == 1)
        {
            collect.emplace_back(lyt);
        }
    }
    return collect;
};

};  // namespace fiction::detail

#endif  // FICTION_EXGS_HPP
