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
std::pair<uint32_t, std::unordered_map<double, charge_distribution_surface<Lyt>>> metastable_layouts(charge_distribution_surface<Lyt>& lyt)
    //void metastable_layouts(charge_distribution_surface<Lyt>& lyt)

{

    auto                                                         t_start = std::chrono::high_resolution_clock::now();

    std::unordered_map<double, charge_distribution_surface<Lyt>> collect{};

    while (lyt.get_charge_index().first <= lyt.get_max_charge_index() - 1)
    {
        lyt.local_potential();
        lyt.system_energy();
        lyt.validity_check();

        if (lyt.get_validity())
        {
            charge_distribution_surface<Lyt> lyt_new{lyt};
            collect.insert(std::pair(lyt_new.get_charge_index().first, lyt_new));
        }

        lyt.increase_charge_index();
    }

    lyt.local_potential();
    lyt.system_energy();

    lyt.validity_check();

    if (lyt.get_validity())
    {
       charge_distribution_surface<Lyt> lyt_new{lyt};
       collect.insert(std::pair(lyt_new.get_charge_index().first, lyt_new));
    }

    auto t_end          = std::chrono::high_resolution_clock::now();
    auto elapsed        = t_end - t_start;
    auto diff_first     = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();

return std::pair(diff_first,collect);

};

};  // namespace fiction::detail

#endif  // FICTION_EXGS_HPP
