//
// Created by Jan Drewniok on 18.12.22.
//

#ifndef FICTION_EXGS_HPP
#define FICTION_EXGS_HPP

#include "fiction/technology/charge_distribution_surface.hpp"

namespace fiction::detail

{

// template <typename Lyt>
// class ground_state_search : public charge_distribution_surface<Lyt>
//{
//   public:
//     explicit ground_state_search(const charge_distribution_surface<Lyt> &lyt) {
//         lyt.chargeconf_to_index(), lyt.index_to_chargeconf(), lyt.initialize_sidb_distance_matrix(),
//         lyt.initialize_sidb_potential_matrix(), max_charge_index = static_cast<uint64_t>(std::pow(3,
//         lyt.num_charges())-1);
//     };

template <typename Lyt>
std::vector<std::pair<charge_distribution_surface<Lyt>, double>> GS(charge_distribution_surface<Lyt>& lyt)
{
    uint64_t                               max_charge_index{};
    std::vector<std::pair<charge_distribution_surface<Lyt>, double>> collect{};

    lyt.initialize_charge_state(sidb_charge_state::NEGATIVE);
    lyt.initialize_sidb_distance_matrix();
    lyt.initialize_sidb_potential_matrix();
    lyt.local_potential();
    lyt.system_energy();
    lyt.chargeconf_to_index();
    lyt.index_to_chargeconf();
    max_charge_index = static_cast<uint64_t>(std::pow(3,lyt.num_charges())-1);


   //for (int i = 0; i<50; i++)
    while (lyt.get_charge_index().first <= (max_charge_index-1))
    {
        lyt.local_potential();
        lyt.system_energy();
        lyt.chargeconf_to_index();
        lyt.index_to_chargeconf();
        //std::cout << lyt.get_charge_index().first << std::endl;
        //std::cout << lyt.get_system_energy() << std::endl;
        lyt.validity_check();
        //std::cout << lyt.get_validity() << std::endl;
        if (lyt.get_validity()==0)
        {
//            std::cout << lyt.get_validity() << std::endl;
//            std::cout << lyt.get_system_energy() << std::endl;
//            std::cout << lyt.get_charge_index().first << std::endl;
            collect.push_back(std::make_pair(lyt,lyt.get_system_energy()));
            //best_layout.insert(charge_distribution_surface<Lyt>(lyt));
            lyt.increase_charge_index();
        }
        else{
        lyt.increase_charge_index();
        }
    }
    return collect;
};

//template <typename T>
//GS(const charge_distribution_surface<T>&) -> GS<T>;

};  // namespace fiction::detail



#endif  // FICTION_EXGS_HPP
