//
// Created by Jan Drewniok on 13.12.22.
//

#ifndef FICTION_VALIDITY_CHECK_HPP
#define FICTION_VALIDITY_CHECK_HPP

#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/algorithms/simulation_sidb/local_potential.hpp>
#include <fiction/algorithms/simulation_sidb/potential_matrix.hpp>

namespace fiction
{
template <typename Lyt, typename Potential = double>
bool validity_check(const charge_distribution_surface<Lyt> & lyt, const local_pot<Lyt, Potential> &loc_pot, const potential_matrix<charge_distribution_surface<Lyt>, Potential> &pot_mat, double mu = simulation_parameter{}.mu)
{
    bool valid=false;
    for (auto& it:loc_pot)
    {
        valid = (((lyt.get_charge_state(it.first) == sidb_charge_state::NEGATIVE) &&
                  ((it.second + mu) < constants::POP_STABILITY_ERR)) ||
                 ((lyt.get_charge_state(it.first) == sidb_charge_state::POSITIVE) &&
                  ((it.second + simulation_parameter{}.mu_p) > constants::POP_STABILITY_ERR)) ||
                 ((lyt.get_charge_state(it.first) == sidb_charge_state::NEUTRAL) &&
                  ((it.second + mu) > constants::POP_STABILITY_ERR) &&
                  (it.second + simulation_parameter{}.mu_p) < constants::POP_STABILITY_ERR));

        if (!valid)
        {
            return false;
        }
    }

    auto hopDel = [&lyt, &loc_pot, &pot_mat](const cell<Lyt> &c1,const cell<Lyt> &c2)
    {
        //std::vector<int> charge_sign_saved = chargesign;
        //float E_ori = this->system_energy_vec(charge_sign_saved);

        int dn_i = (lyt.get_charge_state(c1) == sidb_charge_state::NEGATIVE) ? 1 : -1;
        int dn_j = -dn_i;

        //        return v_local[i] * dn_i + v_local[j] * dn_j +
        //               v_ij(i, j) * ((chargesign[i] + dn_i) * (chargesign[j] + dn_j) - chargesign[i] * chargesign[j]);

        return loc_pot.at(c1)*dn_i + loc_pot.at(c1)*dn_j - pot_mat.at(std::make_pair(c1,c2))*1;
    };

    for (auto &it:loc_pot)
    {
        // do nothing with DB+
        if (lyt.get_charge_state(it.first) == sidb_charge_state::POSITIVE) {
            continue;
}

        for (auto &it_second:loc_pot)
        {

            // attempt hops from more negative charge states to more positive ones
            auto E_del = hopDel(it.first, it_second.first);
            if ((transform_to_sign(lyt.get_charge_state(it.first)) > transform_to_sign(lyt.get_charge_state(it.first))) && (E_del < -constants::POP_STABILITY_ERR))
            {
                return false;
            }
        }
    }
    return true;

}

} // namespace fiction

#endif  // FICTION_VALIDITY_CHECK_HPP
