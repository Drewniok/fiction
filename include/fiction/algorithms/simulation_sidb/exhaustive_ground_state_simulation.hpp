//
// Created by Jan Drewniok on 18.12.22.
//

#ifndef FICTION_EXHAUSTIVE_GROUND_STATE_SIMULATION_HPP
#define FICTION_EXHAUSTIVE_GROUND_STATE_SIMULATION_HPP

#include "fiction/algorithms/simulation_sidb/energy_distribution.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"

#include <mockturtle/utils/stopwatch.hpp>

#include <iostream>

namespace fiction
{

template <typename Lyt>
struct exgs_stats
{
    mockturtle::stopwatch<>::duration             time_total{0};
    std::vector<charge_distribution_surface<Lyt>> valid_lyts{};

    void report(std::ostream& out = std::cout) const
    {
        out << fmt::format("total time  = {:.2f} secs\n", mockturtle::to_seconds(time_total));
        if (!valid_lyts.empty())
        {
            for (auto [energy, count] : get_statistics<Lyt>(valid_lyts))
            {
                out << fmt::format("energy: {} | occurance: {} \n", energy, count);
            }
            out << fmt::format("the ground state energy is  = {:.4f} \n", minimum_energy(valid_lyts));
        }
        else
        {
            std::cout << "no state found | if two state simulation is used, continue with three state" << std::endl;
        }

        out << fmt::format("{} phyiscally valid charge states were found \n", valid_lyts.size());
        std::cout << "_____________________________________________________ \n";
    }
};

/**
 *  All metastable and physically valid charge distribution layouts are computed, stored in a vector and returned.
 *
 * @tparam Lyt cell-level layout.
 * @param lyt charge distribution layout.
 * @return a vector of different charge distribution layouts, all of which satisfy the validity test.
 */
template <typename Lyt>
void exgs(charge_distribution_surface<Lyt>& lyt, exgs_stats<Lyt>& ps,
          const physical_params& phys_params = physical_params{})
{
    mockturtle::stopwatch stop{ps.time_total};
    ps = exgs_stats<Lyt>{};
    lyt.set_physical_parameters(phys_params);
    lyt.set_charge_states(sidb_charge_state::NEGATIVE);
    lyt.chargeconf_to_index();

    while (lyt.get_charge_index().first <= lyt.get_max_charge_index() - 1)
    {
        lyt.local_potential();
        lyt.system_energy();
        lyt.validity_check();

        if (lyt.get_validity())
        {
            charge_distribution_surface<Lyt> lyt_new{lyt};
            ps.valid_lyts.push_back(lyt_new);
        }

        lyt.increase_charge_index();
    }

    lyt.local_potential();
    lyt.system_energy();
    lyt.validity_check();

    if (lyt.get_validity())
    {
        charge_distribution_surface<Lyt> lyt_new{lyt};
        ps.valid_lyts.push_back(lyt_new);
    }
}
}  // namespace fiction

#endif  // FICTION_EXHAUSTIVE_GROUND_STATE_SIMULATION_HPP
