//
// Created by Jan Drewniok on 11.01.23.
//

#ifndef FICTION_EFFACCSIM_HPP
#define FICTION_EFFACCSIM_HPP

#include "fiction/algorithms/network_transformation/fanout_substitution.hpp"
#include "fiction/algorithms/simulation_sidb/get_energy_dist.hpp"
#include "fiction/algorithms/simulation_sidb/minimum_energy.hpp"
#include "fiction/io/print_layout.hpp"
#include "fiction/layouts/clocking_scheme.hpp"
#include "fiction/networks/views/edge_color_view.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/traits.hpp"
#include "fiction/utils/name_utils.hpp"
#include "fiction/utils/network_utils.hpp"
#include "fiction/utils/placement_utils.hpp"

#include <fmt/format.h>
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/node_map.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <mockturtle/views/topo_view.hpp>

#include <algorithm>
#include <optional>
#include <set>
#include <vector>

#if (PROGRESS_BARS)
#include <mockturtle/utils/progress_bar.hpp>
#endif

namespace fiction
{

template <typename Lyt>
struct quicksim_stats
{
    double                                        time_total{0};
    std::vector<charge_distribution_surface<Lyt>> valid_lyts{};

    void report(std::ostream& out = std::cout)
    {
        out << fmt::format("[i] total time  = {:.2f} millisecs\n", time_total / 1000);
          if (!get_statistics<Lyt>(valid_lyts).empty())
        {
            for (auto [energy, count] : get_statistics<Lyt>(valid_lyts))
            {
                out << fmt::format("[i] the lowest state energy is  = {:.4f} \n", minimum_energy(valid_lyts));
                out << fmt::format("energy: {} | occurance: {} \n", energy, count);
            }
        }
        else
        {
            std::cout << "no state found" << std::endl;
        }

        std::cout << "_____________________________________________________ \n";
    }
};

template <typename Lyt>
void quicksim(charge_distribution_surface<Lyt>& lyt, quicksim_stats<Lyt>& ps,
              const physical_params& phys_params = physical_params{}, const uint64_t iteration_steps = 80,
              const double alpha = 0.7)
{

    // instantiate the layout
    lyt.set_physical_parameters(phys_params);

    // measure run time
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<charge_distribution_surface<Lyt>> result{};

    lyt.set_charge_states(sidb_charge_state::NEUTRAL);
    lyt.local_potential();
    lyt.system_energy();
    lyt.validity_check();

    if (lyt.get_validity())
    {
        charge_distribution_surface<Lyt> lyt_new{lyt};
        ps.valid_lyts.push_back(lyt_new);
    }

    float best_energy = MAXFLOAT;
    auto  bound       = static_cast<int>(0.6 * lyt.num_cells());
    for (uint64_t z = 0u; z < iteration_steps; z++)
    {
        for (int i = 0u; i < bound; i++)
        {
            std::vector<int> index_start = {i};
            lyt.set_charge_states(sidb_charge_state::NEUTRAL);
            lyt.assign_charge_state_index(i, sidb_charge_state::NEGATIVE);
            lyt.local_potential();
            lyt.system_energy();

            for (int num = 0; num < lyt.num_cells() / 1.8; num++)
            {
                lyt.adjacent_search(alpha, index_start);
                lyt.validity_check();

                if (lyt.get_validity() && (lyt.get_system_energy() <= best_energy))
                {
                    charge_distribution_surface<Lyt> lyt_new{lyt};
                    ps.valid_lyts.push_back(lyt_new);
                }
            }
        }
    }
    auto end      = std::chrono::high_resolution_clock::now();
    ps.time_total = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

}  // namespace fiction

#endif  // FICTION_EFFACCSIM_HPP
