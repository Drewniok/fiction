//
// Created by Jan Drewniok on 11.01.23.
//

#ifndef FICTION_EFFACCSIM_HPP
#define FICTION_EFFACCSIM_HPP

#include "fiction/algorithms/network_transformation/fanout_substitution.hpp"
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
struct speedysim_stats
{
    mockturtle::stopwatch<>::duration time_total{0};

    std::vector<charge_distribution_surface<Lyt>> valid_lyts{};
    uint64_t                                      num_dbs{0ull};
    std::vector<charge_distribution_surface<Lyt>> gs_lyt{};


    void report(std::ostream& out = std::cout) const
    {
        out << fmt::format("[i] total time  = {:.2f} secs\n", mockturtle::to_seconds(time_total));
        out << fmt::format("[i] num. sidbs = {} \n", num_dbs);
        out << fmt::format("[i] degeneracy of ground state = {} \n", gs_lyt.size());
    }
};

namespace detail
    {

    template <typename Lyt>
    class speedysim_impl
    {
      public:
        speedysim_impl(Lyt &lyt, const physical_params& phys_params, const uint64_t steps = 80, const double al = 0.7) :
                spsim{}, physparam{phys_params}, iteration_steps{steps}, alpha{al}, layout{lyt}
        {}


        void run()
        {
            // measure run time
            mockturtle::stopwatch stop{spsim.time_total};

            // instantiate the layout
            charge_distribution_surface chargelyt{layout,physparam};

            std::vector<charge_distribution_surface<Lyt>> result{};

            chargelyt.set_charge_states(sidb_charge_state::NEUTRAL);
            chargelyt.local_potential();
            chargelyt.system_energy();
            chargelyt.validity_check();

                if (chargelyt.get_validity())
                {
                    charge_distribution_surface<Lyt> lyt_new{chargelyt};
                    valid_layouts.push_back(lyt_new);
                }

                float best_energy = MAXFLOAT;
                auto  bound       = static_cast<int>(0.6 * chargelyt.num_cells());
                for (int z = 0; z < chargelyt.num_cells(); z++)
                {
                    for (int i = 0u; i < bound; i++)
                    {
                        std::vector<int> index_start = {i};
                        chargelyt.set_charge_states(sidb_charge_state::NEUTRAL);
                        chargelyt.assign_charge_state_index(i, sidb_charge_state::NEGATIVE);
                        chargelyt.local_potential();
                        chargelyt.system_energy();

                        for (int num = 0; num < chargelyt.num_cells() / 2 + 4; num++)
                        {
                            chargelyt.next_N(this->alpha, index_start);
                            chargelyt.validity_check();

                            if (chargelyt.get_validity() && (chargelyt.get_system_energy() <= best_energy))
                            {
                                charge_distribution_surface<Lyt> lyt_new{chargelyt};
                                valid_layouts.push_back(lyt_new);
                            }
                        }
                    }
                }
        }

        std::vector<charge_distribution_surface<Lyt>> get_valid_layouts()
        {
            return this->valid_layouts;
        }

      private:

        speedysim_stats &spsim;
        physical_params &physparam;
        uint64_t iteration_steps{};
        double alpha{};
        Lyt layout;
        std::vector<charge_distribution_surface<Lyt>> valid_layouts{};

    };

    }  // namespace detail



    template <typename Lyt, typename Ntk>
    Lyt speedsim(Lyt &lyt, speedysim_stats* pst = nullptr, const physical_params& physical_params_default = physical_params{}, const uint64_t iteration_steps = 70, const double alpha = 0.7)
    {

        detail::speedysim_impl<Lyt> p{lyt, physical_params_default, iteration_steps, alpha};

        p.run();

//        if (pst)
//        {
//            *pst = st;
//        }

        return result;
    }

}  // namespace fiction

#endif  // FICTION_EFFACCSIM_HPP
