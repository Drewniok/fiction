//
// Created by Jan Drewniok on 25.08.25.
//

#include "fiction/algorithms/simulation/sidb/quicksim.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/io/read_sqd_layout.hpp"
#include "fiction/types.hpp"

using namespace fiction;

int main()  // NOLINT
{
    const auto                       cell = read_sqd_layout<sidb_100_cell_clk_lyt_siqad>("and.sqd");
    const sidb_simulation_parameters sim_params{2, -0.32, 5.6, 5.0};
    const quicksim_params            qs_params{sim_params, 300, 0.6};
    const auto                       results = quicksim(cell, qs_params);
    const auto                       gs      = results.value().groundstates();
}
