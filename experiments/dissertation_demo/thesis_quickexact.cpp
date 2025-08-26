//
// Created by Jan Drewniok on 25.08.25.
//

#include "fiction/algorithms/iter/bdl_input_iterator.hpp"
#include "fiction/algorithms/simulation/sidb/quickexact.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/io/read_sqd_layout.hpp"
#include "fiction/layouts/coordinates.hpp"
#include "fiction/types.hpp"

using namespace fiction;

int main()  // NOLINT
{
    const auto                              cell = read_sqd_layout<sidb_100_cell_clk_lyt_siqad>("and.sqd");
    const sidb_simulation_parameters        sim_params{2, -0.32};
    const quickexact_params<siqad::coord_t> qe_params{sim_params};
    const auto                              results = quickexact(cell, qe_params);
    const auto                              gs      = results.groundstates();
}
