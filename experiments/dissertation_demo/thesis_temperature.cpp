//
// Created by Jan Drewniok on 25.08.25.
//

#include "fiction/algorithms/iter/bdl_input_iterator.hpp"
#include "fiction/algorithms/simulation/sidb/critical_temperature.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/io/read_sqd_layout.hpp"
#include "fiction/types.hpp"

#include <vector>

using namespace fiction;

int main()  // NOLINT
{
    const auto cell = read_sqd_layout<sidb_100_cell_clk_lyt_siqad>(
        fmt::format("{}/dissertation_demo/{}", EXPERIMENTS_PATH, "and.sqd"));

    const sidb_simulation_parameters  sim_params{2, -0.32};
    const critical_temperature_params ct_params{{sim_params}};
    critical_temperature_stats        ct_stats{};
    const auto ct = critical_temperature_gate_based(cell, std::vector<tt>{create_and_tt()}, ct_params, &ct_stats);
    std::cout << "Critical Temperature: " << ct << " K" << std::endl;
}
