//
// Created by Jan Drewniok on 25.08.25.
//

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

    sidb_simulation_parameters sim_params{};

    sim_params.base      = 2;
    sim_params.mu_minus  = -0.32;
    sim_params.epsilon_r = 5.6;
    sim_params.lambda_tf = 5.0;

    critical_temperature_params ct_params{};
    ct_params.operational_params.simulation_parameters = sim_params;

    critical_temperature_stats ct_stats{};
    const auto ct = critical_temperature_gate_based(cell, std::vector<tt>{create_and_tt()}, ct_params, &ct_stats);
    std::cout << "Critical Temperature: " << ct << " K" << std::endl;
}
