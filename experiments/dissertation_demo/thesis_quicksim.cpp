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
    using lyt_typ = sidb_100_cell_clk_lyt_siqad;

    const auto cell = read_sqd_layout<lyt_typ>(fmt::format("{}/dissertation_demo/{}", EXPERIMENTS_PATH, "and.sqd"));

    sidb_simulation_parameters sim_params{};
    sim_params.base      = 2;
    sim_params.mu_minus  = -0.32;
    sim_params.epsilon_r = 5.6;
    sim_params.lambda_tf = 5.0;
    sim_params.mu_minus  = -0.32;

    quicksim_params qs_params{};
    qs_params.simulation_parameters = sim_params;
    qs_params.iteration_steps       = 300;
    qs_params.alpha                 = 0.6;

    const auto result_qs = quicksim(cell, qs_params);
    const auto gs        = result_qs.value().groundstates();
}
