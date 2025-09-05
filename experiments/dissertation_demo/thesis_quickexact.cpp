//
// Created by Jan Drewniok on 25.08.25.
//

#include "fiction/algorithms/simulation/sidb/quickexact.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/io/read_sqd_layout.hpp"
#include "fiction/layouts/coordinates.hpp"
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

    quickexact_params<siqad::coord_t> qe_params{};
    qe_params.simulation_parameters = sim_params;

    const auto result_qe = quickexact(cell, qe_params);
    const auto gs_qe     = result_qe.groundstates();
}
