//
// Created by Jan Drewniok on 25.08.25.
//

#include "fiction/algorithms/physical_design/design_sidb_gates.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_engine.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/io/read_sqd_layout.hpp"
#include "fiction/layouts/coordinates.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/truth_table_utils.hpp"

#include <vector>

using namespace fiction;

int main()  // NOLINT
{
    using lyt_typ = sidb_100_cell_clk_lyt;

    const auto skeleton =
        read_sqd_layout<lyt_typ>(fmt::format("{}/dissertation_demo/{}", EXPERIMENTS_PATH, "skeleton_2i1o.sqd"));

    sidb_simulation_parameters sim_params{};
    sim_params.base      = 2;
    sim_params.mu_minus  = -0.32;
    sim_params.epsilon_r = 5.6;
    sim_params.lambda_tf = 5.0;
    sim_params.mu_minus  = -0.32;

    design_sidb_gates_params<offset::ucoord_t> design_params{};
    design_params.operational_params.simulation_parameters = sim_params;

    design_params.operational_params.sim_engine   = sidb_simulation_engine::QUICKEXACT;
    design_params.operational_params.op_condition = is_operational_params::operational_condition::REJECT_KINKS;
    design_params.design_mode =
        design_sidb_gates_params<offset::ucoord_t>::design_sidb_gates_mode::AUTOMATIC_EXHAUSTIVE_GATE_DESIGNER;
    design_params.termination_cond =
        design_sidb_gates_params<offset::ucoord_t>::termination_condition::ALL_COMBINATIONS_ENUMERATED;

    design_params.canvas                 = {{17, 13}, {21, 24}};
    design_params.number_of_canvas_sidbs = 3;

    const auto standard_cells = design_sidb_gates(skeleton, std::vector<tt>{create_and_tt()}, design_params);
    std::cout << "Number of found standard cells: " << standard_cells.size() << std::endl;
}
