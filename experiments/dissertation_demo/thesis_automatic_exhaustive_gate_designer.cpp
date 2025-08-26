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
    const auto                       skeleton = read_sqd_layout<sidb_100_cell_clk_lyt_siqad>("skeleton_2i1o.sqd");
    const sidb_simulation_parameters sim_params{2, -0.32};
    design_sidb_gates_params<siqad::coord_t> params{
        is_operational_params{sim_params, sidb_simulation_engine::QUICKEXACT}};
    params.number_of_canvas_sidbs = 4;
    params.canvas                 = {{15, 8, 0}, {23, 14, 0}};
    params.termination_cond =
        design_sidb_gates_params<siqad::coord_t>::termination_condition::ALL_COMBINATIONS_ENUMERATED;
    params.design_mode =
        design_sidb_gates_params<siqad::coord_t>::design_sidb_gates_mode::AUTOMATIC_EXHAUSTIVE_GATE_DESIGNER;
    const auto standard_cells = design_sidb_gates(skeleton, std::vector<tt>{create_and_tt()}, params);
}
