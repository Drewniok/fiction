//
// Created by Jan Drewniok on 25.08.25.
//

#include "fiction/algorithms/iter/bdl_input_iterator.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/algorithms/simulation/sidb/operational_domain.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_engine.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/io/read_sqd_layout.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/truth_table_utils.hpp"

#include <vector>

using namespace fiction;

int main()  // NOLINT
{
    const auto                       cell = read_sqd_layout<sidb_100_cell_clk_lyt_siqad>("and.sqd");
    const sidb_simulation_parameters sim_params{2, -0.32};
    operational_domain_params        op_params{
        is_operational_params{sim_params, sidb_simulation_engine::QUICKEXACT, bdl_input_iterator_params{},
                              is_operational_params::operational_condition::TOLERATE_KINKS}};
    op_params.sweep_dimensions = {{sweep_parameter::EPSILON_R, 1.0, 10.0, 0.05},
                                  {sweep_parameter::LAMBDA_TF, 1.0, 10.0, 0.05}};
    operational_domain_stats op_stats{};

    const auto op_gs = operational_domain_grid_search(cell, std::vector<tt>{create_and_tt()}, op_params, &op_stats);
    const auto op_rs =
        operational_domain_random_sampling(cell, std::vector<tt>{create_and_tt()}, 2500, op_params, &op_stats);
    const auto op_ff = operational_domain_flood_fill(cell, std::vector<tt>{create_and_tt()}, 250, op_params, &op_stats);
    const auto op_ct =
        operational_domain_contour_tracing(cell, std::vector<tt>{create_and_tt()}, 100, op_params, &op_stats);
}
