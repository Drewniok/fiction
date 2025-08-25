//
// Created by Jan Drewniok on 25.08.25.
//

#include "fiction/algorithms/iter/bdl_input_iterator.hpp"
#include "fiction/algorithms/physical_design/design_sidb_gates.hpp"
#include "fiction/algorithms/simulation/sidb/critical_temperature.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/algorithms/simulation/sidb/operational_domain.hpp"
#include "fiction/algorithms/simulation/sidb/quickexact.hpp"
#include "fiction/algorithms/simulation/sidb/quicksim.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_engine.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/io/read_sqd_layout.hpp"
#include "fiction/layouts/coordinates.hpp"
#include "fiction/types.hpp"

#include <vector>

using namespace fiction;

void run_quickexact_script()  // NOLINT
{
    const auto                              cell = read_sqd_layout<sidb_100_cell_clk_lyt_siqad>("and.sqd");
    const sidb_simulation_parameters        sim_params{2, -0.32};
    const quickexact_params<siqad::coord_t> qe_params{sim_params};
    const auto                              results = quickexact(cell, qe_params);
    const auto                              gs      = results.groundstates();
};

void run_quicksim_script()  // NOLINT
{
    const auto                       layout = read_sqd_layout<sidb_100_cell_clk_lyt_siqad>("and.sqd");
    const sidb_simulation_parameters sim_params{2, -0.32, 5.6, 5.0};
    const quicksim_params            qs_params{sim_params, 300, 0.6};
    const auto                       results = quicksim(layout, qs_params);
    const auto                       gs      = results.value().groundstates();
};

void run_exhaustive_gate_design_script()
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
};

void run_quickcell_script()
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
};

void run_critical_temperature_script()
{
    const auto                        cell = read_sqd_layout<sidb_100_cell_clk_lyt_siqad>("and.sqd");
    const sidb_simulation_parameters  sim_params{2, -0.32};
    const critical_temperature_params ct_params{{sim_params}};
    critical_temperature_stats        ct_stats{};
    const auto                        ct =
        critical_temperature_gate_based(cell, std::vector<tt>{create_and_tt()}, ct_params, &ct_stats);  // ct ≈ 19.8 K
};

void run_operational_domain_script()
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
};

int main()  // NOLINT
{
    run_quickexact_script();
    run_quicksim_script();
    run_exhaustive_gate_design_script();
    run_quickcell_script();
}
