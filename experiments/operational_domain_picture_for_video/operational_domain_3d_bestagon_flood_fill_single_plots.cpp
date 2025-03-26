//
// Created by marcel on 22.07.24.
//

#include "fiction_experiments.hpp"  // experiment class

#include <fiction/algorithms/simulation/sidb/operational_domain.hpp>      // operational domain computation algorithms
#include <fiction/algorithms/simulation/sidb/sidb_simulation_engine.hpp>  // SiDB simulation engines
#include <fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp>  // SiDB simulation parameters
#include <fiction/io/read_sqd_layout.hpp>                                     // reader for SiDB layouts
#include <fiction/io/write_operational_domain.hpp>                            // writer for operational domains
#include <fiction/types.hpp>                    // pre-defined types suitable for the FCN domain
#include <fiction/utils/truth_table_utils.hpp>  // truth tables helper functions

#include <fmt/format.h>                    // output formatting
#include <mockturtle/utils/stopwatch.hpp>  // stopwatch for time measurement

#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace fiction;

int main()  // NOLINT
{
    // simulation parameters
    sidb_simulation_parameters sim_params{};
    sim_params.base     = 2;
    sim_params.mu_minus = -0.32;

    // operational domain parameters
    operational_domain_params op_domain_params{};
    op_domain_params.operational_params.simulation_parameters = sim_params;
    op_domain_params.operational_params.sim_engine            = sidb_simulation_engine::QUICKEXACT;
    op_domain_params.sweep_dimensions         = {{sweep_parameter::EPSILON_R}, {sweep_parameter::LAMBDA_TF}};
    op_domain_params.sweep_dimensions[0].min  = 1.0;
    op_domain_params.sweep_dimensions[0].max  = 10.0;
    op_domain_params.sweep_dimensions[0].step = 0.1;
    op_domain_params.sweep_dimensions[1].min  = 1.0;
    op_domain_params.sweep_dimensions[1].max  = 10.0;
    op_domain_params.sweep_dimensions[1].step = 0.1;

    static const std::string folder = fmt::format("{}sidb_gate_libraries/bestagon_gates/", EXPERIMENTS_PATH);

    static const std::array<std::pair<std::string, std::vector<tt>>, 1> gates = {
        std::make_pair("and", std::vector<tt>{create_and_tt()})};

    for (const auto& [gate, truth_table] : gates)
    {
        auto lyt = read_sqd_layout<sidb_100_cell_clk_lyt_siqad>(fmt::format("{}/{}.sqd", folder, gate), gate);

        // compute the operational domains
        const auto op_domain_gs = operational_domain_flood_fill(lyt, truth_table, 100, op_domain_params);

        return EXIT_SUCCESS;
    }
}
