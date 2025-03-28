//
// Created by Jan Drewniok on 30.11.24.
//

#include "fiction/algorithms/simulation/sidb/defect_clearance.hpp"
#include "fiction/algorithms/simulation/sidb/defect_influence.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/io/read_sqd_layout.hpp"
#include "fiction/io/write_defect_influence_domain.hpp"
#include "fiction/io/write_sqd_layout.hpp"
#include "fiction/technology/sidb_defects.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/truth_table_utils.hpp"
#include "fiction_experiments.hpp"

#include <fmt/format.h>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

using namespace fiction;

int main()  // NOLINT
{

    static const std::string bestagon_folder = fmt::format("{}sidb_gate_libraries/bestagon_gates/", EXPERIMENTS_PATH);
    static const std::string plot_folder     = fmt::format("{}quicktrace/plots/", EXPERIMENTS_PATH);

    static const std::array<std::pair<std::string, std::vector<tt>>, 1> gates = {
        std::make_pair("xnor", std::vector<tt>{create_xnor_tt()})};

    is_operational_params is_op_params{sidb_simulation_parameters{2, -0.32, 5.6, 5.0}};
    is_op_params.sim_engine = sidb_simulation_engine::QUICKEXACT;

    // for this experiment we use a stray SiDB defect
    const auto stray_db = fiction::sidb_defect{fiction::sidb_defect_type::DB, -1, 4.1, 1.8};
    // const auto si_vacancy = fiction::sidb_defect{fiction::sidb_defect_type::SI_VACANCY, -1, 10.6, 5.9};

    defect_influence_params<fiction::cell<sidb_100_cell_clk_lyt_cube>> params{};
    params.additional_scanning_area = {50, 50};
    params.defect                   = stray_db;
    params.operational_params       = is_op_params;

    for (const auto& [gate, truth_table] : gates)
    {
        // Create a folder where the defect influence plots are stored
        std::filesystem::create_directories(plot_folder);

        // Create a folder where the defect influence CSV files of each gate are stored
        const auto plot_folder_for_given_gate = fmt::format("{}{}/", plot_folder, gate);
        std::filesystem::create_directories(plot_folder_for_given_gate);

        // read the Bestagon SiDB layout
        const auto layout = read_sqd_layout<sidb_100_cell_clk_lyt_cube>(fmt::format("{}{}.sqd", bestagon_folder, gate));

        // Write the SQD layout
        write_sqd_layout(layout, fmt::format("{}/{}.sqd", plot_folder_for_given_gate, gate));

        // quicktrace
        defect_influence_stats quicktrace_stats{};
        const auto             defect_inf_quicktrace =
            defect_influence_quicktrace(layout, truth_table, 8, params, &quicktrace_stats);

        write_defect_influence_domain<sidb_100_cell_clk_lyt_cube>(
            defect_inf_quicktrace, fmt::format("{}{}_quicktrace.csv", plot_folder_for_given_gate, gate));
    }

    return EXIT_SUCCESS;
}
