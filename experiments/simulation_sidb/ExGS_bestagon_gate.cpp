//
// Created by Jan Drewniok on 18.01.23.
//

#include <fiction/algorithms/simulation_sidb/TTS.hpp>
#include <fiction/algorithms/simulation_sidb/exhaustive_ground_state_simulation.hpp>
#include <fiction/io/read_sqd_layout.hpp>  // reader for SiDB layouts including surface scan data
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/types.hpp>  // pre-defined types suitable for the FCN domain

#include <array>
#include <filesystem>

using namespace fiction;

int main() // NOLINT
{

    const std::string folder = fmt::format("{}/bestagon_gates/and/", EXPERIMENTS_PATH);

        for (const auto& file : std::filesystem::directory_iterator(folder))
        {
            const auto& benchmark = file.path();

            const auto                                         lyt = read_sqd_layout<sidb_cell_clk_lyt_siq>(benchmark.string());
            const physical_params                              params{2, -0.32};
            charge_distribution_surface<sidb_cell_clk_lyt_siq> chargelyt{lyt};

            exgs_stats<sidb_cell_clk_lyt_siq> stats{};

            exgs<sidb_cell_clk_lyt_siq>(chargelyt, stats, params);
            stats.report();
        }

    return EXIT_SUCCESS;
}
