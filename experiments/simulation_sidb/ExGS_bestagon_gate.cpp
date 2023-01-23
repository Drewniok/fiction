//
// Created by Jan Drewniok on 18.01.23.
//

#include <fiction/algorithms/simulation_sidb/ExGS.hpp>
#include <fiction/algorithms/simulation_sidb/TTS.hpp>
#include <fiction/io/read_sqd_layout.hpp>  // reader for SiDB layouts including surface scan data
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/types.hpp>                    // pre-defined types suitable for the FCN domain
#include <fmt/format.h>  // output formatting
#include <filesystem>
#include <iostream>

using namespace fiction;

int main() // NOLINT
{

    std::vector<std::string> folders = {
        "../../experiments/bestagon/layouts/gates/and/",
    };

    for (const auto& folder : folders)
    {
        for (const auto& file : std::filesystem::directory_iterator(folder))
        {
            std::string benchmark = file.path();

            const auto                                         lyt = read_sqd_layout<sidb_cell_clk_lyt_siq>(benchmark);
            const physical_params                              params{2, -0.32};
            charge_distribution_surface<sidb_cell_clk_lyt_siq> chargelyt{lyt};

            exgs_stats<sidb_cell_clk_lyt_siq> stats{};

            exgs<sidb_cell_clk_lyt_siq>(chargelyt, stats, params);
            stats.report();
        }
    }
    return 0;
}
