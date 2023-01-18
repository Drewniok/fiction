//
// Created by Jan Drewniok on 18.01.23.
//

#include "../fiction_experiments.hpp"

#include <fiction/algorithms/simulation_sidb/ExGS.hpp>
#include <fiction/algorithms/simulation_sidb/TTS.hpp>
#include <fiction/io/read_sqd_layout.hpp>  // reader for SiDB layouts including surface scan data
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/technology/sidb_defects.hpp>  // SiDB defect classes
#include <fiction/types.hpp>                    // pre-defined types suitable for the FCN domain

#include <fmt/format.h>  // output formatting

using namespace fiction;

int main()
{

    std::vector<std::string> folders = {
        "../../experiments/bestagon/layouts/gates/and/",
    };

    for (const auto& folder : folders)
    {
        for (const auto& file : std::filesystem::directory_iterator(folder))
        {
            if (std::filesystem::is_regular_file(file))
            {

                const auto& benchmark = file.path();

                const auto            lyt = read_sqd_layout<sidb_cell_clk_lyt_siq>(benchmark);
                const physical_params params{2, -0.32};
                charge_distribution_surface<sidb_cell_clk_lyt_siq> chargelyt{lyt};

                exgs_stats<sidb_cell_clk_lyt_siq> stats{};

                exgs<sidb_cell_clk_lyt_siq>(chargelyt, stats, params);
                stats.report();
            }
        }
    }
    return 0;
}