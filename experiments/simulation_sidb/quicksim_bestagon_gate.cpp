//
// Created by Jan Drewniok on 01.01.23.
//

#include <fiction/algorithms/simulation_sidb/quicksim.hpp>
#include <fiction/io/read_sqd_layout.hpp>  // reader for SiDB layouts including surface scan data
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/types.hpp>                    // pre-defined types suitable for the FCN domain
#include <fmt/format.h>  // output formatting
#include <filesystem>

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
            auto benchmark = file.path();

            const auto                     lyt = read_sqd_layout<sidb_cell_clk_lyt_siq>(benchmark);
            const fiction::physical_params params{2, -0.32};
            charge_distribution_surface<sidb_cell_clk_lyt_siq> chargelyt{lyt};

            quicksim_stats<sidb_cell_clk_lyt_siq> collect{};

            quicksim<sidb_cell_clk_lyt_siq>(chargelyt, collect, params);
            collect.report();
        }
    }
    return 0;
}
