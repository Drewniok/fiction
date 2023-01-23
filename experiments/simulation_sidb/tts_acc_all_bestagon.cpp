//
// Created by Jan Drewniok 01.01.23
//

#include <fiction/algorithms/simulation_sidb/ExGS.hpp>
#include <fiction/algorithms/simulation_sidb/TTS.hpp>
#include <fiction/io/read_sqd_layout.hpp>  // reader for SiDB layouts including surface scan data
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/types.hpp>                    // pre-defined types suitable for the FCN domain
#include <fmt/format.h>  // output formatting
#include <cstdint>
#include <string>
#include <filesystem>
#include <../experiments/fiction_experiments.hpp>

using namespace fiction;

int main() // NOLINT
{
    experiments::experiment<std::string, double, double, double, double, std::string> simulation_exp{
        "benchmark",
        "gates",
        "single runtime exact (in sec.)",
        "simulation accuracy (in %)",
        "TTS (in sec.)",
        "single runtime of quicksim (in sec.)",
        "SiDB dots"};

    double                sum_sr       = 0u;
    double                sum_tts      = 0u;
    double                sum_acc      = 0u;
    double                sum_sr_quick = 0u;
    std::vector<uint64_t> db_num{};
    uint64_t              benchmark_counter = 0u;

    std::vector<std::string> folders = {
        "../../experiments/bestagon/layouts/gates/and/",       //"../../experiments/bestagon/layouts/gates/cx/",
//        "../../experiments/bestagon/layouts/gates/fo2/",       "../../experiments/bestagon/layouts/gates/ha/",
//        "../../experiments/bestagon/layouts/gates/hourglass/", "../../experiments/bestagon/layouts/gates/inv/",
//        "../../experiments/bestagon/layouts/gates/nand/",      "../../experiments/bestagon/layouts/gates/nor/",
//        "../../experiments/bestagon/layouts/gates/or/",        "../../experiments/bestagon/layouts/gates/wire/",
//        "../../experiments/bestagon/layouts/gates/xnor/",      "../../experiments/bestagon/layouts/gates/xor/",
    };

    for (const auto& folder : folders)
    {
        for (const auto& file : std::filesystem::directory_iterator(folder))
        {
            benchmark_counter += 1;
            const auto& benchmark = file.path();
            std::cout << benchmark << std::endl;

            const auto lyt = read_sqd_layout<sidb_cell_clk_lyt_siq>(benchmark, "bestagon_tts_acc_result");

            const physical_params                              params{2, -0.32};
            charge_distribution_surface<sidb_cell_clk_lyt_siq> chargelyt{lyt};
            exgs_stats<sidb_cell_clk_lyt_siq>                  exgs_stats{};
            exgs<sidb_cell_clk_lyt_siq>(chargelyt, exgs_stats, params);

            tts_stats tts_stat{};
            sim_acc_tts<sidb_cell_clk_lyt_siq>(chargelyt, tts_stat, exgs_stats);

            simulation_exp(benchmark, mockturtle::to_seconds(exgs_stats.time_total), tts_stat.acc, tts_stat.tts,
                           tts_stat.mean_single_runtime, std::to_string(lyt.num_cells()));
            db_num.push_back(lyt.num_cells());
            sum_sr += mockturtle::to_seconds(exgs_stats.time_total);
            sum_sr_quick += tts_stat.mean_single_runtime;
            sum_acc += tts_stat.acc;
            sum_tts += tts_stat.tts;
        }
    }

    auto min_db_num    = std::min_element(db_num.begin(), db_num.end());
    auto max_db_num    = std::max_element(db_num.begin(), db_num.end());
    auto mean_acc      = sum_acc / static_cast<double>(benchmark_counter);
    auto mean_sr_quick = sum_sr_quick / static_cast<double>(benchmark_counter);

    simulation_exp("sum", sum_sr, mean_acc, sum_tts, mean_sr_quick,
                   std::to_string(*min_db_num) + " -- " + std::to_string(*max_db_num));
    simulation_exp.save();
    simulation_exp.table();
    return EXIT_SUCCESS;
}
