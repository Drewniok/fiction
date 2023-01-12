//
// Created by marcel on 31.08.22.
//

#include "../fiction_experiments.hpp"

#include <fiction/algorithms/simulation_sidb/ExGS.hpp>
#include <fiction/algorithms/simulation_sidb/new_approach.hpp>
#include <fiction/algorithms/simulation_sidb/TTS.hpp>
#include <fiction/io/read_sqd_layout.hpp>  // reader for SiDB layouts including surface scan data
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/technology/sidb_defects.hpp>  // SiDB defect classes
#include <fiction/types.hpp>                    // pre-defined types suitable for the FCN domain
#include "fiction/algorithms/simulation_sidb/temperature.hpp"

#include <fmt/format.h>                 // output formatting

#include <cstdint>
#include <string>

int main()
{
    experiments::experiment<std::string, double, double, double, std::string, double> simulation_exp{
        "benchmark",    "gates",     "single runtime exact (in ms.)", "simulation accuracy (in %)",
        "TTS (in ms.)", "SiDB dots", "critical temperature [K]"};

    double                sum_sr  = 0u;
    double                sum_tts = 0u;
    double                sum_acc = 0u;
    std::vector<uint64_t> db_num{};
    uint64_t              benchmark_counter = 0u;
    uint64_t              sum_ct = 0u;

    std::vector<std::string> folders = {
        "../../experiments/bestagon/layouts/defect_robust/",
    };

    for (const auto& folder : folders)
    {
        for (const auto& file : std::filesystem::directory_iterator(folder))
        {
            if (std::filesystem::is_regular_file(file))
            {
                benchmark_counter += 1;
                const auto& benchmark = file.path();
                std::cout << benchmark << std::endl;

                const auto lyt = fiction::read_sqd_layout<fiction::sidb_cell_clk_lyt_siq>(benchmark);

                const fiction::physical_params params{2};

                fiction::charge_distribution_surface charge_layout{lyt, params};

                auto [runtime, exactlyt] = fiction::detail::exgs(charge_layout);
                //                auto runtime = 1;
                //                auto tts = 1;
                //                auto acc = 1;

                auto ct = fiction::critical_temp<fiction::sidb_cell_clk_lyt_siq>(exactlyt, 0.99);

                // auto result = fiction::detail::faccusim(charge_layout);

                // auto ct =
                // fiction::critical_temp<fiction::sidb_cell_clk_lyt_siq>(fiction::detail::faccusim(charge_layout, 1000,
                // 0.5));

                auto [acc, tts] =
                    fiction::sim_acc_tts<fiction::sidb_cell_clk_lyt_siq>(charge_layout, exactlyt, 100, 80);

                simulation_exp(benchmark, runtime, acc, tts, std::to_string(lyt.num_cells()), static_cast<double>(ct));
                db_num.push_back(lyt.num_cells());
                sum_sr += runtime;
                sum_acc += acc;
                sum_tts += tts;
                sum_ct += ct;
            }
        }
    }

    auto min_db_num = std::min_element(db_num.begin(), db_num.end());
    auto max_db_num = std::max_element(db_num.begin(), db_num.end());
    auto mean_acc   = sum_acc / static_cast<double>(benchmark_counter);
    auto mean_ct    = static_cast<double>(sum_ct) / static_cast<double>(benchmark_counter);

    simulation_exp("sum", sum_sr, mean_acc, sum_tts,
                   std::to_string(*min_db_num) + " -- " + std::to_string(*max_db_num), mean_ct);
    simulation_exp.save();
    simulation_exp.table();
    return EXIT_SUCCESS;
}
