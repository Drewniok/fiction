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

#include <fmt/format.h>                 // output formatting
#include <mockturtle/networks/xag.hpp>  // XOR-AND-inverter graphs

#include <cstdint>
#include <string>

int main()
{
    experiments::experiment<std::string, uint64_t, float, uint64_t, std::string> defect_exp{
        "benchmark",          "gates",    "single runtime exact (in ms.)", "simulation accuracy (in %)",
        "TTS (in ms.)", "SiDB dots"};

    static constexpr const uint64_t bench_select =
        fiction_experiments::all & ~fiction_experiments::fontes18 & ~fiction_experiments::trindade16 & ~fiction_experiments::xorgate;

    auto number_bench = fiction_experiments::all_benchmarks(bench_select).size();

    uint32_t              sum_sr  = 0u;
    uint64_t              sum_tts = 0u;
    float                 sum_acc = 0u;
    std::vector<uint64_t> db_num{};

    for (const auto& benchmark : fiction_experiments::all_benchmarks(bench_select))
    {
        const auto lyt = fiction::read_sqd_layout<fiction::sidb_cell_clk_lyt_siq>(benchmark);
        fmt::print("[i] processing {}\n", benchmark);
        mockturtle::xag_network xag{};

        const fiction::simulation_params     params{2};
        fiction::charge_distribution_surface charge_layout{lyt, params};

        //auto [runtime, exactlyt] = fiction::detail::metastable_layouts(charge_layout);

        auto test = fiction::detail::Sim(charge_layout);

        //auto [acc, tts] = fiction::sim_acc_tts<fiction::sidb_cell_clk_lyt_siq>(charge_layout, exactlyt, 1, 80);
        //auto [acc, tts] = fiction::sim_acc_tts<fiction::sidb_cell_clk_lyt_siq>(charge_layout, exactlyt, 1, 80);

//        defect_exp(benchmark, runtime, acc, tts, std::to_string(lyt.num_cells()));
//        db_num.push_back(lyt.num_cells());
//        sum_sr += runtime;
//        sum_acc += acc;
//        sum_tts += tts;
    }

    auto min_db_num = std::min_element(db_num.begin(), db_num.end());
    auto max_db_num = std::max_element(db_num.begin(), db_num.end());
    auto mean_acc   = sum_acc / static_cast<float>(number_bench);

    defect_exp("sum", sum_sr, mean_acc, sum_tts, std::to_string(*min_db_num) + " -- " + std::to_string(*max_db_num));
    defect_exp.save();
    defect_exp.table();
    return EXIT_SUCCESS;
}
