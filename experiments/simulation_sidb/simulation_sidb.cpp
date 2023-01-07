//
// Created by marcel on 31.08.22.
//


#include "../fiction_experiments.hpp"

#include <fiction/io/read_sqd_layout.hpp>                     // reader for SiDB layouts including surface scan data
#include <fiction/technology/sidb_defects.hpp>                // SiDB defect classes
#include <fiction/types.hpp>                                  // pre-defined types suitable for the FCN domain
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/algorithms/simulation_sidb/ExGS.hpp>
#include <fiction/algorithms/simulation_sidb/TTS.hpp>

#include <fmt/format.h>                                        // output formatting
#include <mockturtle/algorithms/miter.hpp>                     // miter structure
#include <mockturtle/networks/xag.hpp>                         // XOR-AND-inverter graphs

#include <cstdint>
#include <sstream>


int main()
{
    experiments::experiment<std::string, uint64_t, float,uint64_t, uint64_t>
        defect_exp{"benchmark", "gates", "single runtime exact (in millisec.)", "simulation accuracy (in %)", "TTS (in millisec.)", "SiDB dots"};


//    static constexpr const uint64_t bench_select = fiction_experiments::all & ~fiction_experiments::fontes18 &
//                                                   ~fiction_experiments::trindade16 & ~fiction_experiments::xor_00
//                                                   & ~fiction_experiments::xor_01
//                                                   & ~fiction_experiments::xor_11 & ~fiction_experiments::xor_wo;

    static constexpr const uint64_t bench_select = fiction_experiments::all & ~fiction_experiments::fontes18 &
                                                   ~fiction_experiments::trindade16;

    std::cout << fiction_experiments::all_benchmarks(bench_select).size() << std::endl;

    for (const auto& benchmark : fiction_experiments::all_benchmarks(bench_select))
    {
        const auto lyt = fiction::read_sqd_layout<fiction::sidb_cell_clk_lyt_siq>(benchmark);
        fmt::print("[i] processing {}\n", benchmark);
        mockturtle::xag_network xag{};

        const fiction::simulation_params params{2};
        fiction::charge_distribution_surface charge_layout{lyt, params};

//
        auto [runtime, exactlyt] =
            fiction::detail::metastable_layouts(charge_layout);

//
//            fiction::detail::metastable_layouts(charge_layout);

        auto [acc, tts] = fiction::sim_acc_tts<fiction::sidb_cell_clk_lyt_siq>(charge_layout, exactlyt, 100, 80);

            defect_exp(benchmark, runtime, 0, 0, lyt.num_cells());

    }
    defect_exp.save();
    defect_exp.table();
    return EXIT_SUCCESS;
}
