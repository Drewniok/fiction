//
// Created by marcel on 31.08.22.
//


#include "../fiction_experiments.hpp"

#include <fiction/algorithms/physical_design/apply_gate_library.hpp>  // layout conversion to cell-level
#include <fiction/algorithms/properties/critical_path_length_and_throughput.hpp>  // critical path and throughput calculations
#include <fiction/io/read_sidb_surface_defects.hpp>                               // reader for simulated SiDB surfaces
#include <fiction/io/read_sqd_layout.hpp>                     // reader for SiDB layouts including surface scan data
#include <fiction/io/write_sqd_layout.hpp>                    // writer for SiQAD files (physical simulation)
#include <fiction/technology/area.hpp>                        // area requirement calculations
#include <fiction/technology/cell_technologies.hpp>           // cell implementations
#include <fiction/technology/sidb_bestagon_library.hpp>       // a pre-defined SiDB gate library
#include <fiction/technology/sidb_defects.hpp>                // SiDB defect classes
#include <fiction/technology/sidb_surface.hpp>                // H-Si(100) 2x1 surface model
#include <fiction/technology/sidb_surface_analysis.hpp>       // SiDB surface analysis
#include <fiction/technology/technology_mapping_library.hpp>  // pre-defined gate types for technology mapping
#include <fiction/types.hpp>                                  // pre-defined types suitable for the FCN domain
#include <fiction/technology/charge_distribution_surface.hpp>
#include <fiction/algorithms/simulation_sidb/ExGS.hpp>
#include <fiction/algorithms/simulation_sidb/new_approach.hpp>
#include <fiction/algorithms/simulation_sidb/TTS.hpp>

#include <fmt/format.h>                                        // output formatting
#include <lorina/genlib.hpp>                                   // Genlib file parsing
#include <lorina/lorina.hpp>                                   // Verilog/BLIF/AIGER/... file parsing
#include <mockturtle/algorithms/cut_rewriting.hpp>             // logic optimization with cut rewriting
#include <mockturtle/algorithms/equivalence_checking.hpp>      // equivalence checking
#include <mockturtle/algorithms/mapper.hpp>                    // Technology mapping on the logic level
#include <mockturtle/algorithms/miter.hpp>                     // miter structure
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>  // NPN databases for cut rewriting of XAGs and AIGs
#include <mockturtle/io/genlib_reader.hpp>                     // call-backs to read Genlib files into gate libraries
#include <mockturtle/io/verilog_reader.hpp>                    // call-backs to read Verilog files into networks
#include <mockturtle/networks/klut.hpp>                        // kLUT network
#include <mockturtle/networks/xag.hpp>                         // XOR-AND-inverter graphs
#include <mockturtle/utils/tech_library.hpp>                   // technology library utils
#include <mockturtle/views/depth_view.hpp>                     // to determine network levels

#include <cstdint>
#include <set>
#include <sstream>
#include <string>
#include <vector>

int main()
{

    experiments::experiment<std::string, uint64_t, float,uint64_t, uint64_t>
        defect_exp{"benchmark", "gates", "single runtime exact (in millisec.)", "simulation accuracy (in %)", "TTS (in millisec.)", "SiDB dots"};


    static constexpr const uint64_t bench_select = fiction_experiments::all & ~fiction_experiments::fontes18 &
                                                   ~fiction_experiments::trindade16 & ~fiction_experiments::xor_00
                                                   & ~fiction_experiments::xor_01
                                                   & ~fiction_experiments::xor_10
                                                   & ~fiction_experiments::xor_11;

    std::cout << fiction_experiments::all_benchmarks(bench_select).size() << std::endl;

    for (const auto& benchmark : fiction_experiments::all_benchmarks(bench_select))
    {
        const auto lyt = fiction::read_sqd_layout<fiction::sidb_cell_clk_lyt_siq>(benchmark);
        fmt::print("[i] processing {}\n", benchmark);
        mockturtle::xag_network xag{};

        const fiction::simulation_params params{2};
        fiction::charge_distribution_surface charge_layout{lyt, params};

        auto [runtime, exactlyt] =
            fiction::detail::metastable_layouts(charge_layout);

        auto [acc, tts] = fiction::sim_acc_tts<fiction::sidb_cell_clk_lyt_siq>(charge_layout, exactlyt, 100, 80);

            defect_exp(benchmark, runtime, acc, tts, lyt.num_cells());

    }
    defect_exp.save();
    defect_exp.table();
    return 0;
}
