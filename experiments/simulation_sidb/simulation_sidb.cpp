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

    using gate_lyt = fiction::hex_even_row_gate_clk_lyt;
    using cell_lyt = fiction::sidb_cell_clk_lyt;

    //static const std::string layouts_folder = fmt::format("{}bestagon/layouts/gates/hex_21_inputsdbp_xor_v1", EXPERIMENTS_PATH);
    //std::cout << layouts_folder;
    //std::cout << EXPERIMENTS_PATH;
    // 740 x 545 dimers = 740 x 1090 DB positions = 12 x 31 Bestagon tiles
    // static const std::string surface_data_path =
    // fmt::format("{}/defect_aware_physical_design/full_scan_area/defects_full70.xml", EXPERIMENTS_PATH);
    // 830 x 326 dimers = 830 x 652 DB positions = 13 x 18 Bestagon tiles
    //    static const std::string surface_data_path =
    //        fmt::format("{}/defect_aware_physical_design/full_scan_area/defects_full56_Oct.xml", EXPERIMENTS_PATH);

    experiments::experiment<std::string, uint32_t, uint32_t>
        defect_exp{"benchmark", "gates", "simulation accuracy (in %)", "TTS (in millisec.)"};

    // instantiate a technology mapping library
    //    std::stringstream library_stream{};
    //    library_stream << fiction::GATE_ZERO << fiction::GATE_ONE << fiction::GATE_BUF << fiction::GATE_INV
    //                   << fiction::GATE_AND2 << fiction::GATE_NAND2 << fiction::GATE_OR2 << fiction::GATE_NOR2
    //                   << fiction::GATE_XOR2 << fiction::GATE_XNOR2;
    //
    //    std::vector<mockturtle::gate> gates{};
    //
    //    // parameters for technology mapping
    //    mockturtle::map_params map_params{};
    //
    //    [[maybe_unused]] const auto read_genlib_result =
    //        lorina::read_genlib(library_stream, mockturtle::genlib_reader{gates});
    //    assert(read_genlib_result == lorina::return_code::success);
    //    mockturtle::tech_library<2> gate_lib{gates};
    //
    //    // parameterize the H-Si(100) 2x1 surface to ignore certain defect types
    //    fiction::sidb_surface_params surface_params{std::set<fiction::sidb_defect_type>{fiction::sidb_defect_type::DB}};
    //
    //    //    fiction::sidb_surface<cell_lyt> surface_lattice{surface_params};
    //
    //    // read surface scan lattice data
    //    const auto surface_lattice = fiction::read_sidb_surface_defects<cell_lyt>(
    //        "../../experiments/defect_aware_physical_design/py_test_surface.txt", "py_test_surface");
    //    //    fiction::read_sqd_layout(surface_lattice, surface_data_path);
    //
    //    const auto lattice_tiling = gate_lyt{{11, 30}};  // our surface data is 12 x 31 Bestagon tiles
    //    //    const auto lattice_tiling = gate_lyt{{12, 17}};  // our surface data is 13 x 18 Bestagon tiles
    //    const auto black_list =
    //        fiction::sidb_surface_analysis<fiction::sidb_bestagon_library>(lattice_tiling, surface_lattice);
    //
    //    // parameters for SMT-based physical design
    //    fiction::exact_physical_design_params<gate_lyt> exact_params{};
    //    exact_params.scheme        = fiction::ptr<gate_lyt>(fiction::row_clocking<gate_lyt>(fiction::num_clks::FOUR));
    //    exact_params.crossings     = true;
    //    exact_params.border_io     = false;
    //    exact_params.desynchronize = false;
    //    exact_params.upper_bound_x = 11;      // 12 x 31 tiles
    //    exact_params.upper_bound_y = 30;      // 12 x 31 tiles
    //                                          //    exact_params.upper_bound_x = 12;         // 13 x 18 tiles
    //                                          //    exact_params.upper_bound_y = 17;         // 13 x 18 tiles
    //    exact_params.timeout    = 3'600'000;  // 1h in ms
    //    exact_params.black_list = black_list;
    //    fiction::exact_physical_design_stats exact_stats{};
    //
    constexpr const uint64_t bench_select = fiction_experiments::cal;
    std::cout << fiction_experiments::all_benchmarks(bench_select).size() << std::endl;

    for (const auto& benchmark : fiction_experiments::all_benchmarks(bench_select))
    {
        std::cout << benchmark;
        const auto lyt = fiction::read_sqd_layout<fiction::sidb_cell_clk_lyt_siq>(benchmark);
        fmt::print("[i] processing {}\n", benchmark);
        mockturtle::xag_network xag{};

        const fiction::simulation_params params{5.6, 5.0 * 1E-9, -0.32, 3.84 * 1E-10, 7.68 * 1E-10, 2.25 * 1E-10, 2};
        fiction::charge_distribution_surface charge_layout{lyt, params};
        std::unordered_map<double, fiction::charge_distribution_surface<fiction::sidb_cell_clk_lyt_siq>> output_exact =
            fiction::detail::metastable_layouts(charge_layout);

        std::unordered_map<double, fiction::charge_distribution_surface<fiction::sidb_cell_clk_lyt_siq>> output_AP =
            fiction::detail::Sim<fiction::sidb_cell_clk_lyt_siq>(charge_layout, 20, 0.8);

        auto acc = fiction::sim_acc_tts<fiction::sidb_cell_clk_lyt_siq>(charge_layout, output_exact, 100);
        std::cout << acc.first;
        //        [[maybe_unused]] const auto read_verilog_result =
        //            lorina::read_verilog(fiction_experiments::benchmark_path(benchmark), mockturtle::verilog_reader(xag));
        //        assert(read_verilog_result == lorina::return_code::success);

        //        // compute depth
        //        mockturtle::depth_view depth_xag{xag};
        //
        //        // rewrite network cuts using the given re-synthesis function
        //        const auto cut_xag = mockturtle::cut_rewriting(xag, resynthesis_function, cut_params);
        //        // compute depth
        //        mockturtle::depth_view depth_cut_xag{cut_xag};
        //
        //        // perform technology mapping
        //        const auto mapped_network = mockturtle::map(cut_xag, gate_lib, map_params);
        //        // compute depth
        //        mockturtle::depth_view depth_mapped_network{mapped_network};
        //
        //        // perform layout generation with an SMT-based exact algorithm
        //        const auto gate_level_layout = fiction::exact<gate_lyt>(mapped_network, exact_params, &exact_stats);

        //       if (output.size()>0)
        //      if (gate_level_layout.has_value())

            // check equivalence
            //            const auto miter = mockturtle::miter<mockturtle::klut_network>(mapped_network, *gate_level_layout); const auto eq    = mockturtle::equivalence_checking(*miter);
            //            assert(eq.has_value());
            //
            //            // compute critical path and throughput
            //            fiction::critical_path_length_and_throughput_stats cp_tp_stats{};
            //            fiction::critical_path_length_and_throughput(*gate_level_layout, &cp_tp_stats);

            // apply gate library
            //            const auto cell_level_layout =
            //                fiction::apply_gate_library<cell_lyt, fiction::sidb_bestagon_library>(*gate_level_layout);

            // compute area
            //            fiction::area_stats                            area_stats{};
            //            fiction::area_params<fiction::sidb_technology> area_ps{};
            //            fiction::area(cell_level_layout, area_ps, &area_stats);
            //
            //            // write a SiQAD simulation file
            //            fiction::write_sqd_layout(cell_level_layout, fmt::format("{}/{}.sqd", layouts_folder, benchmark));
            //            // TODO copy defects to cell layout before writing

            // log results
            defect_exp(benchmark, acc.first, acc.second);
            //        }
            //        else  // no layout was obtained
            //        {
            //            // log results
            //            defect_exp(benchmark, xag.num_pis(), xag.num_pos(), xag.num_gates(), depth_xag.depth(), cut_xag.num_gates(),
            //                       depth_cut_xag.depth(), mapped_network.num_gates(), depth_mapped_network.depth(), 0, 0, 0, 0, 0, 0, 0, mockturtle::to_seconds(exact_stats.time_total), false, 0, 0);
            //        }

            defect_exp.save();
            defect_exp.table();
            //}



    }
    return 0;
}
