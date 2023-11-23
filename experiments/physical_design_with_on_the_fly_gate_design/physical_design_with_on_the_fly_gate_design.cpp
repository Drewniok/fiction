//
// Created by marcel on 16.11.21.
//

#if (FICTION_Z3_SOLVER)

#include "fiction_experiments.hpp"

#include <fiction/algorithms/iter/aspect_ratio_iterator.hpp>
#include <fiction/algorithms/network_transformation/technology_mapping.hpp>  // technology mapping
#include <fiction/algorithms/physical_design/apply_gate_library.hpp>         // layout conversion to cell-level
#include <fiction/algorithms/physical_design/exact.hpp>                      // SMT-based physical design of FCN layouts
#include <fiction/algorithms/physical_design/orthogonal.hpp>
#include <fiction/algorithms/properties/critical_path_length_and_throughput.hpp>  // critical path and throughput calculations
#include <fiction/io/read_sidb_surface_defects.hpp>                               // reader for simulated SiDB surfaces
#include <fiction/io/write_sqd_layout.hpp>  // writer for SiQAD files (physical simulation)
#include <fiction/layouts/coordinates.hpp>
#include <fiction/networks/technology_network.hpp>   // technology-mapped network type
#include <fiction/technology/area.hpp>               // area requirement calculations
#include <fiction/technology/cell_technologies.hpp>  // cell implementations
#include <fiction/technology/sidb_defects.hpp>
#include <fiction/technology/sidb_on_the_fly_gate_library.hpp>  // a dynamic SiDB gate library
#include <fiction/technology/sidb_skeleton_bestagon_library.hpp>
#include <fiction/technology/sidb_surface.hpp>  // SiDB surface with support for atomic defects
#include <fiction/traits.hpp>
#include <fiction/types.hpp>  // pre-defined types suitable for the FCN domain
#include <fiction/utils/truth_table_utils.hpp>

#include <fmt/format.h>                                        // output formatting
#include <lorina/lorina.hpp>                                   // Verilog/BLIF/AIGER/... file parsing
#include <mockturtle/algorithms/cut_rewriting.hpp>             // logic optimization with cut rewriting
#include <mockturtle/algorithms/equivalence_checking.hpp>      // equivalence checking
#include <mockturtle/algorithms/miter.hpp>                     // miter structure
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>  // NPN databases for cut rewriting of XAGs and AIGs
#include <mockturtle/io/verilog_reader.hpp>                    // call-backs to read Verilog files into networks
#include <mockturtle/networks/xag.hpp>                         // XOR-AND-inverter graphs
#include <mockturtle/views/depth_view.hpp>                     // to determine network levels

#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <unordered_set>
#include <utility>

// This script conducts defect-aware placement and routing with defect-aware on-the-fly SiDB gate design. Thereby, SDB
// circuits can be designed in the presence of atomic defects.

int main()  // NOLINT
{
    using gate_lyt = fiction::hex_even_row_gate_clk_lyt;
    using cell_lyt = fiction::sidb_cell_clk_lyt;

    static const std::string layouts_folder =
        fmt::format("{}/physical_design_with_on_the_fly_gate_design/layouts", EXPERIMENTS_PATH);

    const std::vector<double> defect_concentrations = {1};

    for (const auto& concentration : defect_concentrations)
    {
        std::cout << fmt::format("--------------- defect concentration: {} % --------------- ", concentration)
                  << std::endl;
        auto surface_lattice_initial = fiction::read_sidb_surface_defects<cell_lyt>(
            fmt::format("../../experiments/defect_aware_physical_design/{}_percent_with_charged_surface.txt",
                        concentration),
            "py_test_surface");

        fiction::sidb_surface<cell_lyt> surface_lattice{
            {std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max()}};

        surface_lattice_initial.foreach_sidb_defect(
            [&surface_lattice](const auto& cd)
            {
                if (cd.second.type == fiction::sidb_defect_type::DB)
                {
                    surface_lattice.assign_sidb_defect(
                        cd.first, fiction::sidb_defect{fiction::sidb_defect_type::DB, -1, 4.1, 1.8});
                }
                else if (cd.second.type == fiction::sidb_defect_type::SI_VACANCY)
                {
                    surface_lattice.assign_sidb_defect(
                        cd.first, fiction::sidb_defect{fiction::sidb_defect_type::SI_VACANCY, -1, 10.6, 5.9});
                }
                else
                {
                    surface_lattice.assign_sidb_defect(cd.first, cd.second);
                }
            });

        const auto lattice_tiling = gate_lyt{{11, 30}};

        experiments::experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint64_t,
                                uint64_t, uint64_t, uint32_t, uint32_t, uint64_t, uint64_t, double, bool, uint64_t,
                                double>
            bestagon_exp{"bestagon",
                         "benchmark",
                         "inputs",
                         "outputs",
                         "initial nodes",
                         "initial depth",
                         "optimized nodes",
                         "optimized depth",
                         "layout width (in tiles)",
                         "layout height (in tiles)",
                         "layout area (in tiles)",
                         "gates",
                         "wires",
                         "critical path",
                         "throughput",
                         "runtime (in sec)",
                         "equivalent",
                         "SiDB dots",
                         "layout area in nm²"};

        // parameters for SMT-based physical design
        fiction::exact_physical_design_params<gate_lyt> exact_params{};
        exact_params.scheme        = fiction::ptr<gate_lyt>(fiction::row_clocking<gate_lyt>(fiction::num_clks::FOUR));
        exact_params.crossings     = false;
        exact_params.border_io     = false;
        exact_params.desynchronize = true;
        exact_params.upper_bound_x = 11;         // 12 x 31 tiles
        exact_params.upper_bound_y = 30;         // 12 x 31 tiles
        exact_params.timeout       = 3'600'000;  // 1h in ms

        constexpr const uint64_t bench_select =
            fiction_experiments::all & ~fiction_experiments::parity & ~fiction_experiments::two_bit_add_maj &
            ~fiction_experiments::b1_r2 & ~fiction_experiments::clpl & ~fiction_experiments::iscas85 &
            ~fiction_experiments::epfl & ~fiction_experiments::half_adder & ~fiction_experiments::full_adder &
            ~fiction_experiments::one_bit_add_aoig & ~fiction_experiments::one_bit_add_maj;

        //        constexpr const uint64_t bench_select =
        //            fiction_experiments::t;

        for (const auto& benchmark : fiction_experiments::all_benchmarks(bench_select))
        {
            fmt::print("[attempts] processing {}\n", benchmark);
            mockturtle::xag_network xag{};

            const auto read_verilog_result =
                lorina::read_verilog(fiction_experiments::benchmark_path(benchmark), mockturtle::verilog_reader(xag));
            assert(read_verilog_result == lorina::return_code::success);

            // compute depth
            const mockturtle::depth_view depth_xag{xag};

            const fiction::technology_mapping_params tech_map_params = fiction::all_2_input_functions();

            // parameters for cut rewriting
            mockturtle::cut_rewriting_params cut_params{};
            cut_params.cut_enumeration_ps.cut_size = 4;

            const mockturtle::xag_npn_resynthesis<
                mockturtle::xag_network,                    // the input network type
                mockturtle::xag_network,                    // the database network type
                mockturtle::xag_npn_db_kind::xag_complete>  // the kind of database to use
                resynthesis_function{};

            // rewrite network cuts using the given re-synthesis function
            const auto cut_xag = mockturtle::cut_rewriting(xag, resynthesis_function, cut_params);

            // perform technology mapping
            const auto mapped_network = fiction::technology_mapping(cut_xag, tech_map_params);
            // compute depth
            const mockturtle::depth_view depth_mapped_network{mapped_network};

            std::optional<gate_lyt> gate_level_layout = std::nullopt;
            cell_lyt                cell_level_layout{{}, "fail"};

            auto black_list = fiction::sidb_surface_analysis<fiction::sidb_skeleton_bestagon_library>(
                lattice_tiling, surface_lattice, std::make_pair(0, 0));

            uint64_t attempts = 0;

            mockturtle::stopwatch<>::duration time_counter{};

            bool gate_design_failed = true;

            {
                const mockturtle::stopwatch stop{time_counter};

                while (!gate_level_layout.has_value() || gate_design_failed)
                {
                    exact_params.black_list = black_list;
                    fiction::exact_physical_design_stats exact_stats{};
                    if (!gate_level_layout.has_value() && attempts > 0)
                    {
                        break;
                    }
                    std::cout << black_list.size() << std::endl;
                    gate_level_layout = fiction::exact<gate_lyt>(mapped_network, exact_params, &exact_stats);
                    if (gate_level_layout.has_value())
                    {
                        try
                        {
                            const fiction::design_sidb_gates_params<fiction::cube::coord_t> design_gate_params{
                                fiction::sidb_simulation_parameters{2, -0.32},
                                fiction::design_sidb_gates_params<
                                    fiction::cube::coord_t>::design_sidb_gates_mode::EXHAUSTIVE,
                                {{24, 17}, {34, 28}},
                                3,
                                fiction::sidb_simulation_engine::QUICKEXACT,
                                1};

                            auto parameter_gate_library = fiction::sidb_on_the_fly_gate_library_params<cell_lyt>{
                                surface_lattice, design_gate_params};

                            cell_level_layout = fiction::apply_parameterized_gate_library<
                                cell_lyt, fiction::sidb_on_the_fly_gate_library, gate_lyt,
                                fiction::sidb_on_the_fly_gate_library_params<cell_lyt>>(*gate_level_layout,
                                                                                        parameter_gate_library);
                            gate_design_failed = false;
                        }
                        catch (const fiction::gate_design_exception<fiction::tt, gate_lyt>& e)
                        {
                            gate_design_failed = true;
                            black_list[e.which_tile()][e.which_truth_table()].push_back(e.which_port_list());
                        }

                        catch (const std::exception& e)
                        {
                            std::cerr << "Caught std::exception: " << e.what() << std::endl;
                        }
                    }
                    attempts++;
                }
            }

            if (!gate_level_layout.has_value())
            {
                continue;
            }

            // check equivalence
            const auto miter = mockturtle::miter<fiction::technology_network>(mapped_network, *gate_level_layout);
            const auto eq    = mockturtle::equivalence_checking(*miter);
            assert(eq.has_value());

            // compute critical path and throughput
            fiction::critical_path_length_and_throughput_stats cp_tp_stats{};
            fiction::critical_path_length_and_throughput(*gate_level_layout, &cp_tp_stats);

            // apply dynamic gate library

            // compute area
            fiction::area_stats                            area_stats{};
            fiction::area_params<fiction::sidb_technology> area_ps{};
            fiction::area(cell_level_layout, area_ps, &area_stats);
            fiction::sidb_surface<cell_lyt> defect_surface{cell_level_layout};

            surface_lattice.foreach_sidb_defect([&defect_surface](const auto& defect)
                                                { defect_surface.assign_sidb_defect(defect.first, defect.second); });
            // write a SiQAD simulation file
            fiction::write_sqd_layout(defect_surface, fmt::format("{}/{}_{}_percent_after_big_change.sqd",
                                                                  layouts_folder, benchmark, concentration));

            // log results
            bestagon_exp(benchmark, xag.num_pis(), xag.num_pos(), xag.num_gates(), depth_xag.depth(),
                         mapped_network.num_gates(), depth_mapped_network.depth(), gate_level_layout->x() + 1,
                         gate_level_layout->y() + 1, (gate_level_layout->x() + 1) * (gate_level_layout->y() + 1),
                         gate_level_layout->num_gates(), gate_level_layout->num_wires(),
                         cp_tp_stats.critical_path_length, cp_tp_stats.throughput, mockturtle::to_seconds(time_counter),
                         *eq, cell_level_layout.num_cells(), area_stats.area);
            bestagon_exp.save();
            bestagon_exp.table();
        }
    }

    return EXIT_SUCCESS;
}

#else  // FICTION_Z3_SOLVER

#include <cstdlib>
#include <iostream>

int main()  // NOLINT
{
    std::cerr << "[e] Z3 solver is not available, please install Z3 and recompile the code" << std::endl;

    return EXIT_FAILURE;
}

#endif  // FICTION_Z3_SOLVER
