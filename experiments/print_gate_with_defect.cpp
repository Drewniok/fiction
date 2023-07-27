//
// Created by Jan Drewniok on 04.05.23.
//

#include "fiction/algorithms/simulation/sidb/quickexact.hpp"
#include "fiction/io/print_layout.hpp"
#include "fiction/io/read_sqd_layout.hpp"
#include "fiction/technology/sidb_defects.hpp"
#include "fiction/technology/sidb_surface.hpp"
#include "fiction/types.hpp"

#include <string>

using namespace fiction;

int main()  // NOLINT
{

    const quickexact_params<sidb_surface<sidb_cell_clk_lyt_siqad>> sim_params{sidb_simulation_parameters{2, -0.32}};

    auto lyt = read_sqd_layout<sidb_surface<sidb_cell_clk_lyt_siqad>>(
        "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/gate_with_defect/and/layout_found_11.sqd");
    // std::cout << lyt.num_defects() << std::endl;
    // std::cout << lyt.num_cells() << std::endl;
    //    lyt.foreach_sidb_defect(
    //        [&lyt](const auto& defect) {
    //            lyt.assign_sidb_defect(defect.first, sidb_defect{sidb_defect_type::UNKNOWN, -1, 5, 5.6});
    //        });

    sidb_surface<sidb_cell_clk_lyt_siqad> layout_with_all{};

    lyt.foreach_cell([&layout_with_all](const auto& cell)
                     { layout_with_all.assign_cell_type(cell, sidb_cell_clk_lyt_siqad::technology::NORMAL); });

    lyt.foreach_sidb_defect(
        [&layout_with_all](const auto& defect) {
            layout_with_all.assign_sidb_defect(defect.first, sidb_defect{sidb_defect_type::UNKNOWN, -1, 5, 5.6});
        });

    //    const sidb_defect defect{sidb_defect_type::DB, -1, 5, 5.6};
    //
    //    lyt.assign_sidb_defect({20, 8, 0}, defect);
    //    const auto sidb_defects = lyt.get_sidb_defect({20, 8, 0});

    const auto simulation_results = quickexact(layout_with_all, sim_params);

    auto lowest_energy = round_to_n_decimal_places(minimum_energy(simulation_results.charge_distributions), 6);
    charge_distribution_surface<sidb_surface<sidb_cell_clk_lyt_siqad>> lyt_copy{};
    for (const auto& lyt_loop : simulation_results.charge_distributions)
    {
        if (round_to_n_decimal_places(lyt_loop.get_system_energy(), 6) == lowest_energy)
        {
            lyt_copy = lyt_loop;
        }
    }

    const auto sidbs = lyt_copy.get_sidb_order();

    sidb_cell_clk_lyt_siqad new_layout{{50, 20}};

    for (const auto& sidb : sidbs)
    {
        new_layout.assign_cell_type(sidb, sidb_cell_clk_lyt_siqad::technology::NORMAL);
    }

    sidb_surface<sidb_cell_clk_lyt_siqad> defect{new_layout};

    charge_distribution_surface<sidb_surface<sidb_cell_clk_lyt_siqad>> new_charge_layout{defect};

    for (const auto& sidb : sidbs)
    {
        new_charge_layout.assign_charge_state(sidb, lyt_copy.get_charge_state(sidb));
    }

    lyt.foreach_sidb_defect(
        [&new_charge_layout](const auto& defect) {
            new_charge_layout.assign_sidb_defect(defect.first, sidb_defect{sidb_defect_type::UNKNOWN, -1, 5, 5.6});
        });

    print_charge_layout(std::cout, new_charge_layout, false, true, true);

    //    auto lyt_01 = read_sqd_layout<sidb_surface<sidb_cell_clk_lyt_siqad>>(
    //        "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/gate_with_defect/and/layout_found_01.sqd");
    //    const auto [avoidance_distance_01, defect_position_01] = maximal_defect_influence_distance(lyt_01,
    //    sim_params); std::cout << avoidance_distance_01 << std::endl; std::cout << "x: " <<
    //    std::to_string(defect_position_01.x); std::cout << " | y: " << std::to_string(defect_position_01.y); std::cout
    //    << " | z: " << std::to_string(defect_position_01.z) << std::endl;
    //
    //    auto lyt_10 = read_sqd_layout<sidb_surface<sidb_cell_clk_lyt_siqad>>(
    //        "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/gate_with_defect/and/layout_found_10.sqd");
    //    const auto [avoidance_distance_10, defect_position_10] = maximal_defect_influence_distance(lyt_10,
    //    sim_params); std::cout << avoidance_distance_10 << std::endl; std::cout << "x: " <<
    //    std::to_string(defect_position_10.x); std::cout << " | y: " << std::to_string(defect_position_10.y); std::cout
    //    << " | z: " << std::to_string(defect_position_10.z) << std::endl;
    //
    //    auto lyt_11 = read_sqd_layout<sidb_surface<sidb_cell_clk_lyt_siqad>>(
    //        "/Users/jandrewniok/CLionProjects/fiction_fork/experiments/gate_with_defect/and/layout_found_11.sqd");
    //    const auto [avoidance_distance_11, defect_position_11] = maximal_defect_influence_distance(lyt_11,
    //    sim_params); std::cout << avoidance_distance_11 << std::endl; std::cout << "x: " <<
    //    std::to_string(defect_position_11.x); std::cout << " | y: " << std::to_string(defect_position_11.y); std::cout
    //    << " | z: " << std::to_string(defect_position_11.z) << std::endl;

    return EXIT_SUCCESS;
}
