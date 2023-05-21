//
// Created by Jan Drewniok on 21.04.23.
//

#ifndef FICTION_LAYOUT_SIMULATION_HPP
#define FICTION_LAYOUT_SIMULATION_HPP

#include "fiction/algorithms/simulation/sidb/energy_distribution.hpp"
#include "fiction/algorithms/simulation/sidb/exhaustive_ground_state_simulation.hpp"
#include "fiction/algorithms/simulation/sidb/minimum_energy.hpp"
#include "fiction/algorithms/simulation/sidb/quicksim.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_result.hpp"
#include "fiction/io/write_sqd_layout.hpp"
#include "fiction/io/write_sqd_sim_result.hpp"
#include "fiction/layouts/bounding_box.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/technology/sidb_defects.hpp"

#include <fmt/format.h>
#include <mockturtle/utils/stopwatch.hpp>

#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace fiction

{

namespace detail
{

template <typename Lyt>
struct layout_sim_stats
{
    mockturtle::stopwatch<>::duration                          time_total{0};
    std::vector<std::vector<charge_distribution_surface<Lyt>>> all_valid_lyts{};

    void report(std::ostream& out = std::cout) const
    {
        out << fmt::format("total time  = {:.2f} secs\n", mockturtle::to_seconds(time_total));
        if (!all_valid_lyts.empty())
        {
            for (const auto& [energy, count] : energy_distribution<Lyt>(all_valid_lyts))
            {
                out << fmt::format("energy: {} | occurance: {} \n", energy, count);
            }
            out << fmt::format("the ground state energy is  = {:.4f} \n", minimum_energy(all_valid_lyts));
        }
        else
        {
            std::cout << "no state found | if two state simulation is used, continue with three state" << std::endl;
        }
        out << fmt::format("{} physically valid charge states were found \n", all_valid_lyts.size());
        std::cout << "_____________________________________________________ \n";
    }
};

template <typename Lyt>
class layout_simulation_impl
{
  public:
    layout_simulation_impl(Lyt& lyt, const sidb_simulation_parameters& params, layout_sim_stats<Lyt>& st) :
            layout{lyt},
            parameter{params},
            statistic{st}
    {
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
        static_assert(has_siqad_coord_v<Lyt>, "Lyt is not based on SiQAD coordinates");
        //        this->init();
        this->init_top_down();
    }

    void charge_distribution_to_index()
    {
        uint64_t counter  = 0;
        border_cell_index = 0;
        for (uint64_t i = 0; i < border_cells.size(); i++)
        {
            border_cell_index += static_cast<uint64_t>((charge_state_to_sign(border_cell_charge[i]) + 1) *
                                                       std::pow(2, this->num_cells() - 1 - counter));
            counter += 1;
        }
    }

    void index_to_charge_distribution() noexcept
    {

        auto       charge_quot = border_cell_index;
        const auto base        = 2;
        const auto num_charges = border_cells.size() - 1;
        auto       counter     = num_charges;

        if (charge_quot > 0)
        {
            while (charge_quot > 0)
            {
                const auto    charge_quot_int = static_cast<int64_t>(charge_quot);
                const auto    base_int        = static_cast<int64_t>(base);
                const int64_t quotient_int    = charge_quot_int / base_int;
                const int64_t remainder_int   = charge_quot_int % base_int;
                charge_quot                   = static_cast<uint64_t>(quotient_int);

                border_cell_charge[counter] = sign_to_charge_state(static_cast<int8_t>(remainder_int - 1));
                counter -= 1;
            }
        }
        else
        {
            for (auto i = 0u; i < border_cell_charge.size(); i++)
            {
                border_cell_charge[i] = sign_to_charge_state(static_cast<int8_t>(-1));
            }
        }
    }

    //    bool init()
    //    {
    //        //        const auto min_coordinate = wbb.get_min();
    //        //        std::cout << min_coordinate.x << std::endl;
    //        std::vector<cube::coord_t> all_cells = {};
    //        // obtain all cells in the surfaces and order them by their position to achieve a reproducible output
    //        layout.foreach_cell([this, &all_cells](const cell<Lyt>& c)
    //                            { all_cells.push_back(fiction::siqad::to_fiction_coord<cube::coord_t>(c)); });
    //
    //        int32_t x_min = 100000000;
    //        int32_t x_max = 0;
    //        int32_t y_min = 100000000;
    //        for (const auto& cells : all_cells)
    //        {
    //            if (cells.x < x_min)
    //            {
    //                x_min = cells.x;
    //            }
    //            if (cells.x > x_min)
    //            {
    //                x_min = cells.x;
    //            }
    //            if (cells.y < y_min)
    //            {
    //                y_min = cells.y;
    //            }
    //        }
    //        cells       = {};
    //        total_cells = {};
    //        start_cell  = cube::coord_t{x_min, y_min};
    //        std::cout << "finished" << std::endl;
    //    }

    bool init_top_down()
    {
        //        const auto min_coordinate = wbb.get_min();
        //        std::cout << min_coordinate.x << std::endl;
        std::vector<cube::coord_t> all_cells = {};
        // obtain all cells in the surfaces and order them by their position to achieve a reproducible output
        layout.foreach_cell([this, &all_cells](const cell<Lyt>& c)
                            { all_cells.push_back(fiction::siqad::to_fiction_coord<cube::coord_t>(c)); });

        std::sort(all_cells.begin(), all_cells.end());

        start_cell       = all_cells.front();
        start_cell.x     = start_cell.x - 60;
        left_corner_cell = start_cell;
        lowest_cell      = all_cells.back();
        int32_t x_max    = 0;
        for (const auto& cells : all_cells)
        {
            if (cells.x > x_max)
            {
                x_max = cells.x;
            }
        }
        rightest_cell = {x_max, 0};
        cells         = {};
        total_cells   = {};
        std::cout << "finished" << std::endl;
    }

    bool layout_generation_hexagon_real()
    {
        const uint64_t dim_x   = 26;
        const uint64_t dim_y   = 32;
        const uint64_t hex_dim = 60;

        uint64_t allowed   = start_cell.x + dim_x;
        uint64_t allowed_y = start_cell.y + dim_y;

        cells              = {};
        border_cells       = {};
        border_cell_charge = {};

        layout.foreach_cell(
            [&allowed, &allowed_y, this](const cell<Lyt>& c)
            {
                if (std::find(cells.begin(), cells.end(), c) == cells.end())
                {
                    auto cell_conv = fiction::siqad::to_fiction_coord<cube::coord_t>(c);
                    if (cell_conv.y <= allowed_y && cell_conv.y >= start_cell.y && cell_conv.x <= allowed &&
                        cell_conv.x >= start_cell.x)
                    {
                        cells.insert(c);
                        total_cells.insert(c);
                    }
                }
            });

        if (!cells.empty())
        {
            region_counter += 1;
            std::cout << "part cells: " << std::to_string(cells.size()) << std::endl;
            std::cout << "all collected cells: " << std::to_string(total_cells.size()) << std::endl;
            std::cout << allowed << std::endl;
            std::cout << allowed_y << std::endl;

            layout.foreach_cell(
                [this](const cell<Lyt>& c)
                {
                    if (std::find(cells.begin(), cells.end(), c) == cells.end())
                    {
                        uint64_t counter = 0;
                        for (auto it = cells.begin(); it != cells.end(); it++)
                        {
                            if (sidb_nanometer_distance<Lyt>(layout, *it, c, parameter) < 6)
                            {
                                counter += 1;
                            }
                        }
                        if (counter != 0)
                        {
                            border_cells.push_back(c);
                            border_cell_charge.push_back(sidb_charge_state::NEUTRAL);
                        }
                    }
                });

            border_cell_region.insert(std::make_pair(region_counter, border_cells));

            border_cell_max_charge_index = std::pow(2, border_cells.size()) - 1;
        }
        if (allowed < rightest_cell.x)
        {
            start_cell.x = start_cell.x + hex_dim;
        }
        else
        {
            region_col_counter += 1;
            if (region_col_counter % 2 == 0)
            {
                start_cell.x = left_corner_cell.x;
                start_cell.y = left_corner_cell.y + 34 * region_col_counter;
            }
            else
            {
                start_cell.x = left_corner_cell.x + (hex_dim) / 2;
                start_cell.y = left_corner_cell.y + 34 * region_col_counter;
            }
        }
        return EXIT_SUCCESS;
    }

    //    bool layout_generation_hexagon()
    //    {
    //        cells           = {};
    //        uint64_t length = 0;
    //
    //        uint64_t allowed   = start_cell.x + 26;
    //        uint64_t allowed_y = start_cell.y + 32;
    //
    //        border_cells       = {};
    //        border_cell_charge = {};
    //
    //        //        cells.insert(start_cell);
    //        //        total_cells.insert(start_cell);
    //
    //        //        if (start_cell == siqad::to_fiction_coord<cube::coord_t>(left_corner_cell))
    //        //        {
    //        //            cells.insert(siqad::to_siqad_coord(start_cell));
    //        //            total_cells.insert(siqad::to_siqad_coord(start_cell));
    //        //        }
    //
    //        layout.foreach_cell(
    //            [&allowed, &allowed_y, this](const cell<Lyt>& c)
    //            {
    //                if (std::find(cells.begin(), cells.end(), c) == cells.end())
    //                {
    //                    auto cell_conv = fiction::siqad::to_fiction_coord<cube::coord_t>(c);
    //                    if (cell_conv.y <= allowed_y && cell_conv.y >= start_cell.y)
    //                    {
    //                        cells.insert(c);
    //                        total_cells.insert(c);
    //                    }
    //                }
    //            });
    //
    //        //            std::cout << allowed << std::endl;
    //        //            std::cout << allowed_y << std::endl;
    //        //            std::cout << "all collected cells: " << std::to_string(total_cells.size()) << std::endl;
    //        //            std::cout << "all collected cells: " << std::to_string(cells.size()) << std::endl;
    //
    //        std::cout << "part cells: " << std::to_string(cells.size()) << std::endl;
    //        std::cout << "all collected cells: " << std::to_string(total_cells.size()) << std::endl;
    //        std::cout << allowed << std::endl;
    //        std::cout << allowed_y << std::endl;
    //
    //        layout.foreach_cell(
    //            [this](const cell<Lyt>& c)
    //            {
    //                if (std::find(cells.begin(), cells.end(), c) == cells.end())
    //                {
    //                    uint64_t counter = 0;
    //                    for (auto it = cells.begin(); it != cells.end(); it++)
    //                    {
    //                        if (sidb_nanometer_distance<Lyt>(layout, *it, c, parameter) < 6)
    //                        {
    //                            //                            if (c.x > siqad::to_siqad_coord(start_cell).x && c.y >
    //                            //                            siqad::to_siqad_coord(start_cell).y)
    //                            //                            {
    //                            counter += 1;
    //                            //}
    //                        }
    //                    }
    //                    if (counter != 0)
    //                    {
    //                        border_cells.push_back(c);
    //                        border_cell_charge.push_back(sidb_charge_state::NEUTRAL);
    //                    }
    //                }
    //            });
    //
    //        border_cell_max_charge_index = std::pow(2, border_cells.size()) - 1;
    //        start_cell.x                 = start_cell.x + length;
    //        start_cell.y                 = start_cell.y + length;
    //    }

    //    bool layout_generation_hexagon()
    //    {
    //        cells           = {};
    //        uint64_t length = 0;
    //
    //        uint64_t allowed   = start_cell.x + 26;
    //        uint64_t allowed_y = start_cell.y + 32;
    //
    //        border_cells       = {};
    //        border_cell_charge = {};
    //
    //        //        cells.insert(start_cell);
    //        //        total_cells.insert(start_cell);
    //
    //        //        if (start_cell == siqad::to_fiction_coord<cube::coord_t>(left_corner_cell))
    //        //        {
    //        //            cells.insert(siqad::to_siqad_coord(start_cell));
    //        //            total_cells.insert(siqad::to_siqad_coord(start_cell));
    //        //        }
    //
    //            layout.foreach_cell(
    //                [&allowed, &allowed_y, this](const cell<Lyt>& c)
    //                {
    //                    if (std::find(cells.begin(), cells.end(), c) == cells.end())
    //                    {
    //                        auto cell_conv = fiction::siqad::to_fiction_coord<cube::coord_t>(c);
    //                        if (cell_conv.y <= allowed_y && cell_conv.y >= start_cell.y)
    //                        {
    //                            cells.insert(c);
    //                            total_cells.insert(c);
    //                        }
    //                    }
    //                });
    //
    //            //            std::cout << allowed << std::endl;
    //            //            std::cout << allowed_y << std::endl;
    //            //            std::cout << "all collected cells: " << std::to_string(total_cells.size()) << std::endl;
    //            //            std::cout << "all collected cells: " << std::to_string(cells.size()) << std::endl;
    //
    //        std::cout << "part cells: " << std::to_string(cells.size()) << std::endl;
    //        std::cout << "all collected cells: " << std::to_string(total_cells.size()) << std::endl;
    //        std::cout << allowed << std::endl;
    //        std::cout << allowed_y << std::endl;
    //
    //        layout.foreach_cell(
    //            [this](const cell<Lyt>& c)
    //            {
    //                if (std::find(cells.begin(), cells.end(), c) == cells.end())
    //                {
    //                    uint64_t counter = 0;
    //                    for (auto it = cells.begin(); it != cells.end(); it++)
    //                    {
    //                        if (sidb_nanometer_distance<Lyt>(layout, *it, c, parameter) < 6)
    //                        {
    //                            //                            if (c.x > siqad::to_siqad_coord(start_cell).x && c.y >
    //                            //                            siqad::to_siqad_coord(start_cell).y)
    //                            //                            {
    //                            counter += 1;
    //                            //}
    //                        }
    //                    }
    //                    if (counter != 0)
    //                    {
    //                        border_cells.push_back(c);
    //                        border_cell_charge.push_back(sidb_charge_state::NEUTRAL);
    //                    }
    //                }
    //            });
    //
    //        border_cell_max_charge_index = std::pow(2, border_cells.size()) - 1;
    //        start_cell.x                 = start_cell.x + length;
    //        start_cell.y                 = start_cell.y + length;
    //    }

    //    bool layout_generation_top_down()
    //    {
    //        cells           = {};
    //        uint64_t length = 0;
    //
    //        uint64_t allowed   = start_cell.x;
    //        uint64_t allowed_y = start_cell.y + length;
    //
    //        border_cells       = {};
    //        border_cell_charge = {};
    //
    //        //        cells.insert(start_cell);
    //        //        total_cells.insert(start_cell);
    //
    //        //        if (start_cell == siqad::to_fiction_coord<cube::coord_t>(left_corner_cell))
    //        //        {
    //        //            cells.insert(siqad::to_siqad_coord(start_cell));
    //        //            total_cells.insert(siqad::to_siqad_coord(start_cell));
    //        //        }
    //
    //        while (total_cells.size() < layout.num_cells() && cells.size() < 40)
    //        {
    //            layout.foreach_cell(
    //                [&allowed, &allowed_y, this](const cell<Lyt>& c)
    //                {
    //                    if (std::find(cells.begin(), cells.end(), c) == cells.end())
    //                    {
    //                        auto cell_conv = fiction::siqad::to_fiction_coord<cube::coord_t>(c);
    //                        if (cell_conv.y <= allowed_y && cell_conv.y >= start_cell.y)
    //                        {
    //                            cells.insert(c);
    //                            total_cells.insert(c);
    //                        }
    //                        //                        if (cell_conv.x == start_cell.x && cell_conv.y <= allowed_y &&
    //                        //                        cell_conv.y > start_cell.y)
    //                        //                        {
    //                        //                            cells.insert(c);
    //                        //                            total_cells.insert(c);
    //                        //                        }
    //                    }
    //                });
    //            length += 1;
    //            allowed   = start_cell.x + length;
    //            allowed_y = start_cell.y + length;
    //            //            std::cout << allowed << std::endl;
    //            //            std::cout << allowed_y << std::endl;
    //            //            std::cout << "all collected cells: " << std::to_string(total_cells.size()) << std::endl;
    //            //            std::cout << "all collected cells: " << std::to_string(cells.size()) << std::endl;
    //        }
    //        std::cout << "part cells: " << std::to_string(cells.size()) << std::endl;
    //        std::cout << "all collected cells: " << std::to_string(total_cells.size()) << std::endl;
    //        std::cout << allowed << std::endl;
    //        std::cout << allowed_y << std::endl;
    //
    //        layout.foreach_cell(
    //            [this](const cell<Lyt>& c)
    //            {
    //                if (std::find(cells.begin(), cells.end(), c) == cells.end())
    //                {
    //                    uint64_t counter = 0;
    //                    for (auto it = cells.begin(); it != cells.end(); it++)
    //                    {
    //                        if (sidb_nanometer_distance<Lyt>(layout, *it, c, parameter) < 6)
    //                        {
    //                            //                            if (c.x > siqad::to_siqad_coord(start_cell).x && c.y >
    //                            //                            siqad::to_siqad_coord(start_cell).y)
    //                            //                            {
    //                            counter += 1;
    //                            //}
    //                        }
    //                    }
    //                    if (counter != 0)
    //                    {
    //                        border_cells.push_back(c);
    //                        border_cell_charge.push_back(sidb_charge_state::NEUTRAL);
    //                    }
    //                }
    //            });
    //
    //        border_cell_max_charge_index = std::pow(2, border_cells.size()) - 1;
    //        start_cell.x                 = start_cell.x + length;
    //        start_cell.y                 = start_cell.y + length;
    //    }

    //    bool layout_generation()
    //    {
    //        cells           = {};
    //        uint64_t length = 0;
    //
    //        uint64_t allowed   = start_cell.x + length;
    //        uint64_t allowed_y = start_cell.y + length;
    //
    //        border_cells       = {};
    //        border_cell_charge = {};
    //
    //        //        cells.insert(start_cell);
    //        //        total_cells.insert(start_cell);
    //
    //        //        if (start_cell == siqad::to_fiction_coord<cube::coord_t>(left_corner_cell))
    //        //        {
    //        //            cells.insert(siqad::to_siqad_coord(start_cell));
    //        //            total_cells.insert(siqad::to_siqad_coord(start_cell));
    //        //        }
    //
    //        while (total_cells.size() < layout.num_cells() && cells.size() < 27)
    //        {
    //            layout.foreach_cell(
    //                [&allowed, &allowed_y, this](const cell<Lyt>& c)
    //                {
    //                    if (std::find(cells.begin(), cells.end(), c) == cells.end())
    //                    {
    //                        auto cell_conv = fiction::siqad::to_fiction_coord<cube::coord_t>(c);
    //                        if (cell_conv.x < start_cell.x && cell_conv.y <= allowed_y && cell_conv.y >= start_cell.y)
    //                        {
    //                            cells.insert(c);
    //                            total_cells.insert(c);
    //                        }
    //                        if (cell_conv.x >= start_cell.x && cell_conv.y <= allowed_y && cell_conv.x <= allowed)
    //                        {
    //                            cells.insert(c);
    //                            total_cells.insert(c);
    //                        }
    //                        //                        if (cell_conv.x == start_cell.x && cell_conv.y <= allowed_y &&
    //                        //                        cell_conv.y > start_cell.y)
    //                        //                        {
    //                        //                            cells.insert(c);
    //                        //                            total_cells.insert(c);
    //                        //                        }
    //                    }
    //                });
    //            length += 1;
    //            allowed   = start_cell.x + length;
    //            allowed_y = start_cell.y + length;
    //            //            std::cout << allowed << std::endl;
    //            //            std::cout << allowed_y << std::endl;
    //            //            std::cout << "all collected cells: " << std::to_string(total_cells.size()) << std::endl;
    //            //            std::cout << "all collected cells: " << std::to_string(cells.size()) << std::endl;
    //        }
    //        std::cout << "part cells: " << std::to_string(cells.size()) << std::endl;
    //        std::cout << "all collected cells: " << std::to_string(total_cells.size()) << std::endl;
    //        std::cout << allowed << std::endl;
    //        std::cout << allowed_y << std::endl;
    //
    //        layout.foreach_cell(
    //            [this](const cell<Lyt>& c)
    //            {
    //                if (std::find(cells.begin(), cells.end(), c) == cells.end())
    //                {
    //                    uint64_t counter = 0;
    //                    for (auto it = cells.begin(); it != cells.end(); it++)
    //                    {
    //                        if (sidb_nanometer_distance<Lyt>(layout, *it, c, parameter) < 3)
    //                        {
    //                            //                            if (c.x > siqad::to_siqad_coord(start_cell).x && c.y >
    //                            //                            siqad::to_siqad_coord(start_cell).y)
    //                            //                            {
    //                            counter += 1;
    //                            //}
    //                        }
    //                    }
    //                    if (counter != 0)
    //                    {
    //                        border_cells.push_back(c);
    //                        border_cell_charge.push_back(sidb_charge_state::NEUTRAL);
    //                    }
    //                }
    //            });
    //
    //        border_cell_max_charge_index = std::pow(2, border_cells.size()) - 1;
    //        start_cell.x                 = start_cell.x + length;
    //        start_cell.y                 = start_cell.y + length;
    //    }

    uint64_t charge_distribution_external_to_index(std::vector<typename Lyt::cell>&        defect_vector,
                                                   const charge_distribution_surface<Lyt>& lyt) const
    {
        std::sort(defect_vector.begin(), defect_vector.end());
        auto sidbs = lyt.get_sidb_order();
        std::sort(sidbs.begin(), sidbs.end());
        uint64_t counter = 0;
        uint64_t index   = 0;
        for (const auto& cell : defect_vector)
        {
            if (std::find(sidbs.begin(), sidbs.end(), cell) != sidbs.end())
            {
                index += static_cast<uint64_t>((charge_state_to_sign(lyt.get_charge_state(cell)) + 1) *
                                               std::pow(2, defect_vector.size() - 1 - counter));
                counter += 1;
            }
        }
        return index;
    }

    bool defect_map_update()
    {
        defect_cell   = {};
        defect_charge = {};
        for (uint i = 0; i < border_cells.size(); i++)
        {
            defect_cell.push_back(border_cells[i]);
            defect_charge.push_back(charge_state_to_sign(border_cell_charge[i]));
        }
    }

    //    bool run_simulation()
    //    {
    //        uint64_t counter = 0;
    //        while (total_cells.size() < layout.num_cells())
    //        {
    ////            this->layout_generation();
    //            this->layout_generation_hexagon();
    //
    //            Lyt lyt{};
    //
    //            for (const auto& cell : cells)
    //            {
    //                lyt.assign_cell_type(cell, Lyt::cell_type::NORMAL);
    //            }
    //
    //            layout_num = lyt.num_cells();
    //
    //            std::cout << "border cell max index: " << std::to_string(border_cell_max_charge_index) << std::endl;
    //            std::cout << "border cells: " << std::to_string(border_cells.size()) << std::endl;
    //            std::set<uint64_t> charge_index{};
    //
    //            border_cell_index    = 0;
    //            auto lyts_collection = std::vector<charge_distribution_surface<Lyt>>{};
    //            std::vector<std::unordered_map<typename Lyt::cell, const sidb_defect>> all_defect_confs{};
    //            while (border_cell_index <= border_cell_max_charge_index)
    //            {
    //                this->index_to_charge_distribution();
    //
    //                this->defect_map_update();
    //
    //                //                std::vector<sidb_charge_state> defect_charges{};
    //                //                for (const auto& charge : defect_charge)
    //                //                {
    //                //                    defect_charges.push_back(sign_to_charge_state(charge));
    //                //                }
    //
    //                //    std::cout << charge_configuration_to_string(defect_charges) << std::endl;
    //                // std::cout << border_cell_index << std::endl;
    //                std::unordered_map<typename Lyt::cell, const sidb_defect> defect{};
    //                for (uint64_t i = 0; i < defect_cell.size(); i++)
    //                {
    //                    defect.insert({defect_cell[i],
    //                                   sidb_defect{sidb_defect_type::UNKNOWN,
    //                                   static_cast<double>(defect_charge[i])}});
    //                }
    //                border_cell_index += 1;
    //                all_defect_confs.push_back(defect);
    //            }
    //            std::cout << "defect confs:" << std::to_string(all_defect_confs.size()) << std::endl;
    //            // std::cout << "defects: " << defect.size() << std::endl;
    //            exgs_stats<Lyt> exgs_stats{};
    //            // std::vector<uint64_t> charge_index_in{};
    //            //exhaustive_ground_state_simulation(lyt, parameter, &exgs_stats, all_defect_confs);
    //            for (const auto& lyt_loop : exgs_stats.valid_lyts)
    //            {
    //                charge_index.insert(lyt_loop.get_charge_index().first);
    //                // charge_index_in.emplace_back(charge_distribution_external_to_index(defect_cell, lyt_loop));
    //            }
    //            // charge_index_innen.emplace_back(charge_index_in);
    //            all_charge_lyts.emplace_back(exgs_stats.valid_lyts);
    //            for (const auto& solution_lyts : exgs_stats.valid_lyts)
    //            {
    //                lyts_collection.emplace_back(solution_lyts);
    //            }
    //            region_num.emplace_back(counter);
    //
    //            all_defect_charges.emplace_back(defect_charge);
    //            // border_cell_index += 1;
    //
    //            std::vector<charge_distribution_surface<Lyt>> unique_lyts{};
    //            for (const auto& index : charge_index)
    //            {
    //                for (const auto& lyt : lyts_collection)
    //                {
    //                    if (lyt.get_charge_index().first == index)
    //                    {
    //                        unique_lyts.push_back(lyt);
    //                        break;
    //                    }
    //                }
    //            }
    //            lyts_of_regions.emplace_back(unique_lyts);
    //            all_defect_cells.emplace_back(defect_cell);
    //
    //            std::cout << "number valid lyts: " << charge_index.size() << std::endl;
    //
    //            //                        if (charge_index.size() == 1)
    //            //                        {
    //            //                            write_sqd_layout(lyt, "/Users/jandrewniok/Desktop/investi/" +
    //            //                            std::to_string(total_cells.size()));
    //            //                        }
    //
    //            if (border_cell_max_charge_index == 0)
    //            {
    //                fiction::exgs_stats<Lyt> exgs_stats_second{};
    //                exhaustive_ground_state_simulation(lyt, parameter, &exgs_stats_second);
    //                all_charge_lyts.emplace_back(exgs_stats_second.valid_lyts);
    //                region_num.emplace_back(counter);
    //            }
    //            // std::cout << all_charge_lyts.size() << std::endl;
    //            counter += 1;
    //        }
    //
    //        return true;
    //    };

    bool run_simulation_hexagon()
    {
        while (total_cells.size() < layout.num_cells())
        {
            uint64_t counter = 0;

            this->layout_generation_hexagon_real();
            if (!cells.empty())
            {
                Lyt lyt{};

                for (const auto& cell : cells)
                {
                    lyt.assign_cell_type(cell, Lyt::cell_type::NORMAL);
                }

                all_layouts.emplace_back(lyt);

                layout_num = lyt.num_cells();

                std::cout << "border cell max index: " << std::to_string(border_cell_max_charge_index) << std::endl;
                std::cout << "border cells: " << std::to_string(border_cells.size()) << std::endl;
                std::cout << "total cells: " << std::to_string(layout.num_cells()) << std::endl;
                std::set<uint64_t> charge_index{};

                border_cell_index    = 0;
                auto lyts_collection = std::vector<charge_distribution_surface<Lyt>>{};
                std::vector<std::unordered_map<typename Lyt::cell, const sidb_defect>> all_defect_confs{};
                std::vector<std::vector<int8_t>>                                       all_charge_confs{};
                while (border_cell_index <= border_cell_max_charge_index)
                {
                    this->index_to_charge_distribution();

                    this->defect_map_update();

                    auto num_pairs = defect_cell.size() / 2;
                    if (defect_cell.size() % 2 != 0)
                    {
                        std::cout << "not" << std::endl;
                    }

                    uint64_t charge_counter = 0;
                    for (const auto& charge_sign : defect_charge)
                    {
                        if (charge_sign == -1)
                        {
                            charge_counter += 1;
                        }
                    }
                    if (charge_counter != num_pairs)
                    {
                        border_cell_index += 1;
                        continue;
                    }

                    typename Lyt::cell cells1{};
                    typename Lyt::cell cells2{};
                    uint64_t           counter_neg = 0;
                    for (auto i = 0u; i < defect_cell.size(); i++)
                    {
                        if (defect_charge[i] == -1 && counter_neg == 0)
                        {
                            cells1 = defect_cell[i];
                            counter_neg += 1;
                        }
                        else if (defect_charge[i] == -1 && counter_neg == 1)
                        {
                            cells2 = defect_cell[i];
                        }
                    }

                    if (sidb_nanometer_distance<Lyt>(layout, cells1, cells2, parameter) < 1.5)
                    {
                        border_cell_index += 1;
                        continue;
                    }

                    std::vector<int8_t> charges{};
                    for (auto i = 0u; i < defect_charge.size(); i++)
                    {
                        charges.push_back(defect_charge[i]);
                    }

                    //    std::cout << charge_configuration_to_string(defect_charges) << std::endl;
                    // std::cout << border_cell_index << std::endl;
                    std::unordered_map<typename Lyt::cell, const sidb_defect> defect{};
                    std::unordered_map<typename Lyt::cell, int8_t>            charge_collection{};
                    for (uint64_t i = 0; i < defect_cell.size(); i++)
                    {
                        defect.insert({defect_cell[i],
                                       sidb_defect{sidb_defect_type::UNKNOWN, static_cast<double>(defect_charge[i])}});
                        // charge_collection.insert({defect_cell[i], static_cast<int8_t>(defect_charge[i])});
                    }
                    border_cell_index += 1;
                    all_defect_confs.push_back(defect);
                    all_charge_confs.push_back(charges);
                }

                std::cout << "defect confs:" << std::to_string(all_defect_confs.size()) << std::endl;
                exgs_stats<Lyt> exgs_stats{};
                exhaustive_ground_state_simulation(lyt, parameter, &exgs_stats, all_defect_confs);

                for (auto i = 0u; i < exgs_stats.valid_lyts.size(); i++)
                {
                    charge_index.insert(exgs_stats.valid_lyts[i].get_charge_index().first);
                    // border_cells_and_charge.push_back(all_charge_confs);
                    //  charge_index_in.emplace_back(charge_distribution_external_to_index(defect_cell, lyt_loop));
                }
                // charge_index_innen.emplace_back(charge_index_in);
                all_charge_lyts.emplace_back(exgs_stats.valid_lyts);
                region_num.emplace_back(counter);

                // all_defect_charges.emplace_back(defect_charge);
                //  border_cell_index += 1;

                std::vector<charge_distribution_surface<Lyt>> unique_lyts{};
                std::vector<uint64_t>                         unique_defect_confs_index{};
                for (const auto& index : charge_index)
                {
                    for (auto i = 0u; i < exgs_stats.valid_lyts.size(); i++)
                    {
                        if (exgs_stats.valid_lyts[i].get_charge_index().first == index)
                        {
                            unique_lyts.push_back(exgs_stats.valid_lyts[i]);
                            unique_defect_confs_index.push_back(exgs_stats.defect_iter_num_valid_lyts[i].second);
                            break;
                        }
                    }
                }
                lyts_of_regions.emplace_back(unique_lyts);

                all_defect_cells.emplace_back(defect_cell);

                std::vector<std::vector<int8_t>> unique_charge_confs{};
                for (auto i = 0u; i < unique_defect_confs_index.size(); i++)
                {
                    unique_charge_confs.push_back(all_charge_confs[unique_defect_confs_index[i]]);
                }
                all_defect_charges.push_back(unique_charge_confs);

                std::cout << "number valid lyts: " << charge_index.size() << std::endl;

                //                        if (charge_index.size() == 1)
                //                        {
                //                            write_sqd_layout(lyt, "/Users/jandrewniok/Desktop/investi/" +
                //                            std::to_string(total_cells.size()));
                //                        }

                if (border_cell_max_charge_index == 0 && lyt.num_cells() != 0)
                {
                    fiction::exgs_stats<Lyt> exgs_stats_second{};
                    exhaustive_ground_state_simulation(lyt, parameter, &exgs_stats_second);
                    all_charge_lyts.emplace_back(exgs_stats_second.valid_lyts);
                    region_num.emplace_back(counter);
                }
                // std::cout << all_charge_lyts.size() << std::endl;
                counter += 1;
            }
        }
        return true;
    };

    void finding_nn()
    {
        std::vector<std::set<uint64_t>> all_pairs{};
        for (auto i = 0u; i < all_defect_cells.size(); i++)
        {
            std::set<uint64_t> layout_numbers{};
            for (const auto& cell : all_defect_cells[i])
            {
                for (auto j = 0u; j < all_layouts.size(); j++)
                {
                    uint64_t counter = 0;
                    if (j != i)
                    {
                        all_layouts[j].foreach_cell(
                            [this, &cell, &counter](const auto& c)
                            {
                                if (cell == c)
                                {
                                    counter += 1;
                                }
                            });
                    }
                    if (counter != 0)
                    {
                        layout_numbers.insert(j);
                    }
                }
            }
            all_pairs.push_back(layout_numbers);
        }
        all_neighbor_pairs = all_pairs;
    }

    void combining_four()
    {
        auto compareFunc =
            [](const charge_distribution_surface<Lyt>& lyt1, const charge_distribution_surface<Lyt>& lyt2)
        { return lyt1.get_system_energy() < lyt2.get_system_energy(); };

        std::cout << "combining starts" << std::endl;
        std::cout << lyts_of_regions.size() << std::endl;
        auto lyt_one   = lyts_of_regions[0];
        auto lyt_two   = lyts_of_regions[1];
        auto lyt_three = lyts_of_regions[2];
        auto lyt_four  = lyts_of_regions[3];

        std::sort(lyt_one.begin(), lyt_one.end(), compareFunc);
        std::sort(lyt_two.begin(), lyt_two.end(), compareFunc);
        std::sort(lyt_three.begin(), lyt_three.end(), compareFunc);
        std::sort(lyt_four.begin(), lyt_four.end(), compareFunc);

        charge_distribution_surface<Lyt> charge_lyt{layout};
        uint64_t                         counter = 0;
        std::vector<double>              valid_energies{};
        double                           energy_threas = 1000;
        for (const auto& lyts_one : lyt_one)
        {
            for (const auto& lyts_two : lyt_two)
            {
                for (const auto& lyts_three : lyt_three)
                {
                    for (const auto& lyts_four : lyt_four)
                    {
                        {
                            lyts_one.foreach_cell(
                                [this, &charge_lyt, &lyts_one](const auto& c1)
                                { charge_lyt.assign_charge_state(c1, lyts_one.get_charge_state(c1), false); });
                            lyts_two.foreach_cell(
                                [this, &charge_lyt, &lyts_two](const auto& c1)
                                { charge_lyt.assign_charge_state(c1, lyts_two.get_charge_state(c1), false); });
                            lyts_three.foreach_cell(
                                [this, &charge_lyt, &lyts_three](const auto& c1)
                                { charge_lyt.assign_charge_state(c1, lyts_three.get_charge_state(c1), false); });
                            lyts_four.foreach_cell(
                                [this, &charge_lyt, &lyts_four](const auto& c1)
                                { charge_lyt.assign_charge_state(c1, lyts_four.get_charge_state(c1), false); });
                            charge_lyt.update_after_charge_change();
                            if (charge_lyt.is_physically_valid())
                            {
                                if (charge_lyt.get_system_energy() < energy_threas)
                                {
                                    std::vector<charge_distribution_surface<Lyt>> lyts{};
                                    std::cout << charge_lyt.get_system_energy() << std::endl;

                                    sidb_simulation_result<Lyt> sim_result{};
                                    sim_result.algorithm_name = "ExGS";
                                    charge_distribution_surface<Lyt> charge_lyt_copy{charge_lyt};
                                    lyts.emplace_back(charge_lyt_copy);
                                    sim_result.charge_distributions = lyts;
                                    energy_threas                   = charge_lyt.get_system_energy();
                                    write_sqd_sim_result<Lyt>(sim_result, "/Users/jandrewniok/CLionProjects/"
                                                                          "fiction_fork/experiments/result.xml");
                                }
                            }
                            counter += 1;

                            // std::cout << counter << std::endl;
                        }
                    }
                }
            }
        }
    }

    void combining_three()
    {
        uint64_t counter_lyts = 1;
        for (const auto& lyts_region : lyts_of_regions)
        {
            counter_lyts *= lyts_region.size();
        }
        std::cout << "number enumerations: " << std::to_string(counter_lyts) << std::endl;

        auto compareFunc =
            [](const charge_distribution_surface<Lyt>& lyt1, const charge_distribution_surface<Lyt>& lyt2)
        { return lyt1.get_system_energy() < lyt2.get_system_energy(); };

        std::cout << "combining starts" << std::endl;
        std::cout << lyts_of_regions.size() << std::endl;
        auto lyt_one   = lyts_of_regions[0];
        auto lyt_two   = lyts_of_regions[1];
        auto lyt_three = lyts_of_regions[2];
        auto lyt_four  = lyts_of_regions[3];

        std::sort(lyt_one.begin(), lyt_one.end(), compareFunc);
        std::sort(lyt_two.begin(), lyt_two.end(), compareFunc);
        std::sort(lyt_three.begin(), lyt_three.end(), compareFunc);
        std::sort(lyt_four.begin(), lyt_four.end(), compareFunc);

        charge_distribution_surface<Lyt> charge_lyt{layout};
        uint64_t                         counter = 0;
        std::vector<double>              valid_energies{};
        double                           energy_threas = 1000;
        for (const auto& lyts_one : lyt_one)
        {
            for (const auto& lyts_two : lyt_two)
            {
                for (const auto& lyts_three : lyt_three)
                {
                    for (const auto& lyts_four : lyt_four)
                    {
                        lyts_one.foreach_cell(
                            [this, &charge_lyt, &lyts_one](const auto& c1)
                            { charge_lyt.assign_charge_state(c1, lyts_one.get_charge_state(c1), false); });
                        lyts_two.foreach_cell(
                            [this, &charge_lyt, &lyts_two](const auto& c1)
                            { charge_lyt.assign_charge_state(c1, lyts_two.get_charge_state(c1), false); });
                        lyts_three.foreach_cell(
                            [this, &charge_lyt, &lyts_three](const auto& c1)
                            { charge_lyt.assign_charge_state(c1, lyts_three.get_charge_state(c1), false); });
                        lyts_four.foreach_cell(
                            [this, &charge_lyt, &lyts_four](const auto& c1)
                            { charge_lyt.assign_charge_state(c1, lyts_four.get_charge_state(c1), false); });

                        charge_lyt.update_after_charge_change();
                        if (charge_lyt.is_physically_valid())
                        {
                            if (charge_lyt.get_system_energy() < energy_threas)
                            {
                                std::vector<charge_distribution_surface<Lyt>> lyts{};
                                std::cout << charge_lyt.get_system_energy() << std::endl;

                                sidb_simulation_result<Lyt> sim_result{};
                                sim_result.algorithm_name = "ExGS";
                                charge_distribution_surface<Lyt> charge_lyt_copy{charge_lyt};
                                lyts.emplace_back(charge_lyt_copy);
                                sim_result.charge_distributions = lyts;
                                energy_threas                   = charge_lyt.get_system_energy();
                                write_sqd_sim_result<Lyt>(sim_result, "/Users/jandrewniok/CLionProjects/"
                                                                      "fiction_fork/experiments/result.xml");
                            }
                        }
                        counter += 1;

                        if (counter % 10000000 == 0)
                        {
                            std::cout << counter << std::endl;
                        }
                    }
                }
            }
        }

        // std::cout << *std::min_element(valid_energies.begin(), valid_energies.end()) << std::endl;
    }

    void combining_seven()
    {
        auto compareFunc =
            [](const charge_distribution_surface<Lyt>& lyt1, const charge_distribution_surface<Lyt>& lyt2)
        { return lyt1.get_system_energy() < lyt2.get_system_energy(); };

        std::cout << "combining starts: " << std::to_string(lyts_of_regions.size()) << std::endl;
        uint64_t counter_lyts = 1;
        for (const auto& lyts_region : lyts_of_regions)
        {
            counter_lyts *= lyts_region.size();
        }
        std::cout << "number enumerations: " << std::to_string(counter_lyts) << std::endl;

        auto lyt_one   = lyts_of_regions[0];
        auto lyt_two   = lyts_of_regions[1];
        auto lyt_three = lyts_of_regions[2];
        auto lyt_four  = lyts_of_regions[3];
        auto lyt_five  = lyts_of_regions[4];
        auto lyt_six   = lyts_of_regions[5];
        auto lyt_seven = lyts_of_regions[6];

        //                std::sort(lyt_one.begin(), lyt_one.end(), compareFunc);
        //                std::sort(lyt_two.begin(), lyt_two.end(), compareFunc);
        //                std::sort(lyt_three.begin(), lyt_three.end(), compareFunc);
        //                std::sort(lyt_four.begin(), lyt_four.end(), compareFunc);
        //                std::sort(lyt_five.begin(), lyt_five.end(), compareFunc);
        //                std::sort(lyt_six.begin(), lyt_six.end(), compareFunc);
        //                std::sort(lyt_seven.begin(), lyt_seven.end(), compareFunc);
        //
        //                int                                           number = 100;
        //                std::vector<charge_distribution_surface<Lyt>> lyt_ones(
        //                    lyt_one.begin(), lyt_one.begin() + std::min(number,
        //                    static_cast<int>(lyt_one.size())));
        //                std::vector<charge_distribution_surface<Lyt>> lyt_twos(
        //                    lyt_two.begin(), lyt_two.begin() + std::min(number,
        //                    static_cast<int>(lyt_two.size())));
        //                std::vector<charge_distribution_surface<Lyt>> lyt_threes(
        //                    lyt_three.begin(), lyt_three.begin() + std::min(number,
        //                    static_cast<int>(lyt_three.size())));
        //                std::vector<charge_distribution_surface<Lyt>> lyt_fours(
        //                    lyt_four.begin(), lyt_four.begin() + std::min(number,
        //                    static_cast<int>(lyt_four.size())));
        //
        //                std::vector<charge_distribution_surface<Lyt>> lyt_fives(
        //                    lyt_five.begin(), lyt_five.begin() + std::min(number,
        //                    static_cast<int>(lyt_five.size())));
        //                std::vector<charge_distribution_surface<Lyt>> lyt_sixs(
        //                    lyt_six.begin(), lyt_six.begin() + std::min(number,
        //                    static_cast<int>(lyt_six.size())));
        //                std::vector<charge_distribution_surface<Lyt>> lyt_sevens(
        //                    lyt_seven.begin(), lyt_seven.begin() + std::min(number,
        //                    static_cast<int>(lyt_seven.size())));

        charge_distribution_surface<Lyt> charge_lyt{layout};
        uint64_t                         counter = 0;
        std::vector<double>              valid_energies{};
        double                           energy_threas = 1000;
        for (auto i = 0u; i < lyt_one.size(); i++)
        {
            for (const auto& lyts_two : lyt_two)
            {
                for (auto j = 0u; j < lyt_three.size(); j++)
                {
                    for (const auto& lyts_four : lyt_four)
                    {
                        //                                                uint64_t counter_unmatched = 0;
                        //                                                for (const auto& neighbor_cell :
                        //                                                all_defect_cells[0])
                        //                                                {
                        //                                                    lyts_four.foreach_cell(
                        //                                                        [&counter_unmatched, &neighbor_cell,
                        //                                                        &lyts_four, &i, this](const auto& c1)
                        //                                                        {
                        //                                                            if (c1 == neighbor_cell)
                        //                                                            {
                        //                                                                if
                        //                                                                (charge_state_to_sign(lyts_four.get_charge_state(c1))
                        //                                                                !=
                        //                                                                    get_charge_state_defect(0,
                        //                                                                    i, c1))
                        //                                                                {
                        //                                                                    counter_unmatched += 1;
                        //                                                                }
                        //                                                            }
                        //                                                        });
                        //                                                }
                        //                                                if (counter_unmatched != 0)
                        //                                                {
                        //                                                    continue;
                        //                                                }

                        for (const auto& lyts_five : lyt_five)
                        {
                            uint64_t counter_unmatched = 0;
                            for (const auto& neighbor_cell : all_defect_cells[2])
                            {
                                lyts_five.foreach_cell(
                                    [&counter_unmatched, &neighbor_cell, &lyts_five, &j, this](const auto& c1)
                                    {
                                        if (c1 == neighbor_cell)
                                        {
                                            if (charge_state_to_sign(lyts_five.get_charge_state(c1)) !=
                                                get_charge_state_defect(2, j, c1))
                                            {
                                                counter_unmatched += 1;
                                            }
                                        }
                                    });
                            }
                            if (counter_unmatched != 0)
                            {
                                continue;
                            }

                            for (const auto& lyts_six : lyt_six)
                            {
                                for (const auto& lyts_seven : lyt_seven)
                                {
                                    lyt_one[i].foreach_cell(
                                        [this, &charge_lyt, &lyt_one, &i](const auto& c1) {
                                            charge_lyt.assign_charge_state(c1, lyt_one[i].get_charge_state(c1), false);
                                        });
                                    lyts_two.foreach_cell(
                                        [this, &charge_lyt, &lyts_two](const auto& c1)
                                        { charge_lyt.assign_charge_state(c1, lyts_two.get_charge_state(c1), false); });
                                    lyt_three[j].foreach_cell(
                                        [this, &charge_lyt, &lyt_three, &j](const auto& c1) {
                                            charge_lyt.assign_charge_state(c1, lyt_three[j].get_charge_state(c1),
                                                                           false);
                                        });
                                    lyts_four.foreach_cell(
                                        [this, &charge_lyt, &lyts_four](const auto& c1)
                                        { charge_lyt.assign_charge_state(c1, lyts_four.get_charge_state(c1), false); });
                                    lyts_five.foreach_cell(
                                        [this, &charge_lyt, &lyts_five](const auto& c1)
                                        { charge_lyt.assign_charge_state(c1, lyts_five.get_charge_state(c1), false); });
                                    lyts_six.foreach_cell(
                                        [this, &charge_lyt, &lyts_six](const auto& c1)
                                        { charge_lyt.assign_charge_state(c1, lyts_six.get_charge_state(c1), false); });
                                    lyts_seven.foreach_cell(
                                        [this, &charge_lyt, &lyts_seven](const auto& c1) {
                                            charge_lyt.assign_charge_state(c1, lyts_seven.get_charge_state(c1), false);
                                        });
                                    charge_lyt.update_after_charge_change();
                                    if (charge_lyt.is_physically_valid())
                                    {
                                        if (charge_lyt.get_system_energy() < energy_threas)
                                        {
                                            std::vector<charge_distribution_surface<Lyt>> lyts{};
                                            std::cout << charge_lyt.get_system_energy() << std::endl;

                                            sidb_simulation_result<Lyt> sim_result{};
                                            sim_result.algorithm_name = "ExGS";
                                            charge_distribution_surface<Lyt> charge_lyt_copy{charge_lyt};
                                            lyts.emplace_back(charge_lyt_copy);
                                            sim_result.charge_distributions = lyts;
                                            energy_threas                   = charge_lyt.get_system_energy();
                                            write_sqd_sim_result<Lyt>(sim_result,
                                                                      "/Users/jandrewniok/CLionProjects/"
                                                                      "fiction_fork/experiments/result.xml");
                                        }
                                    }
                                    counter += 1;
                                    if (counter % 1000000 == 0)
                                    {
                                        std::cout << counter << std::endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void combining_all_12()
    {
        auto compareFunc =
            [](const charge_distribution_surface<Lyt>& lyt1, const charge_distribution_surface<Lyt>& lyt2)
        { return lyt1.get_system_energy() < lyt2.get_system_energy(); };

        std::cout << "combining starts: " << std::to_string(lyts_of_regions.size()) << std::endl;
        uint64_t counter_lyts = 1;
        for (const auto& lyts_region : lyts_of_regions)
        {
            counter_lyts *= lyts_region.size();
        }
        std::cout << "number enumerations: " << std::to_string(counter_lyts) << std::endl;

        auto lyt_one   = lyts_of_regions[0];
        auto lyt_two   = lyts_of_regions[1];
        auto lyt_three = lyts_of_regions[2];
        auto lyt_four  = lyts_of_regions[3];
        auto lyt_five  = lyts_of_regions[4];
        auto lyt_six   = lyts_of_regions[5];
        auto lyt_seven = lyts_of_regions[6];
        auto lyt_eight = lyts_of_regions[7];
        auto lyt_nine  = lyts_of_regions[8];
        auto lyt_ten   = lyts_of_regions[9];
        auto lyt_11    = lyts_of_regions[10];
        auto lyt_12    = lyts_of_regions[11];

        //        std::sort(lyt_one.begin(), lyt_one.end(), compareFunc);
        //        std::sort(lyt_two.begin(), lyt_two.end(), compareFunc);
        //        std::sort(lyt_three.begin(), lyt_three.end(), compareFunc);
        //        std::sort(lyt_four.begin(), lyt_four.end(), compareFunc);
        //        std::sort(lyt_five.begin(), lyt_five.end(), compareFunc);
        //        std::sort(lyt_six.begin(), lyt_six.end(), compareFunc);
        //        std::sort(lyt_seven.begin(), lyt_seven.end(), compareFunc);
        //        std::sort(lyt_eight.begin(), lyt_eight.end(), compareFunc);
        //        std::sort(lyt_nine.begin(), lyt_nine.end(), compareFunc);
        //        std::sort(lyt_ten.begin(), lyt_ten.end(), compareFunc);
        //        std::sort(lyt_11.begin(), lyt_11.end(), compareFunc);
        //        std::sort(lyt_12.begin(), lyt_12.end(), compareFunc);
        //        std::sort(lyt_13.begin(), lyt_13.end(), compareFunc);
        //
        //        int                                           number = 1000;
        //        std::vector<charge_distribution_surface<Lyt>> lyt_ones(
        //            lyt_one.begin(), lyt_one.begin() + std::min(number, static_cast<int>(lyt_one.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_twos(
        //            lyt_two.begin(), lyt_two.begin() + std::min(number, static_cast<int>(lyt_two.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_threes(
        //            lyt_three.begin(), lyt_three.begin() + std::min(number, static_cast<int>(lyt_three.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_fours(
        //            lyt_four.begin(), lyt_four.begin() + std::min(number, static_cast<int>(lyt_four.size())));
        //
        //        std::vector<charge_distribution_surface<Lyt>> lyt_fives(
        //            lyt_five.begin(), lyt_five.begin() + std::min(number, static_cast<int>(lyt_five.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_sixs(
        //            lyt_six.begin(), lyt_six.begin() + std::min(number, static_cast<int>(lyt_six.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_sevens(
        //            lyt_seven.begin(), lyt_seven.begin() + std::min(number, static_cast<int>(lyt_seven.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_eights(
        //            lyt_eight.begin(), lyt_eight.begin() + std::min(number, static_cast<int>(lyt_eight.size())));
        //
        //        std::vector<charge_distribution_surface<Lyt>> lyt_nines(
        //            lyt_nine.begin(), lyt_nine.begin() + std::min(number, static_cast<int>(lyt_nine.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_tens(
        //            lyt_ten.begin(), lyt_ten.begin() + std::min(number, static_cast<int>(lyt_ten.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_11s(
        //            lyt_11.begin(), lyt_11.begin() + std::min(number, static_cast<int>(lyt_11.size())));
        //
        //        std::vector<charge_distribution_surface<Lyt>> lyt_12s(
        //            lyt_12.begin(), lyt_12.begin() + std::min(number, static_cast<int>(lyt_12.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_13s(
        //            lyt_13.begin(), lyt_13.begin() + std::min(number, static_cast<int>(lyt_13.size())));

        charge_distribution_surface<Lyt> charge_lyt{layout};
        uint64_t                         counter = 0;
        std::vector<double>              valid_energies{};
        double                           energy_threas = 1000;
        for (const auto& lyts_one : lyt_one)
        {
            for (const auto& lyts_two : lyt_two)
            {
                for (const auto& lyts_three : lyt_three)
                {
                    for (const auto& lyts_four : lyt_four)
                    {
                        for (const auto& lyts_five : lyt_five)
                        {
                            for (const auto& lyts_six : lyt_six)
                            {
                                for (const auto& lyts_seven : lyt_seven)
                                {
                                    for (const auto& lyts_eight : lyt_eight)
                                    {
                                        for (const auto& lyts_nine : lyt_nine)
                                        {
                                            for (const auto& lyts_ten : lyt_ten)
                                            {
                                                for (const auto& lyts_11 : lyt_11)
                                                {
                                                    for (const auto& lyts_12 : lyt_12)
                                                    {

                                                        lyts_one.foreach_cell(
                                                            [this, &charge_lyt, &lyts_one](const auto& c1) {
                                                                charge_lyt.assign_charge_state(
                                                                    c1, lyts_one.get_charge_state(c1), false);
                                                            });
                                                        lyts_two.foreach_cell(
                                                            [this, &charge_lyt, &lyts_two](const auto& c1) {
                                                                charge_lyt.assign_charge_state(
                                                                    c1, lyts_two.get_charge_state(c1), false);
                                                            });
                                                        lyts_three.foreach_cell(
                                                            [this, &charge_lyt, &lyts_three](const auto& c1) {
                                                                charge_lyt.assign_charge_state(
                                                                    c1, lyts_three.get_charge_state(c1), false);
                                                            });
                                                        lyts_four.foreach_cell(
                                                            [this, &charge_lyt, &lyts_four](const auto& c1) {
                                                                charge_lyt.assign_charge_state(
                                                                    c1, lyts_four.get_charge_state(c1), false);
                                                            });
                                                        lyts_five.foreach_cell(
                                                            [this, &charge_lyt, &lyts_five](const auto& c1) {
                                                                charge_lyt.assign_charge_state(
                                                                    c1, lyts_five.get_charge_state(c1), false);
                                                            });
                                                        lyts_six.foreach_cell(
                                                            [this, &charge_lyt, &lyts_six](const auto& c1) {
                                                                charge_lyt.assign_charge_state(
                                                                    c1, lyts_six.get_charge_state(c1), false);
                                                            });
                                                        lyts_seven.foreach_cell(
                                                            [this, &charge_lyt, &lyts_seven](const auto& c1) {
                                                                charge_lyt.assign_charge_state(
                                                                    c1, lyts_seven.get_charge_state(c1), false);
                                                            });
                                                        lyts_eight.foreach_cell(
                                                            [this, &charge_lyt, &lyts_eight](const auto& c1) {
                                                                charge_lyt.assign_charge_state(
                                                                    c1, lyts_eight.get_charge_state(c1), false);
                                                            });
                                                        lyts_nine.foreach_cell(
                                                            [this, &charge_lyt, &lyts_nine](const auto& c1) {
                                                                charge_lyt.assign_charge_state(
                                                                    c1, lyts_nine.get_charge_state(c1), false);
                                                            });
                                                        lyts_ten.foreach_cell(
                                                            [this, &charge_lyt, &lyts_ten](const auto& c1) {
                                                                charge_lyt.assign_charge_state(
                                                                    c1, lyts_ten.get_charge_state(c1), false);
                                                            });
                                                        lyts_11.foreach_cell(
                                                            [this, &charge_lyt, &lyts_11](const auto& c1) {
                                                                charge_lyt.assign_charge_state(
                                                                    c1, lyts_11.get_charge_state(c1), false);
                                                            });
                                                        lyts_12.foreach_cell(
                                                            [this, &charge_lyt, &lyts_12](const auto& c1) {
                                                                charge_lyt.assign_charge_state(
                                                                    c1, lyts_12.get_charge_state(c1), false);
                                                            });
                                                        charge_lyt.update_after_charge_change();
                                                        if (charge_lyt.is_physically_valid())
                                                        {
                                                            if (charge_lyt.get_system_energy() < energy_threas)
                                                            {
                                                                std::vector<charge_distribution_surface<Lyt>> lyts{};
                                                                std::cout << charge_lyt.get_system_energy()
                                                                          << std::endl;

                                                                sidb_simulation_result<Lyt> sim_result{};
                                                                sim_result.algorithm_name = "ExGS";
                                                                charge_distribution_surface<Lyt> charge_lyt_copy{
                                                                    charge_lyt};
                                                                lyts.emplace_back(charge_lyt_copy);
                                                                sim_result.charge_distributions = lyts;
                                                                energy_threas = charge_lyt.get_system_energy();
                                                                write_sqd_sim_result<Lyt>(
                                                                    sim_result, "/Users/jandrewniok/CLionProjects/"
                                                                                "fiction_fork/experiments/"
                                                                                "result.xml");
                                                            }
                                                        }
                                                        counter += 1;
                                                        if (counter % 1000000 == 0)
                                                        {
                                                            std::cout << counter << std::endl;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void combining_all_16_old()
    {
        auto compareFunc =
            [](const charge_distribution_surface<Lyt>& lyt1, const charge_distribution_surface<Lyt>& lyt2)
        { return lyt1.get_system_energy() < lyt2.get_system_energy(); };

        std::cout << "combining starts: " << std::to_string(lyts_of_regions.size()) << std::endl;
        uint64_t counter_lyts = 1;
        for (const auto& lyts_region : lyts_of_regions)
        {
            counter_lyts *= lyts_region.size();
        }
        std::cout << "number enumerations: " << std::to_string(counter_lyts) << std::endl;

        auto lyt_one   = lyts_of_regions[0];
        auto lyt_two   = lyts_of_regions[1];
        auto lyt_three = lyts_of_regions[2];
        auto lyt_four  = lyts_of_regions[3];
        auto lyt_five  = lyts_of_regions[4];
        auto lyt_six   = lyts_of_regions[5];
        auto lyt_seven = lyts_of_regions[6];
        auto lyt_eight = lyts_of_regions[7];
        auto lyt_nine  = lyts_of_regions[8];
        auto lyt_ten   = lyts_of_regions[9];
        auto lyt_11    = lyts_of_regions[10];
        auto lyt_12    = lyts_of_regions[11];
        auto lyt_13    = lyts_of_regions[12];
        auto lyt_14    = lyts_of_regions[13];
        auto lyt_15    = lyts_of_regions[14];
        auto lyt_16    = lyts_of_regions[15];

        //        std::sort(lyt_one.begin(), lyt_one.end(), compareFunc);
        //        std::sort(lyt_two.begin(), lyt_two.end(), compareFunc);
        //        std::sort(lyt_three.begin(), lyt_three.end(), compareFunc);
        //        std::sort(lyt_four.begin(), lyt_four.end(), compareFunc);
        //        std::sort(lyt_five.begin(), lyt_five.end(), compareFunc);
        //        std::sort(lyt_six.begin(), lyt_six.end(), compareFunc);
        //        std::sort(lyt_seven.begin(), lyt_seven.end(), compareFunc);
        //        std::sort(lyt_eight.begin(), lyt_eight.end(), compareFunc);
        //        std::sort(lyt_nine.begin(), lyt_nine.end(), compareFunc);
        //        std::sort(lyt_ten.begin(), lyt_ten.end(), compareFunc);
        //        std::sort(lyt_11.begin(), lyt_11.end(), compareFunc);
        //        std::sort(lyt_12.begin(), lyt_12.end(), compareFunc);
        //        std::sort(lyt_13.begin(), lyt_13.end(), compareFunc);
        //
        //        int                                           number = 1000;
        //        std::vector<charge_distribution_surface<Lyt>> lyt_ones(
        //            lyt_one.begin(), lyt_one.begin() + std::min(number, static_cast<int>(lyt_one.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_twos(
        //            lyt_two.begin(), lyt_two.begin() + std::min(number, static_cast<int>(lyt_two.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_threes(
        //            lyt_three.begin(), lyt_three.begin() + std::min(number, static_cast<int>(lyt_three.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_fours(
        //            lyt_four.begin(), lyt_four.begin() + std::min(number, static_cast<int>(lyt_four.size())));
        //
        //        std::vector<charge_distribution_surface<Lyt>> lyt_fives(
        //            lyt_five.begin(), lyt_five.begin() + std::min(number, static_cast<int>(lyt_five.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_sixs(
        //            lyt_six.begin(), lyt_six.begin() + std::min(number, static_cast<int>(lyt_six.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_sevens(
        //            lyt_seven.begin(), lyt_seven.begin() + std::min(number, static_cast<int>(lyt_seven.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_eights(
        //            lyt_eight.begin(), lyt_eight.begin() + std::min(number, static_cast<int>(lyt_eight.size())));
        //
        //        std::vector<charge_distribution_surface<Lyt>> lyt_nines(
        //            lyt_nine.begin(), lyt_nine.begin() + std::min(number, static_cast<int>(lyt_nine.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_tens(
        //            lyt_ten.begin(), lyt_ten.begin() + std::min(number, static_cast<int>(lyt_ten.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_11s(
        //            lyt_11.begin(), lyt_11.begin() + std::min(number, static_cast<int>(lyt_11.size())));
        //
        //        std::vector<charge_distribution_surface<Lyt>> lyt_12s(
        //            lyt_12.begin(), lyt_12.begin() + std::min(number, static_cast<int>(lyt_12.size())));
        //        std::vector<charge_distribution_surface<Lyt>> lyt_13s(
        //            lyt_13.begin(), lyt_13.begin() + std::min(number, static_cast<int>(lyt_13.size())));

        charge_distribution_surface<Lyt> charge_lyt{layout};
        uint64_t                         counter = 0;
        std::vector<double>              valid_energies{};
        double                           energy_threas = 1000;
        for (const auto& lyts_one : lyt_one)
        {
            for (const auto& lyts_two : lyt_two)
            {
                for (const auto& lyts_three : lyt_three)
                {
                    for (const auto& lyts_four : lyt_four)
                    {
                        for (const auto& lyts_five : lyt_five)
                        {
                            for (const auto& lyts_six : lyt_six)
                            {
                                for (const auto& lyts_seven : lyt_seven)
                                {
                                    for (const auto& lyts_eight : lyt_eight)
                                    {
                                        for (const auto& lyts_nine : lyt_nine)
                                        {
                                            for (const auto& lyts_ten : lyt_ten)
                                            {
                                                for (const auto& lyts_11 : lyt_11)
                                                {
                                                    for (const auto& lyts_12 : lyt_12)
                                                    {
                                                        for (const auto& lyts_13 : lyt_13)
                                                        {
                                                            for (const auto& lyts_14 : lyt_14)
                                                            {
                                                                for (const auto& lyts_15 : lyt_15)
                                                                {
                                                                    for (const auto& lyts_16 : lyt_16)
                                                                    {
                                                                        lyts_one.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_one](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_one.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_two.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_two](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_two.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_three.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_three](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_three.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_four.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_four](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_four.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_five.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_five](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_five.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_six.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_six](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_six.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_seven.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_seven](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_seven.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_eight.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_eight](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_eight.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_nine.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_nine](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_nine.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_ten.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_ten](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_ten.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_11.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_11](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_11.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_12.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_12](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_12.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_13.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_13](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_13.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_14.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_14](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_14.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_15.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_15](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_15.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_16.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_16](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_16.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        charge_lyt.update_after_charge_change();
                                                                        if (charge_lyt.is_physically_valid())
                                                                        {
                                                                            if (charge_lyt.get_system_energy() <
                                                                                energy_threas)
                                                                            {
                                                                                std::vector<
                                                                                    charge_distribution_surface<Lyt>>
                                                                                    lyts{};
                                                                                std::cout
                                                                                    << charge_lyt.get_system_energy()
                                                                                    << std::endl;

                                                                                sidb_simulation_result<Lyt>
                                                                                    sim_result{};
                                                                                sim_result.algorithm_name = "ExGS";
                                                                                charge_distribution_surface<Lyt>
                                                                                    charge_lyt_copy{charge_lyt};
                                                                                lyts.emplace_back(charge_lyt_copy);
                                                                                sim_result.charge_distributions = lyts;
                                                                                energy_threas =
                                                                                    charge_lyt.get_system_energy();
                                                                                write_sqd_sim_result<Lyt>(
                                                                                    sim_result,
                                                                                    "/Users/jandrewniok/CLionProjects/"
                                                                                    "fiction_fork/experiments/"
                                                                                    "result.xml");
                                                                            }
                                                                        }
                                                                        counter += 1;
                                                                        if (counter % 1000000 == 0)
                                                                        {
                                                                            std::cout << counter << std::endl;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void combining_all_13()
    {
        auto compareFunc =
            [](const charge_distribution_surface<Lyt>& lyt1, const charge_distribution_surface<Lyt>& lyt2)
        { return lyt1.get_system_energy() < lyt2.get_system_energy(); };

        std::cout << "combining starts: " << std::to_string(lyts_of_regions.size()) << std::endl;
        uint64_t counter_lyts = 1;
        for (const auto& lyts_region : lyts_of_regions)
        {
            counter_lyts *= lyts_region.size();
        }
        std::cout << "number enumerations: " << std::to_string(counter_lyts) << std::endl;

        auto lyt_one   = lyts_of_regions[0];
        auto lyt_two   = lyts_of_regions[1];
        auto lyt_three = lyts_of_regions[2];
        auto lyt_four  = lyts_of_regions[3];
        auto lyt_five  = lyts_of_regions[4];
        auto lyt_six   = lyts_of_regions[5];
        auto lyt_seven = lyts_of_regions[6];
        auto lyt_eight = lyts_of_regions[7];
        auto lyt_nine  = lyts_of_regions[8];
        auto lyt_ten   = lyts_of_regions[9];
        auto lyt_11    = lyts_of_regions[10];
        auto lyt_12    = lyts_of_regions[11];
        auto lyt_13    = lyts_of_regions[12];

        std::sort(lyt_one.begin(), lyt_one.end(), compareFunc);
        std::sort(lyt_two.begin(), lyt_two.end(), compareFunc);
        std::sort(lyt_three.begin(), lyt_three.end(), compareFunc);
        std::sort(lyt_four.begin(), lyt_four.end(), compareFunc);
        std::sort(lyt_five.begin(), lyt_five.end(), compareFunc);
        std::sort(lyt_six.begin(), lyt_six.end(), compareFunc);
        std::sort(lyt_seven.begin(), lyt_seven.end(), compareFunc);
        std::sort(lyt_eight.begin(), lyt_eight.end(), compareFunc);
        std::sort(lyt_nine.begin(), lyt_nine.end(), compareFunc);
        std::sort(lyt_ten.begin(), lyt_ten.end(), compareFunc);
        std::sort(lyt_11.begin(), lyt_11.end(), compareFunc);
        std::sort(lyt_12.begin(), lyt_12.end(), compareFunc);
        std::sort(lyt_13.begin(), lyt_13.end(), compareFunc);

        int                                           number = 1000;
        std::vector<charge_distribution_surface<Lyt>> lyt_ones(
            lyt_one.begin(), lyt_one.begin() + std::min(number, static_cast<int>(lyt_one.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_twos(
            lyt_two.begin(), lyt_two.begin() + std::min(number, static_cast<int>(lyt_two.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_threes(
            lyt_three.begin(), lyt_three.begin() + std::min(number, static_cast<int>(lyt_three.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_fours(
            lyt_four.begin(), lyt_four.begin() + std::min(number, static_cast<int>(lyt_four.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_fives(
            lyt_five.begin(), lyt_five.begin() + std::min(number, static_cast<int>(lyt_five.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_sixs(
            lyt_six.begin(), lyt_six.begin() + std::min(number, static_cast<int>(lyt_six.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_sevens(
            lyt_seven.begin(), lyt_seven.begin() + std::min(number, static_cast<int>(lyt_seven.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_eights(
            lyt_eight.begin(), lyt_eight.begin() + std::min(number, static_cast<int>(lyt_eight.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_nines(
            lyt_nine.begin(), lyt_nine.begin() + std::min(number, static_cast<int>(lyt_nine.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_tens(
            lyt_ten.begin(), lyt_ten.begin() + std::min(number, static_cast<int>(lyt_ten.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_11s(
            lyt_11.begin(), lyt_11.begin() + std::min(number, static_cast<int>(lyt_11.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_12s(
            lyt_12.begin(), lyt_12.begin() + std::min(number, static_cast<int>(lyt_12.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_13s(
            lyt_13.begin(), lyt_13.begin() + std::min(number, static_cast<int>(lyt_13.size())));

        charge_distribution_surface<Lyt> charge_lyt{layout};
        uint64_t                         counter = 0;
        std::vector<double>              valid_energies{};
        double                           energy_threas = 1000;
        for (const auto& lyts_one : lyt_ones)
        {
            for (const auto& lyts_two : lyt_twos)
            {
                for (const auto& lyts_three : lyt_threes)
                {
                    for (const auto& lyts_four : lyt_fours)
                    {
                        for (const auto& lyts_five : lyt_fives)
                        {
                            for (const auto& lyts_six : lyt_sixs)
                            {
                                for (const auto& lyts_seven : lyt_sevens)
                                {
                                    for (const auto& lyts_eight : lyt_eights)
                                    {
                                        for (const auto& lyts_nine : lyt_nines)
                                        {
                                            for (const auto& lyts_ten : lyt_tens)
                                            {
                                                for (const auto& lyts_11 : lyt_11s)
                                                {
                                                    for (const auto& lyts_12 : lyt_12s)
                                                    {
                                                        for (const auto& lyts_13 : lyt_13s)
                                                        {
                                                            lyts_one.foreach_cell(
                                                                [this, &charge_lyt, &lyts_one](const auto& c1) {
                                                                    charge_lyt.assign_charge_state(
                                                                        c1, lyts_one.get_charge_state(c1), false);
                                                                });
                                                            lyts_two.foreach_cell(
                                                                [this, &charge_lyt, &lyts_two](const auto& c1) {
                                                                    charge_lyt.assign_charge_state(
                                                                        c1, lyts_two.get_charge_state(c1), false);
                                                                });
                                                            lyts_three.foreach_cell(
                                                                [this, &charge_lyt, &lyts_three](const auto& c1) {
                                                                    charge_lyt.assign_charge_state(
                                                                        c1, lyts_three.get_charge_state(c1), false);
                                                                });
                                                            lyts_four.foreach_cell(
                                                                [this, &charge_lyt, &lyts_four](const auto& c1) {
                                                                    charge_lyt.assign_charge_state(
                                                                        c1, lyts_four.get_charge_state(c1), false);
                                                                });
                                                            lyts_five.foreach_cell(
                                                                [this, &charge_lyt, &lyts_five](const auto& c1) {
                                                                    charge_lyt.assign_charge_state(
                                                                        c1, lyts_five.get_charge_state(c1), false);
                                                                });
                                                            lyts_six.foreach_cell(
                                                                [this, &charge_lyt, &lyts_six](const auto& c1) {
                                                                    charge_lyt.assign_charge_state(
                                                                        c1, lyts_six.get_charge_state(c1), false);
                                                                });
                                                            lyts_seven.foreach_cell(
                                                                [this, &charge_lyt, &lyts_seven](const auto& c1) {
                                                                    charge_lyt.assign_charge_state(
                                                                        c1, lyts_seven.get_charge_state(c1), false);
                                                                });
                                                            lyts_eight.foreach_cell(
                                                                [this, &charge_lyt, &lyts_eight](const auto& c1) {
                                                                    charge_lyt.assign_charge_state(
                                                                        c1, lyts_eight.get_charge_state(c1), false);
                                                                });
                                                            lyts_nine.foreach_cell(
                                                                [this, &charge_lyt, &lyts_nine](const auto& c1) {
                                                                    charge_lyt.assign_charge_state(
                                                                        c1, lyts_nine.get_charge_state(c1), false);
                                                                });
                                                            lyts_ten.foreach_cell(
                                                                [this, &charge_lyt, &lyts_ten](const auto& c1) {
                                                                    charge_lyt.assign_charge_state(
                                                                        c1, lyts_ten.get_charge_state(c1), false);
                                                                });
                                                            lyts_11.foreach_cell(
                                                                [this, &charge_lyt, &lyts_11](const auto& c1) {
                                                                    charge_lyt.assign_charge_state(
                                                                        c1, lyts_11.get_charge_state(c1), false);
                                                                });
                                                            lyts_12.foreach_cell(
                                                                [this, &charge_lyt, &lyts_12](const auto& c1) {
                                                                    charge_lyt.assign_charge_state(
                                                                        c1, lyts_12.get_charge_state(c1), false);
                                                                });
                                                            lyts_13.foreach_cell(
                                                                [this, &charge_lyt, &lyts_13](const auto& c1) {
                                                                    charge_lyt.assign_charge_state(
                                                                        c1, lyts_13.get_charge_state(c1), false);
                                                                });
                                                            charge_lyt.update_after_charge_change();
                                                            if (charge_lyt.is_physically_valid())
                                                            {
                                                                if (charge_lyt.get_system_energy() < energy_threas)
                                                                {
                                                                    std::vector<charge_distribution_surface<Lyt>>
                                                                        lyts{};
                                                                    std::cout << charge_lyt.get_system_energy()
                                                                              << std::endl;

                                                                    sidb_simulation_result<Lyt> sim_result{};
                                                                    sim_result.algorithm_name = "ExGS";
                                                                    charge_distribution_surface<Lyt> charge_lyt_copy{
                                                                        charge_lyt};
                                                                    lyts.emplace_back(charge_lyt_copy);
                                                                    sim_result.charge_distributions = lyts;
                                                                    energy_threas = charge_lyt.get_system_energy();
                                                                    write_sqd_sim_result<Lyt>(
                                                                        sim_result, "/Users/jandrewniok/CLionProjects/"
                                                                                    "fiction_fork/experiments/"
                                                                                    "result.xml");
                                                                }
                                                            }
                                                            counter += 1;
                                                            if (counter % 1000000 == 0)
                                                            {
                                                                std::cout << counter << std::endl;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void combining_all_11()
    {
        auto compareFunc =
            [](const charge_distribution_surface<Lyt>& lyt1, const charge_distribution_surface<Lyt>& lyt2)
        { return lyt1.get_system_energy() < lyt2.get_system_energy(); };

        std::cout << "combining starts: " << std::to_string(lyts_of_regions.size()) << std::endl;
        uint64_t counter_lyts = 1;
        for (const auto& lyts_region : lyts_of_regions)
        {
            counter_lyts *= lyts_region.size();
        }
        std::cout << "number enumerations: " << std::to_string(counter_lyts) << std::endl;

        auto lyt_one   = lyts_of_regions[0];
        auto lyt_two   = lyts_of_regions[1];
        auto lyt_three = lyts_of_regions[2];
        auto lyt_four  = lyts_of_regions[3];
        auto lyt_five  = lyts_of_regions[4];
        auto lyt_six   = lyts_of_regions[5];
        auto lyt_seven = lyts_of_regions[6];
        auto lyt_eight = lyts_of_regions[7];
        auto lyt_nine  = lyts_of_regions[8];
        auto lyt_ten   = lyts_of_regions[9];
        auto lyt_11    = lyts_of_regions[10];

        std::sort(lyt_one.begin(), lyt_one.end(), compareFunc);
        std::sort(lyt_two.begin(), lyt_two.end(), compareFunc);
        std::sort(lyt_three.begin(), lyt_three.end(), compareFunc);
        std::sort(lyt_four.begin(), lyt_four.end(), compareFunc);
        std::sort(lyt_five.begin(), lyt_five.end(), compareFunc);
        std::sort(lyt_six.begin(), lyt_six.end(), compareFunc);
        std::sort(lyt_seven.begin(), lyt_seven.end(), compareFunc);
        std::sort(lyt_eight.begin(), lyt_eight.end(), compareFunc);
        std::sort(lyt_nine.begin(), lyt_nine.end(), compareFunc);
        std::sort(lyt_ten.begin(), lyt_ten.end(), compareFunc);
        std::sort(lyt_11.begin(), lyt_11.end(), compareFunc);

        int                                           number = 1000;
        std::vector<charge_distribution_surface<Lyt>> lyt_ones(
            lyt_one.begin(), lyt_one.begin() + std::min(number, static_cast<int>(lyt_one.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_twos(
            lyt_two.begin(), lyt_two.begin() + std::min(number, static_cast<int>(lyt_two.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_threes(
            lyt_three.begin(), lyt_three.begin() + std::min(number, static_cast<int>(lyt_three.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_fours(
            lyt_four.begin(), lyt_four.begin() + std::min(number, static_cast<int>(lyt_four.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_fives(
            lyt_five.begin(), lyt_five.begin() + std::min(number, static_cast<int>(lyt_five.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_sixs(
            lyt_six.begin(), lyt_six.begin() + std::min(number, static_cast<int>(lyt_six.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_sevens(
            lyt_seven.begin(), lyt_seven.begin() + std::min(number, static_cast<int>(lyt_seven.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_eights(
            lyt_eight.begin(), lyt_eight.begin() + std::min(number, static_cast<int>(lyt_eight.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_nines(
            lyt_nine.begin(), lyt_nine.begin() + std::min(number, static_cast<int>(lyt_nine.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_tens(
            lyt_ten.begin(), lyt_ten.begin() + std::min(number, static_cast<int>(lyt_ten.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_11s(
            lyt_11.begin(), lyt_11.begin() + std::min(number, static_cast<int>(lyt_11.size())));

        charge_distribution_surface<Lyt> charge_lyt{layout};
        uint64_t                         counter = 0;
        std::vector<double>              valid_energies{};
        double                           energy_threas = 1000;
        for (const auto& lyts_one : lyt_ones)
        {
            for (const auto& lyts_two : lyt_twos)
            {
                for (const auto& lyts_three : lyt_threes)
                {
                    for (const auto& lyts_four : lyt_fours)
                    {
                        for (const auto& lyts_five : lyt_fives)
                        {
                            for (const auto& lyts_six : lyt_sixs)
                            {
                                for (const auto& lyts_seven : lyt_sevens)
                                {
                                    for (const auto& lyts_eight : lyt_eights)
                                    {
                                        for (const auto& lyts_nine : lyt_nines)
                                        {
                                            for (const auto& lyts_ten : lyt_tens)
                                            {
                                                for (const auto& lyts_11 : lyt_11s)
                                                {
                                                    lyts_one.foreach_cell(
                                                        [this, &charge_lyt, &lyts_one](const auto& c1) {
                                                            charge_lyt.assign_charge_state(
                                                                c1, lyts_one.get_charge_state(c1), false);
                                                        });
                                                    lyts_two.foreach_cell(
                                                        [this, &charge_lyt, &lyts_two](const auto& c1) {
                                                            charge_lyt.assign_charge_state(
                                                                c1, lyts_two.get_charge_state(c1), false);
                                                        });
                                                    lyts_three.foreach_cell(
                                                        [this, &charge_lyt, &lyts_three](const auto& c1) {
                                                            charge_lyt.assign_charge_state(
                                                                c1, lyts_three.get_charge_state(c1), false);
                                                        });
                                                    lyts_four.foreach_cell(
                                                        [this, &charge_lyt, &lyts_four](const auto& c1) {
                                                            charge_lyt.assign_charge_state(
                                                                c1, lyts_four.get_charge_state(c1), false);
                                                        });
                                                    lyts_five.foreach_cell(
                                                        [this, &charge_lyt, &lyts_five](const auto& c1) {
                                                            charge_lyt.assign_charge_state(
                                                                c1, lyts_five.get_charge_state(c1), false);
                                                        });
                                                    lyts_six.foreach_cell(
                                                        [this, &charge_lyt, &lyts_six](const auto& c1) {
                                                            charge_lyt.assign_charge_state(
                                                                c1, lyts_six.get_charge_state(c1), false);
                                                        });
                                                    lyts_seven.foreach_cell(
                                                        [this, &charge_lyt, &lyts_seven](const auto& c1) {
                                                            charge_lyt.assign_charge_state(
                                                                c1, lyts_seven.get_charge_state(c1), false);
                                                        });
                                                    lyts_eight.foreach_cell(
                                                        [this, &charge_lyt, &lyts_eight](const auto& c1) {
                                                            charge_lyt.assign_charge_state(
                                                                c1, lyts_eight.get_charge_state(c1), false);
                                                        });
                                                    lyts_nine.foreach_cell(
                                                        [this, &charge_lyt, &lyts_nine](const auto& c1) {
                                                            charge_lyt.assign_charge_state(
                                                                c1, lyts_nine.get_charge_state(c1), false);
                                                        });
                                                    lyts_ten.foreach_cell(
                                                        [this, &charge_lyt, &lyts_ten](const auto& c1) {
                                                            charge_lyt.assign_charge_state(
                                                                c1, lyts_ten.get_charge_state(c1), false);
                                                        });
                                                    lyts_11.foreach_cell(
                                                        [this, &charge_lyt, &lyts_11](const auto& c1) {
                                                            charge_lyt.assign_charge_state(
                                                                c1, lyts_11.get_charge_state(c1), false);
                                                        });
                                                    charge_lyt.update_after_charge_change();
                                                    if (charge_lyt.is_physically_valid())
                                                    {
                                                        std::cout << charge_lyt.get_system_energy() << std::endl;
                                                        if (charge_lyt.get_system_energy() < energy_threas)
                                                        {
                                                            std::vector<charge_distribution_surface<Lyt>> lyts{};
                                                            std::cout << charge_lyt.get_system_energy() << std::endl;

                                                            sidb_simulation_result<Lyt> sim_result{};
                                                            sim_result.algorithm_name = "ExGS";
                                                            charge_distribution_surface<Lyt> charge_lyt_copy{
                                                                charge_lyt};
                                                            lyts.emplace_back(charge_lyt_copy);
                                                            sim_result.charge_distributions = lyts;
                                                            energy_threas = charge_lyt.get_system_energy();
                                                            write_sqd_sim_result<Lyt>(
                                                                sim_result, "/Users/jandrewniok/CLionProjects/"
                                                                            "fiction_fork/experiments/"
                                                                            "result.xml");
                                                        }
                                                    }
                                                    counter += 1;
                                                    if (counter % 1000000 == 0)
                                                    {
                                                        std::cout << counter << std::endl;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void combining_all_less()
    {
        auto compareFunc =
            [](const charge_distribution_surface<Lyt>& lyt1, const charge_distribution_surface<Lyt>& lyt2)
        { return lyt1.get_system_energy() < lyt2.get_system_energy(); };

        std::cout << "combining starts: " << std::to_string(lyts_of_regions.size()) << std::endl;
        uint64_t counter_lyts = 1;
        for (const auto& lyts_region : lyts_of_regions)
        {
            counter_lyts *= lyts_region.size();
        }
        std::cout << "number enumerations: " << std::to_string(counter_lyts) << std::endl;

        auto lyt_one   = lyts_of_regions[0];
        auto lyt_two   = lyts_of_regions[1];
        auto lyt_three = lyts_of_regions[2];
        auto lyt_four  = lyts_of_regions[3];
        auto lyt_five  = lyts_of_regions[4];
        auto lyt_six   = lyts_of_regions[5];
        auto lyt_seven = lyts_of_regions[6];
        auto lyt_eight = lyts_of_regions[7];
        auto lyt_nine  = lyts_of_regions[8];
        auto lyt_ten   = lyts_of_regions[9];
        auto lyt_11    = lyts_of_regions[10];
        auto lyt_12    = lyts_of_regions[11];
        auto lyt_13    = lyts_of_regions[12];
        auto lyt_14    = lyts_of_regions[13];
        auto lyt_15    = lyts_of_regions[14];
        auto lyt_16    = lyts_of_regions[15];

        std::sort(lyt_one.begin(), lyt_one.end(), compareFunc);
        std::sort(lyt_two.begin(), lyt_two.end(), compareFunc);
        std::sort(lyt_three.begin(), lyt_three.end(), compareFunc);
        std::sort(lyt_four.begin(), lyt_four.end(), compareFunc);
        std::sort(lyt_five.begin(), lyt_five.end(), compareFunc);
        std::sort(lyt_six.begin(), lyt_six.end(), compareFunc);
        std::sort(lyt_seven.begin(), lyt_seven.end(), compareFunc);
        std::sort(lyt_eight.begin(), lyt_eight.end(), compareFunc);
        std::sort(lyt_nine.begin(), lyt_nine.end(), compareFunc);
        std::sort(lyt_ten.begin(), lyt_ten.end(), compareFunc);
        std::sort(lyt_11.begin(), lyt_11.end(), compareFunc);
        std::sort(lyt_12.begin(), lyt_12.end(), compareFunc);
        std::sort(lyt_13.begin(), lyt_13.end(), compareFunc);
        std::sort(lyt_14.begin(), lyt_14.end(), compareFunc);
        std::sort(lyt_15.begin(), lyt_15.end(), compareFunc);
        std::sort(lyt_16.begin(), lyt_16.end(), compareFunc);

        int                                           number = 5;
        std::vector<charge_distribution_surface<Lyt>> lyt_ones(
            lyt_one.begin(), lyt_one.begin() + std::min(number, static_cast<int>(lyt_one.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_twos(
            lyt_two.begin(), lyt_two.begin() + std::min(number, static_cast<int>(lyt_two.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_threes(
            lyt_three.begin(), lyt_three.begin() + std::min(number, static_cast<int>(lyt_three.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_fours(
            lyt_four.begin(), lyt_four.begin() + std::min(number, static_cast<int>(lyt_four.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_fives(
            lyt_five.begin(), lyt_five.begin() + std::min(number, static_cast<int>(lyt_five.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_sixs(
            lyt_six.begin(), lyt_six.begin() + std::min(number, static_cast<int>(lyt_six.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_sevens(
            lyt_seven.begin(), lyt_seven.begin() + std::min(number, static_cast<int>(lyt_seven.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_eights(
            lyt_eight.begin(), lyt_eight.begin() + std::min(number, static_cast<int>(lyt_eight.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_nines(
            lyt_nine.begin(), lyt_nine.begin() + std::min(number, static_cast<int>(lyt_nine.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_tens(
            lyt_ten.begin(), lyt_ten.begin() + std::min(number, static_cast<int>(lyt_ten.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_11s(
            lyt_11.begin(), lyt_11.begin() + std::min(number, static_cast<int>(lyt_11.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_12s(
            lyt_12.begin(), lyt_12.begin() + std::min(number, static_cast<int>(lyt_12.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_13s(
            lyt_13.begin(), lyt_13.begin() + std::min(number, static_cast<int>(lyt_13.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_14s(
            lyt_14.begin(), lyt_14.begin() + std::min(number, static_cast<int>(lyt_14.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_15s(
            lyt_15.begin(), lyt_15.begin() + std::min(number, static_cast<int>(lyt_15.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_16s(
            lyt_16.begin(), lyt_16.begin() + std::min(number, static_cast<int>(lyt_16.size())));

        charge_distribution_surface<Lyt> charge_lyt{layout};
        uint64_t                         counter = 0;
        std::vector<double>              valid_energies{};
        double                           energy_threas = 1000;
        for (auto i = 0u; i < lyt_one.size(); i++)
        {
            for (auto j = 0u; j < lyt_two.size(); j++)
            {
                for (auto three = 0u; three < lyt_three.size(); three++)
                {
                    for (auto four = 0u; four < lyt_four.size(); four++)
                    {
                        for (auto five = 0u; five < lyt_five.size(); five++)
                        {
                            for (const auto& lyts_six : lyt_six)
                            {
                                uint64_t counter_unmatched_one = 0;
                                for (const auto& neighbor_cell : all_defect_cells[0])
                                {
                                    lyts_six.foreach_cell(
                                        [&counter_unmatched_one, &neighbor_cell, &lyts_six, &i, this](const auto& c1)
                                        {
                                            if (c1 == neighbor_cell)
                                            {
                                                if (charge_state_to_sign(lyts_six.get_charge_state(c1)) !=
                                                    get_charge_state_defect(0, i, c1))
                                                {
                                                    counter_unmatched_one += 1;
                                                }
                                            }
                                        });
                                }
                                if (counter_unmatched_one != 0)
                                {
                                    continue;
                                }
                                //

                                uint64_t counter_unmatched_two = 0;
                                for (const auto& neighbor_cell : all_defect_cells[1])
                                {
                                    lyts_six.foreach_cell(
                                        [&counter_unmatched_two, &neighbor_cell, &lyts_six, &j, this](const auto& c1)
                                        {
                                            if (c1 == neighbor_cell)
                                            {
                                                if (charge_state_to_sign(lyts_six.get_charge_state(c1)) !=
                                                    get_charge_state_defect(1, j, c1))
                                                {
                                                    counter_unmatched_two += 1;
                                                }
                                            }
                                        });
                                }
                                if (counter_unmatched_two != 0)
                                {
                                    continue;
                                }

                                for (const auto& lyts_seven : lyt_seven)
                                {
                                    uint64_t counter_three = 0;
                                    for (const auto& neighbor_cell : all_defect_cells[2])
                                    {
                                        lyts_seven.foreach_cell(
                                            [&counter_three, &neighbor_cell, &lyts_seven, &three, this](const auto& c1)
                                            {
                                                if (c1 == neighbor_cell)
                                                {
                                                    if (charge_state_to_sign(lyts_seven.get_charge_state(c1)) !=
                                                        get_charge_state_defect(2, three, c1))
                                                    {
                                                        counter_three += 1;
                                                    }
                                                }
                                            });
                                    }
                                    if (counter_three != 0)
                                    {
                                        continue;
                                    }

                                    for (auto eight = 0u; eight < lyt_eight.size(); eight++)
                                    {
                                        uint64_t counter_unmatched_4 = 0;
                                        for (const auto& neighbor_cell_four : all_defect_cells[3])
                                        {
                                            lyt_eight[eight].foreach_cell(
                                                [&counter_unmatched_4, &neighbor_cell_four, &lyt_eight, &four, &eight,
                                                 this](const auto& c1)
                                                {
                                                    if (c1 == neighbor_cell_four)
                                                    {
                                                        if (charge_state_to_sign(lyt_eight[eight].get_charge_state(
                                                                c1)) != get_charge_state_defect(3, four, c1))
                                                        {
                                                            counter_unmatched_4 += 1;
                                                        }
                                                    }
                                                });
                                        }
                                        if (counter_unmatched_4 != 0)
                                        {
                                            continue;
                                        }

                                        for (auto nine = 0u; nine < lyt_nine.size(); nine++)
                                        {
                                            //                                            uint64_t
                                            //                                            counter_unmatched_five = 0;
                                            //                                            for (const auto&
                                            //                                            neighbor_cell_five :
                                            //                                            all_defect_cells[4])
                                            //                                            {
                                            //                                                lyt_nine[nine].foreach_cell(
                                            //                                                    [&counter_unmatched_five,
                                            //                                                    &neighbor_cell_five,
                                            //                                                    &lyt_nine, &five,
                                            //                                                    &nine, this](const
                                            //                                                    auto& c1)
                                            //                                                    {
                                            //                                                        if (c1 ==
                                            //                                                        neighbor_cell_five)
                                            //                                                        {
                                            //                                                            if
                                            //                                                            (charge_state_to_sign(lyt_nine[nine].get_charge_state(c1))
                                            //                                                            !=
                                            //                                                                get_charge_state_defect(4,
                                            //                                                                five, c1))
                                            //                                                            {
                                            //                                                                counter_unmatched_five
                                            //                                                                += 1;
                                            //                                                            }
                                            //                                                        }
                                            //                                                    });
                                            //                                            }
                                            //                                            if (counter_unmatched_five !=
                                            //                                            0)
                                            //                                            {
                                            //                                                continue;
                                            //                                            }

                                            for (auto t = 0u; t < lyt_ten.size(); t++)
                                            {
                                                for (auto l = 0u; l < lyt_11.size(); l++)
                                                {
                                                    for (const auto& lyts_12 : lyt_12)
                                                    {
                                                        //                                                        uint64_t
                                                        //                                                        counter_unmatched_9
                                                        //                                                        = 0;
                                                        //                                                        for
                                                        //                                                        (const
                                                        //                                                        auto&
                                                        //                                                        neighbor_cell_nine
                                                        //                                                        :
                                                        //                                                        all_defect_cells[8])
                                                        //                                                        {
                                                        //                                                            lyts_12.foreach_cell(
                                                        //                                                                [&counter_unmatched_9, &neighbor_cell_nine, &lyts_12, &nine, this](const auto& c1)
                                                        //                                                                {
                                                        //                                                                    if (c1 == neighbor_cell_nine)
                                                        //                                                                    {
                                                        //                                                                        if (charge_state_to_sign(lyts_12.get_charge_state(c1)) !=
                                                        //                                                                            get_charge_state_defect(8, nine, c1))
                                                        //                                                                        {
                                                        //                                                                            counter_unmatched_9 += 1;
                                                        //                                                                        }
                                                        //                                                                    }
                                                        //                                                                });
                                                        //                                                        }
                                                        //                                                        if
                                                        //                                                        (counter_unmatched_9
                                                        //                                                        != 0)
                                                        //                                                        {
                                                        //                                                            continue;
                                                        //                                                        }
                                                        //
                                                        //                                                        uint64_t
                                                        //                                                        counter_unmatched_8
                                                        //                                                        = 0;
                                                        //                                                        for
                                                        //                                                        (const
                                                        //                                                        auto&
                                                        //                                                        neighbor_cell_eight
                                                        //                                                        :
                                                        //                                                        all_defect_cells[7])
                                                        //                                                        {
                                                        //                                                            lyts_12.foreach_cell(
                                                        //                                                                [&counter_unmatched_8, &neighbor_cell_eight, &lyts_12, &eight, this](const auto& c1)
                                                        //                                                                {
                                                        //                                                                    if (c1 == neighbor_cell_eight)
                                                        //                                                                    {
                                                        //                                                                        if (charge_state_to_sign(lyts_12.get_charge_state(c1)) !=
                                                        //                                                                            get_charge_state_defect(7, eight, c1))
                                                        //                                                                        {
                                                        //                                                                            counter_unmatched_8 += 1;
                                                        //                                                                        }
                                                        //                                                                    }
                                                        //                                                                });
                                                        //                                                        }
                                                        //                                                        if
                                                        //                                                        (counter_unmatched_8
                                                        //                                                        != 0)
                                                        //                                                        {
                                                        //                                                            continue;
                                                        //                                                        }

                                                        for (const auto& lyts_13 : lyt_13)
                                                        {
                                                            //                                                            uint64_t counter_unmatched_eleven = 0;
                                                            //                                                            for (const auto& neighbor_cell : all_defect_cells[10])
                                                            //                                                            {
                                                            //                                                                lyts_13.foreach_cell(
                                                            //                                                                    [&counter_unmatched_eleven, &neighbor_cell, &lyts_13, &l, this](const auto& c1)
                                                            //                                                                    {
                                                            //                                                                        if (c1 == neighbor_cell)
                                                            //                                                                        {
                                                            //                                                                            if (charge_state_to_sign(lyts_13.get_charge_state(c1)) !=
                                                            //                                                                                get_charge_state_defect(10, l, c1))
                                                            //                                                                            {
                                                            //                                                                                counter_unmatched_eleven += 1;
                                                            //                                                                            }
                                                            //                                                                        }
                                                            //                                                                    });
                                                            //                                                            }
                                                            //                                                            if (counter_unmatched_eleven != 0)
                                                            //                                                            {
                                                            //                                                                continue;
                                                            //                                                            }
                                                            //
                                                            //                                                            uint64_t counter_unmatched_13 = 0;
                                                            //                                                            for (const auto& neighbor_cell : all_defect_cells[9])
                                                            //                                                            {
                                                            //                                                                lyts_13.foreach_cell(
                                                            //                                                                    [&counter_unmatched_13, &neighbor_cell, &lyts_13, &t, this](const auto& c1)
                                                            //                                                                    {
                                                            //                                                                        if (c1 == neighbor_cell)
                                                            //                                                                        {
                                                            //                                                                            if (charge_state_to_sign(lyts_13.get_charge_state(c1)) !=
                                                            //                                                                                get_charge_state_defect(9, t, c1))
                                                            //                                                                            {
                                                            //                                                                                counter_unmatched_13 += 1;
                                                            //                                                                            }
                                                            //                                                                        }
                                                            //                                                                    });
                                                            //                                                            }
                                                            //                                                            if (counter_unmatched_13 != 0)
                                                            //                                                            {
                                                            //                                                                continue;
                                                            //                                                            }

                                                            for (const auto& lyts_14 : lyt_14)
                                                            {
                                                                for (const auto& lyts_15 : lyt_15)
                                                                {
                                                                    for (const auto& lyts_16 : lyt_16)
                                                                    {
                                                                        lyt_one[i].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_one,
                                                                             &i](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyt_one[i].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyt_two[j].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_two,
                                                                             &j](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyt_two[j].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyt_three[three].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_three,
                                                                             &three](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1,
                                                                                    lyt_three[three].get_charge_state(
                                                                                        c1),
                                                                                    false);
                                                                            });
                                                                        lyt_four[four].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_four,
                                                                             &four](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1,
                                                                                    lyt_four[four].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyt_five[five].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_five,
                                                                             &five](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1,
                                                                                    lyt_five[five].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_six.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_six](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_six.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_seven.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_seven](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_seven.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyt_eight[eight].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_eight,
                                                                             &eight](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1,
                                                                                    lyt_eight[eight].get_charge_state(
                                                                                        c1),
                                                                                    false);
                                                                            });
                                                                        lyt_nine[nine].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_nine,
                                                                             &nine](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1,
                                                                                    lyt_nine[nine].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyt_ten[t].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_ten,
                                                                             &t](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyt_ten[t].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyt_11[l].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_11,
                                                                             &l](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyt_11[l].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_12.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_12](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_12.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_13.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_13](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_13.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_14.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_14](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_14.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_15.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_15](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_15.get_charge_state(c1),
                                                                                    false);
                                                                            });

                                                                        lyts_16.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_16](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_16.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        charge_lyt.update_after_charge_change();
                                                                        if (charge_lyt.is_physically_valid())
                                                                        {
                                                                            if (charge_lyt.get_system_energy() <
                                                                                energy_threas)
                                                                            {
                                                                                std::vector<
                                                                                    charge_distribution_surface<Lyt>>
                                                                                    lyts{};
                                                                                std::cout
                                                                                    << charge_lyt.get_system_energy()
                                                                                    << std::endl;

                                                                                sidb_simulation_result<Lyt>
                                                                                    sim_result{};
                                                                                sim_result.algorithm_name = "ExGS";
                                                                                charge_distribution_surface<Lyt>
                                                                                    charge_lyt_copy{charge_lyt};
                                                                                lyts.emplace_back(charge_lyt_copy);
                                                                                sim_result.charge_distributions = lyts;
                                                                                energy_threas =
                                                                                    charge_lyt.get_system_energy();
                                                                                write_sqd_sim_result<Lyt>(
                                                                                    sim_result,
                                                                                    "/Users/jandrewniok/"
                                                                                    "CLionProjects/"
                                                                                    "fiction_fork/experiments/"
                                                                                    "result.xml");
                                                                            }
                                                                        }
                                                                        counter += 1;
                                                                        if (counter % 100000 == 0)
                                                                        {
                                                                            std::cout << counter << std::endl;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void combining_all_14()
    {
        auto compareFunc =
            [](const charge_distribution_surface<Lyt>& lyt1, const charge_distribution_surface<Lyt>& lyt2)
        { return lyt1.get_system_energy() < lyt2.get_system_energy(); };

        std::cout << "combining starts: " << std::to_string(lyts_of_regions.size()) << std::endl;
        uint64_t counter_lyts = 1;
        for (const auto& lyts_region : lyts_of_regions)
        {
            counter_lyts *= lyts_region.size();
        }
        std::cout << "number enumerations: " << std::to_string(counter_lyts) << std::endl;

        auto lyt_one   = lyts_of_regions[0];
        auto lyt_two   = lyts_of_regions[1];
        auto lyt_three = lyts_of_regions[2];
        auto lyt_four  = lyts_of_regions[3];
        auto lyt_five  = lyts_of_regions[4];
        auto lyt_six   = lyts_of_regions[5];
        auto lyt_seven = lyts_of_regions[6];
        auto lyt_eight = lyts_of_regions[7];
        auto lyt_nine  = lyts_of_regions[8];
        auto lyt_ten   = lyts_of_regions[9];
        auto lyt_11    = lyts_of_regions[10];
        auto lyt_12    = lyts_of_regions[11];
        auto lyt_13    = lyts_of_regions[12];
        auto lyt_14    = lyts_of_regions[13];

        std::sort(lyt_one.begin(), lyt_one.end(), compareFunc);
        std::sort(lyt_two.begin(), lyt_two.end(), compareFunc);
        std::sort(lyt_three.begin(), lyt_three.end(), compareFunc);
        std::sort(lyt_four.begin(), lyt_four.end(), compareFunc);
        std::sort(lyt_five.begin(), lyt_five.end(), compareFunc);
        std::sort(lyt_six.begin(), lyt_six.end(), compareFunc);
        std::sort(lyt_seven.begin(), lyt_seven.end(), compareFunc);
        std::sort(lyt_eight.begin(), lyt_eight.end(), compareFunc);
        std::sort(lyt_nine.begin(), lyt_nine.end(), compareFunc);
        std::sort(lyt_ten.begin(), lyt_ten.end(), compareFunc);
        std::sort(lyt_11.begin(), lyt_11.end(), compareFunc);
        std::sort(lyt_12.begin(), lyt_12.end(), compareFunc);
        std::sort(lyt_13.begin(), lyt_13.end(), compareFunc);
        std::sort(lyt_14.begin(), lyt_14.end(), compareFunc);

        int                                           number = 5;
        std::vector<charge_distribution_surface<Lyt>> lyt_ones(
            lyt_one.begin(), lyt_one.begin() + std::min(number, static_cast<int>(lyt_one.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_twos(
            lyt_two.begin(), lyt_two.begin() + std::min(number, static_cast<int>(lyt_two.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_threes(
            lyt_three.begin(), lyt_three.begin() + std::min(number, static_cast<int>(lyt_three.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_fours(
            lyt_four.begin(), lyt_four.begin() + std::min(number, static_cast<int>(lyt_four.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_fives(
            lyt_five.begin(), lyt_five.begin() + std::min(number, static_cast<int>(lyt_five.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_sixs(
            lyt_six.begin(), lyt_six.begin() + std::min(number, static_cast<int>(lyt_six.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_sevens(
            lyt_seven.begin(), lyt_seven.begin() + std::min(number, static_cast<int>(lyt_seven.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_eights(
            lyt_eight.begin(), lyt_eight.begin() + std::min(number, static_cast<int>(lyt_eight.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_nines(
            lyt_nine.begin(), lyt_nine.begin() + std::min(number, static_cast<int>(lyt_nine.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_tens(
            lyt_ten.begin(), lyt_ten.begin() + std::min(number, static_cast<int>(lyt_ten.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_11s(
            lyt_11.begin(), lyt_11.begin() + std::min(number, static_cast<int>(lyt_11.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_12s(
            lyt_12.begin(), lyt_12.begin() + std::min(number, static_cast<int>(lyt_12.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_13s(
            lyt_13.begin(), lyt_13.begin() + std::min(number, static_cast<int>(lyt_13.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_14s(
            lyt_14.begin(), lyt_14.begin() + std::min(number, static_cast<int>(lyt_14.size())));

        charge_distribution_surface<Lyt> charge_lyt{layout};
        uint64_t                         counter = 0;
        std::vector<double>              valid_energies{};
        double                           energy_threas = 1000;
        for (auto i = 0u; i < lyt_one.size(); i++)
        {
            for (auto j = 0u; j < lyt_two.size(); j++)
            {
                for (auto three = 0u; three < lyt_three.size(); three++)
                {
                    for (auto four = 0u; four < lyt_four.size(); four++)
                    {
                        for (auto five = 0u; five < lyt_five.size(); five++)
                        {
                            for (const auto& lyts_six : lyt_six)
                            {
                                uint64_t counter_unmatched_one = 0;
                                for (const auto& neighbor_cell : all_defect_cells[0])
                                {
                                    lyts_six.foreach_cell(
                                        [&counter_unmatched_one, &neighbor_cell, &lyts_six, &i, this](const auto& c1)
                                        {
                                            if (c1 == neighbor_cell)
                                            {
                                                if (charge_state_to_sign(lyts_six.get_charge_state(c1)) !=
                                                    get_charge_state_defect(0, i, c1))
                                                {
                                                    counter_unmatched_one += 1;
                                                }
                                            }
                                        });
                                }
                                if (counter_unmatched_one != 0)
                                {
                                    continue;
                                }

                                uint64_t counter_unmatched = 0;
                                for (const auto& neighbor_cell : all_defect_cells[1])
                                {
                                    lyts_six.foreach_cell(
                                        [&counter_unmatched, &neighbor_cell, &lyts_six, &j, this](const auto& c1)
                                        {
                                            if (c1 == neighbor_cell)
                                            {
                                                if (charge_state_to_sign(lyts_six.get_charge_state(c1)) !=
                                                    get_charge_state_defect(1, j, c1))
                                                {
                                                    counter_unmatched += 1;
                                                }
                                            }
                                        });
                                }
                                if (counter_unmatched != 0)
                                {
                                    continue;
                                }

                                for (const auto& lyts_seven : lyt_seven)
                                {
                                    uint64_t counter_three = 0;
                                    for (const auto& neighbor_cell : all_defect_cells[2])
                                    {
                                        lyts_seven.foreach_cell(
                                            [&counter_three, &neighbor_cell, &lyts_seven, &three, this](const auto& c1)
                                            {
                                                if (c1 == neighbor_cell)
                                                {
                                                    if (charge_state_to_sign(lyts_seven.get_charge_state(c1)) !=
                                                        get_charge_state_defect(2, three, c1))
                                                    {
                                                        counter_three += 1;
                                                    }
                                                }
                                            });
                                    }
                                    if (counter_three != 0)
                                    {
                                        continue;
                                    }

                                    for (auto eight = 0u; eight < lyt_eight.size(); eight++)
                                    {
                                        uint64_t counter_unmatched_4 = 0;
                                        for (const auto& neighbor_cell_four : all_defect_cells[3])
                                        {
                                            lyt_eight[eight].foreach_cell(
                                                [&counter_unmatched_4, &neighbor_cell_four, &lyt_eight, &four, &eight,
                                                 this](const auto& c1)
                                                {
                                                    if (c1 == neighbor_cell_four)
                                                    {
                                                        if (charge_state_to_sign(lyt_eight[eight].get_charge_state(
                                                                c1)) != get_charge_state_defect(3, four, c1))
                                                        {
                                                            counter_unmatched_4 += 1;
                                                        }
                                                    }
                                                });
                                        }
                                        if (counter_unmatched_4 != 0)
                                        {
                                            continue;
                                        }

                                        for (auto nine = 0u; nine < lyt_nine.size(); nine++)
                                        {
                                            uint64_t counter_unmatched_five = 0;
                                            for (const auto& neighbor_cell_five : all_defect_cells[4])
                                            {
                                                lyt_nine[nine].foreach_cell(
                                                    [&counter_unmatched_five, &neighbor_cell_five, &lyt_nine, &five,
                                                     &nine, this](const auto& c1)
                                                    {
                                                        if (c1 == neighbor_cell_five)
                                                        {
                                                            if (charge_state_to_sign(lyt_nine[nine].get_charge_state(
                                                                    c1)) != get_charge_state_defect(4, five, c1))
                                                            {
                                                                counter_unmatched_five += 1;
                                                            }
                                                        }
                                                    });
                                            }
                                            if (counter_unmatched_five != 0)
                                            {
                                                continue;
                                            }

                                            for (auto t = 0u; t < lyt_ten.size(); t++)
                                            {
                                                for (auto l = 0u; l < lyt_11.size(); l++)
                                                {
                                                    for (const auto& lyts_12 : lyt_12)
                                                    {
                                                        uint64_t counter_unmatched_9 = 0;
                                                        for (const auto& neighbor_cell_nine : all_defect_cells[8])
                                                        {
                                                            lyts_12.foreach_cell(
                                                                [&counter_unmatched_9, &neighbor_cell_nine, &lyts_12,
                                                                 &nine, this](const auto& c1)
                                                                {
                                                                    if (c1 == neighbor_cell_nine)
                                                                    {
                                                                        if (charge_state_to_sign(
                                                                                lyts_12.get_charge_state(c1)) !=
                                                                            get_charge_state_defect(8, nine, c1))
                                                                        {
                                                                            counter_unmatched_9 += 1;
                                                                        }
                                                                    }
                                                                });
                                                        }
                                                        if (counter_unmatched_9 != 0)
                                                        {
                                                            continue;
                                                        }

                                                        uint64_t counter_unmatched_8 = 0;
                                                        for (const auto& neighbor_cell_eight : all_defect_cells[7])
                                                        {
                                                            lyts_12.foreach_cell(
                                                                [&counter_unmatched_8, &neighbor_cell_eight, &lyts_12,
                                                                 &eight, this](const auto& c1)
                                                                {
                                                                    if (c1 == neighbor_cell_eight)
                                                                    {
                                                                        if (charge_state_to_sign(
                                                                                lyts_12.get_charge_state(c1)) !=
                                                                            get_charge_state_defect(7, eight, c1))
                                                                        {
                                                                            counter_unmatched_8 += 1;
                                                                        }
                                                                    }
                                                                });
                                                        }
                                                        if (counter_unmatched_8 != 0)
                                                        {
                                                            continue;
                                                        }

                                                        for (const auto& lyts_13 : lyt_13)
                                                        {
                                                            uint64_t counter_unmatched_eleven = 0;
                                                            for (const auto& neighbor_cell : all_defect_cells[10])
                                                            {
                                                                lyts_13.foreach_cell(
                                                                    [&counter_unmatched_eleven, &neighbor_cell,
                                                                     &lyts_13, &l, this](const auto& c1)
                                                                    {
                                                                        if (c1 == neighbor_cell)
                                                                        {
                                                                            if (charge_state_to_sign(
                                                                                    lyts_13.get_charge_state(c1)) !=
                                                                                get_charge_state_defect(10, l, c1))
                                                                            {
                                                                                counter_unmatched_eleven += 1;
                                                                            }
                                                                        }
                                                                    });
                                                            }
                                                            if (counter_unmatched_eleven != 0)
                                                            {
                                                                continue;
                                                            }

                                                            uint64_t counter_unmatched_13 = 0;
                                                            for (const auto& neighbor_cell : all_defect_cells[9])
                                                            {
                                                                lyts_13.foreach_cell(
                                                                    [&counter_unmatched_13, &neighbor_cell, &lyts_13,
                                                                     &t, this](const auto& c1)
                                                                    {
                                                                        if (c1 == neighbor_cell)
                                                                        {
                                                                            if (charge_state_to_sign(
                                                                                    lyts_13.get_charge_state(c1)) !=
                                                                                get_charge_state_defect(9, t, c1))
                                                                            {
                                                                                counter_unmatched_13 += 1;
                                                                            }
                                                                        }
                                                                    });
                                                            }
                                                            if (counter_unmatched_13 != 0)
                                                            {
                                                                continue;
                                                            }

                                                            for (const auto& lyts_14 : lyt_14)
                                                            {
                                                                lyt_one[i].foreach_cell(
                                                                    [this, &charge_lyt, &lyt_one, &i](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyt_one[i].get_charge_state(c1), false);
                                                                    });
                                                                lyt_two[j].foreach_cell(
                                                                    [this, &charge_lyt, &lyt_two, &j](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyt_two[j].get_charge_state(c1), false);
                                                                    });
                                                                lyt_three[three].foreach_cell(
                                                                    [this, &charge_lyt, &lyt_three,
                                                                     &three](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyt_three[three].get_charge_state(c1),
                                                                            false);
                                                                    });
                                                                lyt_four[four].foreach_cell(
                                                                    [this, &charge_lyt, &lyt_four,
                                                                     &four](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyt_four[four].get_charge_state(c1),
                                                                            false);
                                                                    });
                                                                lyt_five[five].foreach_cell(
                                                                    [this, &charge_lyt, &lyt_five,
                                                                     &five](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyt_five[five].get_charge_state(c1),
                                                                            false);
                                                                    });
                                                                lyts_six.foreach_cell(
                                                                    [this, &charge_lyt, &lyts_six](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyts_six.get_charge_state(c1), false);
                                                                    });
                                                                lyts_seven.foreach_cell(
                                                                    [this, &charge_lyt, &lyts_seven](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyts_seven.get_charge_state(c1), false);
                                                                    });
                                                                lyt_eight[eight].foreach_cell(
                                                                    [this, &charge_lyt, &lyt_eight,
                                                                     &eight](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyt_eight[eight].get_charge_state(c1),
                                                                            false);
                                                                    });
                                                                lyt_nine[nine].foreach_cell(
                                                                    [this, &charge_lyt, &lyt_nine,
                                                                     &nine](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyt_nine[nine].get_charge_state(c1),
                                                                            false);
                                                                    });
                                                                lyt_ten[t].foreach_cell(
                                                                    [this, &charge_lyt, &lyt_ten, &t](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyt_ten[t].get_charge_state(c1), false);
                                                                    });
                                                                lyt_11[l].foreach_cell(
                                                                    [this, &charge_lyt, &lyt_11, &l](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyt_11[l].get_charge_state(c1), false);
                                                                    });
                                                                lyts_12.foreach_cell(
                                                                    [this, &charge_lyt, &lyts_12](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyts_12.get_charge_state(c1), false);
                                                                    });
                                                                lyts_13.foreach_cell(
                                                                    [this, &charge_lyt, &lyts_13](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyts_13.get_charge_state(c1), false);
                                                                    });
                                                                lyts_14.foreach_cell(
                                                                    [this, &charge_lyt, &lyts_14](const auto& c1) {
                                                                        charge_lyt.assign_charge_state(
                                                                            c1, lyts_14.get_charge_state(c1), false);
                                                                    });
                                                                charge_lyt.update_after_charge_change();
                                                                if (charge_lyt.is_physically_valid())
                                                                {
                                                                    if (charge_lyt.get_system_energy() < energy_threas)
                                                                    {
                                                                        std::vector<charge_distribution_surface<Lyt>>
                                                                            lyts{};
                                                                        std::cout << charge_lyt.get_system_energy()
                                                                                  << std::endl;

                                                                        sidb_simulation_result<Lyt> sim_result{};
                                                                        sim_result.algorithm_name = "ExGS";
                                                                        charge_distribution_surface<Lyt>
                                                                            charge_lyt_copy{charge_lyt};
                                                                        lyts.emplace_back(charge_lyt_copy);
                                                                        sim_result.charge_distributions = lyts;
                                                                        energy_threas = charge_lyt.get_system_energy();
                                                                        write_sqd_sim_result<Lyt>(
                                                                            sim_result, "/Users/jandrewniok/"
                                                                                        "CLionProjects/"
                                                                                        "fiction_fork/experiments/"
                                                                                        "result.xml");
                                                                    }
                                                                }
                                                                counter += 1;
                                                                if (counter % 100000 == 0)
                                                                {
                                                                    std::cout << counter << std::endl;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void combining_all_16()
    {
        auto compareFunc =
            [](const charge_distribution_surface<Lyt>& lyt1, const charge_distribution_surface<Lyt>& lyt2)
        { return lyt1.get_system_energy() < lyt2.get_system_energy(); };

        std::cout << "combining starts: " << std::to_string(lyts_of_regions.size()) << std::endl;
        uint64_t counter_lyts = 1;
        for (const auto& lyts_region : lyts_of_regions)
        {
            counter_lyts *= lyts_region.size();
        }
        std::cout << "number enumerations: " << std::to_string(counter_lyts) << std::endl;

        auto lyt_one   = lyts_of_regions[0];
        auto lyt_two   = lyts_of_regions[1];
        auto lyt_three = lyts_of_regions[2];
        auto lyt_four  = lyts_of_regions[3];
        auto lyt_five  = lyts_of_regions[4];
        auto lyt_six   = lyts_of_regions[5];
        auto lyt_seven = lyts_of_regions[6];
        auto lyt_eight = lyts_of_regions[7];
        auto lyt_nine  = lyts_of_regions[8];
        auto lyt_ten   = lyts_of_regions[9];
        auto lyt_11    = lyts_of_regions[10];
        auto lyt_12    = lyts_of_regions[11];
        auto lyt_13    = lyts_of_regions[12];
        auto lyt_14    = lyts_of_regions[13];
        auto lyt_15    = lyts_of_regions[14];
        auto lyt_16    = lyts_of_regions[15];

        std::sort(lyt_one.begin(), lyt_one.end(), compareFunc);
        std::sort(lyt_two.begin(), lyt_two.end(), compareFunc);
        std::sort(lyt_three.begin(), lyt_three.end(), compareFunc);
        std::sort(lyt_four.begin(), lyt_four.end(), compareFunc);
        std::sort(lyt_five.begin(), lyt_five.end(), compareFunc);
        std::sort(lyt_six.begin(), lyt_six.end(), compareFunc);
        std::sort(lyt_seven.begin(), lyt_seven.end(), compareFunc);
        std::sort(lyt_eight.begin(), lyt_eight.end(), compareFunc);
        std::sort(lyt_nine.begin(), lyt_nine.end(), compareFunc);
        std::sort(lyt_ten.begin(), lyt_ten.end(), compareFunc);
        std::sort(lyt_11.begin(), lyt_11.end(), compareFunc);
        std::sort(lyt_12.begin(), lyt_12.end(), compareFunc);
        std::sort(lyt_13.begin(), lyt_13.end(), compareFunc);
        std::sort(lyt_14.begin(), lyt_14.end(), compareFunc);
        std::sort(lyt_15.begin(), lyt_15.end(), compareFunc);
        std::sort(lyt_16.begin(), lyt_16.end(), compareFunc);

        int                                           number = 5;
        std::vector<charge_distribution_surface<Lyt>> lyt_ones(
            lyt_one.begin(), lyt_one.begin() + std::min(number, static_cast<int>(lyt_one.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_twos(
            lyt_two.begin(), lyt_two.begin() + std::min(number, static_cast<int>(lyt_two.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_threes(
            lyt_three.begin(), lyt_three.begin() + std::min(number, static_cast<int>(lyt_three.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_fours(
            lyt_four.begin(), lyt_four.begin() + std::min(number, static_cast<int>(lyt_four.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_fives(
            lyt_five.begin(), lyt_five.begin() + std::min(number, static_cast<int>(lyt_five.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_sixs(
            lyt_six.begin(), lyt_six.begin() + std::min(number, static_cast<int>(lyt_six.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_sevens(
            lyt_seven.begin(), lyt_seven.begin() + std::min(number, static_cast<int>(lyt_seven.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_eights(
            lyt_eight.begin(), lyt_eight.begin() + std::min(number, static_cast<int>(lyt_eight.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_nines(
            lyt_nine.begin(), lyt_nine.begin() + std::min(number, static_cast<int>(lyt_nine.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_tens(
            lyt_ten.begin(), lyt_ten.begin() + std::min(number, static_cast<int>(lyt_ten.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_11s(
            lyt_11.begin(), lyt_11.begin() + std::min(number, static_cast<int>(lyt_11.size())));

        std::vector<charge_distribution_surface<Lyt>> lyt_12s(
            lyt_12.begin(), lyt_12.begin() + std::min(number, static_cast<int>(lyt_12.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_13s(
            lyt_13.begin(), lyt_13.begin() + std::min(number, static_cast<int>(lyt_13.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_14s(
            lyt_14.begin(), lyt_14.begin() + std::min(number, static_cast<int>(lyt_14.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_15s(
            lyt_15.begin(), lyt_15.begin() + std::min(number, static_cast<int>(lyt_15.size())));
        std::vector<charge_distribution_surface<Lyt>> lyt_16s(
            lyt_16.begin(), lyt_16.begin() + std::min(number, static_cast<int>(lyt_16.size())));

        charge_distribution_surface<Lyt> charge_lyt{layout};
        uint64_t                         counter = 0;
        std::vector<double>              valid_energies{};
        double                           energy_threas = 1000;
        for (auto i = 0u; i < lyt_one.size(); i++)
        {
            for (auto j = 0u; j < lyt_two.size(); j++)
            {
                for (auto three = 0u; three < lyt_three.size(); three++)
                {
                    for (auto four = 0u; four < lyt_four.size(); four++)
                    {
                        for (auto five = 0u; five < lyt_five.size(); five++)
                        {
                            for (const auto& lyts_six : lyt_six)
                            {
                                uint64_t counter_unmatched_one = 0;
                                for (const auto& neighbor_cell : all_defect_cells[0])
                                {
                                    lyts_six.foreach_cell(
                                        [&counter_unmatched_one, &neighbor_cell, &lyts_six, &i, this](const auto& c1)
                                        {
                                            if (c1 == neighbor_cell)
                                            {
                                                if (charge_state_to_sign(lyts_six.get_charge_state(c1)) !=
                                                    get_charge_state_defect(0, i, c1))
                                                {
                                                    counter_unmatched_one += 1;
                                                }
                                            }
                                        });
                                }
                                if (counter_unmatched_one != 0)
                                {
                                    continue;
                                }

                                uint64_t counter_unmatched = 0;
                                for (const auto& neighbor_cell : all_defect_cells[1])
                                {
                                    lyts_six.foreach_cell(
                                        [&counter_unmatched, &neighbor_cell, &lyts_six, &j, this](const auto& c1)
                                        {
                                            if (c1 == neighbor_cell)
                                            {
                                                if (charge_state_to_sign(lyts_six.get_charge_state(c1)) !=
                                                    get_charge_state_defect(1, j, c1))
                                                {
                                                    counter_unmatched += 1;
                                                }
                                            }
                                        });
                                }
                                if (counter_unmatched != 0)
                                {
                                    continue;
                                }

                                for (const auto& lyts_seven : lyt_seven)
                                {
                                    uint64_t counter_three = 0;
                                    for (const auto& neighbor_cell : all_defect_cells[2])
                                    {
                                        lyts_seven.foreach_cell(
                                            [&counter_three, &neighbor_cell, &lyts_seven, &three, this](const auto& c1)
                                            {
                                                if (c1 == neighbor_cell)
                                                {
                                                    if (charge_state_to_sign(lyts_seven.get_charge_state(c1)) !=
                                                        get_charge_state_defect(2, three, c1))
                                                    {
                                                        counter_three += 1;
                                                    }
                                                }
                                            });
                                    }
                                    if (counter_three != 0)
                                    {
                                        continue;
                                    }

                                    for (auto eight = 0u; eight < lyt_eight.size(); eight++)
                                    {
                                        uint64_t counter_unmatched_4 = 0;
                                        for (const auto& neighbor_cell_four : all_defect_cells[3])
                                        {
                                            lyt_eight[eight].foreach_cell(
                                                [&counter_unmatched_4, &neighbor_cell_four, &lyt_eight, &four, &eight,
                                                 this](const auto& c1)
                                                {
                                                    if (c1 == neighbor_cell_four)
                                                    {
                                                        if (charge_state_to_sign(lyt_eight[eight].get_charge_state(
                                                                c1)) != get_charge_state_defect(3, four, c1))
                                                        {
                                                            counter_unmatched_4 += 1;
                                                        }
                                                    }
                                                });
                                        }
                                        if (counter_unmatched_4 != 0)
                                        {
                                            continue;
                                        }

                                        for (auto nine = 0u; nine < lyt_nine.size(); nine++)
                                        {
                                            //                                            uint64_t
                                            //                                            counter_unmatched_five = 0;
                                            //                                            for (const auto&
                                            //                                            neighbor_cell_five :
                                            //                                            all_defect_cells[4])
                                            //                                            {
                                            //                                                lyt_nine[nine].foreach_cell(
                                            //                                                    [&counter_unmatched_five,
                                            //                                                    &neighbor_cell_five,
                                            //                                                    &lyt_nine, &five,
                                            //                                                    &nine, this](const
                                            //                                                    auto& c1)
                                            //                                                    {
                                            //                                                        if (c1 ==
                                            //                                                        neighbor_cell_five)
                                            //                                                        {
                                            //                                                            if
                                            //                                                            (charge_state_to_sign(lyt_nine[nine].get_charge_state(c1))
                                            //                                                            !=
                                            //                                                                get_charge_state_defect(4,
                                            //                                                                five, c1))
                                            //                                                            {
                                            //                                                                counter_unmatched_five
                                            //                                                                += 1;
                                            //                                                            }
                                            //                                                        }
                                            //                                                    });
                                            //                                            }
                                            //                                            if (counter_unmatched_five !=
                                            //                                            0)
                                            //                                            {
                                            //                                                continue;
                                            //                                            }

                                            for (auto t = 0u; t < lyt_ten.size(); t++)
                                            {
                                                for (auto l = 0u; l < lyt_11.size(); l++)
                                                {
                                                    for (const auto& lyts_12 : lyt_12)
                                                    {
                                                        //                                                        uint64_t
                                                        //                                                        counter_unmatched_9
                                                        //                                                        = 0;
                                                        //                                                        for
                                                        //                                                        (const
                                                        //                                                        auto&
                                                        //                                                        neighbor_cell_nine
                                                        //                                                        :
                                                        //                                                        all_defect_cells[8])
                                                        //                                                        {
                                                        //                                                            lyts_12.foreach_cell(
                                                        //                                                                [&counter_unmatched_9, &neighbor_cell_nine, &lyts_12,
                                                        //                                                                 &nine, this](const auto& c1)
                                                        //                                                                {
                                                        //                                                                    if (c1 == neighbor_cell_nine)
                                                        //                                                                    {
                                                        //                                                                        if (charge_state_to_sign(
                                                        //                                                                                lyts_12.get_charge_state(c1)) !=
                                                        //                                                                            get_charge_state_defect(8, nine, c1))
                                                        //                                                                        {
                                                        //                                                                            counter_unmatched_9 += 1;
                                                        //                                                                        }
                                                        //                                                                    }
                                                        //                                                                });
                                                        //                                                        }
                                                        //                                                        if
                                                        //                                                        (counter_unmatched_9
                                                        //                                                        != 0)
                                                        //                                                        {
                                                        //                                                            continue;
                                                        //                                                        }
                                                        //
                                                        //                                                        uint64_t
                                                        //                                                        counter_unmatched_8
                                                        //                                                        = 0;
                                                        //                                                        for
                                                        //                                                        (const
                                                        //                                                        auto&
                                                        //                                                        neighbor_cell_eight
                                                        //                                                        :
                                                        //                                                        all_defect_cells[7])
                                                        //                                                        {
                                                        //                                                            lyts_12.foreach_cell(
                                                        //                                                                [&counter_unmatched_8, &neighbor_cell_eight, &lyts_12,
                                                        //                                                                 &eight, this](const auto& c1)
                                                        //                                                                {
                                                        //                                                                    if (c1 == neighbor_cell_eight)
                                                        //                                                                    {
                                                        //                                                                        if (charge_state_to_sign(
                                                        //                                                                                lyts_12.get_charge_state(c1)) !=
                                                        //                                                                            get_charge_state_defect(7, eight, c1))
                                                        //                                                                        {
                                                        //                                                                            counter_unmatched_8 += 1;
                                                        //                                                                        }
                                                        //                                                                    }
                                                        //                                                                });
                                                        //                                                        }
                                                        //                                                        if
                                                        //                                                        (counter_unmatched_8
                                                        //                                                        != 0)
                                                        //                                                        {
                                                        //                                                            continue;
                                                        //                                                        }

                                                        for (const auto& lyts_13 : lyt_13)
                                                        {
                                                            //                                                            uint64_t counter_unmatched_eleven = 0;
                                                            //                                                            for (const auto& neighbor_cell : all_defect_cells[10])
                                                            //                                                            {
                                                            //                                                                lyts_13.foreach_cell(
                                                            //                                                                    [&counter_unmatched_eleven, &neighbor_cell, &lyts_13, &l, this](const auto& c1)
                                                            //                                                                    {
                                                            //                                                                        if (c1 == neighbor_cell)
                                                            //                                                                        {
                                                            //                                                                            if (charge_state_to_sign(lyts_13.get_charge_state(c1)) !=
                                                            //                                                                                get_charge_state_defect(10, l, c1))
                                                            //                                                                            {
                                                            //                                                                                counter_unmatched_eleven += 1;
                                                            //                                                                            }
                                                            //                                                                        }
                                                            //                                                                    });
                                                            //                                                            }
                                                            //                                                            if (counter_unmatched_eleven != 0)
                                                            //                                                            {
                                                            //                                                                continue;
                                                            //                                                            }

                                                            //                                                            uint64_t counter_unmatched_13 = 0;
                                                            //                                                            for (const auto& neighbor_cell : all_defect_cells[9])
                                                            //                                                            {
                                                            //                                                                lyts_13.foreach_cell(
                                                            //                                                                    [&counter_unmatched_13, &neighbor_cell, &lyts_13, &t, this](const auto& c1)
                                                            //                                                                    {
                                                            //                                                                        if (c1 == neighbor_cell)
                                                            //                                                                        {
                                                            //                                                                            if (charge_state_to_sign(lyts_13.get_charge_state(c1)) !=
                                                            //                                                                                get_charge_state_defect(9, t, c1))
                                                            //                                                                            {
                                                            //                                                                                counter_unmatched_13 += 1;
                                                            //                                                                            }
                                                            //                                                                        }
                                                            //                                                                    });
                                                            //                                                            }
                                                            //                                                            if (counter_unmatched_13 != 0)
                                                            //                                                            {
                                                            //                                                                continue;
                                                            //                                                            }

                                                            for (const auto& lyts_14 : lyt_14)
                                                            {
                                                                for (const auto& lyts_15 : lyt_15)
                                                                {
                                                                    for (const auto& lyts_16 : lyt_16)
                                                                    {
                                                                        lyt_one[i].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_one,
                                                                             &i](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyt_one[i].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyt_two[j].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_two,
                                                                             &j](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyt_two[j].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyt_three[three].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_three,
                                                                             &three](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1,
                                                                                    lyt_three[three].get_charge_state(
                                                                                        c1),
                                                                                    false);
                                                                            });
                                                                        lyt_four[four].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_four,
                                                                             &four](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1,
                                                                                    lyt_four[four].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyt_five[five].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_five,
                                                                             &five](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1,
                                                                                    lyt_five[five].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_six.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_six](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_six.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_seven.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_seven](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_seven.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyt_eight[eight].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_eight,
                                                                             &eight](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1,
                                                                                    lyt_eight[eight].get_charge_state(
                                                                                        c1),
                                                                                    false);
                                                                            });
                                                                        lyt_nine[nine].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_nine,
                                                                             &nine](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1,
                                                                                    lyt_nine[nine].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyt_ten[t].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_ten,
                                                                             &t](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyt_ten[t].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyt_11[l].foreach_cell(
                                                                            [this, &charge_lyt, &lyt_11,
                                                                             &l](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyt_11[l].get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_12.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_12](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_12.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_13.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_13](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_13.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_14.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_14](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_14.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        lyts_15.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_15](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_15.get_charge_state(c1),
                                                                                    false);
                                                                            });

                                                                        lyts_16.foreach_cell(
                                                                            [this, &charge_lyt,
                                                                             &lyts_16](const auto& c1) {
                                                                                charge_lyt.assign_charge_state(
                                                                                    c1, lyts_16.get_charge_state(c1),
                                                                                    false);
                                                                            });
                                                                        charge_lyt.update_after_charge_change();
                                                                        if (charge_lyt.is_physically_valid())
                                                                        {
                                                                            if (charge_lyt.get_system_energy() <
                                                                                energy_threas)
                                                                            {
                                                                                std::vector<
                                                                                    charge_distribution_surface<Lyt>>
                                                                                    lyts{};
                                                                                std::cout
                                                                                    << charge_lyt.get_system_energy()
                                                                                    << std::endl;

                                                                                sidb_simulation_result<Lyt>
                                                                                    sim_result{};
                                                                                sim_result.algorithm_name = "ExGS";
                                                                                charge_distribution_surface<Lyt>
                                                                                    charge_lyt_copy{charge_lyt};
                                                                                lyts.emplace_back(charge_lyt_copy);
                                                                                sim_result.charge_distributions = lyts;
                                                                                energy_threas =
                                                                                    charge_lyt.get_system_energy();
                                                                                write_sqd_sim_result<Lyt>(
                                                                                    sim_result,
                                                                                    "/Users/jandrewniok/"
                                                                                    "CLionProjects/"
                                                                                    "fiction_fork/experiments/"
                                                                                    "result.xml");
                                                                            }
                                                                        }
                                                                        counter += 1;
                                                                        if (counter % 100000 == 0)
                                                                        {
                                                                            std::cout << counter << std::endl;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    int8_t get_charge_state_defect(const uint64_t region, const uint64_t i, const typename Lyt::cell& cell)
    {
        auto neighbor_cell = all_defect_cells[region];
        for (auto l = 0u; l < neighbor_cell.size(); l++)
        {
            if (neighbor_cell[l] == cell)
            {
                return all_defect_charges[region][i][l];
            }
        }
    }

  private:
    Lyt&                                                                     layout;
    sidb_simulation_parameters                                               parameter{};
    layout_sim_stats<Lyt>&                                                   statistic{};
    typename cube::coord_t                                                   start_cell{};
    typename cube::coord_t                                                   rightest_cell{};
    typename cube::coord_t                                                   lowest_cell{};
    typename cube::coord_t                                                   left_corner_cell{};
    std::set<typename Lyt::cell>                                             cells{};
    std::vector<std::set<uint64_t>>                                          all_neighbor_pairs{};
    std::vector<typename Lyt::cell>                                          border_cells{};
    std::vector<sidb_charge_state>                                           border_cell_charge{};
    uint64_t                                                                 border_cell_index{};
    uint64_t                                                                 border_cell_max_charge_index{};
    std::vector<typename Lyt::cell>                                          defect_cell{};
    std::vector<int8_t>                                                      defect_charge{};
    std::set<typename Lyt::cell>                                             total_cells{};
    uint64_t                                                                 layout_num{};
    std::vector<std::vector<std::unordered_map<typename Lyt::cell, int8_t>>> border_cells_and_charge{};
    std::vector<std::vector<charge_distribution_surface<Lyt>>>               all_charge_lyts{};
    std::vector<std::vector<charge_distribution_surface<Lyt>>>               lyts_of_regions{};
    std::vector<std::vector<uint64_t>>                                       charge_index_innen{};
    std::vector<std::vector<typename Lyt::cell>>                             all_defect_cells{};
    std::vector<std::vector<std::vector<int8_t>>>                            all_defect_charges{};
    std::vector<Lyt>                                                         all_layouts{};
    std::vector<uint64_t>                                                    region_num{};
    uint64_t                                                                 region_col_counter{0};
    uint64_t                                                                 region_counter{0};
    std::unordered_map<uint64_t, std::vector<typename Lyt::cell>>            border_cell_region{};
};

template <typename Lyt>
bool layout_simulation(Lyt& lyt, const sidb_simulation_parameters& params = sidb_simulation_parameters{},
                       layout_sim_stats<Lyt>* ps = nullptr)
{

    layout_sim_stats<Lyt> st{};

    detail::layout_simulation_impl<Lyt> p{lyt, params, st};

    auto result = p.run_simulation_hexagon();
    p.finding_nn();
    p.combining_all_16();

    if (ps)
    {
        *ps = st;
    }

    return result;
}

}  // namespace detail
}  // namespace fiction
#endif  // FICTION_LAYOUT_SIMULATION_HPP
