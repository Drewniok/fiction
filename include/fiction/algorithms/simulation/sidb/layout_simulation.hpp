//
// Created by Jan Drewniok on 21.04.23.
//

#ifndef FICTION_LAYOUT_SIMULATION_HPP
#define FICTION_LAYOUT_SIMULATION_HPP

#include "fiction/algorithms/simulation/sidb/energy_distribution.hpp"
#include "fiction/algorithms/simulation/sidb/exhaustive_ground_state_simulation.hpp"
#include "fiction/algorithms/simulation/sidb/minimum_energy.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/io/write_sqd_layout.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/technology/sidb_defects.hpp"

#include <fmt/format.h>
#include <mockturtle/utils/stopwatch.hpp>

#include <iostream>
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
        this->init();
    }

    void charge_distribution_to_index()
    {
        uint64_t counter = 0;
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

    bool init()
    {
        std::vector<typename Lyt::cell> all_cells = {};
        // obtain all cells in the surfaces and order them by their position to achieve a reproducible output
        layout.foreach_cell([this, &all_cells](const cell<Lyt>& c) { all_cells.push_back(c); });
        std::sort(all_cells.begin(), all_cells.end());

        left_corner_cell = all_cells.front();
        cells            = {};
        total_cells      = {};
        start_cell       = fiction::siqad::to_fiction_coord<cube::coord_t>(left_corner_cell);
        std::cout << "finished" << std::endl;
    }

    bool layout_generation()
    {
        uint64_t length    = 1;
        uint64_t allowed   = start_cell.x + length;
        uint64_t allowed_y = start_cell.y + length;

        cells              = {};
        border_cells       = {};
        border_cell_charge = {};

        //        if (start_cell == siqad::to_fiction_coord<cube::coord_t>(left_corner_cell))
        //        {
        //            cells.insert(siqad::to_siqad_coord(start_cell));
        //            total_cells.insert(siqad::to_siqad_coord(start_cell));
        //        }

        while (total_cells.size() < layout.num_cells() && cells.size() < 25)
        {
            layout.foreach_cell(
                [&allowed, &allowed_y, this](const cell<Lyt>& c)
                {
                    if (std::find(cells.begin(), cells.end(), c) == cells.end())
                    {
                        auto cell_conv = fiction::siqad::to_fiction_coord<cube::coord_t>(c);
                        if (cell_conv.x < start_cell.x && cell_conv.y <= allowed_y && cell_conv.y >= start_cell.y)
                        {
                            cells.insert(c);
                            total_cells.insert(c);
                        }
                        if (cell_conv.x > start_cell.x && cell_conv.y <= allowed_y && cell_conv.x <= allowed)
                        {
                            cells.insert(c);
                            total_cells.insert(c);
                        }
                        if (cell_conv.x == start_cell.x && cell_conv.y <= allowed_y && cell_conv.x <= allowed)
                        {
                            cells.insert(c);
                            total_cells.insert(c);
                        }
                    }
                });
            length += 1;
            allowed   = start_cell.x + length;
            allowed_y = start_cell.y + length;
        }
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
                        if (sidb_nanometer_distance<Lyt>(layout, *it, c, parameter) < 5)
                        {
                            //                            if (c.x > siqad::to_siqad_coord(start_cell).x && c.y >
                            //                            siqad::to_siqad_coord(start_cell).y)
                            //                            {
                            counter += 1;
                            //}
                        }
                    }
                    if (counter != 0)
                    {
                        border_cells.push_back(c);
                        border_cell_charge.push_back(sidb_charge_state::NEUTRAL);
                    }
                }
            });

        border_cell_max_charge_index = std::pow(2, border_cells.size()) - 1;
        start_cell.x                 = start_cell.x + length;
        start_cell.y                 = start_cell.y + length;
    }

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
            defect_cell.emplace_back(border_cells[i]);
            defect_charge.push_back(charge_state_to_sign(border_cell_charge[i]));
        }
    }

    bool run_simulation()
    {
        uint64_t counter = 0;
        while (total_cells.size() < layout.num_cells())
        {
            this->layout_generation();

            Lyt lyt{};

            for (const auto& cell : cells)
            {
                lyt.assign_cell_type(cell, Lyt::cell_type::NORMAL);
            }

            layout_num = lyt.num_cells();

            std::cout << "border cell max index: " << std::to_string(border_cell_max_charge_index) << std::endl;
            std::cout << "border cells: " << std::to_string(border_cells.size()) << std::endl;
            std::set<uint64_t> charge_index{};

            border_cell_index = 0;
            //            while (border_cell_index <= border_cell_max_charge_index)
            //            {
            //                if (border_cell_index % 1000 == 0)
            //                {
            //                    std::cout << border_cell_index << std::endl;
            //                }
            while (border_cell_index <= 1)
            {

                this->index_to_charge_distribution();

                this->defect_map_update();

                //                std::vector<sidb_charge_state> defect_charges{};
                //                for (const auto& charge : defect_charge)
                //                {
                //                    defect_charges.push_back(sign_to_charge_state(charge));
                //                }
                //
                //                std::cout << charge_configuration_to_string(defect_charges) << std::endl;

                std::unordered_map<typename Lyt::cell, const sidb_defect> defect{};
                for (uint64_t i = 0; i < defect_cell.size(); i++)
                {
                    defect.insert({defect_cell[i],
                                   sidb_defect{sidb_defect_type::UNKNOWN, static_cast<double>(defect_charge[i])}});
                }
                // std::cout << "defects: " << defect.size() << std::endl;
                exgs_stats<Lyt> exgs_stats{};
                // std::vector<uint64_t> charge_index_in{};
                // exhaustive_ground_state_simulation(lyt, parameter, &exgs_stats, defect);
                for (const auto& lyt_loop : exgs_stats.valid_lyts)
                {
                    charge_index.insert(lyt_loop.get_charge_index().first);
                    // charge_index_in.emplace_back(charge_distribution_external_to_index(defect_cell, lyt_loop));
                }
                // charge_index_innen.emplace_back(charge_index_in);
                all_charge_lyts.emplace_back(exgs_stats.valid_lyts);
                region_num.emplace_back(counter);

                all_defect_charges.emplace_back(defect_charge);
                border_cell_index += 1;
            }
            all_defect_cells.emplace_back(defect_cell);
            //            write_sqd_layout(lyt, "/Users/jandrewniok/Desktop/investi/" +
            //                                      std::to_string(border_cell_index) + ".sqd");

            std::cout << "number valid lyts: " << charge_index.size() << std::endl;

            //                        if (charge_index.size() == 1)
            //                        {
            //                            write_sqd_layout(lyt, "/Users/jandrewniok/Desktop/investi/" +
            //                            std::to_string(total_cells.size()));
            //                        }

            if (border_cell_max_charge_index == 0)
            {
                exgs_stats<Lyt> exgs_stats{};
                exhaustive_ground_state_simulation(lyt, parameter, &exgs_stats);
                all_charge_lyts.emplace_back(exgs_stats.valid_lyts);
                region_num.emplace_back(counter);
            }
            // std::cout << all_charge_lyts.size() << std::endl;
            counter += 1;
        }

        return true;
    };

    //    bool getting_final_layout()
    //    {
    //        charge_index_innen = {};
    //        auto region_num_unique = region_num;
    //        std::sort(region_num_unique.begin(), region_num_unique.end());
    //        auto last = std::unique(region_num_unique.begin(), region_num_unique.end());
    //        region_num_unique.erase(last, region_num_unique.end());
    //
    //        for (uint64_t i = 1; i < (region_num_unique.size()); i++)
    //        {
    //            for (uint64_t j = 0; j < region_num.size(); j++)
    //            {
    //                if (region_num[j] == i)
    //                {
    //                    charge_index_innen.emplace_back(charge_distribution_external_to_index(all_defect_cells[i],
    //                    all_charge_lyts[i][j]));
    //                }
    //            }
    //        }
    //
    //        return true;
    //    }

  private:
    Lyt&                                                       layout;
    sidb_simulation_parameters                                 parameter{};
    layout_sim_stats<Lyt>&                                     statistic{};
    typename cube::coord_t                                     start_cell{};
    typename Lyt::cell                                         left_corner_cell{};
    std::set<typename Lyt::cell>                               cells{};
    std::vector<typename Lyt::cell>                            border_cells{};
    std::vector<sidb_charge_state>                             border_cell_charge{};
    uint64_t                                                   border_cell_index{};
    uint64_t                                                   border_cell_max_charge_index{};
    std::vector<typename Lyt::cell>                            defect_cell{};
    std::vector<int8_t>                                        defect_charge{};
    std::set<typename Lyt::cell>                               total_cells{};
    uint64_t                                                   layout_num{};
    std::vector<std::vector<charge_distribution_surface<Lyt>>> all_charge_lyts{};
    std::vector<std::vector<uint64_t>>                         charge_index_innen{};
    std::vector<std::vector<typename Lyt::cell>>               all_defect_cells{};
    std::vector<std::vector<int8_t>>                           all_defect_charges{};
    std::vector<uint64_t>                                      region_num{};
};

template <typename Lyt>
bool layout_simulation(Lyt& lyt, const sidb_simulation_parameters& params = sidb_simulation_parameters{},
                       layout_sim_stats<Lyt>* ps = nullptr)
{

    layout_sim_stats<Lyt> st{};

    detail::layout_simulation_impl<Lyt> p{lyt, params, st};

    auto result = p.run_simulation();
    // auto result1 = p.getting_final_layout();

    if (ps)
    {
        *ps = st;
    }

    return result;
}

}  // namespace detail
}  // namespace fiction
#endif  // FICTION_LAYOUT_SIMULATION_HPP
