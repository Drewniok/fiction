//
// Created by jan-d on 02.09.2022.
//

#include <fiction/io/read_sqd_layout.hpp>
#include <fiction/technology/area.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/types.hpp>
#include <iostream>

#include <vector>
#include <math.h>
#include <tuple>
#include "simu.h"
#include "simu.cpp"
#include <random>

using namespace fiction;

template<typename T>
bool isEqual(std::vector<T> &first, std::vector<T> &second)
{
    if (first.size() != second.size()) {
        return false;
    }
    return std::is_permutation(first.begin(), first.end(), second.begin());
}

// function is used to sort the second and afterwards the first element of the vector
bool sortPairs(std::vector<unsigned long> const& s1, std::vector<unsigned long> const& s2)
{
    if ( s1[1] != s2[1] )
    {
        return s1[1] < s2[1];
    }
    return s1[0] < s2[0];
}


struct Circuits {
    // lattice dimensions
    std::string test10 = "C:/Users/jan-d/OneDrive/Desktop/test10.sqd";
    std::string xor2 = "C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/experiments/layouts/trindade16/xor2.sqd";
};

struct Circuits circuits;

int main()
{
    const auto lyt = read_sqd_layout<sidb_cell_clk_lyt>(circuits.test10);


    std::vector<std::vector<unsigned long>> location;
    std::vector<cell<sidb_cell_clk_lyt>>    all_cells{};
    all_cells.reserve(lyt.num_cells());

    lyt.foreach_cell([&lyt, &all_cells](const auto& c) { all_cells.push_back(c); });

    // transform coordinates in new coordinates, inspired by SIQAD
    for (int i = 0; i < all_cells.size(); i++)
    {
        auto          X = all_cells[i].x;
        unsigned long Y = ((all_cells[i].y) - ((all_cells[i].y) % 2)) * 0.5;
        unsigned long Z = ((all_cells[i].y) % 2);
        location.push_back({X, Y, Z});
    };

    // sort vector such that the dots are saved from top to bottom from left to right
    std::sort(location.begin(), location.end(), sortPairs);

    // print all SiDB location
    for (auto it = location.begin(); it != location.end(); it++)
    {
        std::cout << "X: " << (*it)[0] << " | "
                  << "Y: " << (*it)[1] << " | "
                  << "Z: " << (*it)[2] << std::endl;
    }

    // example of how Siqad works
    // ____________________________________________
    //    cell<sidb_cell_clk_lyt> c{89, 62};
    //
    //    std::cout << fmt::format("Position {} is {}a primary input", c, lyt.is_pi(c) ? "" : "not ") << std::endl;
    //
    //    area_params<sidb_technology> ps{};
    //    area_stats                   st{};
    //    area(lyt, ps, &st);
    //    std::cout << fmt::format("The layout has an area of {} nm²", st.area) << std::endl;
    // ___________________________________________

    // ---------------------------------------------- EXHAUSTIVE SEARCH -------------------------------------------------
    if (1 == 0)
    {
        std::vector<int> charge;
        std::vector<int> initial_sign(location.size(), -1);
        Energyscr        generalisation(location, initial_sign);
        generalisation.toeuc();        // transform lattice coordinates to euclidean coordinates
        generalisation.chargeconf_to_index(3);  // chargesign-vector (e.g. -1 0 0 -1 ...) is converted to an index | number in () is count state | 3 = DB-,DB0,DB+
        generalisation.distance();     // distance-matrix d(i,j) is calculated, only required once
        generalisation.potentials();   // potential-matrix v(i,j) is calculated, also only required once, since chargesign is stil not included
        // std::cout << "max: " << generalisation.chargemax << std::endl;  number of all possible charge distributions can be output
        float system_energy_min = 100000;
        int   counter           = 0;
        while (generalisation.chargeindex <= generalisation.chargemax)
        {
            counter += 1;
            if ((counter % 100) == 0)
            {
                std::cout << "iteration step: " << counter << std::endl;
                generalisation.get_chargesign();
                std::cout << std::endl;
            }
            generalisation
                .total_energy();  // update of v_local[i] which changes as soon as the charge distribution changes
            if ((generalisation.populationValid() == 1) &&
                (generalisation.system_energy() <=
                 system_energy_min))  // if the charge state configuration is the energetically best so far
                                      // and physically valid, it is saved and printed out
            {
                system_energy_min = generalisation.system_energy();  //
                std::cout << "system_energy_min: " << system_energy_min
                          << std::endl;              // output when new best mark is reached
                charge = generalisation.chargesign;  // save new best chargesign
            }
            generalisation.increase_step();  // index is increased by one
            generalisation.index_to_chargeconf(3);    // index is converted to the chargesign-vector
        }
        std::cout << std::endl
                  << "-------------------------------------------------------------------------------" << std::endl;
        // after all charge configurations are evaluated, the energetically lowest charge which is physically valid at the same time is printed
        std::cout << "smallest energy found: " << system_energy_min << std::endl;
        std::cout << "charge configuration: |";
        for (auto it = charge.begin(); it != charge.end(); it++)
        {
            std::cout << (*it) << " | ";
        };
        std::cout << std::endl
                  << "-------------------------------------------------------------------------------" << std::endl;
    };
    // -------------------------------------------------------------------------------------------------------------------------------

    // ---------------------------------------------- MAX-MIN DIVERSITY APPROACH -------------------------------------------------
    if (1 == 1)
    {
        std::vector<int> initial_sign(location.size(), 0);
        Energyscr        generalisation(location, initial_sign);
        generalisation.toeuc();
        generalisation.distance();
        generalisation.potentials();
        generalisation.total_energy();
        std::cout << "system physically valid [yes/no: 1/0]: " << generalisation.populationValid() << std::endl;
        std::random_device              rd;         // obtain a random number from hardware
        std::mt19937                    gen(rd());  // seed the generator
        std::uniform_int_distribution<> distr(
            0, initial_sign.size() - 1);  // define the range, we want to select certain "start" SiDB randomly

        float            system_energy = 100000;
        std::vector<int> charge_config;

        for (int z = 0; z < 1000; z++)
        {
            for (int l = 0; l < 30; l++)
            {
                int i = distr(gen);  // at the beginning, all SiDBs are set to the charge 0. Now, one SiDB is chosen randomly and set to -1
                initial_sign[i - 1] = 0;
                initial_sign[i]     = -1;
                Energyscr first_try(location, initial_sign);
                // -------------- some attributes do not change when charge sign is changed --------------
                first_try.locationeuc = generalisation.locationeuc;
                first_try.db_r        = generalisation.db_r;
                first_try.v_ij        = generalisation.v_ij;
                // ---------------------------------------------------------------------------------------
                first_try.total_energy();  // However, v_local[i] has to recalculated

                std::vector<int> index_start = {i};  // index_start includes all the indices of "-1"-SiDBs
                std::vector<int> index_start_output;
                do {
                    // based on the max-min diversity problem, all -1 SiDBs are selected iteratively
                    index_start_output = first_try.find_new_best_neighbor_GRC(index_start);
                    index_start        = index_start_output;
                } while (!isEqual(index_start, index_start_output));
                first_try.total_energy();

                // check if new charge configuration with lower energy is found
                if (first_try.populationValid() == 1 && first_try.system_energy() <= system_energy)
                {
                    std::cout << "system energy old: " << first_try.system_energy() << std::endl;
                    first_try.get_chargesign();
                    system_energy = first_try.system_energy();
                    charge_config = first_try.chargesign;
                }
                // populationValid_counter() can let an electron jump from a more negatively charged SiDB to a more positively charged one
                std::tuple<int, std::vector<int>, int> output2 = first_try.populationValid_counter();
                first_try.total_energy(); // v_local[i] are updated

                // check if new generated charge configuration is physically valid and stable
                if (first_try.populationValid() == 1 && first_try.system_energy() < system_energy)
                {
                    std::cout << "system energy new: " << first_try.system_energy() << std::endl;
                    first_try.get_chargesign();
                    system_energy = first_try.system_energy();
                    charge_config = first_try.chargesign;
                };
            }
        }

        std::cout << "________________________________" << std::endl;
        //
        std::cout << "smallest energy: " << system_energy << std::endl;

        for (auto it = charge_config.begin(); it != charge_config.end(); it++)
        {
            std::cout << (*it) << " | ";
        }
        return 0;
    };
}


    // if (1==0)
    //{
    //
    //     std::vector<int> initial_sign(location.size(), 0);
    //     // std::vector<int> initial_sign = {-1,0,-1,0,-1,0,0,-1,0,-1};
    //     Energyscr generalisation(location, initial_sign);
    //     generalisation.toeuc();
    //     generalisation.distance();
    //     generalisation.potentials();
    //     generalisation.total_energy();
    //     std::vector<int> pertuber_vector = generalisation.find_perturber_alternative();
    //
    //     if (print_perturber_location == 1)
    //     {
    //         for (auto it = pertuber_vector.begin(); it != pertuber_vector.end(); it++)
    //         {
    //             std::cout << (*it) << std::endl;
    //         }
    //     };
    //
    //     // i and j are two indeces of SiDBs which are the start. Distance between these two os calculated.
    //     float system_energy = 100000;
    //     std::vector<int> charge_config;
    //
    //     for (int i = 0; i < initial_sign.size(); i++)
    //     {
    //         // initial_sign.size()
    //         // int bound = ((i < (initial_sign.size()-50)) ? (i+49) : initial_sign.size());
    //         for (int j = i + 1; j < initial_sign.size(); j++)
    //         {
    //             if (generalisation.db_r(i, j) > 10 * std::pow(10, -9))
    //             {
    //                 break;
    //             }
    //
    //             Energyscr first_try(location, initial_sign);
    //             first_try.change_chargesign(i, j);
    //
    //             first_try.locationeuc = generalisation.locationeuc;
    //             first_try.db_r        = generalisation.db_r;
    //             first_try.v_ij        = generalisation.v_ij;
    //
    //             first_try.get_distance(i, j);
    //             first_try.get_potential(i, j);
    //             std::vector<int> index_start = {i, j};
    //             // std::vector<int> perturber = first_try.find_perturber_alternative();
    //
    //             // i,j SiDBs and perturber are kept on charge state -1
    //             index_start.insert(index_start.end(), pertuber_vector.begin(), pertuber_vector.end());
    //             sort(index_start.begin(), index_start.end());
    //
    //             // ensure that each index is only on time in vector
    //             index_start.erase(unique(index_start.begin(), index_start.end()), index_start.end());
    //             //std::cout << "j: " << j << std::endl;
    //
    //             std::vector<int> new_index1;
    //             for (int i =0; i<3;i++)
    //             {
    //                 new_index1                 = index_start;
    //                 std::vector<int> new_index = first_try.search_same_distance(index_start);
    //                 index_start                = new_index;
    //
    //                 for (auto it = index_start.begin(); it != index_start.end(); it++)
    //                 {
    //
    //                     std::cout << (*it) << std::endl;
    //                     first_try.chargesign[(*it)] = -1;
    //                 }
    //
    //             }
    //             first_try.get_chargesign();
    //             // as soon as there is no new SiDB which can be switched to -1 due to a suitable distance, the chargesign and the energy of the certain charge configuration is returned
    //
    //             //while (!isEqual(new_index1, index_start));
    //             // first_try.get_chargesign();
    //             first_try.total_energy();
    //             // first_try.system_energy();
    //             // std::cout << "system stable?: " << first_try.populationValid() << " | counter: " << first_try.populationValid_counter().first << std::endl; std::cout << "system stable?: " << first_try.populationValid() << " | counter: " << first_try.populationValid_counter().first << std::endl;
    //
    //             // In cases with less than 6 physically invalide state, -v_local[i] + µ is printed to get a feeling.
    //
    //                 //std::vector<int> stability = first_try.populationValid_counter().second;
    //                 //                for (auto it = stability.begin(); it!=stability.end(); it++)
    //                 //                {
    //                 //                    std::cout << "position: " << (*it) << " | ";
    //                 //
    //                 //                }
    //                 std::cout << "correction: " << first_try.populationValid() << std::endl;
    //                 first_try.get_chargesign();
    //                 std::cout << "system energy: " << first_try.system_energy() << std::endl;
    //                 std::pair<int, std::vector<int>> output = first_try.populationValid_counter();
    //                 first_try.total_energy();
    //                 std::vector<int> stability = output.second;
    //                 std::cout << "Valid mistakes: " << output.first << std::endl;
    //                 std::cout << "system energy: " << first_try.system_energy() << std::endl;
    //                 std::cout << "correction: " << first_try.populationValid() << std::endl;
    //                 first_try.get_chargesign();
    //                 std::cout << "i: " << i << std::endl;
    //                 std::cout << "j: " << j << std::endl;
    //                 Energyscr second_try(location, first_try.chargesign);
    //                 second_try.locationeuc = generalisation.locationeuc;
    //                 second_try.db_r        = generalisation.db_r;
    //                 second_try.v_ij        = generalisation.v_ij;
    //                 second_try.change_chargesign_one(stability);
    //                 second_try.total_energy();
    //                 second_try.get_chargesign();
    //                 std::cout << "correction_new: " << second_try.populationValid() << std::endl;
    //                 std::pair<int, std::vector<int>> output2 = second_try.populationValid_counter();
    //                 std::vector<int> stability_new = output2.second;
    //                 second_try.total_energy();
    //                 std::cout << "Valid mistakes_new: " << output2.first << std::endl;
    //                 std::cout << "correction_new: " << second_try.populationValid() << std::endl;
    //                 std::cout << "Valid mistakes_new: " << output2.first << std::endl;
    //                 std::cout << "system energy_new: " << second_try.system_energy() << std::endl;
    //
    //                 //                for (auto it = stability_new.begin(); it!=stability_new.end(); it++)
    //                 //                {
    //                 //                    std::cout << "position_new: " << (*it) << " | ";
    //                 //                }
    //                 //                std::cout << "________________________________" << std::endl;
    //
    //                 Energyscr third_try(location, second_try.chargesign);
    //                 third_try.locationeuc = generalisation.locationeuc;
    //                 third_try.db_r        = generalisation.db_r;
    //                 third_try.v_ij        = generalisation.v_ij;
    //                 third_try.change_chargesign_one(stability_new);
    //                 third_try.total_energy();
    //                 third_try.get_chargesign();
    //                 std::cout << "correction 3: " << third_try.populationValid() << std::endl;
    //                 std::pair<int, std::vector<int>> output3 = third_try.populationValid_counter();
    //                 std::vector<int> stability_new3 = output3.second;
    //                 third_try.total_energy();
    //                 std::cout << "Valid mistakes 3: " << output3.first << std::endl;
    //                 std::cout << "correction 3: " << third_try.populationValid() << std::endl;
    //                 std::cout << "Valid mistakes_new_third: " << output3.first << std::endl;
    //                 std::cout << "system energy_new_third: " << third_try.system_energy() << std::endl;
    //
    //                 Energyscr fourth_try(location, third_try.chargesign);
    //                 fourth_try.locationeuc = generalisation.locationeuc;
    //                 fourth_try.db_r        = generalisation.db_r;
    //                 fourth_try.v_ij        = generalisation.v_ij;
    //                 fourth_try.change_chargesign_one(stability_new3);
    //                 fourth_try.total_energy();
    //                 fourth_try.get_chargesign();
    //                 std::cout << "correction 4: " << fourth_try.populationValid() << std::endl;
    //                 std::pair<int, std::vector<int>> output4 = fourth_try.populationValid_counter();
    //                 std::vector<int> stability_new4 = output4.second;
    //                 fourth_try.total_energy();
    //                 std::cout << "Valid mistakes 4: " << output4.first << std::endl;
    //                 std::cout << "correction 4: " << fourth_try.populationValid() << std::endl;
    //                 std::cout << "Valid mistakes 4: " << output4.first << std::endl;
    //                 std::cout << "system energ 4: " << fourth_try.system_energy() << std::endl;
    //
    //                 //std::vector<int> stability_new4 = fourth_try.populationValid_counter().second;
    //                 //                for (auto it = stability_new3.begin(); it!=stability_new3.end(); it++)
    //                 //                {
    //                 //                    std::cout << "position_new_third: " << (*it) << " | ";
    //                 //                }
    //                 std::cout << "________________________________" << std::endl;
    //
    //             if ((first_try.populationValid() == 1) && (first_try.system_energy() < system_energy))
    //             {
    //                 system_energy = first_try.system_energy();
    //                 charge_config = first_try.chargesign;
    //             }
    //         }
    //     }
    //
    //     std::cout << "smallest energy: " << system_energy << std::endl;
    //
    //     for (auto it = charge_config.begin(); it != charge_config.end(); it++)
    //     {
    //         std::cout << (*it) << " | ";
    //     }
    //
    // }

