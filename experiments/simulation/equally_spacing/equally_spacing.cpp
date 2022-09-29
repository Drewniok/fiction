//
// Created by jan-d on 02.09.2022.
//

#include "simu.cpp"
#include "simu.h"


#include <iterator>
#include <fiction/io/read_sqd_layout.hpp>
#include <fiction/technology/area.hpp>
#include <fiction/io/write_sqd_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/types.hpp>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <execution>
#include <iostream>
#include <random>
#include <tuple>
#include <vector>
#include <string>
#include <filesystem>

using namespace fiction;
using std::filesystem::directory_iterator;


// function is used to sort the second and afterwards the first element of the vector
bool sortPairs(std::vector<unsigned long> const& s1, std::vector<unsigned long> const& s2)
{
    if (s1[1] != s2[1])
    {
        return s1[1] < s2[1];
    }
    return s1[0] < s2[0];
}

namespace Circuits
{
    // lattice dimensions
    //const std::string test10       = "C:/Users/jan-d/OneDrive/Desktop/test10.sqd";
    const std::string xor2_gate    = "C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/experiments/layouts/trindade16/xor2.sqd";
    const std::string or_gate_10      = "C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/bestagon-gates/2i1o_or/21_hex_inputsdbp_or_v17_10.sqd";
    const std::string or_gate      = "C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/bestagon-gates/2i1o_or/21_hex_inputsdbp_or_v17_10.sqd";
    const std::string and_gate_new = "C:/Users/jan-d/Downloads/21_hex_inputsdbp_and_v7.sqd";
    const std::string inv_diagonal = "C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/bestagon-gates/1i1o_inv_diag/hex_11_inputsdbp_inv_diag_v0_manual.sqd";
    const std::string and_gate = "C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/bestagon-gates/2i1o_and/21_hex_inputsdbp_and_v19.sqd";
    const std::string hourglass = "C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/bestagon-gates/2i2o_hourglass/22_hex_inputsdbp_hourglass_v0.sqd";
    const std::string  test8 = "C:/Users/jan-d/OneDrive/Desktop/test8.sqd";
    const std::string majority = "C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/experiments/layouts/fontes18/majority.sqd";
    const std::string fo2_gate = "C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/bestagon-gates/1i2o_fo2/12_hex_inputsdbp_fo2_v6.sqd";
    const std::string ha = "C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/bestagon-gates/2i2o_ha/22_hex_inputsdbp_ha_v2.sqd";
    const std::string xor_gate = "C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/bestagon-gates/2i2o_xor/hex_21_inputsdbp_xor_v1_11";
    const std::string HA_gate = "C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/experiments/layouts/trindade16/HA.sqd";
    };

    std::string path = "C:/Users/jan-d/OneDrive/Dokumente/PhD/FCN/sqd/";

int main()
{
    int identical = 0;

//for (const auto &file : directory_iterator(path))
//{

for (int circuit_num=0; circuit_num<1;circuit_num++)
{
//    std::vector<std::vector<unsigned long>> location;
//
//    std::random_device                                       dev;
//    std::mt19937                                             rng(dev());
//    std::uniform_int_distribution<std::mt19937::result_type> dist6(0, 20);
//    std::uniform_int_distribution<std::mt19937::result_type> dist10(0, 15);
//    std::uniform_int_distribution<std::mt19937::result_type> dist1(0, 1);
//
//    int number_of_sidb = 26;
//    // generate random layout
//    std::vector<int> numbers_x;
//    for (int i = 0; i < number_of_sidb; i++)  // add 0-99 to the vector
//        numbers_x.push_back(dist6(rng));
//    unsigned seed_x = std::chrono::system_clock::now().time_since_epoch().count();
//    std::shuffle(numbers_x.begin(), numbers_x.end(), std::default_random_engine(seed_x));
//
//    std::vector<int> numbers_z;
//    for (int i = 1; i < number_of_sidb + 1; i++)  // add 0-99 to the vector
//        numbers_z.push_back(dist10(rng));
//    unsigned seed_z = std::chrono::system_clock::now().time_since_epoch().count();
//    std::shuffle(numbers_z.begin(), numbers_z.end(), std::default_random_engine(seed_z));
//
//    std::vector<int> numbers_y;
//    for (int i = 0; i < number_of_sidb; i++)  // add 0-99 to the vector
//        numbers_y.push_back(dist1(rng));
//    unsigned seed_y = std::chrono::system_clock::now().time_since_epoch().count();
//    std::shuffle(numbers_y.begin(), numbers_y.end(), std::default_random_engine(seed_z));
//
//    std::set<std::vector<unsigned long>> collect_set;
//
//    for (int l = 0; l < numbers_y.size(); l++)
//    {
//        collect_set.insert({static_cast<unsigned long>(numbers_x[l]), static_cast<unsigned long>(numbers_z[l]),
//                            static_cast<unsigned long>(numbers_y[l])});
//    };
//
//    for (auto it = collect_set.begin(); it != collect_set.end(); it++)
//    {
//        location.push_back(*it);
//    }




    std::string selection = Circuits::xor2_gate;
   // std::string selection = path;
   //std::string selection = file.path();
//

            const auto lyt = read_sqd_layout<sidb_cell_clk_lyt>(selection);

            std::vector<std::vector<unsigned long>> location;
            std::vector<cell<sidb_cell_clk_lyt>>    all_cells{};
            all_cells.reserve(lyt.num_cells());

            lyt.foreach_cell([&all_cells](const auto& c) { all_cells.push_back(c); });

            // transform coordinates in new coordinates, inspired by SIQAD
            for (const auto& c : all_cells)
            {
                auto          X = c.x;
                unsigned long Y = ((c.y) - ((c.y) % 2)) * 0.5;
                unsigned long Z = ((c.y) % 2);
                location.push_back({X, Y, Z});
            }

//            sidb_cell_clk_lyt location_to_fiction_layout{{10,10}};
//
//            for (const auto& c : location)
//            {
//                location_to_fiction_layout.assign_cell_type(convert_coordinate(c), sidb_technology::cell_type::NORMAL);
//            }
//
//            write_sqd_layout(location_to_fiction_layout, "filename.sqd");

    //    std::for_each(std::execution::par, all_cells.cbegin(), all_cells.cend(),
    //                  [&location](const auto& c)
    //                  {
    //                      auto          X = c.x;
    //                      unsigned long Y = ((c.y) - ((c.y) % 2)) * 0.5;
    //                      unsigned long Z = ((c.y) % 2);
    //                      location.push_back({X, Y, Z});
    //                  });

    std::cout << std::endl
              << fmt::format("There are {} SiDBs in the circuit", location.size()) << std::endl
              << std::endl;

    // sort vector such that the dots are saved from top to bottom from left to right
    std::sort(location.begin(), location.end(), sortPairs);

    // print all SiDB location

//    for (const auto& it : location)
//    {
//        std::cout << "X: " << it[0] << " | "
//                  << "Y: " << it[1] << " | "
//                  << "Z: " << it[2] << std::endl;
//    }

    std::vector<int> charge;
    std::vector<int> initial_sign(location.size(), -1);

    Energyscr        check(location, initial_sign);
    check.toeuc();
    check.distance();
    int warning = 0;
    for (int row = 0; row < location.size(); row++)
    {
        for (int col = 0; col < location.size(); col++)
        {
            if ((check.db_r(row, col) < 6E-10) && (check.db_r(row, col) != 0))
            {
                warning += 1;
            }
        }
    }

    if (warning > 0)
    {
        continue;
    }

    // std::vector<int> v{0,1,2,4,8,16,32,64,128,256,512};
    // int circuit_num = 1;
//   std::ofstream outFile(selection + std::to_string(circuit_num) +  "_location.txt");
//   std::ofstream chargeFile(selection + std::to_string(circuit_num) + "_charge.txt");


    std::ofstream outFile(selection +  "_location.txt");
    std::ofstream chargeFile(selection + "_charge.txt");

    // the important part
    std::cout << selection << std::endl;

    outFile << "x;"
            << "y;"
            << "z;" << std::endl;
    for (int i = 0; i < location.size(); i++)
    {
        for (const auto& e : location[i]) outFile << std::to_string(e) << ";";
        outFile << std::endl;
    };

    for (int i = 0; i < location.size(); i++)
    {
        chargeFile << std::to_string(i) << ";";
    };
    chargeFile << "type;";
    chargeFile << "runtime;";
    chargeFile << "energy" << std::endl;

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

    // ---------------------------------------------- EXHAUSTIVE SEARCH
    // -------------------------------------------------
    auto        t1        = std::chrono::high_resolution_clock::now();
    float system_energy_min = 100000;
    if (!true)
    {
        std::vector<int> charge;
        std::vector<int> initial_sign(location.size(), -1);
        Energyscr        generalisation(location, initial_sign);
        generalisation.toeuc();                 // transform lattice coordinates to euclidean coordinates
        generalisation.chargeconf_to_index(2);  // chargesign-vector (e.g. -1 0 0 -1 ...) is converted to an index |
                                                // number in () is count state | 3 = DB-,DB0,DB+
        generalisation.distance();              // distance-matrix d(i,j) is calculated, only required once
        generalisation.potentials();            // potential-matrix v(i,j) is calculated, also only required once, since
                                                // chargesign is stil not included
        std::cout << "max: " << generalisation.chargemax << std::endl;  // number of all possible charge distributions
        // can be output

        int counter = 0;
        while (generalisation.chargeindex <= generalisation.chargemax)
        {
            counter += 1;
            if ((counter % 1000000) == 0)
            {
                // std::cout << "iteration step: " << counter << std::endl;
                //generalisation.get_chargesign();
                //std::cout << std::endl;
            }
            generalisation
                .total_energy();  // update of v_local[i] which changes as soon as the charge distribution changes
            if (generalisation.populationValid() &&
                (generalisation.system_energy() <=
                 system_energy_min))  // if the charge state configuration is the energetically best so far
                                      // and physically valid, it is saved and printed out
            {
                system_energy_min = generalisation.system_energy();  //
                std::cout << "system_energy_min: " << system_energy_min
                          << std::endl;              // output when new best mark is reached
                charge = generalisation.chargesign;  // save new best chargesign
            }
            generalisation.increase_step();         // index is increased by one
            generalisation.index_to_chargeconf(2);  // index is converted to the chargesign-vector
        }
        std::cout << std::endl
                  << "-------------------------------------------------------------------------------" << std::endl;
        // after all charge configurations are evaluated, the energetically lowest charge which is physically valid at
        // the same time is printed
        std::cout << "smallest energy found: " << system_energy_min << std::endl;
        std::cout << "charge configuration: |";
        for (auto it = charge.begin(); it != charge.end(); it++)
        {
            std::cout << (*it) << " | ";
            chargeFile << std::to_string(*it) << ";";
        }
        chargeFile << "EXGS;";
        std::cout << std::endl
                  << "-------------------------------------------------------------------------------" << std::endl;
        const auto                    t11        = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff_first = t11 - t1;
        std::cout << std::endl
                  << std::endl
                  << fmt::format("It took {} s to simulate the {} circuit", diff_first.count(), selection) << std::endl;
        chargeFile << std::to_string(diff_first.count()) << ";";
        chargeFile << std::to_string(system_energy_min);
    }
    // -------------------------------------------------------------------------------------------------------------------------------

    // ---------------------------------------------- MAX-MIN DIVERSITY APPROACH
    // -------------------------------------------------
    float system_energy = MAX_FLOAT;

    auto t12 = std::chrono::high_resolution_clock::now();
    if (true)
    {
        //std::vector<int> initial_sign(location.size(), 0);
        //and
        //std::vector<int> initial_sign = {-1,-1,0,0,-1,-1,0,0,-1,-1,0,0,0,-1,-1,0,-1,0,-1,0,-1,0,-1};

        //hourglass
        //std::vector<int> initial_sign = {-1,-1,0,0,-1,-1,0,0,0,-1,-1,0,-1,0,0,-1,0,0,-1,-1,0,0,-1,-1,0,0,-1,-1,-1};
        Energyscr        generalisation(location, initial_sign);
        generalisation.toeuc();
        generalisation.distance();
        generalisation.potentials();
        generalisation.total_energy();
        std::cout << std::endl
                  << "system physically valid [yes/no: 1/0] at the beginning: " << generalisation.populationValid()
                  << std::endl;
        std::random_device              rd;         // obtain a random number from hardware
        std::mt19937                    gen(rd());  // seed the generator
        std::uniform_int_distribution<> distr(
            0, initial_sign.size() - 1);  // define the range, we want to select certain "start" SiDB randomly

        std::vector<int> charge_config;

        static const std::vector<bool> iterator_helper(100, false);
        // this loop executes in parallel if the compiler so wishes
        std::vector<std::vector<float>> collection_all;

        // std::for_each(std::execution::par, iterator_helper.cbegin(), iterator_helper.cend(), [&](const auto& b){
        int threashold_num = 1000;
        int count_equal_energy = 0;
        for (int z = 0; z < threashold_num; z++)
        {

            for (int i = 0; i < location.size(); i++)
            {
                // int i = distr(gen);  // at the beginning, all SiDBs are set to the charge 0. Now, one SiDB is chosen
                //  randomly and set to -1
                //std::cout << z << std::endl;
                std::vector<int> initial_sign(location.size(), 0);
                initial_sign[i - 1] = 0;
                initial_sign[i]     = -1;
                Energyscr first_try(location, initial_sign);
                // -------------- some attributes do not change when charge sign is changed --------------
                first_try.locationeuc = generalisation.locationeuc;
                first_try.db_r        = generalisation.db_r;
                first_try.v_ij        = generalisation.v_ij;
                // ---------------------------------------------------------------------------------------
                first_try.total_energy();  // However, v_local[i] has to recalculated

                //            auto first_try = generalisation;
                //            first_try = ;
                std::vector<int> index_start      = {i};  // index_start includes all the indices of "-1"-SiDBs
                auto             index_start_size = index_start.size();
                auto             index_end_size   = index_start.size();

                for (int set_dbs = 0; set_dbs < location.size() - 1; set_dbs++)
                {
                    index_start_size = index_start.size();
                    // based on the max-min diversity problem, all -1 SiDBs are selected iteratively
                    first_try.find_new_best_neighbor_GRC(index_start);
                    // first_try.get_chargesign();
                    index_end_size = index_start.size();
                    // auto index_size = index_start.size();
                    // std::cout <<  "index_size:" << index_start_size << std::endl;
                    // std::cout <<  "index_end_size:" << index_end_size << std::endl;
                    //first_try.get_chargesign();
                    first_try.total_energy();
                    // first_try.get_chargesign();
                    //   check if new charge configuration with lower energy is found

                    if (first_try.populationValid() && first_try.system_energy() <= system_energy)
                    {
                        if(first_try.system_energy() == system_energy)
                        {
                            count_equal_energy +=1;
                        }
                        std::cout << "system energy old: " << first_try.system_energy() << std::endl;
                        //first_try.get_chargesign();
                        system_energy = first_try.system_energy();
                        charge_config = first_try.chargesign;
                        break;
                    }
                    // populationValid_counter() can let an electron jump from a more negatively charged SiDB to a more
                    // positively charged one

//                    std::tuple<int, std::vector<int>, int> output2 = first_try.populationValid_counter();
//                    first_try.total_energy();  // v_local[i] are updated
//                    if ((first_try.populationValid()) && (first_try.system_energy() <= system_energy))
//                    {
//                        if(first_try.system_energy() == system_energy)
//                        {
//                            count_equal_energy +=1;
//                        }
//                        std::cout << "system energy new: " << first_try.system_energy() << std::endl;
//                        first_try.get_chargesign();
//                        system_energy = first_try.system_energy();
//                        charge_config = first_try.chargesign;
//                        break;
//                    };
                    // first_try.get_chargesign();
                    collection_all.push_back(
                        {first_try.system_energy(), static_cast<float>(first_try.populationValid())});

                };

            };
            if (count_equal_energy > 10)
            {break;}
        }
        std::cout << "________________________________" << std::endl;
        //
        std::cout << "smallest energy: " << system_energy << std::endl;
        chargeFile << std::endl;
        for (const auto& it : charge_config)
        {
            std::cout << it << " | ";
            chargeFile << std::to_string(it) << ";";
        }
        chargeFile << "EQ;";


        //        for (const auto& it : collection_all)
        //        {
        //            if (it[1] == 1)
        //            {std::cout << it[0] << " | " << it[1] << std::endl;}
        //        }
    }
    const auto t2     = std::chrono::high_resolution_clock::now();
    const auto ms_int = std::chrono::duration_cast<std::chrono::seconds>(t2 - t12);

    if (system_energy_min == system_energy)
    {
        identical += 1;
    }
    std::chrono::duration<double> diff = t2 - t12;
    chargeFile << std::to_string(diff.count());
    chargeFile << ";" << std::to_string(system_energy);
    std::cout << std::endl
              << std::endl
              << fmt::format("It took {} s to simulate the {} circuit", diff.count(), selection) << std::endl;
    //}

    std::cout << "identical: " << identical << std::endl;
}
return 0;
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
//             // as soon as there is no new SiDB which can be switched to -1 due to a suitable distance, the chargesign
//             and the energy of the certain charge configuration is returned
//
//             //while (!isEqual(new_index1, index_start));
//             // first_try.get_chargesign();
//             first_try.total_energy();
//             // first_try.system_energy();
//             // std::cout << "system stable?: " << first_try.populationValid() << " | counter: " <<
//             first_try.populationValid_counter().first << std::endl; std::cout << "system stable?: " <<
//             first_try.populationValid() << " | counter: " << first_try.populationValid_counter().first << std::endl;
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
