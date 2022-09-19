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
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


using namespace fiction;

template<typename T>
bool isEqual(std::vector<T> &first, std::vector<T> &second)
{
    if (first.size() != second.size()) {
        return false;
    }
    return std::is_permutation(first.begin(), first.end(), second.begin());
}

// function is used to sort the secodn and afterwards the first element of the vector
bool sortPairs(std::vector<unsigned long> const& s1, std::vector<unsigned long> const& s2)
{
    if ( s1[1] != s2[1] )
    {
        return s1[1] < s2[1];
    }
    return s1[0] < s2[0];
}

// set variable if perturber location should be printed on time
int print_perturber_location = 1;

int main()
{
    const auto lyt = read_sqd_layout<sidb_cell_clk_lyt>("C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/bestagon-gates/2i1o_xor/hex_21_inputsdbp_xor_v1.sqd");
    //const auto lyt = read_sqd_layout<sidb_cell_clk_lyt>("C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/experiments/layouts/trindade16/HA.sqd");
  // "C:\Users\jan-d\Downloads\sidb-bestagon-gate-library-xml-fix\sidb-bestagon-gate-library-xml-fix\experiments\layouts\trindade16\xor2.sqd"
   // "C:\Users\jan-d\Downloads\sidb-bestagon-gate-library-xml-fix\sidb-bestagon-gate-library-xml-fix\bestagon-gates\2i1o_xor\hex_21_inputsdbp_xor_v1_11.sqd"
    //"C:\Users\jan-d\Downloads\sidb-bestagon-gate-library-xml-fix\sidb-bestagon-gate-library-xml-fix\bestagon-gates\1i2o_fo2\12_hex_inputsdbp_fo2_v6.sqd"
    //"C:/Users/jan-d/OneDrive/Desktop/test9.sqd"
    //"C:\Users\jan-d\Downloads\sidb-bestagon-gate-library-xml-fix\si
    //const auto lyt = read_sqd_layout<sidb_cell_clk_lyt>("C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/bestagon-gates/1i1o_inv_diag/hex_11_inputsdbp_inv_diag_v0_manual.sqd");
    //"C:\Users\jan-d\Downloads\sidb-bestagon-gate-library-xml-fix\sidb-bestagon-gate-library-xml-fix\bestagon-gates\1i1o_inv_diag\hex_11_inputsdbp_inv_diag_v0_manual.sqd"
    //const auto lyt = read_sqd_layout<sidb_cell_clk_lyt>("C:/Users/jan-d/OneDrive/Desktop/test3.sqd");
    std::vector<std::vector<unsigned long>> location;
    std::vector<cell<sidb_cell_clk_lyt>> all_cells{};
    all_cells.reserve(lyt.num_cells());

    lyt.foreach_cell(
        [&lyt, &all_cells](const auto& c)
        {
            all_cells.push_back(c);
        });

    // transforming coordinates in new coordinated, inspired by SIQAD
    for (int i =0; i<all_cells.size(); i++)
    {
        auto X = all_cells[i].x;
        unsigned long Y = ((all_cells[i].y) - ((all_cells[i].y) % 2))*0.5;
        unsigned long Z = ((all_cells[i].y) % 2);
       location.push_back({X,Y,Z});
    };

    // sort vector such that the dots are saved from top to bottom from left to right
    std::sort(location.begin(), location.end(), sortPairs);

    // print all SiDB location
    for (auto it = location.begin(); it != location.end();it++)
    {
        std::cout << "X: " << (*it)[0] << " | " << "Y: " << (*it)[1] << " | " << "Z: " << (*it)[2] << std::endl;
    }






//    cell<sidb_cell_clk_lyt> c{89, 62};
//
//    std::cout << fmt::format("Position {} is {}a primary input", c, lyt.is_pi(c) ? "" : "not ") << std::endl;
//
//    area_params<sidb_technology> ps{};
//    area_stats                   st{};
//    area(lyt, ps, &st);
//    std::cout << fmt::format("The layout has an area of {} nm²", st.area) << std::endl;

//if (1==0)
//{
//    std::vector<int> charge;
//    std::vector<int> initial_sign(location.size(), -1);
//    // std::vector<int> initial_sign = {-1,0,-1,0,-1,0,0,-1,0,-1};
//    Energyscr generalisation(location, initial_sign);
//    generalisation.toeuc();
//    generalisation.euctoindex(2);
//    generalisation.get_index();
//    generalisation.distance();
//    generalisation.potentials();
//    std::cout << "max: " << generalisation.chargemax << std::endl;
//    float system_energy_min = 100000;
//
//    while (generalisation.chargeindex > 0)
//    {
//
//        // generalisation.indextoeuc(2);
//        generalisation.total_energy();
//        generalisation.get_chargesign();
//        std::cout << "correction: " << generalisation.populationValid() << std::endl;
//        std::cout << "system energy: " << generalisation.system_energy() << std::endl;
//
//        std::cout << "Valid mistakes: " << generalisation.populationValid_counter().first << std::endl;
//        generalisation.total_energy();
//        std::cout << "system energy: " << generalisation.system_energy() << std::endl;
//        std::cout << "correction: " << generalisation.populationValid() << std::endl;
//        generalisation.get_chargesign();
//        if ((1 == generalisation.populationValid()) && (generalisation.system_energy() < system_energy_min))
//        {
//            system_energy_min = generalisation.populationValid()*generalisation.system_energy();
//            charge = generalisation.chargesign;
//        }
//        //std::cout << "correction: " << generalisation.populationValid() << std::endl;
//
//        generalisation.increase_step();
//        //generalisation.get_index();
//        // generalisation.euctoindex(2);
//        generalisation.indextoeuc(2);
//        std::cout << "_______________________________________" << std::endl;
//        // generalisation.get_chargesign();
//
//        // generalisation.euctoindex(2);
//    }
//    std::cout << " system_energy_min: " << system_energy_min << std::endl;
//                    for (auto it = charge.begin(); it!=charge.end(); it++)
//                    {
//                        std::cout << (*it) << " | ";
//                    }
//    std::cout << std::endl << "________________________________" << std::endl;
//
//}




//if (1==0)
//{
//
//    std::vector<int> initial_sign(location.size(), 0);
//    // std::vector<int> initial_sign = {-1,0,-1,0,-1,0,0,-1,0,-1};
//    Energyscr generalisation(location, initial_sign);
//    generalisation.toeuc();
//    generalisation.distance();
//    generalisation.potentials();
//    generalisation.total_energy();
//    std::vector<int> pertuber_vector = generalisation.find_perturber_alternative();
//
//    if (print_perturber_location == 1)
//    {
//        for (auto it = pertuber_vector.begin(); it != pertuber_vector.end(); it++)
//        {
//            std::cout << (*it) << std::endl;
//        }
//    };
//
//    // i and j are two indeces of SiDBs which are the start. Distance between these two os calculated.
//    float system_energy = 100000;
//    std::vector<int> charge_config;
//
//    for (int i = 0; i < initial_sign.size(); i++)
//    {
//        // initial_sign.size()
//        // int bound = ((i < (initial_sign.size()-50)) ? (i+49) : initial_sign.size());
//        for (int j = i + 1; j < initial_sign.size(); j++)
//        {
//            if (generalisation.db_r(i, j) > 10 * std::pow(10, -9))
//            {
//                break;
//            }
//
//            Energyscr first_try(location, initial_sign);
//            first_try.change_chargesign(i, j);
//
//            first_try.locationeuc = generalisation.locationeuc;
//            first_try.db_r        = generalisation.db_r;
//            first_try.v_ij        = generalisation.v_ij;
//
//            first_try.get_distance(i, j);
//            first_try.get_potential(i, j);
//            std::vector<int> index_start = {i, j};
//            // std::vector<int> perturber = first_try.find_perturber_alternative();
//
//            // i,j SiDBs and perturber are kept on charge state -1
//            index_start.insert(index_start.end(), pertuber_vector.begin(), pertuber_vector.end());
//            sort(index_start.begin(), index_start.end());
//
//            // ensure that each index is only on time in vector
//            index_start.erase(unique(index_start.begin(), index_start.end()), index_start.end());
//            //std::cout << "j: " << j << std::endl;
//
//            std::vector<int> new_index1;
//            for (int i =0; i<3;i++)
//            {
//                new_index1                 = index_start;
//                std::vector<int> new_index = first_try.search_same_distance(index_start);
//                index_start                = new_index;
//
//                for (auto it = index_start.begin(); it != index_start.end(); it++)
//                {
//
//                    std::cout << (*it) << std::endl;
//                    first_try.chargesign[(*it)] = -1;
//                }
//
//            }
//            first_try.get_chargesign();
//            // as soon as there is no new SiDB which can be switched to -1 due to a suitable distance, the chargesign and the energy of the certain charge configuration is returned
//
//            //while (!isEqual(new_index1, index_start));
//            // first_try.get_chargesign();
//            first_try.total_energy();
//            // first_try.system_energy();
//            // std::cout << "system stable?: " << first_try.populationValid() << " | counter: " << first_try.populationValid_counter().first << std::endl; std::cout << "system stable?: " << first_try.populationValid() << " | counter: " << first_try.populationValid_counter().first << std::endl;
//
//            // In cases with less than 6 physically invalide state, -v_local[i] + µ is printed to get a feeling.
//
//                //std::vector<int> stability = first_try.populationValid_counter().second;
//                //                for (auto it = stability.begin(); it!=stability.end(); it++)
//                //                {
//                //                    std::cout << "position: " << (*it) << " | ";
//                //
//                //                }
//                std::cout << "correction: " << first_try.populationValid() << std::endl;
//                first_try.get_chargesign();
//                std::cout << "system energy: " << first_try.system_energy() << std::endl;
//                std::pair<int, std::vector<int>> output = first_try.populationValid_counter();
//                first_try.total_energy();
//                std::vector<int> stability = output.second;
//                std::cout << "Valid mistakes: " << output.first << std::endl;
//                std::cout << "system energy: " << first_try.system_energy() << std::endl;
//                std::cout << "correction: " << first_try.populationValid() << std::endl;
//                first_try.get_chargesign();
//                std::cout << "i: " << i << std::endl;
//                std::cout << "j: " << j << std::endl;
//                Energyscr second_try(location, first_try.chargesign);
//                second_try.locationeuc = generalisation.locationeuc;
//                second_try.db_r        = generalisation.db_r;
//                second_try.v_ij        = generalisation.v_ij;
//                second_try.change_chargesign_one(stability);
//                second_try.total_energy();
//                second_try.get_chargesign();
//                std::cout << "correction_new: " << second_try.populationValid() << std::endl;
//                std::pair<int, std::vector<int>> output2 = second_try.populationValid_counter();
//                std::vector<int> stability_new = output2.second;
//                second_try.total_energy();
//                std::cout << "Valid mistakes_new: " << output2.first << std::endl;
//                std::cout << "correction_new: " << second_try.populationValid() << std::endl;
//                std::cout << "Valid mistakes_new: " << output2.first << std::endl;
//                std::cout << "system energy_new: " << second_try.system_energy() << std::endl;
//
//                //                for (auto it = stability_new.begin(); it!=stability_new.end(); it++)
//                //                {
//                //                    std::cout << "position_new: " << (*it) << " | ";
//                //                }
//                //                std::cout << "________________________________" << std::endl;
//
//                Energyscr third_try(location, second_try.chargesign);
//                third_try.locationeuc = generalisation.locationeuc;
//                third_try.db_r        = generalisation.db_r;
//                third_try.v_ij        = generalisation.v_ij;
//                third_try.change_chargesign_one(stability_new);
//                third_try.total_energy();
//                third_try.get_chargesign();
//                std::cout << "correction 3: " << third_try.populationValid() << std::endl;
//                std::pair<int, std::vector<int>> output3 = third_try.populationValid_counter();
//                std::vector<int> stability_new3 = output3.second;
//                third_try.total_energy();
//                std::cout << "Valid mistakes 3: " << output3.first << std::endl;
//                std::cout << "correction 3: " << third_try.populationValid() << std::endl;
//                std::cout << "Valid mistakes_new_third: " << output3.first << std::endl;
//                std::cout << "system energy_new_third: " << third_try.system_energy() << std::endl;
//
//                Energyscr fourth_try(location, third_try.chargesign);
//                fourth_try.locationeuc = generalisation.locationeuc;
//                fourth_try.db_r        = generalisation.db_r;
//                fourth_try.v_ij        = generalisation.v_ij;
//                fourth_try.change_chargesign_one(stability_new3);
//                fourth_try.total_energy();
//                fourth_try.get_chargesign();
//                std::cout << "correction 4: " << fourth_try.populationValid() << std::endl;
//                std::pair<int, std::vector<int>> output4 = fourth_try.populationValid_counter();
//                std::vector<int> stability_new4 = output4.second;
//                fourth_try.total_energy();
//                std::cout << "Valid mistakes 4: " << output4.first << std::endl;
//                std::cout << "correction 4: " << fourth_try.populationValid() << std::endl;
//                std::cout << "Valid mistakes 4: " << output4.first << std::endl;
//                std::cout << "system energ 4: " << fourth_try.system_energy() << std::endl;
//
//                //std::vector<int> stability_new4 = fourth_try.populationValid_counter().second;
//                //                for (auto it = stability_new3.begin(); it!=stability_new3.end(); it++)
//                //                {
//                //                    std::cout << "position_new_third: " << (*it) << " | ";
//                //                }
//                std::cout << "________________________________" << std::endl;
//
//            if ((first_try.populationValid() == 1) && (first_try.system_energy() < system_energy))
//            {
//                system_energy = first_try.system_energy();
//                charge_config = first_try.chargesign;
//            }
//        }
//    }
//
//    std::cout << "smallest energy: " << system_energy << std::endl;
//
//    for (auto it = charge_config.begin(); it != charge_config.end(); it++)
//    {
//        std::cout << (*it) << " | ";
//    }
//
//}


if (1==1)
{

    std::vector<int> initial_sign(location.size(), 0);
    //std::vector<int> initial_sign = {-1,-1,0,0,-1,-1,0,0,-1,-1,0,0,-1,-1,0,0,-1,-1,0,0,-1,-1,0,0,-1,-1,0,0,-1,-1};
    Energyscr generalisation(location, initial_sign);
    generalisation.toeuc();
    generalisation.distance();
    generalisation.potentials();
    generalisation.total_energy();

    // i and j are two indeces of SiDBs which are the start. Distance between these two os calculated.
    float            system_energy = 100000;
    std::vector<int> charge_config;

    for (int z = 0; z < 300; z++)
    {
        std::cout << "z: " << z << std::endl;
        for (int i = 0; i < initial_sign.size(); i++)
        {
            // initial_sign.size()
            // int bound = ((i < (initial_sign.size()-50)) ? (i+49) : initial_sign.size());
            initial_sign[i - 1] = 0;
            initial_sign[i]     = -1;
            //std::cout << "i: " << i << std::endl;
            Energyscr first_try(location, initial_sign);
            first_try.locationeuc = generalisation.locationeuc;
            first_try.db_r        = generalisation.db_r;
            first_try.v_ij        = generalisation.v_ij;
            first_try.total_energy();

            std::vector<int> new_index1;
            std::vector<int> index_start = {i};

            //        std::vector<int> perturber = first_try.find_perturber_alternative();
            //
            //        // i,j SiDBs and perturber are kept on charge state -1
            //        index_start.insert(index_start.end(), perturber.begin(), perturber.end());
            //        sort(index_start.begin(), index_start.end());
            //
            //        // ensure that each index is only on time in vector
            //        index_start.erase(unique(index_start.begin(), index_start.end()), index_start.end());

            for (int k = 0; k < initial_sign.size() - 1; k++)
            // int true_value = 0;
            // while (!isEqual(new_index1,index_start))
            // do
            {
                //std::cout << "i: " << i << std::endl;
                //std::cout << "k: " << k << std::endl;
                new_index1                 = index_start;
                std::vector<int> new_index = first_try.find_new_best_neighbor_GRC(index_start);

                //                            for (auto it = index_start.begin(); it != index_start.end(); it++)
                //                            {
                //                                std::cout << "position_old: " << (*it) << " | ";
                //                            };
                //
                //                            for (auto it = new_index.begin(); it != new_index.end(); it++)
                //                            {
                //                                std::cout << "position_new: " << (*it) << " | ";
                //                            }

                //            if (isEqual(new_index, new_index1))
                //            {
                //                std::cout << "stop" << std::endl;
                //                break;
                //            }

                index_start = new_index;
                // first_try.get_chargesign();
                first_try.total_energy();
                if (first_try.populationValid() == 1 && first_try.system_energy()<system_energy)
                {
                    std::cout << "system energy_old: " << first_try.system_energy() << std::endl;
                    first_try.get_chargesign();
                    system_energy = first_try.system_energy();
                    charge_config = first_try.chargesign;
                }

                std::tuple<int, std::vector<int>, int> output2 = first_try.populationValid_counter();
                // std::cout << "correction_valid: " << std::get<0>(output2) << std::endl;
                // std::cout << "correction_metastable: " << std::get<2>(output2) << std::endl;
                //  std::cout << "correction: " << first_try.populationValid() << std::endl;
                first_try.total_energy();
                if (first_try.populationValid() == 1 && first_try.system_energy()<system_energy)
                {
                    std::cout << "system energy_new: " << first_try.system_energy() << std::endl;
                    first_try.get_chargesign();
                    std::cout << "correction_new: " << first_try.populationValid() << std::endl;
                    system_energy = first_try.system_energy();
                    charge_config = first_try.chargesign;

                    // std::cout << "correction: " << first_try.populationValid() << std::endl;
                    // true_value = first_try.populationValid();
                    // first_try.total_energy();
                    // std::pair<int, std::vector<int>> output2 = first_try.populationValid_counter();
                    // std::cout << "correction: " << output2.first << std::endl;
                    // std::cout << "Valid changed?: " << first_try.populationValid() << std::endl;
                };
            }
        }
        // std::cout << "correction: " << first_try.populationValid() << std::endl;
        // first_try.get_chargesign();
        // std::cout << "system energy: " << first_try.system_energy() << std::endl;
        // first_try.total_energy();

        //            Energyscr second_try(location, first_try.chargesign);
        //            second_try.locationeuc = generalisation.locationeuc;
        //            second_try.db_r        = generalisation.db_r;
        //            second_try.v_ij        = generalisation.v_ij;
        //            second_try.change_chargesign_one(stability);
        //            second_try.total_energy();
        //            second_try.get_chargesign();
        //            std::cout << "correction_new: " << second_try.populationValid() << std::endl;
        //            std::pair<int, std::vector<int>> output2 = second_try.populationValid_counter();
        //            std::vector<int> stability_new = output2.second;
        //            second_try.total_energy();
        //            std::cout << "Valid mistakes_new: " << output2.first << std::endl;
        //            std::cout << "correction_new: " << second_try.populationValid() << std::endl;
        //            std::cout << "Valid mistakes_new: " << output2.first << std::endl;
        //            std::cout << "system energy_new: " << second_try.system_energy() << std::endl;
        //
        //            //                for (auto it = stability_new.begin(); it!=stability_new.end(); it++)
        //            //                {
        //            //                    std::cout << "position_new: " << (*it) << " | ";
        //            //                }
        //            //                std::cout << "________________________________" << std::endl;
        //
        //            Energyscr third_try(location, second_try.chargesign);
        //            third_try.locationeuc = generalisation.locationeuc;
        //            third_try.db_r        = generalisation.db_r;
        //            third_try.v_ij        = generalisation.v_ij;
        //            third_try.change_chargesign_one(stability_new);
        //            third_try.total_energy();
        //            third_try.get_chargesign();
        //            std::cout << "correction 3: " << third_try.populationValid() << std::endl;
        //            std::pair<int, std::vector<int>> output3 = third_try.populationValid_counter();
        //            std::vector<int> stability_new3 = output3.second;
        //            third_try.total_energy();
        //            std::cout << "Valid mistakes 3: " << output3.first << std::endl;
        //            std::cout << "correction 3: " << third_try.populationValid() << std::endl;
        //            std::cout << "Valid mistakes_new_third: " << output3.first << std::endl;
        //            std::cout << "system energy_new_third: " << third_try.system_energy() << std::endl;
        //
        //            Energyscr fourth_try(location, third_try.chargesign);
        //            fourth_try.locationeuc = generalisation.locationeuc;
        //            fourth_try.db_r        = generalisation.db_r;
        //            fourth_try.v_ij        = generalisation.v_ij;
        //            fourth_try.change_chargesign_one(stability_new3);
        //            fourth_try.total_energy();
        //            fourth_try.get_chargesign();
        //            std::cout << "correction 4: " << fourth_try.populationValid() << std::endl;
        //            std::pair<int, std::vector<int>> output4 = fourth_try.populationValid_counter();
        //            std::vector<int> stability_new4 = output4.second;
        //            fourth_try.total_energy();
        //            std::cout << "Valid mistakes 4: " << output4.first << std::endl;
        //            std::cout << "correction 4: " << fourth_try.populationValid() << std::endl;
        //            std::cout << "Valid mistakes 4: " << output4.first << std::endl;
        //            std::cout << "system energ 4: " << fourth_try.system_energy() << std::endl;

        // std::vector<int> stability_new4 = fourth_try.populationValid_counter().second;
        //                 for (auto it = stability_new3.begin(); it!=stability_new3.end(); it++)
        //                 {
        //                     std::cout << "position_new_third: " << (*it) << " | ";
        //                 }


    }
    std::cout << "________________________________" << std::endl;
    //
    std::cout << "smallest energy: " << system_energy << std::endl;

    for (auto it = charge_config.begin(); it != charge_config.end(); it++)
    {
        std::cout << (*it) << " | ";
    }
}
    return 0;
};

