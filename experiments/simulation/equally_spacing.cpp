//
// Created by jan-d on 02.09.2022.
//

#include "fiction/io/read_sqd_layout.hpp"
#include "fiction/io/write_sqd_layout.hpp"
#include "fiction/technology/area.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/types.hpp"
#include "simu.cpp"
#include "simu.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <execution>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <tuple>
#include <vector>

using namespace fiction;
using std::filesystem::directory_iterator;

template<typename T>
float average(std::vector<T> const& v){
    if(v.empty()){
        return 0;
    }

    auto const count = static_cast<T>(v.size());
    return std::reduce(v.begin(), v.end()) / count;
}


template<typename T>
std::vector<T> unique_values( std::vector<T> & input_vec ){

    std::vector<T> uniques( input_vec.size() );
    typename std::vector<T>::iterator it;
    it = std::unique_copy (input_vec.begin(), input_vec.end(), uniques.begin() );

    std::sort( uniques.begin(), it );

    it = std::unique_copy( uniques.begin(), it, uniques.begin() );

    uniques.resize( std::distance(uniques.begin(), it) );

    return uniques;
}

template<typename T>
std::vector<int> count_unique_values( std::vector<T> & input_vec ){

    std::vector<T> uniques = unique_values<T>( input_vec );
    std::vector<int> counts( uniques.size() );

    for(size_t i = 0; i < counts.size(); ++i)
        counts[i] = std::count( input_vec.begin(), input_vec.end(), uniques[i] );

    return counts;
}


template<typename T>
std::ostream & operator << (std::ostream & o, const std::vector<T> & v){

    o << "[ ";
    for (size_t i = 0; i < v.size()-1; i++)
        o << v[i] << ", ";

    o << v[v.size()-1];

    std::cout << " ]" << std::endl;

    return o;
}

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

std::string mux = "/Users/jandrewniok/CLionProjects/fiction/experiments/bestagon/layouts/sam_layouts/mux_mu32.sqd";
std::string and_gate = "/Users/jandrewniok/CLionProjects/fiction/experiments/bestagon/layouts/select_gates/21_hex_inputsdbp_and_v19.sqd";
std::string or_gate = "/Users/jandrewniok/CLionProjects/fiction/experiments/bestagon/layouts/select_gates/22_hex_inputsdbp_hourglass_v0.sqd";

    };

    //std::string path = "C:/Users/jan-d/OneDrive/Dokumente/PhD/FCN/sqd/";
    std::string path = "/Users/jandrewniok/CLionProjects/fiction/experiments/bestagon/layouts/select_gates";
    std::string folder = "/gates_sqd";
    //std::string path = "/Users/jandrewniok/CLionProjects/fiction/experiments/bestagon/layouts/select_gates/wrapper_files/";

int main()
{

        bool EXGS_on   = false;


        std::vector<float> energy_EXGS = {1.01035, 1.01035, 1.07667, 0.92899, 1.07623, 0.9814 , 1.19412,
                                            0.97742, 0.97742, 1.02576, 1.02445, 0.9838 , 1.00219, 0.73757,
                                            0.73921, 0.78406, 0.91751, 0.91751, 0.97982, 0.97419, 1.03774,
                                            1.24894, 1.24894, 1.31221, 1.2935 , 1.36088, 1.01939, 1.01939,
                                            1.06535, 1.08088, 1.12825, 1.39355, 1.39355, 1.44769, 1.44769,
                                            1.50717, 1.03817, 1.03817, 1.106  , 0.90262, 0.90262, 0.95849,
                                            0.74904, 0.74904, 0.81316, 0.58676, 0.58676, 0.65266, 1.16061,
                                            1.16061, 0.95537, 0.94709, 0.99646, 0.93664, 0.93664, 0.96742,
                                            0.99119, 0.7797};

        int loop_count = 0;
        std::ofstream result_file(path + "all_results.txt");
        result_file << "TTS;" << "energy1;"<< "layout;" << "number_sidb" << std::endl;
        for (const auto& file : directory_iterator(path + folder))
        {


            //std::string selection = file.path();
            //std::cout << file.path() << std::endl;
//        for (int circuit_num =0; circuit_num<300;circuit_num++)
//        {
            std::vector<float> smallest_energy;
            std::vector<float> time;



//                 std::vector<std::vector<unsigned long>> location;
//
//                 std::random_device                                       dev;
//                 std::mt19937                                             rng(dev());
//                 std::uniform_int_distribution<std::mt19937::result_type> dist6(0, 60);
//                 std::uniform_int_distribution<std::mt19937::result_type> dist10(0, 2);
//                 std::uniform_int_distribution<std::mt19937::result_type> dist1(0, 1);
//                 std::uniform_int_distribution<std::mt19937::result_type> dist_number_sidb(10, 12);
//
//                 int number_of_sidb = dist_number_sidb(rng);
//                 // generate random layout
//                 std::vector<int> numbers_x;
//                 for (int i = 0; i < number_of_sidb; i++)  // add 0-99 to the vector
//                     numbers_x.push_back(dist6(rng));
//                 unsigned seed_x = std::chrono::system_clock::now().time_since_epoch().count();
//                 std::shuffle(numbers_x.begin(), numbers_x.end(), std::default_random_engine(seed_x));
//
//                 std::vector<int> numbers_z;
//                 for (int i = 1; i < number_of_sidb + 1; i++)  // add 0-99 to the vector
//                     numbers_z.push_back(dist10(rng));
//                 unsigned seed_z = std::chrono::system_clock::now().time_since_epoch().count();
//                 std::shuffle(numbers_z.begin(), numbers_z.end(), std::default_random_engine(seed_z));
//
//                 std::vector<int> numbers_y;
//                 for (int i = 0; i < number_of_sidb; i++)  // add 0-99 to the vector
//                     numbers_y.push_back(dist1(rng));
//                 unsigned seed_y = std::chrono::system_clock::now().time_since_epoch().count();
//                 std::shuffle(numbers_y.begin(), numbers_y.end(), std::default_random_engine(seed_z));
//
//                 std::set<std::vector<unsigned long>> collect_set;
//
//                 for (int l = 0; l < numbers_y.size(); l++)
//                 {
//                     collect_set.insert({static_cast<unsigned long>(numbers_x[l]), static_cast<unsigned long>(numbers_z[l]),
//                                         static_cast<unsigned long>(numbers_y[l])});
//                 };
//
//                 for (auto it = collect_set.begin(); it != collect_set.end(); it++)
//                 {
//                     location.push_back(*it);
//                 }

            std::string selection = file.path();
             //std::string selection = path;
            std::cout << selection << std::endl;

            //
            //

            const auto lyt = read_sqd_layout<sidb_cell_clk_lyt>(selection);

            std::vector<std::vector<unsigned long>> location;
            std::vector<cell<sidb_cell_clk_lyt>>    all_cells{};
            all_cells.reserve(lyt.num_cells());

            lyt.foreach_cell([&all_cells](const auto& c) { all_cells.push_back(c); });

            // transform coordinates in new coordinates, inspired by SIQAD
            for (const auto& c : all_cells)
            {
                auto          X = static_cast<unsigned long>(c.x);
                auto Y = static_cast<unsigned long>(((c.y) - ((c.y) % 2)) * 0.5);
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

            //        std::for_each(std::execution::par, all_cells.cbegin(), all_cells.cend(),
            //                      [&location](const auto& c)
            //                      {
            //                          auto          X = c.x;
            //                          unsigned long Y = ((c.y) - ((c.y) % 2)) * 0.5;
            //                          unsigned long Z = ((c.y) % 2);
            //                          location.push_back({X, Y, Z});
            //                      });

//            std::cout << std::endl
//                      << fmt::format("There are {} SiDBs in the circuit", location.size()) << std::endl
//                      << std::endl;

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

            Energyscr check(location, initial_sign);
            check.toeuc();

            std::string filename = selection.substr(0, selection.find("sqd/"));
            filename = selection.substr(filename.size()+4, selection.find("sqd/"));


//            std::ofstream File_python(path + "/wrapper_files/random/random/loc1/" + filename + "_wrapper.txt");
//            std::ofstream energy_python(path + "/wrapper_files/random/random/energy1/" + filename + "_energy.txt");

            float energy_file_read;
            std::fstream my_file;
            //std::cout << path + "/wrapper_files/energy/" + filename << std::endl;
            my_file.open(path + "/wrapper_files/energy/" + filename + "_energy.txt", std::ios::in);
            if (!my_file) {
                std::cout << "No such file";
            }
            else {
                float energy_file;

                while (1) {
                    my_file >> energy_file;
                    if (my_file.eof())
                        break;

                    //std::cout << ch;
                }
                energy_file_read = energy_file;
            }

            std::cout << energy_file_read << std::endl;



            check.distance();
            // std::vector<int> perturber = check.find_perturber_alternative();
            int warning = 0;
            for (int row = 0; row < location.size(); row++)
            {
                for (int col = 0; col < location.size(); col++)
                {
                    if ((check.db_r(row, col) < 8E-10) && (check.db_r(row, col) != 0))
                    {
                        warning += 1;
                    }
                }
            }

            if (warning > 0)
            {
                continue;
            }


            std::ofstream File_python(path + "/random/loc3/" + "_wrapper.txt");
            std::ofstream energy_python(path + "/random/energy3/" + "_energy.txt");
            check.location_infile(File_python);
            // std::vector<int> v{0,1,2,4,8,16,32,64,128,256,512};

//            int circuit_num    = 1;
//            int number_of_sidb = 1;
//
//            std::ofstream outFile(selection + "/location_random3/" + std::to_string(circuit_num) +
//                                  std::to_string(number_of_sidb) + "_location.txt");
//            std::ofstream chargeFile(selection + "/charge_random3/" + std::to_string(circuit_num) +
//                                     std::to_string(number_of_sidb) + "_charge.txt");

            //
  //          std::ofstream outFile(selection  + "/XX/location.txt");
  //          std::ofstream chargeFile(selection + "/XX/charge.txt");

            // the important part
            //std::cout << selection << std::endl;

//            outFile << "x;"
//                    << "y;"
//                    << "z;" << std::endl;
//            for (int i = 0; i < location.size(); i++)
//            {
//                for (const auto& e : location[i]) outFile << std::to_string(e) << ";";
//                outFile << std::endl;
//            };
//
//            for (int i = 0; i < location.size(); i++)
//            {
//                chargeFile << std::to_string(i) << ";";
//            };


  //          chargeFile << "type;";
  //          chargeFile << "runtime;";
   //         chargeFile << "energy1" << std::endl;
;
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
            auto  t1                = std::chrono::high_resolution_clock::now();
            float system_energy_min = energy_file_read;
            if (EXGS_on)
            {
                std::vector<int> charge;
                std::vector<int> initial_sign(location.size(), -1);
                Energyscr        generalisation(location, initial_sign);
                generalisation.toeuc();                 // transform lattice coordinates to euclidean coordinates
                generalisation.chargeconf_to_index(2);  // chargesign-vector (e.g. -1 0 0 -1 ...) is converted to an index | number in () is count state | 3 = DB-,DB0,DB+
                generalisation.distance();    // distance-matrix d(i,j) is calculated, only required once
                generalisation.potentials();  // potential-matrix v(i,j) is calculated, also only required once, since
                                              // chargesign is stil not included
               // std::cout << "max: " << generalisation.chargemax
                //          << std::endl;  // number of all possible charge distributions
                // can be output

                int counter = 0;
                while (generalisation.chargeindex <= generalisation.chargemax)
                {
                    counter += 1;
                    if ((counter % 1000000) == 0)
                    {
                        // std::cout << "iteration step: " << counter << std::endl;
                        // generalisation.get_chargesign();
                        // std::cout << std::endl;
                    }
                    generalisation.total_energy();  // update of v_local[i] which changes as soon as the charge distribution changes
                    if (generalisation.populationValid() &&
                        (generalisation.system_energy() <=
                         system_energy_min))  // if the charge state configuration is the energetically best so far
                                              // and physically valid, it is saved and printed out
                    {
                        system_energy_min = generalisation.system_energy();  //
                 //       std::cout << "system_energy_min: " << system_energy_min
                 //                 << std::endl;              // output when new best mark is reached
                        charge = generalisation.chargesign;  // save new best chargesign
                    }
                    generalisation.increase_step();         // index is increased by one
                    generalisation.index_to_chargeconf(2);  // index is converted to the chargesign-vector
                }
             //   std::cout << std::endl
//                          << "-------------------------------------------------------------------------------"
//                          << std::endl;
                // after all charge configurations are evaluated, the energetically lowest charge which is physically valid at the same time is printed
        //        std::cout << "smallest energy1 found: " << system_energy_min << std::endl;
        //        std::cout << "charge configuration: |";
                //for (auto it = charge.begin(); it != charge.end(); it++)
                //{
                    //std::cout << (*it) << " | ";
//                  chargeFile << std::to_string(*it) << ";";
                //}
         //       chargeFile << "EXGS;";
        //        std::cout << std::endl
//                          << "-------------------------------------------------------------------------------"
//                          << std::endl;
                const auto                    t11        = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> diff_first = t11 - t1;
//                std::cout << std::endl
//                          << std::endl
//                          << fmt::format("It took {} s to simulate the {} circuit", diff_first.count(), selection)
//                          << std::endl;
  //              chargeFile << std::to_string(diff_first.count()) << ";";
  //              chargeFile << std::to_string(system_energy_min);
            }
            energy_python << std::to_string(system_energy_min) << std::endl;
            // -------------------------------------------------------------------------------------------------------------------------------

            // ---------------------------------------------- MAX-MIN DIVERSITY APPROACH
            // -------------------------------------------------

            int  identical = 0;
            int unidentical = 0;

            Energyscr generalisation(location, initial_sign);
            generalisation.toeuc();
            generalisation.distance();
            generalisation.potentials();
            generalisation.total_energy();
            //std::vector<int> perturber = generalisation.find_perturber_alternative();

            for (int q = 0; q<1000; q++)
            {
                float energy_exact = energy_file_read;
                loop_count += 1;

            float system_energy = MAX_FLOAT;

            auto t12 = std::chrono::high_resolution_clock::now();

            if (true)
            {

                // define the range, we want to select certain "start" SiDB randomly

                std::vector<int> charge_config;

                //static const std::vector<bool> iterator_helper(100, false);
                // this loop executes in parallel if the compiler so wishes
                //std::vector<std::vector<float>> collection_all;

                // std::for_each(std::execution::par, iterator_helper.cbegin(), iterator_helper.cend(), [&](const auto& b){
                int threashold_num     = 10;
                int count_equal_energy = 0;
                std::vector<int> initial_sign(location.size(), 0);
                for (int z = 0; z < threashold_num; z++)
                {

                    // std::cout << z << std::endl;
                    for (int i = 0; i < location.size(); i++)
                    {


                    //int i = 0;  // at the beginning, all SiDBs are set to the charge 0. Now, one SiDB is chosen
                    initial_sign[i] = -1;
                    initial_sign[i-1] = 0;

                        //  randomly and set to -1
                        // std::cout << z << std::endl;
                        Energyscr first_try(location, initial_sign);
                        // -------------- some attributes do not change when charge sign is changed --------------
                        first_try.locationeuc = generalisation.locationeuc;
                        first_try.db_r        = generalisation.db_r;
                        first_try.v_ij        = generalisation.v_ij;
                        first_try.v_local = std::vector<float>(first_try.v_ij.size1(),0);
                            // ---------------------------------------------------------------------------------------
                        //system_energy_local = first_try.total_energy();  // However, v_local[i] has to recalculated
                        //std::vector<int> perturber = first_try.find_perturber_alternative();
                        // perturber.push_back(i);

                        std::vector<int> index_start = {i};  // index_start includes all the indices of "-1"-SiDBs
                        //perturber.push_back(i);
                        first_try.total_energy_EQ(i);
                        auto index_start_size = index_start.size();
                        auto index_end_size   = index_start.size();


                        for (int set_dbs = 0; set_dbs < location.size()/2+5; set_dbs++)
                        {
                            index_start_size = index_start.size();
                            // based on the max-min diversity problem, all -1 SiDBs are selected iteratively
                            //float save_energy = first_try.potential_energy;
                            //std::vector<float> save_v_local = first_try.v_local;

                            first_try.find_new_best_neighbor_GRC(index_start);
                            index_end_size = index_start.size();

                            //first_try.total_energy();
                            // first_try.get_chargesign();
                            //    check if new charge configuration with lower energy1 is found
                            if (first_try.populationValid() && first_try.potential_energy <= system_energy)
                            {
                                if (first_try.potential_energy == system_energy)
                                {
                                    count_equal_energy += 1;
                                }
                                //std::cout << "system energy1 old: " << first_try.system_energy() << std::endl;
                                //first_try.get_chargesign();
                                system_energy = first_try.potential_energy;
                                charge_config = first_try.chargesign;

                                //system_energy = first_try.total_energy;
                                break;
                            }

                        };
                    };
                    if (count_equal_energy > 5)
                    {
                        break;
                    }
                }

                //std::cout << "________________________________" << std::endl;
                //
               // std::cout << "smallest energy1: " << system_energy << std::endl;
                smallest_energy.push_back(std::round((system_energy)*100000)/100000);
   //             chargeFile << std::endl;
   //             chargeFile << "EQ;";


            const auto t2     = std::chrono::high_resolution_clock::now();
            const auto ms_int = std::chrono::duration_cast<std::chrono::seconds>(t2 - t12);

            if ((system_energy_min == system_energy) && (EXGS_on))
            {
                identical += 1;
            }
            else if ((system_energy != MAX_FLOAT) && (system_energy_min != system_energy) && (EXGS_on))
            {
                unidentical += 1;
            }

            else if ((!EXGS_on) && (std::abs(system_energy - energy_exact) < 0.00001))
            {
                identical += 1;
            }

            else if ((!EXGS_on) && (std::abs(system_energy - energy_exact) > 0.00001))
            {
                unidentical += 1;
            }
            std::chrono::duration<double> diff = t2 - t12;
  //          chargeFile << std::to_string(diff.count());
            if (system_energy != MAX_FLOAT)
            {
  //              chargeFile << ";" << std::to_string(system_energy);
            }
//            std::cout << std::endl
//                      << std::endl
//                      << fmt::format("It took {} s to simulate the {} circuit", diff.count(), selection) << std::endl;


            //std::cout << "identical: " << identical << std::endl;
            //std::cout << "unidentical: " << unidentical << std::endl;
            time.push_back(diff.count());
        }

        }
        std::cout << "identical: " << identical << std::endl;

    auto const a = average(time);
 //   std::cout << "average time: " << a << "\n";
    std::vector<float> v_uniques_int = unique_values<float>( smallest_energy );
    std::cout << "Vector of unique values: " << v_uniques_int;

    std::vector<int> v_counts_int = count_unique_values<float>( smallest_energy );
    float write_energy = v_uniques_int[0];
    std::cout << "Vector of unique counts: " << v_counts_int;
    float TTS;
    if ((v_uniques_int[0]==MAX_FLOAT) && (v_uniques_int.size() == 1))
    {
        TTS = 50;
        write_energy = 1000000;
        std::cout << "WARNING" << std::endl;
    }

    else if (identical == 0)
    {
        write_energy = 1000000;
        TTS = 50;
    }

    else if (v_counts_int[0] == 1000)
    { TTS = a;}

    else
    {    TTS = a * log(1.0-0.997) / log(1.0-static_cast<float>(v_counts_int[0])/1000.0);
    }

    std::cout << "\n\n";

    result_file << std::to_string(TTS) << ";" << std::to_string(v_uniques_int[0]) << ";" << filename << ";" << std::to_string(location.size()) << std::endl;
    //result_file << std::to_string(TTS) << ";" << std::to_string(write_energy) << ";" << std::to_string(circuit_num) << ";" << std::to_string(location.size()) << std::endl;


    }
return 0;
}
