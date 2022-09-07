//
// Created by jan-d on 02.09.2022.
//

#include <fiction/io/read_sqd_layout.hpp>
#include <fiction/technology/area.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/types.hpp>
//#include <fiction/algorithms/simulation/...>
#include <iostream>

#include <vector>
#include <math.h>
#include <tuple>
#include "simu.h"
#include "simu.cpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


#include <src/bqp.cpp>
#include <src/tabu_search.cpp>
#include <src/utils.cpp>
#include <string>

#include <fmt/format.h>

#include <iostream>
#include <set>

using namespace fiction;


template<typename T>
bool isEqual(std::vector<T> &first, std::vector<T> &second)
{
    if (first.size() != second.size()) {
        return false;
    }

    return std::is_permutation(first.begin(), first.end(), second.begin());
}

bool sortPairs(std::vector<unsigned long> const& s1, std::vector<unsigned long> const& s2)
{
    // If the values of the first column are not equal,
    // just use them to order s1 and s2.
    if ( s1[1] != s2[1] )
    {
        return s1[1] < s2[1];
    }
    // If the values of the first column are equal,
    // use the values of the second column to order s1 and s2.
    return s1[0] < s2[0];
}

bool sortPairsx(std::vector<unsigned long> &x, std::vector<unsigned long>  &y)
{
    return x[0] < y[0];
};


int main()
{
    const auto lyt = read_sqd_layout<sidb_cell_clk_lyt>("C:/Users/jan-d/Downloads/sidb-bestagon-gate-library-xml-fix/sidb-bestagon-gate-library-xml-fix/bestagon-gates/2i1o_and/21_hex_inputsdbp_and_v19_01.sqd");

    //std::cout << fmt::format("Read a layout of {} x {} size", lyt.x() + 1, lyt.y() + 1) << std::endl;

    std::vector<std::vector<unsigned long>> location;

    std::vector<cell<sidb_cell_clk_lyt>> all_cells{};
    all_cells.reserve(lyt.num_cells());

    lyt.foreach_cell(
        [&lyt, &all_cells](const auto& c)
        {
             //collect.push_back({c.x(),c.y(),c.z()});
//            const auto adjacent_cells = lyt.adjacent_coordinates<std::set<cell<sidb_cell_clk_lyt>>>(c);
//            std::cout << fmt::format("DB at position {} with neighbors {}", c, fmt::join(adjacent_cells, ", "))
//                      << std::endl;
            //const auto adjacent_cells = lyt.adjacent_coordinates<std::set<cell<sidb_cell_clk_lyt>>>(c);
             //std::cout << c << std::endl;
            //std::cout << fmt::format("DB at position {}", c) << std::endl;
            all_cells.push_back(c);
        });

    //std::cout << all_cells.size() << std::endl;

    for (int i  =0; i<all_cells.size(); i++)
    {
        auto X = all_cells[i].x;
        auto Y = ((all_cells[i].y) - ((all_cells[i].y) % 2))/2;
        auto Z = ((all_cells[i].y) % 2);
        //std::cout << "X: " << X << " | " << "Y: " << Y << " | " << "Z: " << Z << std::endl << std::endl;
      location.push_back({X,Y,Z});
      //auto test = all_cells[i].x;
      //std::cout << typeof(test) << std::endl;
    };


    //std::sort(location.begin(), location.end(), sortPairs);
    std::sort(location.begin(), location.end(), sortPairs);


    for (auto it = location.begin(); it != location.end();it++)
    {
        std::cout << "X: " << (*it)[0] << " | " << "Y: " << (*it)[1] << " | " << "Z: " << (*it)[2] << std::endl;
        //std::cout << (*it)[0] << std::endl;
    }





//    cell<sidb_cell_clk_lyt> c{89, 62};
//
//    std::cout << fmt::format("Position {} is {}a primary input", c, lyt.is_pi(c) ? "" : "not ") << std::endl;
//
//    area_params<sidb_technology> ps{};
//    area_stats                   st{};
//    area(lyt, ps, &st);
//    std::cout << fmt::format("The layout has an area of {} nm²", st.area) << std::endl;

    // simulate_t_dependence(lyt);
    std::vector<int> initial_sign(location.size(),0);
    Energyscr first_try(location,initial_sign);
    first_try.toeuc();;
    first_try.distance();
    first_try.potentials();
    first_try.total_energy();
    std::vector<int> pertuber_vector = first_try.find_perturber();
    //pertuber_vector.push_back(130);
    //pertuber_vector.push_back(132);

    for (auto it = pertuber_vector.begin(); it != pertuber_vector.end();it++)
    {
        std::cout << (*it) << std::endl;
        //std::cout << (*it)[0] << std::endl;
    }


    for (int i = 0; i<initial_sign.size();i++)
    {
        for (int j = i+1; j<initial_sign.size();j++)
        {
            Energyscr first_try(location,initial_sign);
            //first_try.find_perturber();
            //first_try.get_chargesign();
            first_try.change_chargesign(i,j);

            first_try.toeuc();
            first_try.distance();
            first_try.potentials();
            first_try.get_distance(i,j);
            std::vector<int> index_start = {i,j};
            std::vector<int> perturber = first_try.find_perturber();
            //perturber.push_back(130);
            //pertuber_vector.push_back(132);
            //first_try.get_chargesign();
            index_start.insert(index_start.end(), perturber.begin(),perturber.end());
            sort( index_start.begin(), index_start.end() );
            index_start.erase( unique( index_start.begin(), index_start.end() ), index_start.end() );

            std::vector<int> new_index1;

            do {
                new_index1 = index_start;
                std::vector<int> new_index = first_try.search_same_distance(index_start);
                index_start = new_index;
            }

            while (!isEqual(new_index1,index_start));

            first_try.get_chargesign();
            first_try.total_energy();
            first_try.system_energy();

            std::cout << "system stable?: " << first_try.populationValid() << " | counter: " << first_try.populationValid_counter().first << std::endl;
            if (first_try.populationValid_counter().first < 5)
            {
                std::vector<float> stability = first_try.populationValid_counter().second;
                for (auto it = stability.begin(); it!=stability.end(); it++)
                {
                    std::cout << (*it) << std::endl;
                }
            }
        }

    }

    return 0;
};

