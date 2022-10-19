//
// Created by jan-d on 02.07.2022.
//
#ifndef UNTITLED_SIMU_H
#define UNTITLED_SIMU_H

#include <iostream>
#include <limits>
#include <tuple>
#include <vector>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using DBeucCo = std::vector<float>;
using FPMat   = boost::numeric::ublas::matrix<float>;

inline constexpr float MAX_FLOAT = std::numeric_limits<float>::infinity();

class config
{

  public:
    // Constructor
    config(const std::vector<std::vector<unsigned long>>& input1, const std::vector<int>& input2) :
            locationind(input1),
            chargesign(input2){};

    void toeuc();  // siqad specific coordinates are used to calculate the real euclidean positions of our SiDBs
    void get_locationeuc();  // print euclidean coordinates of each single placed SiDB
    void location_infile(std::ofstream &file);
    void chargeconf_to_index(const int& base);  // charge configuration vector (-1,0,1,-1,...) is assigned to a unique index
    void index_to_chargeconf(const int& base);  // index to charge configuration vector
    void get_chargesign() const;                      // // print the charge configuration vector
    void get_index() const;                           // print the charge index
    void increase_step();                       // charge index is increased by one
    void distance();                            // distance matrix is calculated
    void potentials();
    void identify_outsider();

    std::vector<std::vector<unsigned long>> locationind;
    std::vector<DBeucCo>                    locationeuc;  // euclidean coordinates of our SiDB
    std::vector<int>                        chargesign;   // charge states of all SiDBs are selected
    unsigned int                                     chargemax{};
    FPMat                                   db_r;  // distance matrix
    float                                   distance_value{};
    float                                   potential_value{};
    std::vector<float>                      v_local;  // local potential at each SiDB loc
    unsigned int                            chargeindex = 0;
    FPMat                                   v_ij;
    float                                   potential_energy=0;
    std::vector<int>                        outsider;

};

class Energyscr : public config
{
  public:
    // Constructor
    Energyscr(const std::vector<std::vector<unsigned long>>& input1, const std::vector<int>& input2) :
            config(input1, input2){};

    void                                   total_energy() ;
    float                                   total_energy_EQ(int &index);
    float                                  system_energy() const;
    float                                  system_energy_vec(std::vector<int> &n_in) const;
    void                                    get_distance(int& i, int& j);
    void                                    get_potential(int& i, int& j);
    bool                                   populationValid() const;
    std::tuple<int, std::vector<int>, int> populationValid_counter();
    void                                  find_new_best_neighbor_GRC(std::vector<int>& index_db);
    bool                                   sum_distance(std::vector<int>& input, int& old, int& now);
    std::vector<int>                   find_perturber();
    std::vector<int>                   find_perturber_alternative();

};

#endif  // UNTITLED_SIMU_H
