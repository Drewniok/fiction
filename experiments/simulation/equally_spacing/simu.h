//
// Created by jan-d on 02.07.2022.
//
#include <iostream>
#ifndef UNTITLED_SIMU_H
#define UNTITLED_SIMU_H
#include <tuple>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using DBeucCo = std::vector<float>;
using FPMat = boost::numeric::ublas::matrix<float>;


class config {

    public:

      // Constructor
    config(const std::vector<std::vector<unsigned long>> &input1, const std::vector<int> &input2) : locationind(input1), chargesign(input2)  {};

    // Destructor
    ~config(){};


    void toeuc(); // siqad specific coordinates are used to calculate the real euclidean positions of our SiDBs
    void get_locationeuc(); // print euclidean coordinates of each single placed SiDB
    void chargeconf_to_index(const int &base); // charge configuration vector (-1,0,1,-1,...) is assigned to a unique index
    void index_to_chargeconf(const int &base); // index to charge configuration vector
    void get_chargesign(); // // print the charge configuration vector
    void get_index(); // print the charge index
    void increase_step(); // charge index is increased by one
    void distance(); // distance matrix is calculated
    void potentials();

    std::vector<std::vector<unsigned long>> locationind;
    std::vector<DBeucCo> locationeuc; // euclidean coordinates of our SiDB
    std::vector<int> chargesign; // charge states of all SiDBs are selected
    int chargemax;
    FPMat db_r; // distance matrix
    float distance_value;
    float potential_value;
    std::vector<float> v_local; // local potential at each SiDB location
    int chargeindex=0;
    FPMat v_ij;
};

class Energyscr : public config {

public:

    // Constructor
    Energyscr(const std::vector<std::vector<unsigned long>> &input1, const std::vector<int> &input2) : config(input1, input2) {};

    // Destructor
    ~Energyscr(){};


    void total_energy();
    float system_energy();
    bool populationValid() const;
    std::tuple<int, std::vector<int>, int> populationValid_counter();
    std::vector<int> find_new_best_neighbor_GRC(std::vector<int> &index_db);
    bool sum_distance(std::vector<int> &input, int &old, int &now);


    void change_chargesign(int &i, int &j);
    void change_chargesign_one(std::vector<int> &vector);
    void get_distance(int &i,int &j);
    void get_potential(int &i,int &j);
    std::vector<float> laplace(std::vector<float> &vector_input, const float &f);
    std::vector<int> find_new_best_neighbor(std::vector<int> &index_db);
    std::vector<int> search_same_distance(std::vector<int> &index_db);
    std::vector<int> search_same_distance_new(std::vector<int> &index_db);
    std::vector<int> search_same_potential(std::vector<int> &index_db);
    std::vector<int> find_perturber();
    std::vector<int> find_perturber_alternative();
    void set_charge();
    std::pair<float, std::vector<int>> shortestPath(int &u, int &v, int &k);

};



#endif //UNTITLED_SIMU_H
