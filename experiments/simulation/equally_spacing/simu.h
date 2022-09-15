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

using DBindexCo = std::vector<int>;
using DBeucCo = std::vector<float>;
using FPMat = boost::numeric::ublas::matrix<float>;


class config {

    public:

    config(const std::vector<std::vector<unsigned long>> &input1, const std::vector<int> &input2) : locationind(input1), chargesign(input2)  {};


    void toeuc();


    void get_locationeuc();
    void euctoindex(const int &base);
    void indextoeuc(const int &base);
    void get_chargesign();

    void get_index();
    void increase_step();
    void print();
    void calculate_V_d();
    void distance();


    std::vector<std::vector<unsigned long>> locationind;
    std::vector<DBeucCo> locationeuc;
    std::vector<int> chargesign;
    int chargemax;
    FPMat db_r;
    float distance_value;
    float potential_value;
    std::vector<float> v_local;

    int chargeindex=0;

    private:


};

class Energyscr : public config {
public:
    Energyscr(const std::vector<std::vector<unsigned long>> &input1, const std::vector<int> &input2) : config(input1, input2) {};

    void potentials();
    void change_chargesign(int &i, int &j);
    void change_chargesign_one(std::vector<int> &vector);
    void get_distance(int &i,int &j);
    void get_potential(int &i,int &j);
    std::vector<float> laplace(std::vector<float> &vector_input, const float &f);
    //auto gradient(std::vector<float> &vector_input);
    std::pair<std::vector<float>, int> step(std::vector<float> &vector_input);
    std::vector<int> search_same_distance(std::vector<int> &index_db);
    std::vector<int> search_same_distance_new(std::vector<int> &index_db);
    std::vector<int> search_same_potential(std::vector<int> &index_db);
    float system_energy();
    void total_energy();
    std::vector<int> find_perturber();
    std::vector<int> find_perturber_alternative();
    bool populationValid() const;
    std::pair<int, std::vector<int>> populationValid_counter();
    void set_charge();
    std::pair<float, std::vector<int>> shortestPath(int &u, int &v, int &k);

    FPMat v_ij;
private:




};

class EnergyNoscr : public config {
public:
    EnergyNoscr(const std::vector<std::vector<unsigned long>> &input1, const std::vector<int> &input2) : config(input1, input2) {};
    void potentials();
    void total_energy();

private:

    FPMat v_ij;

};


#endif //UNTITLED_SIMU_H
