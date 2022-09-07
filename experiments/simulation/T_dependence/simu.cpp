//
// Created by jan-d on 02.07.2022.
//
#include <math.h>
#include "simu.h"
#include <iostream>
#include <cmath>
#include <math.h>
#include <string>
#include <climits>
#include <numeric>

float INF = 1000000000.0;

struct Params {

    // Gitterkonstanten
    float a = 3.84 * std::pow(10,-10);
    float b = 2.25 * std::pow(10,-10);
    float c = 7.68 * std::pow(10,-10);

    //

    float epsilon = 8.854 * std::pow(10,-12);
    float epsilon_screen = 5.6;
    float k = 1 / (4*3.141592653 * epsilon * epsilon_screen);
    float e = 1.602 * std::pow(10,-19);
    float mu = -0.32;
    float mu_p = mu - 0.59;
    float tf = 5 *  std::pow(10,-9);
    const float POP_STABILITY_ERR = 1E-6;

};

struct Params params;


void config::toeuc() {
    std::vector<DBeucCo> m;
    for (unsigned int i = 0; i < locationind.size(); ++i) {

           m.push_back({(locationind[i][0]) * params.a, locationind[i][1] * params.c + locationind[i][2] * params.b});

    }
    locationeuc = m;
};

void config::get_locationeuc() {
    for (unsigned int i = 0; i<locationeuc.size(); ++i)
    {

       std::cout << locationeuc[i][0] << " | " << locationeuc[i][1] << std::endl;

    }
};

void config::euctoindex(const int &base) {
    chargeindex = 0;
    for (unsigned int i=0; i<chargesign.size(); i++) {
        chargeindex += (chargesign[i]) * pow(base, chargesign.size() - i - 1);
    }
};

void config::indextoeuc(const int &base)
{
    int n_sidb = locationind.size();
    for (int i=0; i<n_sidb; i++) {
        chargesign.push_back(-1);
        //std::cout << chargesign[i] << std::endl;
    };
    //std::cout << chargesign.size() << std::endl;
    int i = n_sidb - 1; // i is the number of the n_dbs
    while (chargeindex > 0) {
        div_t d; // Structure to represent the result value of an integral division performed by function div.
        d = div(chargeindex, base);
        chargeindex = d.quot;
        chargesign[i--] = d.rem;
    }

};



void config::get_chargesign() {
    for (unsigned int i = 0; i<chargesign.size(); i++) {
        //std::cout << i << std::endl;
        if (i == 0) {
            std::cout << "chargeconfig: " << chargesign[i] << "|" ;
        }
        else if (i != chargesign.size()-1){
            std::cout << chargesign[i] << "|";
        }
        else{
            std::cout << chargesign[i] << "|" << std::endl;
        }
    }
};





void config::get_index() {
    std::cout << "Index: " << chargeindex << std::endl;
};

void config::increase_step() {
    chargeindex += 1;
};

float energy(FPMat &m, unsigned int i, unsigned int j)
{
    float e = params.k * 1/m(i,j) * params.e; // eV
    return e;
};

float energy_screened(FPMat &m, unsigned int i, unsigned int j)
{
    float e = params.k * 1/m(i,j)

            * std::exp(-m(i,j)/params.tf) * params.e; // eV
    return e;
};

void config::distance(){
    db_r.resize(locationeuc.size(),locationeuc.size());
    for (unsigned int i = 0; i < db_r.size1(); i++)
    {
        db_r(i,i)=0;
        for (unsigned int j = i+1; j < db_r.size2(); j++)
        {
            db_r(i,j) = std::sqrt(std::pow(locationeuc[i][0] - locationeuc[j][0],2) + std::pow(locationeuc[i][1] - locationeuc[j][1],2));
            //std::cout << db_r(i,j) << std::endl;
            db_r(j,i) = db_r(i,j);
        }

    }
   // std::cout << db_r << std::endl;
};

void Energyscr::get_distance(int &i, int &j)
{
    distance_value = db_r(i,j);
    //std::cout << db_r(i,j) << std::endl;
};

std::vector<int> Energyscr::search_same_distance(std::vector<int> &index_db)
{

        for (int l =0; l<db_r.size1();l++) {
            //std::cout << "k: " << k << std::endl;
            //std::cout << "l: " << l << std::endl;
            //std::cout << "distance: " << distance_value << std::endl;
            //std::cout << "db_r(i,k): " << db_r(i, k) << std::endl;
            //std::cout << "db_r(j,k): " << db_r(j, k) << std::endl;
            float left = 0.7;
            float right = 1.30;
            int number_correct = 0;
            //std::vector<int> possible_pos;
            for (int n = 0; n < index_db.size(); n++) {
                //std::cout << "db_r(l,n): " << db_r(l, n) << std::endl;
//            if ((left * distance_value <= db_r(i, k)) && (db_r(i, k) <= right * distance_value) &&
//                (left * distance_value <= db_r(j, k)) && (db_r(j, k) <=  right * distance_value) && (left * distance_value <= db_r(i, l)) && (db_r(i, l) <=  right * distance_value) &&
//                (left * distance_value <= db_r(j, l)) && (db_r(j, l) <=  right * distance_value) &&
//                (left * distance_value <= db_r(k, l)) && (db_r(k, l) <=  right * distance_value)) {
                //std::cout << "l: " << l << std::endl;
                //std::cout << "n: " << index_db[n] << std::endl;
                //std::cout << "db_r(l,n): " << db_r(l, index_db[n]) << std::endl;
                //std::cout << "distance_value: " << distance_value << std::endl;
                if (((left * distance_value <= db_r(l, index_db[n])) && (db_r(l, index_db[n]) <= right * distance_value) ) || (db_r(l, index_db[n]) >= left * distance_value)) {
                    number_correct += 1;
                }
            }
            //std::cout << "number correct: " << number_correct << std::endl;
            //std::cout << "index_db_size: " << index_db.size() << std::endl << std::endl;

            if (number_correct ==index_db.size()) {
                chargesign[l] = -1;
                index_db.push_back(l);
                //possible_pos.push_back(l);
            }

//                std::cout << "new SiDB | k: " << k << std::endl;
//                std::cout << "new SiDB | l: " << l << std::endl;
//                chargesign[k] = -1;
//                chargesign[l] = -1;
                //get_chargesign();
                //system_energy();
                //potentials();
                //total_energy();
                //std::cout << "stability: " << populationValid() << std::endl;
            }
        return index_db;
    };


std::vector<int> Energyscr::find_perturber(){
    std::vector<int> perturber_collect;
    for (int i = 0; i<db_r.size1(); i++)
    {
        int c = 0;
        float threas = 1.7 * std::pow(10,-9);
        for (int j = 0; j<db_r.size1(); j++)
        {
            if (db_r(i,j)>threas && db_r(i,j)!=0)
            {
                c+=1;
            }
        }
        //std::cout << "c: " << c << std::endl;
        if (c == db_r.size1()-1)
        {
            chargesign[i] = -1;
            perturber_collect.push_back(i);
        }
    }
    return perturber_collect;
};

void Energyscr::set_charge() {
    float min = 10.0;
    //std::cout << min << std::endl;
    int index = 0;
    for (int i = 0; i < v_local.size(); i++) {
        //std::cout << chargesign[i] << std::endl;
        if ((-v_local[i] < min) && (chargesign[i] != -1)) {
            min = -v_local[i];
            //std::cout << i << std::endl;
            //std::cout << v_local[i] << std::endl;
            index = i;
        }
    }
    //std::cout << min << std::endl;
    if (chargesign[index]>-1)
    {
        chargesign[index] = chargesign[index] -1;
    }
};

void Energyscr::change_chargesign(int &i,int &j)
{
    chargesign[i]=-1;
    chargesign[j]=-1;
};

bool Energyscr::populationValid() const
{
    // Check whether v_local at each site meets population validity constraints
    // Note that v_local components have flipped signs from E_sys

    bool valid;
    int collect_invalid=0;
    const float &zero_equiv = params.POP_STABILITY_ERR;
    for (int i=0; i<chargesign.size(); i++) {
        valid = ((chargesign[i] == -1 && -v_local[i] + params.mu < -zero_equiv)   // DB- condition
                 || (chargesign[i] == 1  && -v_local[i] + params.mu_p > zero_equiv)  // DB+ condition
                 || (chargesign[i] == 0  && -v_local[i] + params.mu > zero_equiv
                     && -v_local[i] + params.mu_p < -zero_equiv));
        collect_invalid += valid;
        if (!valid) {
            return false;
        }
    }
    //std::cout << collect_invalid << std::endl;
    return true;
}

std::pair<int, std::vector<float>> Energyscr::populationValid_counter()
{
    // Check whether v_local at each site meets population validity constraints
    // Note that v_local components have flipped signs from E_sys

    bool valid;
    int collect_invalid=0;
    std::vector<float> data_collect;
    const float &zero_equiv = params.POP_STABILITY_ERR;
    for (int i=0; i<chargesign.size(); i++) {
        valid = ((chargesign[i] == -1 && -v_local[i] + params.mu < -zero_equiv)   // DB- condition
                 || (chargesign[i] == 1  && -v_local[i] + params.mu_p > zero_equiv)  // DB+ condition
                 || (chargesign[i] == 0  && -v_local[i] + params.mu > zero_equiv
                     && -v_local[i] + params.mu_p < -zero_equiv));
        if (valid!=1)
        {
            data_collect.push_back(-v_local[i] + params.mu);
            //chargesign[i] =
            collect_invalid += 1;}
    }

    //std::cout << collect_invalid << std::endl;
    return std::make_pair(collect_invalid, data_collect);
}

void Energyscr::potentials() {
    v_ij.resize(locationeuc.size(), locationeuc.size());
    for (unsigned int i = 0; i < v_ij.size1(); i++) {
        v_ij(i,i)=0;

        for (unsigned int j = i+1; j < v_ij.size2(); j++) {
            //std::cout << "i: " << i << " | j: " << j << std::endl;
            v_ij(i, j) = energy_screened(db_r, i, j);
            v_ij(j, i) =  v_ij(i, j);
        }
    }
   //std::cout << v_ij << std::endl;
};



void Energyscr::total_energy(){
    std::vector<float> m(v_ij.size1(),0);
   for (unsigned int i = 0; i < v_ij.size1(); i++)
   {
       float collect = 0;
       for (unsigned int j = 0; j < v_ij.size2(); j++)
       {
           //std::cout << "i: " << i << " | j: " << j << std::endl;
           //std::cout << "chargesign["<<j<<"] :" << chargesign[j] << std::endl;
           collect += v_ij(j,i) * chargesign[j];
       }
       m[i] = collect;
       //std::cout << "potential at " << i << ": " << collect << std::endl;
   }
   v_local = m;
};

void Energyscr::system_energy(){
    float collect_all=0;
    for (unsigned int i = 0; i < v_ij.size1(); i++)
    {

        for (unsigned int j = i+1; j < v_ij.size2(); j++)
        {
            //std::cout << "i: " << i << " | j: " << j << std::endl;
            //std::cout << "chargesign["<<j<<"] :" << chargesign[j] << std::endl;
            collect_all += v_ij(j,i) * chargesign[j]*chargesign[i];
        }
        //std::cout << "potential at " << i << ": " << collect << std::endl;
    }
    std::cout << "system energy:" << collect_all << std::endl;
};

//std::pair<float, std::vector<int>> Energyscr::shortestPath(int &u, int &v, int &k) {
//    // Table to be filled up using DP. The value sp[i][j][e] will store
//    // weight of the shortest path from i to j with exactly k edges
//    int size = v_ij.size1();
//    float sp[size][size][k + 1];
//    std::vector<int> path[size][size][k + 1];
//
//    std::vector<int> sampling(size);
//    std::iota(sampling.begin(), sampling.end(),0);
//
//
//        for (int i = 0; i < size; i++)  // for source
//        {
//            //std::cout << "i: " << i << std::endl;
//            for (int j = 0; j < size; j++) // for destination
//            {
//                //std::cout << "j: " << j << std::endl;
//
//
//                //path[i][j] = ;
//                //std::fill(path[i][j].begin(),path[i][j].end(),0);
//                std::vector<int> path_cer;
//                for (int e = 0; e <= size; e++) {
//                    path[i][j][e] = std::vector<int>{0};
//                    sp[i][j][e] = INF;
//                    //std::string path_cer = "0";
//                    // from base cases
//                    if (e == 0 && i == j) {
//                        sp[i][j][e] = 0;
//                        //path[i][j][e].push_back(i);
//                    }
//
//                    if (e == 1 && v_ij(i, j) != INF) {
//                        sp[i][j][e] = v_ij(i, j);
//                        path[i][j][e] = std::vector<int>{i};
//                        sampling.erase(std::remove(sampling.begin(), sampling.end(),i), sampling.end());
//                        sampling.erase(std::remove(sampling.begin(), sampling.end(),j), sampling.end());
//                    }
//
//                    //go to adjacent only when number of edges is more than 1
//                    std::vector<float> collect;
//                    std::vector<int> collect_index;
//                    if (e > 1) {
//                        //path_cer.push_back(i);
//
//
//
//                        for (auto it = sampling.begin(); it != sampling.end(); it++) {
//                            //for (int a = i+1; a < size; a++) {
//                            //std::cout << *it << std::endl;
//                            // There should be an edge from i to a and a
//                            // should not be same as either i or j
//                            if ((v_ij(i, *it) != INF) && (i != *it) &&
//                                (j != *it) && (sp[*it][j][e - 1] != INF) && (sp[*it][j][e - 1] >= 0)) {
//
////                                std::cout << "a: " << v_ij(i, *it) << std::endl;
////                                std::cout << "b: " << sp[*it][j][e-1] << std::endl;
////                                std::cout << "i: " << i << std::endl;
////                                std::cout << "j: " << j << std::endl;
////                                std::cout << "e: " << e << std::endl;
////                                std::cout << *it << std::endl;
//                                //std::cout << "Summe: " << v_ij(i, *it) + sp[*it][j][e-1] << std::endl;
//                                collect.push_back(v_ij(i, *it) + sp[*it][j][e - 1]);
//                                collect_index.push_back(*it);
//
//
////                                sp[i][j][e] = std::fmin(sp[i][j][e], v_ij(i, a) +
////                                                                     sp[a][j][e-1]);
//
//
//
//                            }
//                            else
//                                continue;
//
//                        }
//
//                        auto result = std::min_element(collect.begin(), collect.end());
//                        int min_index = std::distance(collect.begin(), result);
//                        //std::cout << *result << std::endl;
//                        int value = collect_index[min_index];
//                        //std::cout << "e: " << e << std::endl;
//                        //std::cout << value << std::endl;
//                        //std::cout << collect[min_index] << std::endl;
//                        path_cer.push_back(value);
//                        //sampling.erase(std::remove(sampling.begin(), sampling.end(),value), sampling.end());
//                        //std::cout << "j: " << j << " | i: " << i << " | "<< e << " | " << collect_index[min_index] << std::endl;
////                         for (auto it = path_cer.begin(); it!= path_cer.end(); it++)
////                        {
////                            std::cout << *it << std::endl;
////                            std::cout << e << std::endl;
////                        }
//
//                        //path_cer.push_back(j);
//
//                        //path[i][j][e-1].insert(path[i][j][e-1].end(),path_cer.begin(), path_cer.end());
//                        //path[i][j][e] = collect_index;
//                        path[i][j][e] = path_cer;
//                        sp[i][j][e] = collect[min_index];
//                        //std::cout << path[i][j][e][0] << std::endl;
//                    }
//                    else
//                        continue;
//
//                }
//
//
//            }
//            std::cout << path[i][1][2][0] << std::endl;
//        };
//
//
////std::cout << sizeof(path) << std::endl;
////            for (auto it = path[0][5][2].begin(); it!= path[0][5][2].end(); it++)
////            {
////                std::cout << *it << std::endl;
////                std::cout << 1 << std::endl;
////                //std::cout <<  << std::endl;
////            }
//
////    for (auto i = path_cer.begin(); i!= path_cer.end(); i++)
////    {
////        std::cout << *i;
////    }
//        return std::make_pair(sp[u][v][k], path[u][v][k]);
//};
//
//
//template <size_t N>
//int Sum(const int (&intArray)[N])
//{
//    cout << "Array size in function: " << N << endl;
//    return std::accumulate(std::begin(intArray), std::end(intArray), 0);
//}
//
//bool check (int &dim, int &size, const int &m[][][][], int &i, int &j, int &loopind) {
//
//    return include;
//};

std::pair<float, std::vector<int>> Energyscr::shortestPath(int &u, int &v, int &k) {

#include <iostream>
#include <climits>

#include <algorithm>

    using namespace std;
    int size = v_ij.size1();
    //std::cout << size << std::endl;
// Define number of vertices in the graph and infinite value

#define INF INT_MAX
    float sp[size][size][k+1];
    std::vector<int> path[size][size][k+1];
// A Dynamic programming based function to find the shortest path from
// u to v with exactly k edges.
        // Table to be filled up using DP. The value sp[i][j][e] will store
        // weight of the shortest path from i to j with exactly k edges


        // Loop for number of edges from 0 to k
        for (int e = 1; e <= k; e++)
        {
            //std::cout << "e: " << e << std::endl;
            for (int i = 0; i < size; i++)  // for source
            {
                for (int j = 0; j < size; j++) // for destination
                {
                    // initialize value
                    sp[i][j][e] = INF;

                    // from base cases
                    if (e == 0 && i == j)
                        sp[i][j][e] = 0;
                    if (e == 1 && v_ij(i,j) != INF) {
                        sp[i][j][e] = v_ij(i, j);
                        path[i][j][e].push_back(i);
                        path[i][j][e].push_back(j);
                    }
                    //go to adjacent only when number of edges is more than 1
                    if (e > 1)
                    {
                        std::vector<int> index;
                        std::vector<float> collect;
                        for (int a = 0; a < size; a++)
                        {
                            // There should be an edge from i to a and a
                            // should not be same as either i or j

                            bool inor;
                            int variable=0;
                            inor = std::find(path[i][j][e-1].begin(), path[i][j][e-1].end(), a) !=  path[i][j][e-1].end();
                            if (inor == true)
                            {
                                variable +=1;
                            }

//                            std::cout << "e: " << e << std::endl;
//                            std::cout << "variable: " << variable << std::endl;
//                            std::cout << "a: " << a << std::endl;


                                //std::cout << "a_new: " << a << std::endl;
                                //std::cout << a << std::endl;
                                if (v_ij(i, j) != INF && i != a &&
                                    j != a && sp[a][j][e - 1] != INF && variable ==0) {
                                     {
                                        index.push_back(a);
                                        collect.push_back(v_ij(i, a) + sp[a][j][e - 1]);

                                    }

                                    //sp[i][j][e] = min(sp[i][j][e], v_ij(i, a) + sp[a][j][e - 1]);
                                }
                        }



                        auto result = std::min_element(collect.begin(), collect.end());
                        //for (auto it = result.begin(); it!=result.end(); it++)
                        //{
                           //std::cout << *result << std::endl;
                        //}
                        int min_index = std::distance(collect.begin(), result);
                        //std::cout << min_index << std::endl;
                        int value = index[min_index];
                        //path[i][j][e] = path[i][j][e-1];
                        for (int z = 0; z < path[i][j][e-1].size();z++)
                        {
                            //std::cout << z << std::endl;
                            //std::cout << path[i][j][e-1][z] << std::endl;
                            path[i][j][e].push_back(path[i][j][e-1][z]);
                        }
                        path[i][j][e].push_back(value);

                        sp[i][j][e] = collect[min_index];

                    }
                };
//                std::cout << "i: " << i << std::endl;
//                std::cout << sp[i][1][e] << std::endl;
            };
//            std::cout << "e: " << e << std::endl;
//            std::cout << sp[2][8][e] << std::endl;
//            std::cout << "e: " << e << std::endl;
//            std::cout << sp[2][8][e] << std::endl;
            //std::cout << path[2][8][e][1] << std::endl;
        };
        //std::cout << sp[1][6][4] << std::endl;
        return std::make_pair(sp[u][v][k],path[u][v][k]);
    //return std::make_pair(sp[u][v][k],path[u][v][k]);
    };







