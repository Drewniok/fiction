//
// Created by jan-d on 02.07.2022.
//
#include "simu.h"

#include <cmath>
#include <iostream>



namespace Params
{
    // lattice dimensions
    static constexpr float a = 3.84f * 1E-10;
    static constexpr float b = 2.25f * 1E-10;
    static constexpr float c = 7.68f * 1E-10;

    // physical parameters
    static constexpr float epsilon        = 8.854f * 1E-12;
    static constexpr float epsilon_screen = 5.6;
    static constexpr float k              = 1.0f / (4 * 3.141592653 * epsilon * epsilon_screen);
    static constexpr float e              = 1.602f * 1E-19;
    static constexpr float mu             = -0.32;
    static constexpr float mu_p           = mu - 0.59f;
    static constexpr float tf             = 5.0f * 1E-9;

    static constexpr float POP_STABILITY_ERR = 1E-6;

};

// siqad specific lattice coordinates are converted in euclidean coordinates
void config::toeuc()
{
    locationeuc.clear();
    locationeuc.reserve(locationind.size());

    for (const auto& i : locationind)
    {
        locationeuc.push_back({(static_cast<float>(i[0])) * Params::a, (static_cast<float>(i[1])) * Params::c + (static_cast<float>(i[2])) * Params::b});
    }
};

// print euclidean coordinates of each single placed SiDB
void config::get_locationeuc()
{
    for (const auto& i : locationeuc)
    {
        std::cout << i[0] << " | " << i[1] << std::endl;
    }
};

// charge configuration vector (-1,0,1,-1,...) is transformed to an unique index
void config::chargeconf_to_index(const int& base)
{
    chargeindex = 0;
    chargemax   = static_cast<unsigned int>(std::pow(base, chargesign.size()) - 1);
    for (unsigned int i = 0; i < chargesign.size(); i++)
    {
        chargeindex += static_cast<unsigned int>((chargesign[i] + 1) * pow(base, chargesign.size() - i - 1));
    }
};

// index is converted back to the charge config. vector
void config::index_to_chargeconf(const int& base)
{
    int n_sidb       = locationind.size();
    int i            = n_sidb - 1;  // i is the number of the n_dbs
    int chargeindex1 = chargeindex;
    while (chargeindex1 > 0)
    {
        div_t d;  // Structure to represent the result value of an integral division performed by function div.
        d               = div(chargeindex1, base);
        chargeindex1    = d.quot;
        chargesign[i--] = d.rem - 1;
    }
};

// print the charge config. vector
void config::get_chargesign()
{
    for (unsigned int i = 0; i < chargesign.size(); i++)
    {
        if (i == 0)
        {
            std::cout << "charge configuration: |" << chargesign[i] << "|";
        }
        else if (i != chargesign.size() - 1)
        {
            std::cout << chargesign[i] << "|";
        }
        else
        {
            std::cout << chargesign[i] << "|" << std::endl;
        }
    }
};

// print the charge index
void config::get_index()
{
    std::cout << "Index: " << chargeindex << std::endl;
};

// charge index is increased by one
void config::increase_step()
{
    chargeindex += 1;
};

// potential calculation
float energy_screened(FPMat& m, unsigned int i, unsigned int j)
{
    float potential = Params::k * 1 / m(i, j)

                      * std::exp(-m(i, j) / Params::tf) * Params::e;  // eV
    return potential;
};

// distance matrix is calculated
void config::distance()
{
    db_r.resize(locationeuc.size(), locationeuc.size());
    for (unsigned int i = 0; i < db_r.size1(); i++)
    {
        db_r(i, i) = 0;
        for (unsigned int j = i + 1; j < db_r.size2(); j++)
        {
            db_r(i, j) = std::sqrt(std::pow(locationeuc[i][0] - locationeuc[j][0], 2) +
                                   std::pow(locationeuc[i][1] - locationeuc[j][1], 2));
            // symmetry is introduced
            db_r(j, i) = db_r(i, j);
        }
    }
};

// define distance between two SiDBs as separate variable
void Energyscr::get_distance(int& i, int& j)
{
    distance_value = db_r(i, j);
};

// define potential between two SiDBs as separate variable
void Energyscr::get_potential(int& i, int& j)
{
    potential_value = v_ij(i, j);
};

// check if the summed distance between the current SiDB (now) and the other available ones is larger than between the
// old one and the others.
bool Energyscr::sum_distance(std::vector<int>& input, int& old, int& now)
{
    float sum_old = 0;
    float sum_now = 0;
    for (int j = 0; j < input.size(); j++)
    {
        sum_old += db_r(j, old);
        sum_now += db_r(j, now);
    }
    return sum_now > sum_old;
};

std::vector<int> Energyscr::search_same_potential(std::vector<int>& index_db)
{
    float            min_v     = 100000;
    float            min_far_v = 10000;
    int              index     = -1;
    int              index_inc = -1;
    int              count     = 0;
    std::vector<int> possible_pos_correct;
    std::vector<int> possible_pos_incorrect;
    std::vector<int> possible_pos_index;
    std::vector<int> possible_pos_index_inc;

    for (int l = 0; l < db_r.size1(); l++)
    {
        float left             = 0.9;
        float right            = 1.1;
        int   number_correct   = 0;
        int   number_incorrect = 0;
        float pot_sum          = 0;
        float pot_sum_incor    = 0;
        for (int n = 0; n < index_db.size(); n++)
        {
            if (((left * potential_value <= v_ij(l, index_db[n])) &&
                 (v_ij(l, index_db[n]) <= right * potential_value)) &&
                (v_ij(l, index_db[n]) != 0))
            {
                number_correct += 1;
                pot_sum += v_ij(l, index_db[n]) / index_db.size();
            }
            if (((v_ij(l, index_db[n]) < left * potential_value)) && (v_ij(l, index_db[n]) != 0))
            {
                number_incorrect += 1;
                pot_sum_incor += v_ij(l, index_db[n]) / index_db.size();
            };
        };
        if ((number_correct == index_db.size()) && std::abs(potential_value - pot_sum) < min_v)
        {
            count += 1;
            index = l;
            min_v = std::abs(potential_value - pot_sum);
        }

        else if ((number_correct != index_db.size()) && (number_incorrect == index_db.size()) &&
                 std::abs(potential_value - pot_sum_incor) < min_far_v)
        {
            index_inc = l;
            min_far_v = std::abs(potential_value - pot_sum_incor);
        }
    }
    if (count > 0 && index != -1)
    {
        chargesign[index] = -1;
        index_db.push_back(index);
    }
    else if (index_inc != -1)
    {
        chargesign[index_inc] = -1;
        index_db.push_back(index_inc);
    }
    return index_db;
};

std::vector<int> Energyscr::find_new_best_neighbor(std::vector<int>& index_db)
{
    float first_distance;
    float max_value;
    int   index_l = -1;
    int   count   = 0;
    for (int l = 0; l < db_r.size1(); l++)
    {
        if ((std::find(index_db.begin(), index_db.end(), l) != index_db.end()))
        {
            continue;
        }
        else
        {
            count += 1;
        }
        float distance_min = db_r(l, index_db[0]);
        for (int n = 0; n < index_db.size(); n++)
        {
            if (db_r(l, index_db[n]) < distance_min)
            {
                distance_min = db_r(l, index_db[n]);
            }
        }

        if (count == 1)
        {
            max_value = distance_min;
            index_l   = l;
        }

        if (count > 1)
        {
            if (distance_min > max_value)
            {
                max_value = distance_min;
                index_l   = l;
            }
            else if (distance_min == max_value)
            {
                if (sum_distance(index_db, index_l, l))
                {
                    max_value = distance_min;
                    index_l   = l;
                }
            }
        }
    }

    if (index_l > -1)
    {
        chargesign[index_l] = -1;
        index_db.push_back(index_l);
    }

    return index_db;
};

// analogy to the max-min diversity problem (algorithm is inspired from https://doi.org/10.1016/j.cor.2008.05.011)
void Energyscr::find_new_best_neighbor_GRC(
    std::vector<int>& index_db)  // input is a vector with the indices of all already SiDBs with -1
{
    float max_value;
    int   index_l = -1;
    int   count   = 0;
    for (int l = 0; l < db_r.size1(); l++)
    {

        if ((std::find(index_db.begin(), index_db.end(), l) !=
             index_db.end()))  // take no l-index which is already occupied
        {
            continue;
        }
        else
        {
            count += 1;
        }


        float distance_min = db_r(l, index_db[0]);
        for (int n : index_db)
        {
            if (db_r(l, n) < distance_min)
            {
                distance_min = db_r(l, n);
            }
        }

        if (count == 1)
        {
            max_value = distance_min;
            index_l   = l;
        }

        if (count > 1)
        {
            if (distance_min > max_value)
            {
                max_value = distance_min;
                index_l   = l;
            }
            else if (distance_min == max_value)
            {
                if (sum_distance(index_db, index_l, l))
                {
                    max_value = distance_min;
                    index_l   = l;
                }
            }
        }
    }
//std::cout << "max value" << index_l << std::endl;
    std::vector<int> random;
    int              count1 = 0;
    for (int l = 0; l < db_r.size1(); l++)
    {
        if ((std::find(index_db.begin(), index_db.end(), l) != index_db.end()))
        {
            continue;
        }

        else
        {
            count1 += 1;
        }

        float distance_min = db_r(l, index_db[0]);
        for (int n : index_db)
        {
            if (db_r(l, n) < distance_min)
            {
                distance_min = db_r(l, n);
                // std::cout << distance_min << std::endl;
            }
        };


        if (count1 > 0)
        {
            if (distance_min >= 0.8 * max_value)
            {
                random.push_back(l);
            }
        }
    };

//    for (auto it = random.begin(); it<random.end(); it++)
//    {
//        std::cout << *it << std::endl;
//    }
    int random_index = rand() % random.size();
    if ((random_index > -1) && (random.size()!=0))
    {
        int random_element         = random[random_index];
        //std::cout << "random element: " << random_element << std::end;
        chargesign[random_element] = -1;
        index_db.push_back(random_element);
    }
}

std::vector<int> Energyscr::search_same_distance_new(std::vector<int>& index_db)
{
    float            min_v     = 100000;
    float            min_far_v = 10000;
    int              index     = -1;
    int              index_inc = -1;
    int              count     = 0;
    std::vector<int> possible_pos_correct;
    std::vector<int> possible_pos_incorrect;
    std::vector<int> possible_pos_index;
    std::vector<int> possible_pos_index_inc;

    for (int l = 0; l < db_r.size1(); l++)
    {

        float left             = 0.7;
        float right            = 1.3;
        int   number_correct   = 0;
        int   number_incorrect = 0;
        float pot_sum          = 0;
        float pot_sum_incor    = 0;
        for (int n = 0; n < index_db.size(); n++)
        {

            if (((left * distance_value <= db_r(l, index_db[n])) && (db_r(l, index_db[n]) <= right * distance_value)) &&
                (db_r(l, index_db[n]) != 0))
            {
                number_correct += 1;
                pot_sum += db_r(l, index_db[n]) / index_db.size();
            }
            if (((v_ij(l, index_db[n]) < left * distance_value)) && (db_r(l, index_db[n]) != 0))
            {
                number_incorrect += 1;
                pot_sum_incor += db_r(l, index_db[n]) / index_db.size();
            };
        };

        if ((number_correct == index_db.size()) && std::abs(distance_value - pot_sum) < min_v)
        {
            count += 1;
            index = l;
            min_v = std::abs(distance_value - pot_sum);
        }

        else if ((number_correct != index_db.size()) && (number_incorrect == index_db.size()) &&
                 std::abs(distance_value - pot_sum_incor) < min_far_v)
        {
            index_inc = l;
            min_far_v = std::abs(distance_value - pot_sum_incor);
        }
    }

    if (count > 0 && index != -1)
    {
        chargesign[index] = -1;
        index_db.push_back(index);
    }

    else if (index_inc != -1)
    {
        chargesign[index_inc] = -1;
        index_db.push_back(index_inc);
    }

    return index_db;
};

std::vector<int> Energyscr::search_same_distance(std::vector<int>& index_db)
{
    for (int l = 0; l < db_r.size1(); l++)
    {
        float left           = 0.7;
        float right          = 1.3;
        int   number_correct = 0;
        for (int n : index_db)
        {

            if (((left * distance_value <= db_r(l, n)) && (db_r(l, n) <= right * distance_value)) ||
                (db_r(l, n) >= left * distance_value))
            {
                number_correct += 1;
            }
        }

        if (number_correct == index_db.size())
        {
            chargesign[l] = -1;
            index_db.push_back(l);
        }
    }
    return index_db;
};

// find perturbers in a circuit

std::vector<int> Energyscr::find_perturber()
{
    std::vector<int> perturber_collect;
    for (int i = 0; i < db_r.size1(); i++)
    {
        int   c      = 0;
        float threas = 3.3 * std::pow(10, -9);
        for (int j = 0; j < db_r.size1(); j++)
        {
            if (db_r(i, j) > threas && db_r(i, j) != 0)
            {
                c += 1;
            }
        }
        if (c == db_r.size1() - 1)
        {
            chargesign[i] = -1;
            perturber_collect.push_back(i);
        }
    }
    return perturber_collect;
};

std::vector<int> Energyscr::find_perturber_alternative()
{
    std::vector<int> perturber_collect;
    for (int i = 0; i < db_r.size1(); i++)
    {
        int                right   = 0;
        int                left    = 0;
        int                top     = 0;
        int                bottom  = 0;
        int                c       = 0;
        int                d       = 0;
        float              threas  = 6 * 1E-9;
        float              threas1 = 1.4 * 1E-9;
        std::vector<float> collect_index;
        for (int j = 0; j < db_r.size1(); j++)
        {
            if ((db_r(i, j) < threas) && (db_r(i, j) != 0))
            {
                d += 1;
                if (db_r(i, j) > threas1)
                {
                    c += 1;
                }

                if ((locationeuc[i][0] > locationeuc[j][0]))
                {
                    left += 1;
                }

                if ((locationeuc[i][0] < locationeuc[j][0]))
                {
                    right += 1;
                }

                if ((locationeuc[i][1] < locationeuc[j][1]))
                {
                    top += 1;
                }

                if (locationeuc[i][1] > locationeuc[j][1])
                {
                    bottom += 1;
                }
            }
        }

        if (((c == top) || (c == bottom) || (c == right) || (c == left)) && (c != 0) && (c == d))
        {
            chargesign[i] = -1;
            perturber_collect.push_back(i);
        }
    }
    return perturber_collect;
};

void Energyscr::set_charge()
{
    float min   = 10.0;
    int   index = 0;
    for (int i = 0; i < v_local.size(); i++)
    {
        // std::cout << chargesign[i] << std::endl;
        if ((-v_local[i] < min) && (chargesign[i] != -1))
        {
            min   = -v_local[i];
            index = i;
        }
    }
    if (chargesign[index] > -1)
    {
        chargesign[index] = chargesign[index] - 1;
    }
};

void Energyscr::change_chargesign(int& i, int& j)
{
    chargesign[i] = -1;
    chargesign[j] = -1;
};

void Energyscr::change_chargesign_one(std::vector<int>& vector)
{
    for (auto it = vector.begin(); it != vector.end(); it++)
    {
        chargesign[*it] = (chargesign[*it] == -1 ? 0 : -1);
        // chargesign[*it] = -1;
    }
};

bool Energyscr::populationValid() const
{
    // Check whether v_local at each site meets population validity constraints
    // this includes the checking if the charge state is valid and if no hop can be performed which reduces the energy
    // even further
    bool         valid;
    int          collect_invalid = 0;
    const float& zero_equiv      = Params::POP_STABILITY_ERR;
    for (int i = 0; i < chargesign.size(); i++)
    {
        valid = ((chargesign[i] == -1 && -v_local[i] + Params::mu < zero_equiv)       // DB- condition
                 || (chargesign[i] == 1 && -v_local[i] + Params::mu_p > -zero_equiv)  // DB+ condition
                 || (chargesign[i] == 0 && -v_local[i] + Params::mu > -zero_equiv &&
                     -v_local[i] + Params::mu_p < zero_equiv));
        collect_invalid += valid;
        if (!valid)
        {
            return false;
        }
    }

    auto hopDel = [this](const int& i, const int& j)
    {
        int dn_i = (chargesign[i] == -1) ? 1 : -1;
        int dn_j = -dn_i;
//        return v_local[i] * dn_i + v_local[j] * dn_j +
//               v_ij(i, j) * ((chargesign[i] + dn_i) * (chargesign[j] + dn_j) - chargesign[i] * chargesign[j]);

        // this is from Siqad. However, it is quite likely wrong
        return v_local[i]*dn_i + v_local[j]*dn_j - v_ij(i,j);
    };

    for (unsigned int i = 0; i < chargesign.size(); i++)
    {
        // do nothing with DB+
        if (chargesign[i] == 1)
            continue;

        for (unsigned int j = 0; j < chargesign.size(); j++)
        {

            // attempt hops from more negative charge states to more positive ones
            float E_del = hopDel(i, j);
            if ((chargesign[j] > chargesign[i]) && (E_del < -zero_equiv))
            {
                return false;
            }
        }
    }
    return true;
}

std::tuple<int, std::vector<int>, int> Energyscr::populationValid_counter()
{

    bool             valid;
    int              collect_invalid  = -1;
    int              collect_unstable = -1;
    std::vector<int> data_collect;
    const float&     zero_equiv = Params::POP_STABILITY_ERR;
    for (int i = 0; i < chargesign.size(); i++)
    {
        valid = ((chargesign[i] == -1 && -v_local[i] + Params::mu < zero_equiv)       // DB- condition
                 || (chargesign[i] == 1 && -v_local[i] + Params::mu_p > -zero_equiv)  // DB+ condition
                 || (chargesign[i] == 0 && -v_local[i] + Params::mu > -zero_equiv &&
                     -v_local[i] + Params::mu_p < zero_equiv));
        if (valid != 1)
        {
            data_collect.push_back(i);
            collect_invalid += 1;
            continue;
        }

        auto hopDel = [this](const int& i, const int& j)
        {
            int dn_i = (chargesign[i] == -1) ? 1 : -1;
            int dn_j = -dn_i;
//            return v_local[i] * dn_i + v_local[j] * dn_j +
//                   v_ij(i, j) * ((chargesign[i] + dn_i) * (chargesign[j] + dn_j) - chargesign[i] * chargesign[j]);
            return v_local[i]*dn_i + v_local[j]*dn_j - v_ij(i,j);
        };

        for (unsigned int i = 0; i < chargesign.size(); i++)
        {
            // do nothing with DB+
            if (chargesign[i] == 1)
                continue;

            for (unsigned int j = 0; j < chargesign.size(); j++)
            {

                // attempt hops from more negative charge states to more positive ones
                float E_del = hopDel(i, j);
                if ((chargesign[j] > chargesign[i]) && (E_del < -zero_equiv))
                {
                    // data_collect.push_back(i);
                    //chargesign[i] = 0;
                    chargesign[j] = -1;
                    collect_unstable += 1;
                }
            }
        }
    }
    return std::make_tuple(collect_invalid, data_collect, collect_unstable);
}

void config::potentials()
{
    v_ij.resize(locationeuc.size(), locationeuc.size());
    for (unsigned int i = 0; i < v_ij.size1(); i++)
    {
        v_ij(i, i) = 0;
        for (unsigned int j = i + 1; j < v_ij.size2(); j++)
        {

            if (db_r(i, j) > 50 * std::pow(10, -9))
            {
                v_ij(i, j) = 0;
            }
            else
            {
                v_ij(i, j) = energy_screened(db_r, i, j);
            }
            v_ij(j, i) = v_ij(i, j);
        }
    }
};

// calculate v_local for each SiDB location for one specific charge distribution
void Energyscr::total_energy()
{
    std::vector<float> m(v_ij.size1(), 0);
    for (unsigned int i = 0; i < v_ij.size1(); i++)
    {
        float collect = 0;
        for (unsigned int j = 0; j < v_ij.size2(); j++)
        {
            collect += v_ij(j, i) * static_cast<float>(chargesign[j]);
        }
        m[i] = collect;
    }
    v_local = m;
};

float Energyscr::system_energy()
{
    float collect_all = 0;
    for (unsigned int i = 0; i < v_ij.size1(); i++)
    {

        for (unsigned int j = i + 1; j < v_ij.size2(); j++)
        {
            collect_all += v_ij(j, i) * static_cast<float>(chargesign[j] * chargesign[i]);
        }
    }
    return collect_all;
};

// std::pair<float, std::vector<int>> Energyscr::shortestPath(int &u, int &v, int &k) {
//     // Table to be filled up using DP. The value sp[i][j][e] will store
//     // weight of the shortest path from i to j with exactly k edges
//     int size = v_ij.size1();
//     float sp[size][size][k + 1];
//     std::vector<int> path[size][size][k + 1];
//
//     std::vector<int> sampling(size);
//     std::iota(sampling.begin(), sampling.end(),0);
//
//
//         for (int i = 0; i < size; i++)  // for source
//         {
//             //std::cout << "i: " << i << std::endl;
//             for (int j = 0; j < size; j++) // for destination
//             {
//                 //std::cout << "j: " << j << std::endl;
//
//
//                 //path[i][j] = ;
//                 //std::fill(path[i][j].begin(),path[i][j].end(),0);
//                 std::vector<int> path_cer;
//                 for (int e = 0; e <= size; e++) {
//                     path[i][j][e] = std::vector<int>{0};
//                     sp[i][j][e] = MAX_FLOAT;
//                     //std::string path_cer = "0";
//                     // from base cases
//                     if (e == 0 && i == j) {
//                         sp[i][j][e] = 0;
//                         //path[i][j][e].push_back(i);
//                     }
//
//                     if (e == 1 && v_ij(i, j) != MAX_FLOAT) {
//                         sp[i][j][e] = v_ij(i, j);
//                         path[i][j][e] = std::vector<int>{i};
//                         sampling.erase(std::remove(sampling.begin(), sampling.end(),i), sampling.end());
//                         sampling.erase(std::remove(sampling.begin(), sampling.end(),j), sampling.end());
//                     }
//
//                     //go to adjacent only when number of edges is more than 1
//                     std::vector<float> collect;
//                     std::vector<int> collect_index;
//                     if (e > 1) {
//                         //path_cer.push_back(i);
//
//
//
//                         for (auto it = sampling.begin(); it != sampling.end(); it++) {
//                             //for (int a = i+1; a < size; a++) {
//                             //std::cout << *it << std::endl;
//                             // There should be an edge from i to a and a
//                             // should not be same as either i or j
//                             if ((v_ij(i, *it) != MAX_FLOAT) && (i != *it) &&
//                                 (j != *it) && (sp[*it][j][e - 1] != MAX_FLOAT) && (sp[*it][j][e - 1] >= 0)) {
//
//                                 std::cout << "a: " << v_ij(i, *it) << std::endl;
//                                 std::cout << "b: " << sp[*it][j][e-1] << std::endl;
//                                 std::cout << "i: " << i << std::endl;
//                                 std::cout << "j: " << j << std::endl;
//                                 std::cout << "e: " << e << std::endl;
//                                 std::cout << *it << std::endl;
//                                 //std::cout << "Summe: " << v_ij(i, *it) + sp[*it][j][e-1] << std::endl;
//                                 collect.push_back(v_ij(i, *it) + sp[*it][j][e - 1]);
//                                 collect_index.push_back(*it);
//
//
//                                 sp[i][j][e] = std::fmin(sp[i][j][e], v_ij(i, a) +
//                                                                      sp[a][j][e-1]);
//
//
//
//                             }
//                             else
//                                 continue;
//
//                         }
//
//                         auto result = std::min_element(collect.begin(), collect.end());
//                         int min_index = std::distance(collect.begin(), result);
//                         //std::cout << *result << std::endl;
//                         int value = collect_index[min_index];
//                         //std::cout << "e: " << e << std::endl;
//                         //std::cout << value << std::endl;
//                         //std::cout << collect[min_index] << std::endl;
//                         path_cer.push_back(value);
//                         //sampling.erase(std::remove(sampling.begin(), sampling.end(),value), sampling.end());
//                         //std::cout << "j: " << j << " | i: " << i << " | "<< e << " | " << collect_index[min_index]
//                         << std::endl;
//                          for (auto it = path_cer.begin(); it!= path_cer.end(); it++)
//                         {
//                             std::cout << *it << std::endl;
//                             std::cout << e << std::endl;
//                         }
//
//                         //path_cer.push_back(j);
//
//                         //path[i][j][e-1].insert(path[i][j][e-1].end(),path_cer.begin(), path_cer.end());
//                         //path[i][j][e] = collect_index;
//                         path[i][j][e] = path_cer;
//                         sp[i][j][e] = collect[min_index];
//                         //std::cout << path[i][j][e][0] << std::endl;
//                     }
//                     else
//                         continue;
//
//                 }
//
//
//             }
//             std::cout << path[i][1][2][0] << std::endl;
//         };
//
//
// std::cout << sizeof(path) << std::endl;
//             for (auto it = path[0][5][2].begin(); it!= path[0][5][2].end(); it++)
//             {
//                 std::cout << *it << std::endl;
//                 std::cout << 1 << std::endl;
//                 //std::cout <<  << std::endl;
//             }
//     for (auto i = path_cer.begin(); i!= path_cer.end(); i++)
//     {
//         std::cout << *i;
//     }
//         return std::make_pair(sp[u][v][k], path[u][v][k]);
// };
//
//
// template <size_t N>
// int Sum(const int (&intArray)[N])
//{
//     cout << "Array size in function: " << N << endl;
//     return std::accumulate(std::begin(intArray), std::end(intArray), 0);
// }
//
// bool check (int &dim, int &size, const int &m[][][][], int &i, int &j, int &loopind) {
//
//     return include;
// };

std::pair<float, std::vector<int>> Energyscr::shortestPath(int& u, int& v, int& k)
{

    // using namespace std;  // das macht man nicht :D
    int size = v_ij.size1();
    // std::cout << size << std::endl;
    // Define number of vertices in the graph and MAX_INTinite value

    float            sp[size][size][k + 1];
    std::vector<int> path[size][size][k + 1];
    // A Dynamic programming based function to find the shortest path from
    // u to v with exactly k edges.
    // Table to be filled up using DP. The value sp[i][j][e] will store
    // weight of the shortest path from i to j with exactly k edges

    // Loop for number of edges from 0 to k
    for (int e = 1; e <= k; e++)
    {
        // std::cout << "e: " << e << std::endl;
        for (int i = 0; i < size; i++)  // for source
        {
            for (int j = 0; j < size; j++)  // for destination
            {
                // initialize value
                sp[i][j][e] = MAX_FLOAT;

                // from base cases
                if (e == 0 && i == j)
                    sp[i][j][e] = 0;
                if (e == 1 && v_ij(i, j) != MAX_FLOAT)
                {
                    sp[i][j][e] = v_ij(i, j);
                    path[i][j][e].push_back(i);
                    path[i][j][e].push_back(j);
                }
                // go to adjacent only when number of edges is more than 1
                if (e > 1)
                {
                    std::vector<int>   index;
                    std::vector<float> collect;
                    for (int a = 0; a < size; a++)
                    {
                        // There should be an edge from i to a and a
                        // should not be same as either i or j

                        bool inor;
                        int  variable = 0;
                        inor =
                            std::find(path[i][j][e - 1].begin(), path[i][j][e - 1].end(), a) != path[i][j][e - 1].end();
                        if (inor == true)
                        {
                            variable += 1;
                        }

                        //                            std::cout << "e: " << e << std::endl;
                        //                            std::cout << "variable: " << variable << std::endl;
                        //                            std::cout << "a: " << a << std::endl;

                        // std::cout << "a_new: " << a << std::endl;
                        // std::cout << a << std::endl;
                        if (v_ij(i, j) != MAX_FLOAT && i != a && j != a && sp[a][j][e - 1] != MAX_FLOAT && variable == 0)
                        {
                            {
                                index.push_back(a);
                                collect.push_back(v_ij(i, a) + sp[a][j][e - 1]);
                            }

                            // sp[i][j][e] = min(sp[i][j][e], v_ij(i, a) + sp[a][j][e - 1]);
                        }
                    }

                    auto result = std::min_element(collect.begin(), collect.end());
                    // for (auto it = result.begin(); it!=result.end(); it++)
                    //{
                    // std::cout << *result << std::endl;
                    //}
                    int min_index = std::distance(collect.begin(), result);
                    // std::cout << min_index << std::endl;
                    int value = index[min_index];
                    // path[i][j][e] = path[i][j][e-1];
                    for (int z = 0; z < path[i][j][e - 1].size(); z++)
                    {
                        // std::cout << z << std::endl;
                        // std::cout << path[i][j][e-1][z] << std::endl;
                        path[i][j][e].push_back(path[i][j][e - 1][z]);
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
        // std::cout << path[2][8][e][1] << std::endl;
    };
    // std::cout << sp[1][6][4] << std::endl;
    return std::make_pair(sp[u][v][k], path[u][v][k]);
    // return std::make_pair(sp[u][v][k],path[u][v][k]);
};
