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
void config::location_infile(std::ofstream &file)
{
    for (const auto& i : locationeuc)
    {
        file << i[0]*1*std::pow(10,10) << "; " << i[1]*std::pow(10,10) << std::endl;
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
void config::get_chargesign() const
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
void config::get_index() const
{
    std::cout << "Index: " << chargeindex << std::endl;
};

// charge index is increased by one
void config::increase_step()
{
    chargeindex += 1;
};


void config::identify_outsider()
{

    int i_up=0;
    float up=0;
    int i_down=0;
    float down=MAX_FLOAT;
    int i_right=0;
    float right=0;
    int i_left=0;
    float left=MAX_FLOAT;

    for (int i = 0; i < locationeuc.size(); i++)
    {
        if (locationeuc[i][0] < left)
        {
            left= locationeuc[i][0];
            i_left = i;

        }

        if (locationeuc[i][0] > right)
        {
            right= locationeuc[i][0];
            i_left = i;

        }

        if (locationeuc[i][1] > up)
        {
            up =  locationeuc[i][1];
            i_up = i;

        }

        if (locationeuc[i][1] < down)
        {
            down =  locationeuc[i][1];
            i_down = i;

        }
    }

    outsider = {i_left,i_right,i_down,i_up};

    for (int i = 0; i < locationeuc.size(); i++)
    {
        if (locationeuc[i][0] == left && i!=i_left)
        {
            outsider.push_back(i);
        }

        if (locationeuc[i][0] == right && i!=i_right)
        {
            outsider.push_back(i);
        }

        if (locationeuc[i][1] == up && i!=i_up)
        {
            outsider.push_back(i);
        }

        if (locationeuc[i][1] == down && i!=i_down)
        {
            outsider.push_back(i);
        }
    }
//    for (auto it = outsider.begin();it != outsider.end(); it++)
//    {
//        std::cout << *it << std::endl;
//    }

}

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



// analogy to the max-min diversity problem (algorithm is inspired from https://doi.org/10.1016/j.cor.2008.05.011)
void Energyscr::find_new_best_neighbor_GRC(
    std::vector<int>& index_db)  // input is a vector with the indices of all already SiDBs with -1
{
    float max_value =0;
    int   index_l = -1;
    int   count   = 0;
    for (auto it = chargesign.begin(); it!=chargesign.end();it++)
    {
        if ((*it)!=0)
            continue;
        count += 1;
        int l = static_cast<int>(std::distance(chargesign.begin(), it));


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

        else if (count > 1)
        {
            if (distance_min > max_value)
            {
                max_value = distance_min;
                index_l   = l;
            }
            else if (distance_min == max_value)
            {
               if (sum_distance(index_db, index_l, l))
                    //if (1)
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

    for (auto it = chargesign.begin(); it!=chargesign.end();it++)
    {
        if ((*it)!=0)
            continue;

        count1 += 1;
        int l = static_cast<int>(std::distance(chargesign.begin(), it));


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
            const float& zero_equiv      = Params::POP_STABILITY_ERR;
            if (distance_min >= 0.7 * max_value)
            {
                random.push_back(l);
            }
        }
    };

//    for (auto it = random_3.begin(); it<random_3.end(); it++)
//    {
//        std::cout << *it << std::endl;
//    }
    int random_index = rand() % random.size();
    if ((random_index > -1) && (random.size()!=0))
    {
        int random_element         = random[random_index];
        //std::cout << "random_3 element: " << random_element << std::end;
        chargesign[random_element] = -1;
        index_db.push_back(random_element);
        potential_energy += -v_local[random_element];
        for (unsigned int i = 0; i < v_ij.size1(); i++)
            {
//                if ((std::find(index_db.begin(), index_db.end(), i) !=
//                     index_db.end()))  // take no l-index which is already occupied
//                {
//                    continue;
//                }
                v_local[i] += -v_ij(random_element, i);
            }


        };
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


bool Energyscr::populationValid() const
{
    // Check whether v_local at each site meets population validity constraints
    // this includes the checking if the charge state is valid and if no hop can be performed which reduces the energy1
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
        //std::vector<int> charge_sign_saved = chargesign;
        //float E_ori = this->system_energy_vec(charge_sign_saved);

        int dn_i = (chargesign[i] == -1) ? 1 : -1;
        int dn_j = -dn_i;

            //        return v_local[i] * dn_i + v_local[j] * dn_j +
//               v_ij(i, j) * ((chargesign[i] + dn_i) * (chargesign[j] + dn_j) - chargesign[i] * chargesign[j]);

       float E_fast = v_local[i]*dn_i + v_local[j]*dn_j - v_ij(i,j)*1;

            return E_fast;
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
            return v_local[i]*dn_i + v_local[j]*dn_j - v_ij(i,j)*1;
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
                    chargesign[i] = 0;
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

            if (db_r(i, j) > 500000 * std::pow(10, -9))
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

// calculate v_local for each SiDB loc for one specific charge distribution
void Energyscr::total_energy()
{
    std::vector<float> m(v_ij.size1(), 0);
    for (unsigned int i = 0; i < v_ij.size1(); i++)
    {
        float collect = 0;
        for (unsigned int j = 0; j < v_ij.size2(); j++)
        {
            float energy = v_ij(j, i) * static_cast<float>(chargesign[j]);
            collect += energy;
        }
        m[i] = collect;
    }
    v_local = m;
};


// calculate v_local for each SiDB loc for one specific charge distribution
float Energyscr::total_energy_EQ(int &index)
{

    for (unsigned int pos= 0; pos < v_ij.size1(); pos++)
    {
        v_local[pos] = -v_ij(index, pos);
    }
};

//float Energyscr::system_energy_new() const
//{
//    float collect_all = 0;
//    for (unsigned int i = 0; i < v_ij.size1(); i++)
//    {
//
//        for (unsigned int j = i + 1; j < v_ij.size2(); j++)
//        {
//            collect_all += v_ij(j, i) * static_cast<float>(chargesign[j] * chargesign[i]);
//        }
//    }
//    return collect_all;
//};

float Energyscr::system_energy() const
{
    float collect_all=0;
    for (unsigned int i = 0; i < v_ij.size1(); i++)
    {

            collect_all += 0.5*v_local[i] * chargesign[i];
    }

    return collect_all;
};



