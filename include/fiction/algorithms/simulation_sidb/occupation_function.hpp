//
// Created by Jan Drewniok on 13.01.23.
//

#ifndef FICTION_OCCUPATION_FUNCTION_HPP
#define FICTION_OCCUPATION_FUNCTION_HPP

#include <map>
#include <cmath>

double occu_prop(std::map<double, uint64_t>& energy_deg, const double &temp)
{
    double part_func = 0;

    for (const auto& it : energy_deg)
    {
        part_func += static_cast<double>(it.second) * std::exp(-it.first * 12000 / temp);
    }

    auto it = energy_deg.begin();

    return static_cast<double>(it->second) * std::exp(-it->first * 12000 / temp) / part_func;

};


#endif  // FICTION_OCCUPATION_FUNCTION_HPP
