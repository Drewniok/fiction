//
// Created by Jan Drewniok on 13.05.24.
//

#ifndef FICTION_DETECT_BDL_WIRES_HPP
#define FICTION_DETECT_BDL_WIRES_HPP

#include "fiction/algorithms/simulation/sidb/detect_bdl_pairs.hpp"

#include <vector>
#include <optional>
#include <vector>
#include <set>

namespace fiction
{


template <typename Lyt>
using bdl_wire = std::vector<bdl_pair<Lyt>>;

enum class bdl_wire_direction
{
    TOP_DOWN,
    DOWN_TOP
};

template <typename Lyt>
std::vector<bdl_wire_direction> determine_wire_direction(const std::vector<bdl_pair<Lyt>> &input_pairs, const std::vector<bdl_wire<Lyt>> &bdl_wires)
{
    std::vector<bdl_wire_direction> directions{};
    for (const auto &bdl : input_pairs)
    {
        auto bdl_copy = bdl;
        bdl_copy.type = sidb_technology::cell_type::NORMAL;
        for (const auto &wire : bdl_wires)
        {
            auto it = std::find(wire.begin(), wire.end(), bdl_copy);
            if (it != wire.end())
            {
                bool down_to_top = false;
                for (const auto &wire_bdl : wire)
                {
                    if (bdl_copy.lower.y > wire_bdl.lower.y)
                    {
                        down_to_top = true;
                        break;
                    }
                }
                if (down_to_top)
                {
                    directions.push_back(bdl_wire_direction::DOWN_TOP);
                    break;
                }
                else
                {
                    directions.push_back(bdl_wire_direction::TOP_DOWN);
                    break;
                }
            }
        }
    }
    return directions;
}

template <typename Lyt>
std::optional<bdl_pair<Lyt>> bdl_parter_upper(const bdl_pair<Lyt>& given_bdl, const std::set<bdl_pair<Lyt>>& bdl_pairs)
{
    for (const auto& bdl : bdl_pairs)
    {
        if (sidb_nm_distance<Lyt>(Lyt{}, given_bdl.upper, bdl.lower) < 4 && given_bdl != bdl &&
            given_bdl.upper.y >= bdl.lower.y)
        {
            return std::optional<bdl_pair<Lyt>>(bdl);
        }
    }
    return std::nullopt;
};

template <typename Lyt>
std::optional<bdl_pair<Lyt>> bdl_parter_lower(const bdl_pair<Lyt>& given_bdl, const std::set<bdl_pair<Lyt>>& bdl_pairs)
{
    for (const auto& bdl : bdl_pairs)
    {
        if (sidb_nm_distance<Lyt>(Lyt{}, given_bdl.lower, bdl.upper) < 4 && given_bdl != bdl &&
            given_bdl.lower.y <= bdl.upper.y)
        {
            return std::optional<bdl_pair<Lyt>>(bdl);
        }
    }
    return std::nullopt;
};

enum class WIRE
{
    ALL,
    INPUT_W_INPUT_BDL,
    INPUT_WO_INPUT_BDL,
    OUTPUT
};


template <typename Lyt>
[[nodiscard]] std::vector<bdl_wire<Lyt>> detect_bdl_wires(const Lyt& layout, const WIRE wire_selection = WIRE::ALL)
{
    // collecting all BDL pairs of the layout, including input and output BDL pairs.
    std::set<bdl_pair<Lyt>> bdl_pairs{};
    const auto              normal_bdls = detect_bdl_pairs(layout, sidb_technology::cell_type::NORMAL);

    auto input_bdls  = detect_bdl_pairs(layout, sidb_technology::cell_type::INPUT);
    auto output_bdls = detect_bdl_pairs(layout, sidb_technology::cell_type::OUTPUT);

    for (const auto& bdl : normal_bdls)
    {
        bdl_pairs.insert(bdl);
    }

    // adding the output BDL pairs to the set and changing the BDL type.
    for (auto& bdl : output_bdls)
    {
        bdl.type = sidb_technology::cell_type::NORMAL;
        bdl_pairs.insert(bdl);
    }

    // adding the input BDL pairs to the set and changing the BDL type.
    for (auto& bdl : input_bdls)
    {
        bdl.type = sidb_technology::cell_type::NORMAL;
        bdl_pairs.insert(bdl);
    }

    std::vector<bdl_wire<Lyt>> wires{};
    wires.reserve(input_bdls.size() + output_bdls.size() + normal_bdls.size());

    while (!bdl_pairs.empty())
    {
        bdl_wire<Lyt> chain{};

        bool       neighbor_bdl_found = true;
        auto       front_bdl_pair     = *bdl_pairs.begin();
        const auto start_pair         = front_bdl_pair;

        chain.push_back(front_bdl_pair);
        bdl_pairs.erase(front_bdl_pair);

        while (neighbor_bdl_found)
        {
            const auto lower_partner = bdl_parter_lower<Lyt>(front_bdl_pair, bdl_pairs);
            if (lower_partner.has_value())
            {
                chain.push_back(lower_partner.value());
                bdl_pairs.erase(lower_partner.value());
                front_bdl_pair = lower_partner.value();
            }
            else
            {
                front_bdl_pair = start_pair;
                if (bdl_parter_upper<Lyt>(front_bdl_pair, bdl_pairs).has_value())
                {
                    chain.push_back(lower_partner.value());
                    bdl_pairs.erase(lower_partner.value());
                }
                else
                {
                    neighbor_bdl_found = false;
                    wires.push_back(chain);
                }
            }
        }
    }

    if (wire_selection == WIRE::INPUT_WO_INPUT_BDL)
    {
        std::vector<bdl_wire<Lyt>> in_wires{};
        in_wires.reserve(input_bdls.size());

        for (auto& bdl : input_bdls)
        {
            bdl.type = sidb_technology::cell_type::NORMAL;
            for (auto& chain : wires)
            {
                auto it = std::find(chain.begin(), chain.end(), bdl);
                if (it != chain.end())
                {
                    chain.erase(it);
                    in_wires.push_back(chain);
                    break;
                }
            }
        }
        return in_wires;
    }

    if (wire_selection == WIRE::INPUT_W_INPUT_BDL)
    {
        std::vector<bdl_wire<Lyt>> in_wires{};
        in_wires.reserve(input_bdls.size());

        for (auto& bdl : input_bdls)
        {
            //bdl.type = sidb_technology::cell_type::INPUT;
            for (auto& chain : wires)
            {
                auto it = std::find(chain.begin(), chain.end(), bdl);
                if (it != chain.end())
                {
                    in_wires.push_back(chain);
                    break;
                }
            }
        }
        return in_wires;
    }

    if (wire_selection == WIRE::OUTPUT)
    {
        std::vector<bdl_wire<Lyt>> out_wires{};
        for (auto& bdl : output_bdls)
        {
            bdl.type = sidb_technology::cell_type::NORMAL;
            for (auto& chain : wires)
            {
                auto it = std::find(chain.begin(), chain.end(), bdl);
                if (it != chain.end())
                {
                    out_wires.push_back(chain);
                    break;
                }
            }
        }
        return out_wires;
    }

    std::vector<bdl_wire<Lyt>> in_and_out_wires{};
    for (auto& chain : wires)
    {
        bool input_pair_found = false;
        for (auto& bdl : input_bdls)
        {
            bdl.type = sidb_technology::cell_type::NORMAL;
            if (auto it = std::find(chain.begin(), chain.end(), bdl); it != chain.end())
            {
                input_pair_found = true;
                chain.erase(it);
                in_and_out_wires.push_back(chain);
                break;
            }
        }
        if (!input_pair_found)
        {
            in_and_out_wires.push_back(chain);
        }
    }

    return in_and_out_wires;
}


}

#endif  // FICTION_DETECT_BDL_WIRES_HPP

