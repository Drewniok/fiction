//
// Created by Jan Drewniok on 01.12.22.
//

#ifndef FICTION_EXGS_HPP
#define FICTION_EXGS_HPP

#include "fiction/traits.hpp"
#include "fiction/types.hpp"

namespace fiction{

namespace detail{

template <typename Lyt>
class simulation_init
{
  public:
    //
    explicit simulation_init(Lyt &lyt): Lyt(lyt){};


    void to_siqad_coord()
    {
        std::vector<std::vector<unsigned long>> location;
        std::vector<cell<sidb_cell_clk_lyt>>    all_cells{};
        all_cells.reserve(lyt.num_cells());

        lyt.foreach_cell([&all_cells](const auto& c) { all_cells.push_back(c); });

        // transform coordinates in new coordinates, inspired by SIQAD
        for (const auto& c : all_cells)
        {
            auto          x = static_cast<unsigned long>(c.x);
            auto          y = static_cast<unsigned long>(((c.y) - ((c.y) % 2)) * 0.5);
            unsigned long Z = ((c.y) % 2);
            location.push_back({x, y, Z});
        }
    };

  private:
    Lyt &lyt;
};

}

}

#endif  // FICTION_EXGS_HPP
