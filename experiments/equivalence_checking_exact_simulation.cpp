//
// Created by Jan Drewniok on 03.03.24.
//

#include <fiction/algorithms/simulation/sidb/check_equivalence_simulation_result.hpp>
#include <fiction/algorithms/simulation/sidb/exhaustive_ground_state_simulation.hpp>
#include <fiction/algorithms/simulation/sidb/quickexact.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp>
#include <fiction/layouts/coordinates.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/types.hpp>
#include <fiction/utils/layout_utils.hpp>
#include <fiction/utils/math_utils.hpp>

#include <fmt/format.h>

#include <cstdint>

using namespace fiction;

int main()  // NOLINT
{
    const auto all_cells_in_region = all_coordinates_in_spanned_area<offset::ucoord_t>({0, 0}, {10, 10});

    const auto all_distributions =
        determine_all_combinations_of_distributing_k_entities_on_n_positions(4, all_cells_in_region.size());

    const auto params = sidb_simulation_parameters{3, -0.32};

    uint64_t non_equivalence_counter = 0;

    for (const auto& distribution : all_distributions)
    {
        sidb_cell_clk_lyt lyt{};
        for (const auto idx : distribution)
        {
            lyt.assign_cell_type(all_cells_in_region[idx], sidb_technology::cell_type::NORMAL);
        }
        auto result_exgs = exhaustive_ground_state_simulation(lyt, params);
        auto result_quickexact =
            quickexact(lyt, quickexact_params<sidb_cell_clk_lyt>{
                                params, quickexact_params<sidb_cell_clk_lyt>::automatic_base_number_detection::OFF});
        if (!check_equivalence_simulation_result(result_exgs, result_quickexact))
        {
            non_equivalence_counter++;
        }
    }
    std::cout << fmt::format("non equivalent layouts = {}", non_equivalence_counter) << std::endl;
}