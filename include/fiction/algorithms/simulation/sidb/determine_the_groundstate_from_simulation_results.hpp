//
// Created by Jan Drewniok on 12.02.24.
//

#ifndef FICTION_DETERMINE_THE_GROUNDSTATE_FROM_CDS_HPP
#define FICTION_DETERMINE_THE_GROUNDSTATE_FROM_CDS_HPP

#include <fiction/algorithms/simulation/sidb/sidb_simulation_result.hpp>
#include <fiction/technology/charge_distribution_surface.hpp>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <set>
#include <vector>

namespace fiction
{

/**
 * This function calculates the ground state charge distributions from the provided simulation results.
 * The ground state charge distributions are those with energy closest to the minimum energy found in the simulation
 * results.
 *
 * @tparam Lyt The layout type used in the simulation results.
 * @param simulation_results The simulation results containing charge distributions.
 * @return A vector of charge distributions with the minimal energy.
 *
 * @note The size of the return vector is only greater than one if degenerate states exist.
 */
template <typename Lyt>
[[nodiscard]] std::vector<charge_distribution_surface<Lyt>>
determine_the_groundstate_from_simulation_results(const sidb_simulation_result<Lyt>& simulation_results) noexcept
{
    std::vector<charge_distribution_surface<Lyt>> groundstate_charge_distributions{};
    std::set<uint64_t>                            charge_indices{};

    // Find all unique charge indices. This is done because simulation results can have multiple identical charge
    // distributions.
    for (const auto& cds : simulation_results.charge_distributions)
    {
        charge_indices.insert(cds.get_charge_index_and_base().first);
    }

    // Find the minimum energy
    double min_energy = std::numeric_limits<double>::infinity();
    for (const auto& cds : simulation_results.charge_distributions)
    {
        min_energy = std::min(min_energy, cds.get_system_energy());
    }

    for (const auto charge_index : charge_indices)
    {
        auto cds_it = std::find_if(
            simulation_results.charge_distributions.begin(), simulation_results.charge_distributions.end(),
            [&](const auto& cds)
            {
                return cds.get_charge_index_and_base().first == charge_index &&
                       std::abs(cds.get_system_energy() - min_energy) < std::numeric_limits<double>::epsilon();
            });

        if (cds_it != simulation_results.charge_distributions.end())
        {
            groundstate_charge_distributions.push_back(*cds_it);
        }
    }

    return groundstate_charge_distributions;
}

}  // namespace fiction

#endif  // FICTION_DETERMINE_THE_GROUNDSTATE_FROM_CDS_HPP