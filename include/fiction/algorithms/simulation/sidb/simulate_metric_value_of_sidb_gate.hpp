//
// Created by Jan Drewniok on 20.11.23.
//

#ifndef FICTION_SIMULATE_METRIC_VALUE_OF_SIDB_GATE_HPP
#define FICTION_SIMULATE_METRIC_VALUE_OF_SIDB_GATE_HPP

#include "fiction/algorithms/simulation/sidb/assess_physical_population_stability.hpp"
#include "fiction/algorithms/simulation/sidb/critical_temperature.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/algorithms/simulation/sidb/maximum_defect_influence_position_and_distance.hpp"
#include "fiction/algorithms/simulation/sidb/operational_domain.hpp"

#include <vector>

namespace fiction
{

/**
 * Various metrics that can be used to assess the design of a gate. Each metric provides a specific aspect of
 * evaluation, such as critical temperature, population stability, maximum defect influence
 * distance, operational domain, or none (no specific metric applied). The metrics serve as criteria to evaluate the
 * performance and characteristics of the gate, helping to analyze and optimize gate implementations.
 */
enum class gate_metric
{
    /**
     * The critical temperature of the gate.
     */
    CRITICAL_TEMPERATURE,
    /**
     * The population stability of the gate.
     */
    POPULATION_STABILITY,
    /**
     * The defect avoidance distance of the gate.
     */
    MAXIMUM_DEFECT_INFLUENCE_DISTANCE,
    /**
     * The operational domain of the gate.
     */
    OPERATIONAL_DOMAIN,
    /**
     * No metric is applied.
     */
    NONE,
};

struct simulate_metric_value_of_sidb_gate_params
{
    /**
     * All Parameters for physical SiDB simulations.
     */
    sidb_simulation_parameters phys_params{2};
    /**
     * The simulation engine to be used for the operational domain computation.
     */
    sidb_simulation_engine      sim_engine{sidb_simulation_engine::QUICKEXACT};
    critical_temperature_params critical_temp_params{quicksim_params{phys_params, 0, 0.0},
                                                     critical_temperature_params::simulation_engine::EXACT, 0.99, 350};
    /**
     * Parameters used to simulate the population stability.
     */
    assess_physical_population_stability_params population_stability_params{phys_params};
    /**
     * Parameters used to simulate the maximum defect influence distance.
     */
    maximum_defect_influence_distance_params influence_defect_distance_param{
        sidb_defect{sidb_defect_type::UNKNOWN, -1, 10.6, 5.9}, phys_params};
    /**
     * Parameters used to simulate the operational domain area.
     */
    operational_domain_params op_domain_params{
        phys_params, sim_engine, operational_domain::sweep_parameter::EPSILON_R, 1.0,
        10.0,        0.2,        operational_domain::sweep_parameter::LAMBDA_TF, 1.0,
        10.0,        0.2};
};

namespace detail
{

template <typename Lyt, typename TT>
class simulate_metric_value_of_sidb_gate_impl
{
  public:
    simulate_metric_value_of_sidb_gate_impl(const Lyt& lyt, const std::vector<TT>& tt, const gate_metric& met,
                                            const simulate_metric_value_of_sidb_gate_params& params) :
            layout{lyt},
            metric{met},
            truth_table{tt},
            parameter{params},
            bii(bdl_input_iterator<Lyt>{layout, params.op_domain_params.bdl_params})
    {}

    [[nodiscard]] double metric_value_of_sidb_gate() noexcept
    {
        switch (metric)
        {
            case gate_metric::CRITICAL_TEMPERATURE:
            {
                return critical_temperature_gate_based<Lyt>(layout, truth_table, parameter.critical_temp_params);
            }
            case gate_metric::POPULATION_STABILITY:
            {
                return assess_population_stability_of_sidb_gate();
            }
            case gate_metric::OPERATIONAL_DOMAIN:
            {
                operational_domain_stats op_domain_stats{};
                operational_domain_flood_fill(layout, truth_table, 0, parameter.op_domain_params,
                                              operational_domain::parameter_point(parameter.phys_params.epsilon_r,
                                                                                  parameter.phys_params.lambda_tf),
                                              &op_domain_stats);
                return op_domain_stats.percentual_operational_area;
            }
            case gate_metric::MAXIMUM_DEFECT_INFLUENCE_DISTANCE:
            {
                return maximum_defect_influence_position_and_distance_of_sidb_gate();
            }
            case gate_metric::NONE: return 0.0;
        }
        return 0.0;
    }

  private:
    const Lyt&                                      layout{};
    const gate_metric&                              metric{};
    const std::vector<TT>&                          truth_table;
    const simulate_metric_value_of_sidb_gate_params parameter{};
    /**
     * Iterator that iterates over all possible input states.
     */
    bdl_input_iterator<Lyt> bii;

    [[nodiscard]] double assess_population_stability_of_sidb_gate() noexcept
    {
        auto minimal_pop_stability_for_all_inputs = std::numeric_limits<double>::infinity();
        // number of different input combinations
        for (auto i = 0u; i < truth_table.front().num_bits(); ++i, ++bii)
        {
            const auto pop_stability =
                assess_physical_population_stability<Lyt>(layout, parameter.population_stability_params);
            if (!pop_stability.empty())
            {
                const auto stability_for_given_input = pop_stability.front().distance_corresponding_to_potential;
                if (pop_stability.front().distance_corresponding_to_potential < minimal_pop_stability_for_all_inputs)
                {
                    minimal_pop_stability_for_all_inputs = stability_for_given_input;
                }
            }
        }
        return minimal_pop_stability_for_all_inputs;
    }

    [[nodiscard]] double maximum_defect_influence_position_and_distance_of_sidb_gate() noexcept
    {
        double maximum_defect_influence_distance = 0.0;
        // number of different input combinations
        for (auto i = 0u; i < truth_table.front().num_bits(); ++i, ++bii)
        {
            const auto influence_distance =
                maximum_defect_influence_position_and_distance(layout, parameter.influence_defect_distance_param)
                    .second;
            if (influence_distance > maximum_defect_influence_distance)
            {
                maximum_defect_influence_distance = influence_distance;
            }
        }
        return maximum_defect_influence_distance;
    }
};
}  // namespace detail

template <typename Lyt, typename TT>
[[nodiscard]] double simulate_metric_value_of_sidb_gate(const Lyt& lyt, const std::vector<TT>& spec,
                                                        const gate_metric&                               metric,
                                                        const simulate_metric_value_of_sidb_gate_params& params)
{
    detail::simulate_metric_value_of_sidb_gate_impl<Lyt, TT> p{lyt, spec, metric, params};
    return p.metric_value_of_sidb_gate();
}
}  // namespace fiction

#endif  // FICTION_SIMULATE_METRIC_VALUE_OF_SIDB_GATE_HPP
