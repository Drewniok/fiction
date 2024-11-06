//
// Created by Jan Drewniok 18.10.2024
//

#ifndef PYFICTION_DEFECT_OPERATIONAL_DOMAIN_HPP
#define PYFICTION_DEFECT_OPERATIONAL_DOMAIN_HPP

#include "pyfiction/documentation.hpp"
#include "pyfiction/types.hpp"

#include <fiction/algorithms/simulation/sidb/defect_operational_domain.hpp>

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace pyfiction
{

namespace detail
{

template <typename Lyt>
void defect_operational_domain_detail(pybind11::module& m)
{
    using namespace pybind11::literals;

    m.def("defect_operational_domain_grid_search", &fiction::defect_operational_domain_grid_search<Lyt, py_tt>, "lyt"_a,
          "spec"_a, "step_size"_a, "params"_a = fiction::defect_operational_domain_params{}, "stats"_a = nullptr,
          DOC(fiction_defect_operational_domain_grid_search));

    m.def("defect_operational_domain_random_sampling", &fiction::defect_operational_domain_random_sampling<Lyt, py_tt>,
          "lyt"_a, "spec"_a, "samples"_a, "params"_a = fiction::defect_operational_domain_params{}, "stats"_a = nullptr,
          DOC(fiction_defect_operational_domain_random_sampling));

    m.def("defect_operational_domain_quicktrace", &fiction::defect_operational_domain_quicktrace<Lyt, py_tt>, "lyt"_a,
          "spec"_a, "samples"_a, "params"_a = fiction::defect_operational_domain_params{}, "stats"_a = nullptr,
          DOC(fiction_defect_operational_domain_quicktrace));
}

}  // namespace detail

inline void defect_operational_domain(pybind11::module& m)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<fiction::defect_operational_domain<py_sidb_100_lattice_cube>>(m, "defect_operational_domain_100_cube")
        .def(py::init<>());

    py::class_<fiction::defect_operational_domain<py_sidb_111_lattice_cube>>(m, "defect_operational_domain_111_cube")
        .def(py::init<>());

    py::class_<fiction::defect_operational_domain_params>(m, "defect_operational_domain_params",
                                                          DOC(fiction_defect_operational_domain_params))
        .def(py::init<>())
        .def_readwrite("defect_influence_params", &fiction::defect_operational_domain_params::defect_influence_params,
                       DOC(fiction_defect_operational_domain_params_defect_influence_params))
        .def_readwrite("operational_params", &fiction::defect_operational_domain_params::operational_params,
                       DOC(fiction_defect_operational_domain_params_operational_params));

    py::class_<fiction::defect_operational_domain_stats>(m, "defect_operational_domain_stats",
                                                         DOC(fiction_defect_operational_domain_stats))
        .def(py::init<>())
        .def_readonly("time_total", &fiction::defect_operational_domain_stats::time_total,
                      DOC(fiction_defect_operational_domain_stats_duration))
        .def_readonly("num_simulator_invocations", &fiction::defect_operational_domain_stats::num_simulator_invocations,
                      DOC(fiction_defect_operational_domain_stats_num_simulator_invocations))
        .def_readonly("num_evaluated_parameter_combinations",
                      &fiction::defect_operational_domain_stats::num_evaluated_defect_positions,
                      DOC(fiction_defect_operational_domain_stats_num_evaluated_defect_positions))
        .def_readonly("num_operational_defect_positions",
                      &fiction::defect_operational_domain_stats::num_operational_defect_positions,
                      DOC(fiction_defect_operational_domain_stats_num_operational_defect_positions))
        .def_readonly("num_non_operational_defect_positions",
                      &fiction::defect_operational_domain_stats::num_non_operational_defect_positions,
                      DOC(fiction_defect_operational_domain_stats_num_non_operational_defect_positions));

    // NOTE be careful with the order of the following calls! Python will resolve the first matching overload!
    detail::defect_operational_domain_detail<py_sidb_100_lattice_cube>(m);
    detail::defect_operational_domain_detail<py_sidb_111_lattice_cube>(m);
}

}  // namespace pyfiction

#endif  // PYFICTION_DEFECT_OPERATIONAL_DOMAIN_HPP
