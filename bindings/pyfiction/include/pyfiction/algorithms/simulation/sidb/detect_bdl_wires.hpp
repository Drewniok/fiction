//
// Created by Jan Drewniok on 17.07.24.
//

#ifndef PYFICTION_DETECT_BDL_WIRES_HPP
#define PYFICTION_DETECT_BDL_WIRES_HPP

#include "pyfiction/documentation.hpp"
#include "pyfiction/types.hpp"

#include <fiction/algorithms/simulation/sidb/detect_bdl_wires.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace pyfiction
{

namespace detail
{

template <typename Lyt>
void detect_bdl_wires(pybind11::module& m)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def("detect_bdl_wires", &fiction::detect_bdl_wires<Lyt>, "lyt"_a, "params"_a = fiction::detect_bdl_wires_params{},
          "wire_selection"_a = fiction::bdl_wire_selection::ALL);
}

}  // namespace detail

inline void detect_bdl_wires(pybind11::module& m)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::enum_<fiction::bdl_wire_direction>(m, "bdl_wire_direction")
        .value("NORTH_SOUTH", fiction::bdl_wire_direction::NORTH_SOUTH)
        .value("SOUTH_NORTH", fiction::bdl_wire_direction::SOUTH_NORTH)
        .value("NO_DIRECTION", fiction::bdl_wire_direction::NO_DIRECTION);

    py::enum_<fiction::bdl_wire_selection>(m, "bdl_wire_selection")
        .value("ALL", fiction::bdl_wire_selection::ALL)
        .value("INPUT", fiction::bdl_wire_selection::INPUT)
        .value("OUTPUT", fiction::bdl_wire_selection::OUTPUT);

    py::class_<fiction::detect_bdl_wires_params>(m, "detect_bdl_wires_params")
        .def(py::init<>())
        .def_readwrite("threshold_bdl_interdistance", &fiction::detect_bdl_wires_params::threshold_bdl_interdistance)
        .def_readwrite("params_bdl_pairs", &fiction::detect_bdl_wires_params::params_bdl_pairs);

    // NOTE be careful with the order of the following calls! Python will resolve the first matching overload!

    detail::detect_bdl_wires<py_sidb_100_lattice>(m);
    detail::detect_bdl_wires<py_sidb_111_lattice>(m);
}

}  // namespace pyfiction

#endif  // PYFICTION_DETECT_BDL_WIRES_HPP