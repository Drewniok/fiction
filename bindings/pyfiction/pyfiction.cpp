//
// Created by marcel on 02.06.22.
//

#include "pyfiction/algorithms/network_transformation/fanout_substitution.hpp"
#include "pyfiction/algorithms/network_transformation/network_balancing.hpp"
#include "pyfiction/algorithms/path_finding/a_star.hpp"
#include "pyfiction/algorithms/path_finding/distance.hpp"
#include "pyfiction/algorithms/path_finding/enumerate_all_paths.hpp"
#include "pyfiction/algorithms/path_finding/k_shortest_paths.hpp"
#include "pyfiction/algorithms/physical_design/apply_gate_library.hpp"
#include "pyfiction/algorithms/physical_design/color_routing.hpp"
#include "pyfiction/algorithms/physical_design/exact.hpp"
#include "pyfiction/algorithms/physical_design/orthogonal.hpp"
#include "pyfiction/algorithms/properties/critical_path_length_and_throughput.hpp"
#include "pyfiction/algorithms/simulation/logic_simulation.hpp"
#include "pyfiction/algorithms/verification/design_rule_violations.hpp"
#include "pyfiction/algorithms/verification/equivalence_checking.hpp"
#include "pyfiction/io/read_fqca_layout.hpp"
#include "pyfiction/io/write_dot_layout.hpp"
#include "pyfiction/io/write_fqca_layout.hpp"
#include "pyfiction/io/write_qca_layout.hpp"
#include "pyfiction/io/write_qcc_layout.hpp"
#include "pyfiction/io/write_qll_layout.hpp"
#include "pyfiction/io/write_sqd_layout.hpp"
#include "pyfiction/io/write_svg_layout.hpp"
#include "pyfiction/layouts/cartesian_layout.hpp"
#include "pyfiction/layouts/cell_level_layout.hpp"
#include "pyfiction/layouts/clocked_layout.hpp"
#include "pyfiction/layouts/coordinates.hpp"
#include "pyfiction/layouts/gate_level_layout.hpp"
#include "pyfiction/layouts/hexagonal_layout.hpp"
#include "pyfiction/layouts/obstruction_layout.hpp"
#include "pyfiction/networks/logic_network.hpp"
#include "pyfiction/technology/area.hpp"
#include "pyfiction/utils/routing_utils.hpp"

#include <pybind11/pybind11.h>

PYBIND11_MODULE(pyfiction, m)
{
    // docstring
    m.doc() = "Python bindings fiction, a framework for Design Automation for Field-coupled Nanotechnologies";

    /**
     * Layouts
     */
    pyfiction::coordinates(m);
    pyfiction::cartesian_layout(m);
    pyfiction::hexagonal_layout(m);
    pyfiction::clocked_layouts(m);
    pyfiction::gate_level_layouts(m);
    pyfiction::cell_level_layouts(m);
    pyfiction::obstruction_layouts(m);
    /**
     * Networks
     */
    pyfiction::logic_network(m);
    /**
     * Algorithms: Physical Design
     */
    pyfiction::exact(m);
    pyfiction::orthogonal(m);
    // NOTE: currently not functioning because the Python interpreter can only run as a single instance
    // pyfiction::one_pass_synthesis(m);
    pyfiction::apply_gate_library(m);
    pyfiction::color_routing(m);
    /**
     * Algorithms: Network Transformation
     */
    pyfiction::fanout_substitution(m);
    pyfiction::network_balancing(m);
    /**
     * Algorithms: Path Finding
     */
    pyfiction::distance(m);
    pyfiction::a_star(m);
    pyfiction::yen_k_shortest_paths(m);
    pyfiction::enumerate_all_clocking_paths(m);
    /**
     * Algorithms: Properties
     */
    pyfiction::critical_path_length_and_throughput(m);
    /**
     * Algorithms: Simulation
     */
    pyfiction::logic_simulation(m);
    /*
     * Algorithms: Verification
     */
    pyfiction::equivalence_checking(m);
    pyfiction::design_rule_violations(m);
    /**
     * Technology
     */
    pyfiction::area(m);
    /**
     * Input/Output
     */
    pyfiction::write_dot_layout(m);
    pyfiction::write_qca_layout(m);
    pyfiction::write_svg_layout(m);
    pyfiction::write_sqd_layout(m);
    pyfiction::write_qcc_layout(m);
    pyfiction::write_qll_layout(m);
    pyfiction::write_fqca_layout(m);
    pyfiction::read_fqca_layout(m);
    /**
     * Utils
     */
     pyfiction::routing_utils(m);
}