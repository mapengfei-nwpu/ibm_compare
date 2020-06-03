#include <pybind11/stl.h>
#include <vector>
#include <dolfin.h>
using namespace dolfin;

std::array<dolfin::Point, 2> &output(std::array<dolfin::Point, 2> &points) {

    return points;
}

namespace py = pybind11;
PYBIND11_MODULE(test_struct_vector, m) {
    m.doc() = "output what input. ";
    m.def("output", &output, py::return_value_policy::reference);
}
