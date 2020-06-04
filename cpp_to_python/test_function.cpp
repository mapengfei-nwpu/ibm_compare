#include<dolfin.h>
#include<memory>
#include<pybind11/pybind11.h>

using namespace dolfin;
int size(std::shared_ptr<Function> function){
    return function->function_space()->dim();
}

namespace py = pybind11;
PYBIND11_MODULE(test_function, m) {
    m.doc() = "pybind11 test dolfin function"; // optional module docstring
    m.def("size", &size, "Return the size of a function.");
}