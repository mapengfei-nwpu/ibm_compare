#include<pybind11/stl.h>
#include<vector>

namespace py = pybind11;

std::vector<float> vec_add(std::vector<float>& a, std::vector<float>& b) {

    std::vector<float> out;
    assert(a.size() == b.size());
    for (int i = 0; i < a.size(); i++)
    {
        out.push_back(a[i] + b[i]);
    }
    return out;

}

PYBIND11_MODULE(test_vector, m) {
    m.doc() = "This is a simple demo using C++ STL";
    m.def("vec_add", &vec_add, py::return_value_policy::reference);
}
