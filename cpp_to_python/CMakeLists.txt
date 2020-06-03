cmake_minimum_required(VERSION 3.5)
project(test_function)

# Find pybind11
find_package(pybind11 REQUIRED)

# Find dolfin
find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

# Add pybin11 module
pybind11_add_module(test_function test_function.cpp)
pybind11_add_module(test_vector test_vector.cpp)
pybind11_add_module(test_struct_vector test_struct_vector.cpp)
pybind11_add_module(IBMesh IBMesh.cpp)


# Target the library
target_link_libraries(test_function PRIVATE dolfin)
target_link_libraries(test_struct_vector PRIVATE dolfin)
target_link_libraries(IBMesh PRIVATE dolfin)
