cmake_minimum_required(VERSION 3.5)
project(test_function)

# Find pybind11
find_package(pybind11 REQUIRED)

# Find dolfin
find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

# Target libraries

pybind11_add_module(test_function test_function.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE dolfin)