# Set minimum cmake version
cmake_minimum_required(VERSION 3.5)

# Set project name
set(PROJECT_NAME pengfeiIB_project)
project(${PROJECT_NAME})

# Set CMake behavior
cmake_policy(SET CMP0004 NEW)

# Find dolfin
find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

# Print all warning informations
# Compile for valgrind
add_definitions("-Wall -g")

# Executable
add_executable(${PROJECT_NAME} main.cpp)

# Target libraries
target_link_libraries(${PROJECT_NAME} dolfin)

# Test targets
set(test_parameters -np 3 ${MPIEXEC_PARAMS} "./${PROJECT_NAME}")
add_test(NAME ${PROJECT_NAME}_mpi COMMAND "mpirun" ${test_parameters})
add_test(NAME ${PROJECT_NAME}_serial COMMAND ${PROJECT_NAME})
