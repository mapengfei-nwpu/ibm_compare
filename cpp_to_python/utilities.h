
#include <dolfin.h>
#include <dolfin/geometry/SimplexQuadrature.h>
namespace dolfin{

std::vector<std::array<double, 2>> get_global_dof_coordinates(const Function &f);

void get_gauss_rule(const Function &f, std::vector<double> &coordinates, std::vector<double> &values, std::vector<double> &weights);

template <typename T> std::vector<T> my_mpi_gather(std::vector<T> local);



};