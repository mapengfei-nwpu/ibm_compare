#ifndef _UTILITIES_H_
#define _UTILITIES_H_


#include <dolfin.h>
namespace dolfin{

std::vector<std::array<double, 2>> get_global_dof_coordinates(const Function &f);

void get_gauss_rule(const Function &f, std::vector<double> &coordinates, std::vector<double> &values, std::vector<double> &weights);

/// 把每个处理器上的vector收集起来，放在一个vector中。
template <typename T>
std::vector<T> my_mpi_gather(std::vector<T> local)
{
	auto mpi_size = dolfin::MPI::size(MPI_COMM_WORLD);

	/// Collect vector on every processor.
	std::vector<std::vector<T>> mpi_collect(mpi_size);
	dolfin::MPI::all_gather(MPI_COMM_WORLD, local, mpi_collect);
	std::vector<T> global;

	/// Unwrap vector<vector> to vector.
	for (size_t i = 0; i < mpi_collect.size(); i++)
	{
		for (size_t j = 0; j < mpi_collect[i].size(); j++)
		{
			global.push_back(mpi_collect[i][j]);
		}
	}
	return global;
}

};

#endif