#include <dolfin.h>
#include <dolfin/geometry/SimplexQuadrature.h>
using namespace dolfin;

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

/// f 必须属于拉格朗日空间，这样，它的每个点的值就是对应的自由度。
/// dofs are placed on the coordinates of dofs
/// therefore, dofs are assigned on every dof_coordinate.
std::vector<std::array<double, 2>> get_global_dof_coordinates(const Function &f)
{
	/// some shorcut
	auto mesh = f.function_space()->mesh();
	auto mpi_comm = mesh->mpi_comm();
	auto mpi_size = dolfin::MPI::size(MPI_COMM_WORLD);

	/// get local coordinate_dofs
	auto local_dof_coordinates = f.function_space()->tabulate_dof_coordinates();

	/// collect local coordinate_dofs on every process
	auto dof_coordinates_long = my_mpi_gather(local_dof_coordinates);

	/// unwrap it.
	std::vector<std::array<double, 2>> dof_coordinates(dof_coordinates_long.size() / 4);
	for (size_t i = 0; i < dof_coordinates_long.size(); i += 4)
	{
		dof_coordinates[i / 4][0] = dof_coordinates_long[i];
		dof_coordinates[i / 4][1] = dof_coordinates_long[i + 1];
	}

	/// whatch the type of this function return.
	/// it could be changed to "std::vector<std::array<double,3>>" if necessary.
	return dof_coordinates;
}

/// 计算出一个网格中的所有高斯点和对应的权重。
/// 并求出每个高斯点上函数f的值。
void get_gauss_rule(
	const Function &f,
	std::vector<double> &coordinates,
	std::vector<double> &values,
	std::vector<double> &weights)
{
	// Construct Gauss quadrature rules
	SimplexQuadrature gq(2, 9);
	auto mesh = f.function_space()->mesh();

	/// iterate cell of local mesh.
	for (CellIterator cell(*mesh); !cell.end(); ++cell)
	{
		// Create ufc_cell associated with dolfin cell.
		ufc::cell ufc_cell;
		cell->get_cell_data(ufc_cell);

		// Compute quadrature rule for the cell.
		auto qr = gq.compute_quadrature_rule(*cell);
		dolfin_assert(qr.second.size() == qr.first.size() / 2);
		for (size_t i = 0; i < qr.second.size(); i++)
		{
			Array<double> v(2);
			Array<double> x(2);

			/// Call evaluate function
			x[0] = qr.first[2 * i];
			x[1] = qr.first[2 * i + 1];
			f.eval(v, x, *cell, ufc_cell);

			/// push back what we get.
			coordinates.push_back(x[0]);
			coordinates.push_back(x[1]);
			weights.push_back(qr.second[i]);
			values.push_back(v[0]);
			values.push_back(v[1]);
		}
	}
	values = my_mpi_gather(values);
	weights = my_mpi_gather(weights);
	coordinates = my_mpi_gather(coordinates);
}
