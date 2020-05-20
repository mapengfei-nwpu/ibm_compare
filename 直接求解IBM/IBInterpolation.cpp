#include <dolfin.h>
#include <dolfin/geometry/SimplexQuadrature.h>
#include <numeric>
#include "IBMesh.h"
using namespace dolfin;

template <typename T>
std::vector<T> my_mpi_gather(std::vector<T> local)
{

	auto mpi_size = dolfin::MPI::size(MPI_COMM_WORLD);

	/// collect local coordinate_dofs on every process
	std::vector<std::vector<T>> mpi_collect(mpi_size);
	dolfin::MPI::all_gather(MPI_COMM_WORLD, local, mpi_collect);
	std::vector<T> global;

	/// unwrap mpi_dof_coordinates.
	for (size_t i = 0; i < mpi_collect.size(); i++)
	{
		for (size_t j = 0; j < mpi_collect[i].size(); j++)
		{
			global.push_back(mpi_collect[i][j]);
		}
	}
	return global;
}

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

std::vector<double> get_global_dofs(const Function &f)
{
	/// some shorcut
	auto mesh = f.function_space()->mesh();
	auto mpi_comm = mesh->mpi_comm();
	auto mpi_size = dolfin::MPI::size(MPI_COMM_WORLD);

	/// get local values.
	std::vector<double> local_values;
	f.vector()->get_local(local_values);

	/// collect local values on every process.
	auto values = my_mpi_gather(local_values);
	return values;
}

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

class DeltaInterplation
{
public:
	/// information about mesh structure.
	IBMesh &um;
	std::vector<double> side_lengths;

	/// construct function.
	DeltaInterplation(IBMesh &uniform_mesh) : um(uniform_mesh)
	{
		side_lengths = um.side_length();
	}

	void fluid_to_solid(Function &fluid, Function &solid)
	{
		/// Calculate global dof coordinates and dofs.
		auto dof_coordinates = get_global_dof_coordinates(solid);
		std::vector<double> values(dof_coordinates.size() * 2);

		fluid_to_solid_raw(fluid, values, dof_coordinates);
		solid.vector()->set_local(values);
	}

	void fluid_to_solid_raw(
		Function &fluid,
		std::vector<double> &solid_values,
		std::vector<std::array<double, 2>> &solid_coordinates)
	{
		/// iterate every dof coordinate.
		auto value_size = fluid.value_size();
		for (size_t i = 0; i < solid_values.size() / value_size; i++)
		{
			Array<double> x(2, solid_coordinates[i].data());
			Array<double> v(2);
			fluid.eval(v, x);
			for (size_t l = 0; l < value_size; l++)
				solid_values[i * value_size + l] = v[l];
		}
	}

	/// Assign the solid displacement with the velocity of fluid.
	std::vector<double> solid_to_fluid(Function &fluid, Function &solid)
	{
		std::vector<double> solid_dof_coordinates;
		std::vector<double> solid_values;
		std::vector<double> weights;
		get_gauss_rule(solid, solid_dof_coordinates, solid_values, weights);

		size_t temp_size = weights.size();
		std::vector<std::array<double, 2>> gauss_points(temp_size);
		std::vector<std::array<double, 2>> gauss_values(temp_size);

		for (size_t i = 0; i < temp_size; i++)
		{
			gauss_points[i][0] = solid_dof_coordinates[2 * i];
			gauss_points[i][1] = solid_dof_coordinates[2 * i + 1];
			gauss_values[i][0] = solid_values[2 * i];
			gauss_values[i][1] = solid_values[2 * i + 1];
		}
		/// std::cout << temp_size << std::endl;
		return solid_to_fluid_raw(fluid, gauss_values, gauss_points, weights);
	}

	std::vector<double> solid_to_fluid_raw(
		Function &fluid,
		std::vector<std::array<double, 2>> &gauss_values,
		std::vector<std::array<double, 2>> &gauss_points,
		std::vector<double> &weights)
	{
		dolfin_assert(solid_values.size() == gauss_points.size());
		dolfin_assert(solid_values.size() == gauss_weights.size());

		std::vector<double> results(fluid.function_space()->dim());

		for (size_t i = 0; i < gauss_values.size(); i++)
		{
			Point point(2, gauss_points[i].data());
			Cell cell = find_cell_contianing_point(point, fluid);
			auto dab = basis_values_gauss(point, cell, fluid);
			integrate_basis(weights[i], dab, gauss_values[i], results);
		}
		return results;
	}

	std::pair<std::vector<double>, std::vector<std::size_t>> basis_values_gauss(
		const Point &point,
		const Cell &cell,
		const Function &f)
	{
		auto element = f.function_space()->element();

		auto value_size = f.value_size();
		auto space_dimension = element->space_dimension();
		std::vector<double> basis_values(value_size * space_dimension);

		ufc::cell ufc_cell;
		cell.get_cell_data(ufc_cell);

		std::vector<double> coordinate_dofs;
		cell.get_coordinate_dofs(coordinate_dofs);

		element->evaluate_basis_all(
			basis_values.data(),
			point.coordinates(),
			coordinate_dofs.data(),
			ufc_cell.orientation);
		
		auto cell_dofmap = f.function_space()->dofmap()->cell_dofs(cell.index());
		std::vector<size_t> cell_dofmap_vector;
		for (size_t i = 0; i < cell_dofmap.size(); i++)
		{
			cell_dofmap_vector.push_back(cell_dofmap[i]);
		}
		
		return std::make_pair(basis_values, cell_dofmap_vector);
	}
	void integrate_basis(
		double weight,
		std::pair<std::vector<double>, std::vector<std::size_t>> &dab,
		const std::array<double, 2> &solid_values,
		std::vector<double> &results)
	{

		size_t space_dimension = 12;
		size_t value_size = 2;

		
		auto basis_value = dab.first;
		auto cell_dofmap = dab.second;

		auto count = cell_dofmap.size();

		/// Every pair in dab consists of a dof index and a basis value.
		/// dab[index].first and dab[index].second
		for (size_t i = 0; i < count; i++)
		{
			for (size_t j = 0; j < value_size; j++)
			{
				results[cell_dofmap[i]] += weight * solid_values[j] * basis_value[2*i+j];
			}
		}
	}

	Cell find_cell_contianing_point(const Point &point, const Function &f)
	{
		/// must be regular grid.
		auto mesh = f.function_space()->mesh();

		auto index_1 = 2 * um.hash(point);
		auto index_2 = 2 * um.hash(point) + 1;

		Cell cell_1(*mesh, index_1);
		Cell cell_2(*mesh, index_2);

		return cell_1.contains(point) ? cell_1 : cell_2;
	}
};
