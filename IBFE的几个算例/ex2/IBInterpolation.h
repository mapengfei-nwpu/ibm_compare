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

	/// Assign the solid displacement with the velocity of fluid.
	void fluid_to_solid(Function &fluid, Function &solid)
	{
		/// calculate global dof coordinates and dofs.
		auto dof_coordinates = get_global_dof_coordinates(solid);
		std::vector<double> values(dof_coordinates.size() * 2);
		fluid_to_solid_raw(fluid, values, dof_coordinates);

		/// collect all solid_values vectors on every processes and
		/// assign them to a function.
		auto size = solid.vector()->local_size();
		auto offset = solid.vector()->local_range().first;
		std::vector<double> local_values(size);
		for (size_t i = 0; i < size; ++i)
			local_values[i] = values[i + offset];
		solid.vector()->set_local(local_values);

		/// Finalize assembly of tensor.
		/// this step is quite important.
		solid.vector()->apply("insert");
	}

	void fluid_to_solid_raw(Function &fluid, std::vector<double> &solid_values,
							std::vector<std::array<double, 2>> &solid_coordinates)
	{
				/// the meshes of v and um should be the same.
		/// TODO : compare two meshes
		dolfin_assert(fluid.value_size() == solid_values.size() / solid_coordinates.size());

		/// smart shortcut
		auto value_size   = fluid.value_size();

		/// iterate every solid_values coordinate.
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
	void solid_to_fluid(Function &fluid, Function &solid)
	{
		/// calculate global dof coordinates and dofs of solid.
		/// auto solid_dof_coordinates = get_global_dof_coordinates(solid);
		/// auto solid_values = get_global_dofs(solid);
		std::vector<double> solid_dof_coordinates;
		std::vector<double> solid_values;
		std::vector<double> weights;
		get_gauss_rule(solid, solid_dof_coordinates, solid_values, weights);

		/// interpolates solid values into fluid mesh.
		/// the returned fluid_values is the dofs of fluid function.

		auto fluid_values = solid_to_fluid_raw(fluid, solid_values, solid_dof_coordinates, weights);

		/// assemble fluid_values into a function.
		auto local_size = fluid.vector()->local_size();
		auto offset = fluid.vector()->local_range().first;
		std::vector<double> local_values(local_size);
		for (size_t i = 0; i < local_size; ++i)
		{
			local_values[i] = fluid_values[i + offset];
		}

		fluid.vector()->set_local(local_values);
		fluid.vector()->apply("insert");

	}

	std::vector<double> solid_to_fluid_raw(
		Function &fluid,
		std::vector<double> &solid_values,
		std::vector<double> &solid_coordinates,
		std::vector<double> &weights)
	{
		/// smart shortcut
		auto rank = dolfin::MPI::rank(fluid.function_space()->mesh()->mpi_comm());
		auto mesh = fluid.function_space()->mesh();		// pointer to a mesh
		auto dofmap = fluid.function_space()->dofmap(); // pointer to a dofmap
		auto hmax = mesh->hmax();

		/// get the element of function space
		auto element = fluid.function_space()->element();
		auto value_size = fluid.value_size();
		auto global_fluid_size = fluid.function_space()->dim();

		/// Get local to global dofmap
		std::vector<size_t> local_to_global;
		dofmap->tabulate_local_to_global_dofs(local_to_global);

		/// initial local fluid values.
		std::vector<double> local_fluid_values(global_fluid_size);

		/// iterate every solid_values coordinate.
		for (size_t i = 0; i < solid_values.size() / value_size; i++)
		{

			/// get indices of adjacent cells on fluid mesh.
			Point solid_point(solid_coordinates[2 * i], solid_coordinates[2 * i + 1]);
			auto adjacents = um.get_adjacents(solid_point);

			/// iterate adjacent cells and collect element nodes in these cells.
			/// it has nothing to do with cell type.
			std::map<size_t, double> indices_to_delta;
			for (size_t j = 0; j < adjacents.size(); j++)
			{
				/// step 1 : get coordinates of cell dofs
				Cell cell(*mesh, adjacents[j]);
				std::vector<double> coordinate_dofs;
				cell.get_coordinate_dofs(coordinate_dofs);
				boost::multi_array<double, 2> coordinates;
				element->tabulate_dof_coordinates(coordinates, coordinate_dofs, cell);

				/// step 2 : get the dof map
				auto cell_dofmap = dofmap->cell_dofs(cell.index());

				/// iterate node coordinates of the cell.
				for (size_t k = 0; k < cell_dofmap.size() / value_size; k++)
				{
					Point cell_point(coordinates[k][0], coordinates[k][1]);
					double param = delta(solid_point, cell_point, hmax);
					if (cell_dofmap[k] < fluid.vector()->local_size() && param > 0.0)
					{
						indices_to_delta[cell_dofmap[k]] = param;
					}
				}
			}

			/// delta distribution.
			for (auto it = indices_to_delta.begin(); it != indices_to_delta.end(); it++)
			{
				for (size_t l = 0; l < value_size; l++)
				{
					local_fluid_values[local_to_global[it->first + l]] += solid_values[i * value_size + l] * it->second * weights[i];
				}
			}
		}

		//////////////////  TODO : this part can use MPI_reduce directly  //////////////////////
		std::vector<double> fluid_values(global_fluid_size);
		std::vector<std::vector<double>> mpi_collect(dolfin::MPI::size(mesh->mpi_comm()));
		dolfin::MPI::all_gather(mesh->mpi_comm(), local_fluid_values, mpi_collect);
		for (size_t i = 0; i < fluid_values.size(); i++)
		{
			for (size_t j = 0; j < mpi_collect.size(); j++)
			{
				fluid_values[i] += mpi_collect[j][i];
			}
		}
		return fluid_values;
	}

	/////////////////////////////////////////////
	//  thses methods must not be modified!!  ///
	/////////////////////////////////////////////
	double phi(double r)
	{
		r = fabs(r);
		if (r > 2)
			return 0;
		else
			return 0.25 * (1 + cos(FENICS_PI * r * 0.5));
	}
	double delta(Point p0, Point p1, double h = 0.0625)
	{
		double ret = 1.0;
   		h = side_lengths[0]*0.5;
		for (unsigned i = 0; i < 2; ++i)
		{
			double dx = p0.coordinates()[i] - p1.coordinates()[i];
			ret *= 1. / h * phi(dx / h);
		}
		return ret;
	}
	/////////////////////////////////////////////
	//  thses methods must not be modified!!  ///
	/////////////////////////////////////////////
};
