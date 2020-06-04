#include <dolfin.h>
#include <dolfin/geometry/SimplexQuadrature.h>
#include <numeric>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include "IBMesh.h"
#include "utilities.h"
using namespace dolfin;


class IBInterpolation
{
public:
	/// information about mesh structure.
	std::shared_ptr<IBMesh> um;

	/// construct function.
	IBInterpolation(std::shared_ptr<IBMesh> um) : um(um)
	{
		///
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
		return solid_to_fluid_raw(fluid, gauss_values, gauss_points, weights);
	}

	std::vector<double> solid_to_fluid_raw(
		Function &fluid,
		std::vector<std::array<double, 2>> &gauss_values,
		std::vector<std::array<double, 2>> &gauss_points,
		std::vector<double> &weights)
	{
		/// dolfin_assert(solid_values.size() == gauss_points.size());
		/// dolfin_assert(solid_values.size() == gauss_weights.size());

		std::vector<double> results(fluid.function_space()->dim());

		for (size_t i = 0; i < gauss_values.size(); i++)
		{
			Point point(2, gauss_points[i].data());
			Cell cell;
            /// 单线程运行是肯定能找到cell的
            if (um->find_cell(point, cell)){
			    auto dab = basis_values_gauss(point, cell, fluid);
			    integrate_basis(weights[i], dab, gauss_values[i], results);
            }
		}
		return results;
	}

    /// 在找到了高斯点所在的单元的情况下，计算出各个基函数在高斯点处的值
    /// 并且求出各个基函数的全局索引。
    /// TODO : 对于Vector，目前所知，需要调用A.add_local(data, dofs)和
    ///        A.apply("add")来插入数据，这样的话不需要拷贝数据。
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
};

namespace py = pybind11;
PYBIND11_MODULE(IBInterpolation, m)
{
    py::class_<IBInterpolation>(m, "IBInterpolation")
        .def(py::init<std::shared_ptr<IBMesh>>())
		.def("solid_to_fluid", &IBInterpolation::solid_to_fluid)
		.def("fluid_to_solid", &IBInterpolation::fluid_to_solid)
		;
}