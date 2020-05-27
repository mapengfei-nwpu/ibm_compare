#include <dolfin.h>
#include <dolfin/geometry/SimplexQuadrature.h>
#include <numeric>
using namespace dolfin;
/// calculate derivatives at gauss point.
void get_gauss_rule(
    const Function &displace,
    std::vector<std::vector<double>> &coordinates,
    std::vector<std::vector<double>> &values,
    std::vector<double> &weights)
{
    // Construct Gauss quadrature rules
    // dimension 2 and order 9
    SimplexQuadrature gq(2, 9);
    auto mesh = displace.function_space()->mesh();
    auto dofmap = displace.function_space()->dofmap();
    auto element = displace.function_space()->element();
    auto value_size = displace.value_size();
    auto derivative_value_size = value_size * value_size;
    auto space_dimension = element->space_dimension();

    for (CellIterator cell(*mesh); !cell.end(); ++cell)
    {
        ufc::cell ufc_cell;
        cell->get_cell_data(ufc_cell);

        std::vector<double> coordinate_dofs;
        cell->get_coordinate_dofs(coordinate_dofs);

        std::vector<double> basis_derivative_values(value_size * space_dimension * 2);

        /// Compute quadrature rule for the cell.
        /// qr.second and qr.first are the coordinate and weight of gauss point respectively.
        auto qr = gq.compute_quadrature_rule(*cell);
        for (size_t i = 0; i < qr.second.size(); i++)
        {
            std::vector<double> point({qr.first[2 * i], qr.first[2 * i + 1]});
            double weight = qr.second[i];
            element->evaluate_basis_derivatives_all(
                1,
                basis_derivative_values.data(),
                point.data(),
                coordinate_dofs.data(),
                ufc_cell.orientation);
            auto cell_dofmap = dofmap->cell_dofs(cell->index());
            std::vector<double> derivative_value(derivative_value_size);
            for (size_t i = 0; i < cell_dofmap.size(); i++)
            {
                for (size_t j = 0; j < derivative_value_size; j++)
                {
                    derivative_value[j] +=
                        (*displace.vector())[cell_dofmap[i]] *
                        basis_derivative_values[derivative_value_size * i + j];
                }
            }
            values.push_back(derivative_value);
            coordinates.push_back(point);
            weights.push_back(weight);
        }
    }
}

/// 一个函数离散成了三个变量：高斯点，函数值，权重。
/// 现在要对它求积，将它与另一个函数空间的基函数点积，然后积分。

/// 这个函数传入：一个高斯点，高斯点上的函数值，
///        返回：这个单元内非零的基函数在这个点上的值，对应基函数的自由度索引。

void calculate_basis_derivative_values(
    const FunctionSpace &function_space,
    const Cell &cell,
    const std::vector<double> &point,
    std::vector<size_t> &cell_dofmap,
    std::vector<double> &cell_basis_derivatives)
{
    /// std::vector<size_t> &cell_dofmap 和 std::vector<double> &cell_basis_derivatives
    /// 需要在调用这个函数之前分配好空间。
    auto element = function_space.element();

    ufc::cell ufc_cell;
    cell.get_cell_data(ufc_cell);

    std::vector<double> coordinate_dofs;
    cell.get_coordinate_dofs(coordinate_dofs);

    /// store values derivative basis.
    element->evaluate_basis_derivatives_all(
        1,
        cell_basis_derivatives.data(),
        point.data(),
        coordinate_dofs.data(),
        ufc_cell.orientation);

    auto cell_dofmap_eigen = function_space.dofmap()->cell_dofs(cell.index());

    std::vector<size_t> cell_dofmap;

    /// store cell dof_map.
    for (size_t i = 0; i < cell_dofmap_eigen.size(); i++)
    {
        cell_dofmap[i] = cell_dofmap_eigen[i];
    }
}

std::vector<double> assemble(
    const Function &displace,
    const FunctionSpace &function_space)
{

    auto space_dimension = function_space.element->space_dimension();
    std::vector<double> results(function_space.dim());

    std::vector<std::vector<double>> coordinates;
    std::vector<std::vector<double>> values;
    std::vector<double> weights;

    get_gauss_rule(displace, coordinates, values, weights);

    for (size_t i = 0; i < coordinates.size(); i++)
    {
        auto point = coordinates[i];
        auto value = values[i];
        auto weight = weights[i];
        auto cell = find_cell();

        std::vector<size_t> cell_dofmap(space_dimension);
        std::vector<double> cell_basis_derivatives(space_dimension * 4);
        calculate_basis_derivative_values(function_space, cell, point, cell_dofmap, cell_basis_derivatives);

        for (size_t j = 0; j < cell_dofmap.size(); j++)
        {
            double result = 0.0;
            for (size_t k = 0; k < 4; k++)
            {
                result += cell_basis_derivatives[4 * i + k] * value[k];
            }
            results[cell_dofmap[i]] += weight * result;
        }
    }
}

Cell find_cell(const Point &point, const FunctionSpace &function_space, const IBMesh &um)
{
    /// must be regular grid.
    auto mesh = function_space.mesh();

    auto index_1 = 2 * um.hash(point);
    auto index_2 = 2 * um.hash(point) + 1;
    Cell cell_1(*mesh, index_1);
    Cell cell_2(*mesh, index_2);

    return cell_1.contains(point) ? cell_1 : cell_2;
}