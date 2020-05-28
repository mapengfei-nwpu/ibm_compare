#include <dolfin.h>
#include <dolfin/geometry/SimplexQuadrature.h>
#include <numeric>
using namespace dolfin;
/// calculate derivatives at gauss point.
Cell find_cell(const Point &point, const IBMesh &um)
{
    /// must be regular grid.
    auto mesh = um.mesh_ptr;

    auto index_1 = 2 * um.hash(point);
    auto index_2 = 2 * um.hash(point) + 1;
    Cell cell_1(*mesh, index_1);
    Cell cell_2(*mesh, index_2);

    return cell_1.contains(point) ? cell_1 : cell_2;
}

void calculate_values_at_gauss_points(
    const Function &displace,
    std::vector<std::vector<double>> &values)
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
        
        auto cell_dofmap = dofmap->cell_dofs(cell->index());

        /// Compute quadrature rule for the cell.
        /// qr.second and qr.first are the coordinate and weight of gauss point respectively.
        auto qr = gq.compute_quadrature_rule(*cell);
        for (size_t i = 0; i < qr.second.size(); i++)
        {
            std::vector<double> point({qr.first[2 * i], qr.first[2 * i + 1]});
            element->evaluate_basis_derivatives_all(
                1,
                basis_derivative_values.data(),
                point.data(),
                coordinate_dofs.data(),
                ufc_cell.orientation);
            std::vector<double> derivative_value(derivative_value_size);
            for (size_t k = 0; k < cell_dofmap.size(); k++)
            {
                for (size_t j = 0; j < derivative_value_size; j++)
                {
                    derivative_value[j] +=
                        (*displace.vector())[cell_dofmap[k]] *
                        basis_derivative_values[derivative_value_size * k + j];
                }
            }
            for (size_t k = 0; k < 4; k++){
                std::cout<<derivative_value[k]<<std::endl;
            }
            values.push_back(derivative_value);
        }
    }
}

void calculate_gauss_points_and_weights(
    const Function &displace,
    std::vector<std::vector<double>> &points,
    std::vector<double> &weights)
{
    // Construct Gauss quadrature rules
    // dimension 2 and order 9
    SimplexQuadrature gq(2, 9);
    auto mesh = displace.function_space()->mesh();

    for (CellIterator cell(*mesh); !cell.end(); ++cell)
    {
        /// Compute quadrature rule for the cell.
        /// qr.second and qr.first are the coordinate and weight of gauss point respectively.
        auto qr = gq.compute_quadrature_rule(*cell);
        for (size_t i = 0; i < qr.second.size(); i++)
        {
            std::vector<double> point({qr.first[2 * i], qr.first[2 * i + 1]});
            points.push_back(point);
            weights.push_back(qr.second[i]);
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

    for (size_t i = 0; i < cell_dofmap.size(); i++)
    {
        cell_dofmap[i] = cell_dofmap_eigen[i];
    }
}

std::vector<double> source_assemble(
    std::vector<std::vector<double>> &points,
    std::vector<std::vector<double>> &values,
    std::vector<double> &weights,
    const FunctionSpace &function_space,
    const IBMesh &um)
{
    std::vector<double> results(function_space.dim());
    auto space_dimension = function_space.element()->space_dimension();

    for (size_t i = 0; i < points.size(); i++)
    {
        // std::cout<<points.size()<<std::endl;
        auto point = points[i];
        auto value = values[i];
        auto weight = weights[i];
        Point point_temp(point[0],point[1]);
        auto cell = find_cell(point_temp, um);

        std::vector<size_t> cell_dofmap(space_dimension);
        std::vector<double> cell_basis_derivatives(space_dimension * 4);

        calculate_basis_derivative_values(function_space, cell, point, cell_dofmap, cell_basis_derivatives);

        for (size_t j = 0; j < cell_dofmap.size(); j++)
        {
            double result = 0.0;
            for (size_t k = 0; k < 4; k++)
            {
                result += cell_basis_derivatives[4 * j + k] * value[k];
            }
            results[cell_dofmap[j]] += weight * result;
        }
    }
    return results;
}

