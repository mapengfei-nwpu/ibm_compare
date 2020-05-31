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
    std::vector<double> &weights,
    std::vector<std::vector<double>> &points,
    std::vector<std::vector<double>> &values)
{
    // Construct Gauss quadrature rules
    // dimension 2 and order 9
    SimplexQuadrature gq(2, 16);
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
            for (size_t k = 0; k < static_cast<size_t>(cell_dofmap.size()); k++)
            {
                for (size_t j = 0; j < derivative_value_size; j++)
                {
                    derivative_value[j] +=
                        (*displace.vector())[cell_dofmap[k]] *
                        basis_derivative_values[derivative_value_size * k + j];
                }
            }
            

            {
                double a = derivative_value[0];
                double b = derivative_value[1];
                double c = derivative_value[2];
                double d = derivative_value[3];
                
                double det = 1.0/(a*d-b*c);
                derivative_value[0] =  (b*b+d*d)/(det*det);
                derivative_value[1] = -(c*d+a*b)/(det*det);
                derivative_value[2] = -(c*d+a*b)/(det*det);
                derivative_value[3] =  (a*a+c*c)/(det*det);

                /// std::cout << derivative_value[0] << ", " << derivative_value[1] << ", " << derivative_value[2] << ", " << derivative_value[3] <<std::endl;
            }
            weights.push_back(qr.second[i]);
            points.push_back(point);
            values.push_back(derivative_value);
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

    ufc::cell ufc_cell;
    cell.get_cell_data(ufc_cell);

    std::vector<double> coordinate_dofs;
    cell.get_coordinate_dofs(coordinate_dofs);

    /// store values derivative basis.
    auto element = function_space.element();
    element->evaluate_basis_derivatives_all(
        1,
        cell_basis_derivatives.data(),
        point.data(),
        coordinate_dofs.data(),
        ufc_cell.orientation);
    auto cell_dofmap_eigen = function_space.dofmap()->cell_dofs(cell.index());

    for (size_t i = 0; i < cell_dofmap.size(); i++)
    {
        cell_dofmap[i]  = cell_dofmap_eigen[i];
    }
}

std::vector<double> source_assemble(
    const std::vector<double> &weights,
    const std::vector<std::vector<double>> &points,
    const std::vector<std::vector<double>> &values,
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

        Point point_temp(2, point.data());
        auto cell = find_cell(point_temp, um);
        /*
        if (cell.contains(point_temp)) {
            std::cout<< "contain." << std::endl;
        } else {
            std::cout<< "contain.contain.contain.contain.contain.contain.contain.contain.contain.contain.contain.contain.contain.contain.contain.contain.contain.contain.contain." << std::endl;

        }*/
        /// std::cout<<point_temp<<std::endl;
        /// std::vector<double> cell_coordinates;
        /// cell.get_vertex_coordinates(cell_coordinates);
        /// for (size_t j = 0; j < cell_coordinates.size(); j++)
        /// {
        ///    std::cout << cell_coordinates[j] << std::endl;
        /// }
        

        std::vector<size_t> cell_dofmap(space_dimension);
        std::vector<double> cell_basis_derivatives(space_dimension * 4);

        calculate_basis_derivative_values(function_space, cell, point, cell_dofmap, cell_basis_derivatives);

        for (size_t j = 0; j < cell_dofmap.size(); j++)
        {
            for (size_t k = 0; k < 4; k++)
            {
                ///  *cell_basis_derivatives[j*4 + k]
                results[cell_dofmap[j]] += weight* value[k];
            }
        }
    }
    return results;
}

