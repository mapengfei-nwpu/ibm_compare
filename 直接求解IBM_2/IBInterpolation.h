#include <dolfin.h>
#include <dolfin/geometry/SimplexQuadrature.h>
#include <numeric>
using namespace dolfin;

/// calculate derivatives at gauss point.
void get_gauss_rule(
	const Function &f,
	std::vector<std::vector<double>> &coordinates,
	std::vector<std::vector<double>> &values,
	std::vector<double> &weights)
{
	// Construct Gauss quadrature rules
    // dimension 2 and order 9
	SimplexQuadrature gq(2, 9);
	auto mesh = f.function_space()->mesh();
    auto dofmap = f.function_space()->dofmap();
    auto element = f.function_space()->element();
    auto value_size = f.value_size();
    auto derivative_value_size = value_size*value_size;
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
            std::vector<double> point({qr.first[2*i],qr.first[2*i+1]});
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
                        (*f.vector())[cell_dofmap[i]] * 
                        basis_derivative_values[derivative_value_size*i + j];
                }
		    }
            values.push_back(derivative_value);
            coordinates.push_back(point);
            weights.push_back(weight);
		}
	}
}
