#include "IBInterpolation.h"

#include <dolfin.h>
#include "Poisson.h"
using namespace dolfin;

class Source : public Expression
{
public:
    // Constructor
    Source() : Expression(2) {}

    // Evaluate pressure at inflow
    void eval(Array<double> &values, const Array<double> &x) const
    {
        values[0] = 4*x[0]+3*x[1];
        values[1] = 2*x[0]+1*x[1];
    }
};


int main()
{
    // Create chanel mesh
    size_t nnn = 32;
    Point point0(0.0, 0.0);
    Point point1(1.0, 1.0);
    IBMesh ba({point0, point1}, {nnn, nnn});

    auto V = std::make_shared<Poisson::FunctionSpace>(ba.mesh());
    auto u = std::make_shared<Function>(V);

    auto source = std::make_shared<Source>();
    u->interpolate(*source);

    std::vector<std::vector<double>> coordinates;
	std::vector<std::vector<double>> values;
	std::vector<double> weights;

    get_gauss_rule(*u, coordinates, values, weights);

    std::cout << "size: "
              << coordinates.size()
              << values.size()
              << weights.size()
              << std::endl;

    for (size_t i = 0; i < coordinates.size(); i++)
    {
        auto coordinate = coordinates[i];
        auto value = values[i];
        auto weight = weights[i];
        std::cout << "coordinate: " 
                  << coordinate[0] << ", " 
                  << coordinate[1] 
                  << std::endl;
        std::cout << "value:  " << value[0] << ", " 
                                << value[1] << ", "     
                                << value[2] << ", " 
                                << value[3]
                                << std::endl; 
        std::cout << "weight: " << weight
                                << std::endl;

    }
    
}