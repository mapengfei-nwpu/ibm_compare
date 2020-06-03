#include <vector>
#include <set>
#include <stdio.h>
#include <iostream>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/MeshDomains.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/MeshConnectivity.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshEntityIterator.h>
#include <dolfin/mesh/MeshFunction.h>

#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Expression.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/common/Array.h>
#include <dolfin/fem/fem_utils.h>
#include <dolfin/common/types.h>
#include <dolfin/common/Hierarchical.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/FunctionAXPY.h>
#include <dolfin/geometry/BoundingBoxTree.h>

using namespace dolfin;

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
namespace py = pybind11;
PYBIND11_MAKE_OPAQUE(std::vector<double>);


class IBInterpolation
{
private:
    Mesh fluid;
    Mesh solid;
    size_t fluid_size;
    size_t solid_size;
    std::vector<double> sp2fc;
    std::vector<double> coeff;
    void construct_sp2fc();
    void test_sp2fc();
    void construct_coeff();
public:
    IBInterpolation(const Mesh& fluid, const Mesh& solid);
    void fluid_to_solid(Function& fluid, Function& solid);
    void solid_to_fluid(Function& fluid, Function& solid);
    ~IBInterpolation();
};

IBInterpolation::IBInterpolation(const Mesh& fluid, const Mesh& solid):
    fluid(fluid),
    solid(solid),
    solid_size(solid.num_vertices()),
    fluid_size(fluid.num_vertices()),
    sp2fc(solid_size),
    coeff(solid_size*4)
{
    construct_sp2fc();
    test_sp2fc();
    construct_coeff();
}

IBInterpolation::~IBInterpolation()
{
}

void IBInterpolation::construct_sp2fc(){
    std::vector<double> solid_coord = solid.coordinates();
    for (size_t i = 0; i < solid_size; i++)
    {
        //const Cell point_in_which_cell(const Mesh& mesh, const double* x){
        // this ideal is from function "void Function::eval(Array<double>& values, const Array<double>& x);"
        const double x[3] = {solid_coord[i*3],solid_coord[i*3+1],solid_coord[i*3+2]};
        const Point point(3, x);
        // Get index of first cell containing point
        unsigned int id = fluid.bounding_box_tree()->compute_first_entity_collision(point);

        // If not found, use the closest cell
        //if (id == std::numeric_limits<unsigned int>::max())
            // id = fluid.bounding_box_tree()->compute_closest_entity(point).first;
        printf("%d\n",id);
        sp2fc[i] = id;
    }
    
}
void IBInterpolation::test_sp2fc(){
    std::vector<double> solid_coord = solid.coordinates();
    std::vector<double> fluid_coord = fluid.coordinates();
    for (size_t i = 0; i < solid_size; i++)
    {
        Cell cell(fluid, sp2fc[i]);
        for(MeshEntityIterator e(cell,0);!e.end();++e){
            size_t index = e->index();
            printf("fluid cell coordinate:(x:%.8lf,y:%.8lf,z:%.8lf)\n",fluid_coord[index*3+0],fluid_coord[index*3+1],fluid_coord[index*3+2]);
        }
        printf(    "solid coordinate:     (x:%.8lf,y:%.8lf,z:%.8lf)\n",solid_coord[i*3+0],solid_coord[i*3+1],solid_coord[i*3+2]);
    }
}
void IBInterpolation::construct_coeff(){
    std::vector<double> solid_coord = solid.coordinates();
    std::vector<double> fluid_coord = fluid.coordinates();
    for (size_t i = 0; i < solid_size; i++)
    {
        double M[4][4];
        Cell cell(fluid, sp2fc[i]);
        const unsigned int* point = cell.entities(0);
        for(size_t j=0; j<4; j++){
            size_t index = point[j];
            M[0][j]=1.0;
            M[1][j]=fluid_coord[index*3+0];
            M[2][j]=fluid_coord[index*3+1];
            M[3][j]=fluid_coord[index*3+2];
        }
        double P[4] = {1.0,solid_coord[i*3+0],solid_coord[i*3+1],solid_coord[i*3+2]};
        double chushu = M[1][0]*M[2][1]*M[3][2] - M[1][0]*M[2][1]*M[3][3] - M[1][0]*M[2][2]*M[3][1] + M[1][0]*M[2][2]*M[3][3] + M[1][0]*M[2][3]*M[3][1] - M[1][0]*M[2][3]*M[3][2] - M[1][1]*M[2][0]*M[3][2] + M[1][1]*M[2][0]*M[3][3] + M[1][1]*M[2][2]*M[3][0] - M[1][1]*M[2][2]*M[3][3] - M[1][1]*M[2][3]*M[3][0] + M[1][1]*M[2][3]*M[3][2] + M[1][2]*M[2][0]*M[3][1] - M[1][2]*M[2][0]*M[3][3] - M[1][2]*M[2][1]*M[3][0] + M[1][2]*M[2][1]*M[3][3] + M[1][2]*M[2][3]*M[3][0] - M[1][2]*M[2][3]*M[3][1] - M[1][3]*M[2][0]*M[3][1] + M[1][3]*M[2][0]*M[3][2] + M[1][3]*M[2][1]*M[3][0] - M[1][3]*M[2][1]*M[3][2] - M[1][3]*M[2][2]*M[3][0] + M[1][3]*M[2][2]*M[3][1];       
        coeff[i*4]=(-M[1][1]*M[2][2]*M[3][3] + M[1][1]*M[2][3]*M[3][2] + M[1][2]*M[2][1]*M[3][3] - M[1][2]*M[2][3]*M[3][1] - M[1][3]*M[2][1]*M[3][2] + M[1][3]*M[2][2]*M[3][1] + P[1]*(M[2][1]*M[3][2] - M[2][1]*M[3][3] - M[2][2]*M[3][1] + M[2][2]*M[3][3] + M[2][3]*M[3][1] - M[2][3]*M[3][2]) - P[2]*(M[1][1]*M[3][2] - M[1][1]*M[3][3] - M[1][2]*M[3][1] + M[1][2]*M[3][3] + M[1][3]*M[3][1] - M[1][3]*M[3][2]) + P[3]*(M[1][1]*M[2][2] - M[1][1]*M[2][3] - M[1][2]*M[2][1] + M[1][2]*M[2][3] + M[1][3]*M[2][1] - M[1][3]*M[2][2]))/chushu; 
        coeff[i*4+1]=(M[1][0]*M[2][2]*M[3][3] - M[1][0]*M[2][3]*M[3][2] - M[1][2]*M[2][0]*M[3][3] + M[1][2]*M[2][3]*M[3][0] + M[1][3]*M[2][0]*M[3][2] - M[1][3]*M[2][2]*M[3][0] - P[1]*(M[2][0]*M[3][2] - M[2][0]*M[3][3] - M[2][2]*M[3][0] + M[2][2]*M[3][3] + M[2][3]*M[3][0] - M[2][3]*M[3][2]) + P[2]*(M[1][0]*M[3][2] - M[1][0]*M[3][3] - M[1][2]*M[3][0] + M[1][2]*M[3][3] + M[1][3]*M[3][0] - M[1][3]*M[3][2]) - P[3]*(M[1][0]*M[2][2] - M[1][0]*M[2][3] - M[1][2]*M[2][0] + M[1][2]*M[2][3] + M[1][3]*M[2][0] - M[1][3]*M[2][2]))/chushu;
        coeff[i*4+2]=(-M[1][0]*M[2][1]*M[3][3] + M[1][0]*M[2][3]*M[3][1] + M[1][1]*M[2][0]*M[3][3] - M[1][1]*M[2][3]*M[3][0] - M[1][3]*M[2][0]*M[3][1] + M[1][3]*M[2][1]*M[3][0] + P[1]*(M[2][0]*M[3][1] - M[2][0]*M[3][3] - M[2][1]*M[3][0] + M[2][1]*M[3][3] + M[2][3]*M[3][0] - M[2][3]*M[3][1]) - P[2]*(M[1][0]*M[3][1] - M[1][0]*M[3][3] - M[1][1]*M[3][0] + M[1][1]*M[3][3] + M[1][3]*M[3][0] - M[1][3]*M[3][1]) + P[3]*(M[1][0]*M[2][1] - M[1][0]*M[2][3] - M[1][1]*M[2][0] + M[1][1]*M[2][3] + M[1][3]*M[2][0] - M[1][3]*M[2][1]))/chushu;
        coeff[i*4+3]=(M[1][0]*M[2][1]*M[3][2] - M[1][0]*M[2][2]*M[3][1] - M[1][1]*M[2][0]*M[3][2] + M[1][1]*M[2][2]*M[3][0] + M[1][2]*M[2][0]*M[3][1] - M[1][2]*M[2][1]*M[3][0] - P[1]*(M[2][0]*M[3][1] - M[2][0]*M[3][2] - M[2][1]*M[3][0] + M[2][1]*M[3][2] + M[2][2]*M[3][0] - M[2][2]*M[3][1]) + P[2]*(M[1][0]*M[3][1] - M[1][0]*M[3][2] - M[1][1]*M[3][0] + M[1][1]*M[3][2] + M[1][2]*M[3][0] - M[1][2]*M[3][1]) - P[3]*(M[1][0]*M[2][1] - M[1][0]*M[2][2] - M[1][1]*M[2][0] + M[1][1]*M[2][2] + M[1][2]*M[2][0] - M[1][2]*M[2][1]))/chushu;
        printf("chushu:%.8lf\n",chushu);
        printf("coefficient 1:%.8lf\n",coeff[i*4]);
    }
}

void IBInterpolation::fluid_to_solid(Function& Ffluid, Function& Fsolid){

    auto v2d_fluid = vertex_to_dof_map(*Ffluid.function_space());
    auto v2d_solid = vertex_to_dof_map(*Fsolid.function_space());
    std::vector<double> dof_fluid(fluid_size, 0.0);
    std::vector<double> dof_solid(solid_size, 0.0);
    Ffluid.vector()->get_local(dof_fluid);
    for (size_t i = 0; i < solid_size; i++){
        Cell cell(fluid,sp2fc[i]);
        const unsigned int* point = cell.entities(0);
        for (size_t j = 0; j < 4; j++){
            dof_solid[v2d_solid[i]] += dof_fluid[v2d_fluid[point[j]]]*coeff[i*4+j];
            printf(" dof of p1: %.8lf\n", i, dof_fluid[v2d_fluid[point[j]]]);
        }
        printf(" The dof of NO.%d dof_solid %.8lf\n", i, dof_solid[v2d_solid[i]]);
    }
    Fsolid.vector()->set_local(dof_solid);
}
void IBInterpolation::solid_to_fluid(Function& Ffluid, Function& Fsolid){

    auto v2d_fluid = vertex_to_dof_map(*Ffluid.function_space());
    auto v2d_solid = vertex_to_dof_map(*Fsolid.function_space());
    std::vector<double> dof_fluid(fluid_size, 0.0);
    std::vector<double> dof_solid(solid_size, 0.0);
    Fsolid.vector()->get_local(dof_solid);
    for (size_t i = 0; i < solid_size; i++){
        Cell cell(fluid,sp2fc[i]);
        const unsigned int* point = cell.entities(0);
        for (size_t j = 0; j < 4; j++){
            dof_fluid[v2d_fluid[point[j]]] += dof_solid[v2d_solid[i]]*coeff[i*4+j];
            // printf(" dof of p1: %.8lf\n", i, dof_fluid[v2d_fluid[point[j]]]);
        }
        // printf(" The dof of NO.%d dof_solid %.8lf\n", i, dof_solid[v2d_solid[i]]);
    }
    Ffluid.vector()->set_local(dof_fluid);
}

PYBIND11_MODULE(SIGNATURE, m)
{
    py::class_<IBInterpolation>(m, "IBInterpolation")
        .def(py::init<const Mesh &, const Mesh &>())
        .def("fluid_to_solid",&IBInterpolation::fluid_to_solid)
        .def("solid_to_fluid",&IBInterpolation::solid_to_fluid);
}
/*
from fenics import *
from mshr   import *
import numpy as np

cpp_code = open('IBInterpolation.cpp', 'r').read()
testMesh = compile_cpp_code(cpp_code)

# define entire computing area
cylinder = Cylinder(Point(0,0,0),Point(0,0,1),1,1)
sphere = Sphere(Point(0,0,0.5),0.5)
fluid = generate_mesh(cylinder,10)
solid = generate_mesh(sphere,10)

Vf = FunctionSpace(fluid, 'P', 1)
Vs = FunctionSpace(solid, 'P', 1)


u0 = Expression('x[0]+x[1]', degree=1)
uf = interpolate(u0, Vf)
us = Function(Vs)

IB = testMesh.IBInterpolation(fluid, solid)
IB.fluid_to_solid(uf._cpp_object,us._cpp_object)

File("us.pvd")<< us
*/

/*
from fenics import *
from mshr   import *
import numpy as np

cpp_code = open('IBInterpolation.cpp', 'r').read()
testMesh = compile_cpp_code(cpp_code)

# define entire computing area
cylinder = Cylinder(Point(0,0,0),Point(0,0,1),1,1)
sphere = Sphere(Point(0,0,0.5),0.5)
fluid = generate_mesh(cylinder,10)
solid = generate_mesh(sphere,10)

Vf = FunctionSpace(fluid, 'P', 1)
Vs = FunctionSpace(solid, 'P', 1)


u0 = Expression('x[0]+x[1]', degree=1)
us = interpolate(u0, Vs)
uf = Function(Vf)

IB = testMesh.IBInterpolation(fluid, solid)
IB.solid_to_fluid(uf._cpp_object,us._cpp_object)

File("uf.pvd")<< uf
*/