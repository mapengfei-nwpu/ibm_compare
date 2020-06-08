from fenics import *
from IBInterpolation import *
from IBMesh import *

# 测试一个例子
# 注意points中两个点的顺序，现在x,y坐标必须从小到大。
points = [Point(0, 0, 0), Point(1, 1, 0)]
seperations = [2, 2, 0]

regular_mesh = IBMesh(points, seperations)
fluid_mesh = regular_mesh.mesh()
solid_mesh = UnitSquareMesh(2,2)

Vf = VectorFunctionSpace(fluid_mesh, "P", 1)
Vs = VectorFunctionSpace(solid_mesh, "P", 1)

IB = IBInterpolation(regular_mesh)

uf = project(Expression(("x[0]","x[1]"),degree=1), Vf)
us = Function(Vs)

IB.fluid_to_solid(uf._cpp_object, us._cpp_object)
File("us.pvd") << us

uf = Function(Vf)
us = project(Expression(("x[0]","x[1]"),degree=1), Vs)
result = IB.solid_to_fluid(uf._cpp_object, us._cpp_object)

u = TrialFunction(Vf)
v = TestFunction(Vf)
u0 = Function(Vf)

A=assemble(inner(u,v)*dx)
b=assemble(inner(Expression(("x[0]","x[1]"),degree=1),v)*dx)

for i in range(len(b)):
    print(b[i]-result[i])