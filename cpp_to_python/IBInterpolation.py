from fenics import *
from IBInterpolation import *
from IBMesh import *

# 测试一个例子
# 注意points中两个点的顺序，现在x,y坐标必须从小到大。
points = [Point(-1, -1, 0), Point(2, 2, 0)]
seperations = [100, 100, 0]

regular_mesh = IBMesh(points, seperations)
fluid_mesh = regular_mesh.mesh()
solid_mesh = UnitSquareMesh(200,200)

Vf = VectorFunctionSpace(fluid_mesh, "P", 2)
Vs = VectorFunctionSpace(solid_mesh, "P", 2)

IB = IBInterpolation(regular_mesh)

uf = project(Expression(("x[0]","x[1]"),degree=1), Vf)
us = Function(Vs)

IB.fluid_to_solid(uf._cpp_object, us._cpp_object)
File("us.pvd") << us



uf = Function(Vf)
us = project(Expression(("x[0]","x[1]"),degree=1), Vs)
result = IB.solid_to_fluid(uf._cpp_object, us._cpp_object)
uf.vector().set_local(result)
File("uf.pvd") << uf

u = TrialFunction(Vf)
v = TestFunction(Vf)
u0 = Function(Vf)

A=assemble(inner(u,v)*dx)
solve(A, u0.vector(), uf.vector())
File("u0.pvd") << u0

max_err = 0
for i in range(100):
    for j in range(100):
        a = u0(i/100, j/100)
        max_err = max(max_err,abs(a[0]-i/100))
        max_err = max(max_err,abs(a[1]-i/100))

print(max_err)