from fenics import *
from IBInterpolation import *
from IBMesh import *

# 测试一个例子
# 注意points中两个点的顺序，现在x,y坐标必须从小到大。
points = [Point(-2, -2, 0), Point(2, 2, 0)]
seperations = [40, 40, 0]

regular_mesh = IBMesh(points, seperations)
fluid_mesh = regular_mesh.mesh()
solid_mesh = UnitSquareMesh(10,10)

Vf = VectorFunctionSpace(fluid_mesh, "P", 1)
Vs = VectorFunctionSpace(solid_mesh, "P", 1)

IB = IBInterpolation(regular_mesh)

uf = project(Expression(("x[0]*x[0]","x[1]*x[1]"),degree=1), Vf)
us = Function(Vs)

IB.fluid_to_solid(uf._cpp_object, us._cpp_object)
File("us.pvd") << us

n = 100
rd = np.random.rand(2*n)*0.8+0.1
for i in range(n):
    result = us(rd[i*2],rd[i*2+1])
    print(result[0] - rd[i*2]  *rd[i*2])
    print(result[1] - rd[i*2+1]*rd[i*2+1])



us = project(Expression(("x[0]*x[0]","x[1]*x[1]"),degree=1), Vs)
uf = Function(Vf)

IB.solid_to_fluid(uf._cpp_object, us._cpp_object)

for i in range(n):
    result = uf(rd[i*2],rd[i*2+1])
    print(result[0] - rd[i*2]  *rd[i*2])
    print(result[1] - rd[i*2+1]*rd[i*2+1])

