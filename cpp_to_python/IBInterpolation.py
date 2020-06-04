from fenics import *
import IBInterpolation
import IBMesh

# 测试一个例子
# 注意points中两个点的顺序，x,y坐标必须从小到大。
points = [Point(0, 0, 0), Point(1, 1, 0)]
seperations = [10, 10, 0]

regular_mesh = IBMesh.IBMesh(points, seperations)
mesh = regular_mesh.mesh()
cell = Cell(mesh, 0) 
