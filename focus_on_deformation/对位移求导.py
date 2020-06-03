from fenics import *
import matplotlib.pyplot as plt

# 解释了为什么要用

mesh = RectangleMesh(Point(0, 0), Point(1, 1), 10, 10)
V = VectorFunctionSpace(mesh, "P", 1)
T = TensorFunctionSpace(mesh, "P", 1)

position = Expression(("x[0] + x[0]*x[0]", "x[1] + x[1]*x[1]"), degree = 2)
displace = Expression(("x[0]*x[0]", "x[1]*x[1]"), degree = 2)

u0 = Function(V)
u0.interpolate(position)

for coordinate in mesh.coordinates():
    coordinate[0] = coordinate[0] + coordinate[0]*coordinate[0]
    coordinate[1] = coordinate[1] + coordinate[1]*coordinate[1]

u1 = project(grad(u0), T)
vec = u1.vector()[:]
for i in range(441):
    print(vec[4*i],", ", vec[4*i+1], ", ", vec[4*i+2], ", ", vec[4*i+3])


plot(mesh)
plt.show()