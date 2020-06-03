from fenics import *
import matplotlib.pyplot as plt

# 解释了为什么要用

mesh = RectangleMesh(Point(0, 0), Point(1, 1), 10, 10)
V = VectorFunctionSpace(mesh, "P", 1)
T = TensorFunctionSpace(mesh, "P", 1)

position = Expression(("x[0] + x[0]*x[0]", "x[1] + x[1]*x[1]"), degree = 2)
displace = Expression(("x[0]*x[0]", "x[1]*x[1]"), degree = 2)

u0 = Function(V)

import test_function

test_function.size(u0)