from fenics import *
import matplotlib.pyplot as plt

# 解释了为什么要用

mesh = RectangleMesh(Point(0, 0), Point(1, 1), 10, 10)
V = VectorFunctionSpace(mesh, "P", 1)
u0 = Function(V)

import test_function

test_function.size(u0._cpp_object)