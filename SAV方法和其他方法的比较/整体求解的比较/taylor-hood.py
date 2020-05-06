from dolfin import *
from solutions import solutions
import numpy as np

# Load exact solutions
solution = solutions[1]
f = Expression((solution["fx"], solution["fy"]), degree=5, t=0)
p_exact = Expression(solution["p"], degree=5, t=0)
u_exact = Expression((solution["ux"], solution["uy"]), degree=5, t=0)

# Generate mesh
n = 64
mesh = RectangleMesh(Point(0, 0), Point(1, 1), n, n)

# Define function spaces
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)

# Define boundary conditions
bcs = [DirichletBC(W.sub(0), u_exact, "on_boundary"),
       DirichletBC(W.sub(1), p_exact, "on_boundary")]

# Set parameters
dt = 1e-5
T = 0.1
nu = 0.01
k = Constant(dt)
wn = Function(W)
(un, pn) = wn.split(True)    # from now on, wn can't be used.

# Define variational problem
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
F = inner((u-un)/k, v)*dx + inner(grad(un)*un, v)*dx + nu*inner(grad(u), grad(v))*dx - (div(v)*p + q*div(u))*dx - inner(f, v)*dx
a = lhs(F)
L = rhs(F)

# Assemble matrix and vector
A = assemble(a)

t = dt
un.interpolate(u_exact)
while t < T + DOLFIN_EPS:
    # Update boundary condition
    p_exact.t = t
    u_exact.t = t
    f.t = t
    # Assemble right side term
    b = assemble(L)
    # Compute solution
    w_ = Function(W)
    [bc.apply(A, b) for bc in bcs]
    solve(A, w_.vector(), b)
    # Update coefficients
    (u_, p_) = w_.split(True)
    un.assign(u_)
    # Print errors
    print("||u||_2: ", np.sqrt(assemble(inner((u_-u_exact), (u_-u_exact))*dx)))
    print("||p||_2: ", np.sqrt(assemble((p_-p_exact)*(p_-p_exact)*dx)))
