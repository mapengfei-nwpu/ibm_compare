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
bcs1 = [DirichletBC(W.sub(0), u_exact, "on_boundary"),
        DirichletBC(W.sub(1), p_exact, "x[0]<DOLFIN_EPS && x[1]<DOLFIN_EPS", "pointwise")]
bcs2 = [DirichletBC(W.sub(0), (0, 0), "on_boundary"),
        DirichletBC(W.sub(1), 0, "x[0]<DOLFIN_EPS && x[1]<DOLFIN_EPS", "pointwise")]

# Define time step
dt = 0.01
k = Constant(dt)

# Define parameters
nu = 0.01
delta = 0.1
T = 10


# Define test function and trial function
w=TrialFunctions(W)
m=TestFunctions(W)
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
# Define coefficient functions
(un, _) = Function(W).split(True)
# Define variational form
F1 = inner((u - un)/k, v)*dx + nu * inner(grad(u), grad(v)) * dx - div(v)*p*dx + q*div(u)*dx - inner(f, v)*dx
F2 = inner(u/k, v)*dx + inner(grad(un) * un, v)*dx + nu * inner(grad(u), grad(v))*dx - div(v)*p*dx + q*div(u)*dx

a1 = lhs(F1)
A1 = assemble(a1)
A1.apply()
L1 = rhs(F1)

a2 = lhs(F2)
A2 = assemble(a2)
L2 = rhs(F2)

t = dt
while t < T + DOLFIN_EPS:
    # Update boundary condition
    f.t = t
    u_exact.t = t
    p_exact.t = t
    # solve midstep values
    [bc.apply(A1, b1) for bc in bcs1]
    [bc.apply(A2, b2) for bc in bcs2]
    solve(A1, w1.vector(), b1)
    solve(A2, w2.vector(), b2)
    # solve final result
    wn = project(w1+w2,W)
    # update and calculate errors
    (u_, p_) = wn.split(True)
    un.assign(u_)
    print(np.sqrt(assemble(inner(u_-u_exact, u_-u_exact)*dx)))
    print(np.sqrt(assemble(inner(p_-p_exact, p_-p_exact)*dx)))
    t += dt
