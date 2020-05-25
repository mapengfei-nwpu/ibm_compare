# SAV程序，取S=1
from dolfin import *
from solutions import solutions
import numpy as np

# Load exact solutions
solution = solutions[1]
f = Expression((solution["fx"], solution["fy"]), degree=5, t=0)
p_exact = Expression(solution["p"], degree=5, t=0)
u_exact = Expression((solution["ux"], solution["uy"]), degree=5, t=0)
SS = Expression("S", degree=0, S=1)
# Generate mesh
n = 32
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
dt = 1/n/n
k = Constant(dt)

# Define parameters
nu = 0.01
delta = 0.1
T = 0.1


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
L1 = rhs(F1)

a2 = lhs(F2)
A2 = assemble(a2)
L2 = rhs(F2)


def calculate_s(w1,w2,un):
    (u1,p1)=w1.split()
    (u2,p2)=w2.split()
    temp = 0.5*assemble(inner(un, un)*dx) + delta
    r = sqrt(temp)
    a2 = nu*assemble(inner(grad(u2), grad(u2))*dx) + 2/dt*temp
    a1 = 2*nu*assemble(inner(grad(u2), grad(u1))*dx) - 2 * \
        r/dt*sqrt(temp)-assemble(inner(f, u2)*dx)
    a0 = nu*assemble(inner(grad(u1), grad(u1))*dx) - assemble(inner(f, u1)*dx)
    s1 = (-a1+sqrt(a1*a1-4*a2*a0))/2/a2
    s2 = (-a1-sqrt(a1*a1-4*a2*a0))/2/a2
    print(s1, s2)
    if max(s1, s2) > 0 and min(s1, s2) < 0:
        return max(s1, s2)
    if abs(s1-1) > abs(s2-1):
        return s2
    else:
        return s1



t = dt
un.interpolate(u_exact)
w1 = Function(W)
w2 = Function(W)
while t < T + DOLFIN_EPS:
    # Update boundary condition
    f.t = t
    u_exact.t = t
    p_exact.t = t
    # solve midstep values
    b1 = assemble(L1)
    b2 = assemble(L2)
    [bc.apply(A1, b1) for bc in bcs1]
    [bc.apply(A2, b2) for bc in bcs2]
    solve(A1, w1.vector(), b1)
    solve(A2, w2.vector(), b2)
    S = calculate_s(w1, w2, un)
    print(S)
    SS.S = S
    # solve final result
    wn = project(w1+SS*w2,W)
    # update and calculate errors
    (u_, p_) = wn.split(True)
    un.assign(u_)
    t += dt

print(np.sqrt(assemble(inner(u_-u_exact, u_-u_exact)*dx)))
print(np.sqrt(assemble(inner(p_-p_exact, p_-p_exact)*dx)))
