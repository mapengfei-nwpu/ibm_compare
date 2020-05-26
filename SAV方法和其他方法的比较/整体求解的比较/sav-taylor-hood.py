from dolfin import *
from solutions import solutions
import numpy as np

nu = 1
# Load exact solutions
solution = solutions[1]
f = Expression((solution["fx"], solution["fy"]), degree=5, t=0, nu = nu)
p_exact = Expression(solution["p"], degree=5, t=0)
u_exact = Expression((solution["ux"], solution["uy"]), degree=5, t=0)

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
        DirichletBC(W.sub(1), p_exact, "x[0]<DOLFIN_EPS && x[1]<DOLFIN_EPS","pointwise")]
bcs2 = [DirichletBC(W.sub(0), (0, 0), "on_boundary"),
        DirichletBC(W.sub(1), 0, "x[0]<DOLFIN_EPS && x[1]<DOLFIN_EPS","pointwise")]

# Define time step
dt = 0.01/n/n
k = Constant(dt)

# Define parameters
delta = 0.1
T = 0.1


class SAVSolver1:
    def __init__(self):
        # Define test function and trial function
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        # Define coefficient functions
        wn = Function(W)
        (self.un, _) = wn.split(True)
        # Define variational form
        F = inner((u-self.un)/k, v)*dx + nu * inner(grad(u), grad(v)) * \
            dx - div(v)*p*dx + q*div(u)*dx - inner(f, v)*dx
        a = lhs(F)
        self.L = rhs(F)
        self.A = assemble(a)

    def solve(self, un):
        self.w1 = Function(W)
        self.un.assign(un)
        b = assemble(self.L)
        [bc.apply(self.A, b) for bc in bcs1]
        solve(self.A, self.w1.vector(), b)
        return self.w1


class SAVSolver2:
    def __init__(self):
        # Define test function and trial function
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        # Define coefficient functions
        wn = Function(W)
        (self.un, _) = wn.split(True)
        # Define variational form
        F = inner(u/k, v)*dx + inner(grad(self.un)*self.un, v)*dx + \
            nu * inner(grad(u), grad(v))*dx - div(v)*p*dx + q*div(u)*dx
        a = lhs(F)
        self.L = rhs(F)
        self.A = assemble(a)

    def solve(self, un):
        self.w1 = Function(W)
        self.un.assign(un)
        b = assemble(self.L)
        [bc.apply(self.A, b) for bc in bcs2]
        solve(self.A, self.w1.vector(), b)
        return self.w1


class SAVUpdate:
    def __init__(self):
        w = TrialFunctions(W)
        m = TestFunctions(W)
        self.w1 = Function(W)
        self.w2 = Function(W)
        self.wn = Function(W)
        self.S = Expression("S", degree=1, S=1.0)
        a = inner(w[0], m[0])*dx+w[1]*m[1]*dx
        self.A = assemble(a)
        (u1, p1) = self.w1.split()
        (u2, p2) = self.w2.split()
        self.L = inner(u1 + self.S*u2, m[0]) * dx+(p1 + self.S*p2)*m[1]*dx

    def solve(self, w1, w2, S):
        self.w1.assign(w1)
        self.w2.assign(w2)
        self.S.S = S
        b = assemble(self.L)
        solve(self.A, self.wn.vector(), b, "bicgstab", "default")
        return self.wn


def S_solve(u0, u1, u2):
    temp = 0.5*assemble(inner(u0, u0)*dx) + delta
    r = sqrt(temp)
    a2 = nu*assemble(inner(grad(u2), grad(u2))*dx) + 2.0/dt*temp
    a1 = 2*nu*assemble(inner(grad(u2), grad(u1))*dx)  -assemble(inner(f, u2)*dx)- 2 * r/dt*sqrt(temp) #
    a0 = nu*assemble(inner(grad(u1), grad(u1))*dx)  - assemble(inner(f, u1)*dx)
    s1 = (-a1+sqrt(a1*a1-4*a2*a0))/2/a2
    s2 = (-a1-sqrt(a1*a1-4*a2*a0))/2/a2
    print("a0, a1, a2: ", a0, a1, a2)
    print("r: ",r)
    print(s1, s2)
    if max(s1, s2) > 0 and min(s1, s2) < 0:
        return max(s1, s2)
    if abs(s1-1) > abs(s2-1):
        return s2
    return s1


w0 = Function(W)
(u0, _) = w0.split(True)
u0.interpolate(u_exact)

sav_solver_1 = SAVSolver1()
sav_solver_2 = SAVSolver2()
sav_update = SAVUpdate()

t = dt
while t < T + DOLFIN_EPS:
    # Update boundary condition
    f.t = t
    u_exact.t = t
    p_exact.t = t
    # solve midstep values
    w1 = sav_solver_1.solve(u0)
    w2 = sav_solver_2.solve(u0)
    (u1, _) = w1.split(True)
    (u2, _) = w2.split(True)
    S = S_solve(u0, u1, u2)
    # solve final result
    wn = sav_update.solve(w1, w2, S)
    # update and calculate errors
    (un, pn) = wn.split(True)
    u0.assign(un)
    print(np.sqrt(assemble(inner(un-u_exact, un-u_exact)*dx)))
    print(np.sqrt(assemble(inner(pn-p_exact, pn-p_exact)*dx)))
    t += dt
