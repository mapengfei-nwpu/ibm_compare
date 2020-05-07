import matplotlib.pyplot as plt
from solutions import solutions
from dolfin import *
# Set exact results
# f = Expression(("2*cos(pi*x[1])*sin(pi*x[0])*cos(t) + 4*pi*cos(pi*x[0])*sin(pi*x[0]) - 4*pi*cos(pi*x[0])*sin(pi*x[0])*cos(t)*cos(t) + 4*pi*pi*cos(pi*x[1])*sin(pi*x[0])*sin(t) + 2*pi*cos(pi*x[0])*sin(pi*x[1])*cos(t)",
#                 "4*pi*cos(pi*x[1])*sin(pi*x[1]) - 2*cos(pi*x[0])*sin(pi*x[1])*cos(t) - 4*pi*cos(pi*x[1])*sin(pi*x[1])*cos(t)*cos(t) - 4*pi*pi*cos(pi*x[0])*sin(pi*x[1])*sin(t) + 2*pi*cos(pi*x[1])*sin(pi*x[0])*cos(t)"), t=0, degree = 3)
# p_exact = Expression("2*sin(pi*x[1])*sin(pi*x[0])*cos(t)", degree=3, t=0)
# u_exact = Expression(("2*cos(pi*x[1])*sin(pi*x[0])*sin(t)", "-2*sin(pi*x[1])*cos(pi*x[0])*sin(t)"), degree=3, t=0)
solution = solutions[0]
f = Expression((solution["fx"], solution["fy"]), degree=2, t=0)
p_exact = Expression(solution["p"], degree=2, t=0)
u_exact = Expression((solution["ux"], solution["uy"]), degree=2, t=0)

# Set exact results
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"
parameters['krylov_solver']['nonzero_initial_guess'] = True
# Set parameter values
dt = 0.0001
nu = 1
delta = 0.001
k = Constant(dt)
# Define function spaces (P2-P1)
mesh = RectangleMesh(Point(0, 0), Point(1, 1), 40, 40)
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)
uu = TrialFunction(V)
pp = TrialFunction(Q)
vv = TestFunction(V)
qq = TestFunction(Q)
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)
# Define trial and test functions
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
# Define boundary conditions
bcu1 = [DirichletBC(W.sub(0), u_exact, "on_boundary")]
bcu2 = [DirichletBC(W.sub(0), (0, 0), "on_boundary")]
bcp = [DirichletBC(W.sub(1), 0, "near(x[0],0) && near(x[1],-1)","pointwise")]


class SAVSolver1:
    def __init__(self):
        # Create functions
        self.w0 = Function(W)
        self.w1 = Function(W)
        (self.u0, self.p0) = self.w0.split()
        (self.u1, self.p1) = self.w1.split()
        F1 = (1/dt)*inner((u-self.u0),v)*dx + nu*(inner(grad(u), grad(v)))*dx - div(v)*p + q*div(u)*dx - inner(f,v)*dx
        a1 = lhs(F1)
        self.L1 = rhs(F1)
        self.A1 = assemble(a1)

    def solve(self, u0):
        self.u0.assign(u0)
        b1 = assemble(self.L1)
        [bc.apply(self.A1, b1) for bc in bcu1]
        [bc.apply(self.A1, b1) for bc in bcp]
        solve(self.A1, self.w1.vector(), b1)
        # return self.u1, self.p1
        return self.w1.split(True)


class SAVSolver2:
    def __init__(self):
        # Create functions
        self.w0 = Function(W)
        self.w1 = Function(W)
        (self.u0, self.p0) = self.w0.split()
        (self.u1, self.p1) = self.w1.split()
        F1 = (1/dt)*inner(u,v)*dx + inner(grad(self.u0)*self.u0, v)*dx + nu*(inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx
        a1 = lhs(F1)
        self.L1 = rhs(F1)
        self.A1 = assemble(a1)
        
    def solve(self, u0):
        self.u0.assign(u0)
        b1 = assemble(self.L1)
        [bc.apply(self.A1, b1) for bc in bcu1]
        [bc.apply(self.A1, b1) for bc in bcp]
        solve(self.A1, self.w1.vector(), b1)
        # return self.u1, self.p1
        return self.w1.split(True)


def S_solve(u0, u1, u2):
    temp = 0.5*assemble(inner(u0, u0)*dx) + delta
    r = sqrt(temp)
    a2 = nu*assemble(inner(grad(u2), grad(u2))*dx) + 0.5*dt*temp
    a1 = 2*nu*assemble(inner(grad(u2), grad(u1))*dx) - 2*r/dt*sqrt(temp)
    a0 = nu*assemble(inner(grad(u1), grad(u1))*dx)
    s1 = (-a1+sqrt(a1*a1-4*a2*a0))/2/a2
    s2 = (-a1-sqrt(a1*a1-4*a2*a0))/2/a2
    if max(s1, s2) > 0 and min(s1, s2) < 0:
        return max(s1, s2)
    if abs(s1-1) > abs(s2-1):
        return s2
    else:
        return s1


# Create files for storing solution
ufile = File("results/velocity.pvd")
pfile = File("results/pressure.pvd")

class SAVUpdate:
    def __init__(self):
        self.u1 = Function(V)
        self.p1 = Function(Q)
        self.u2 = Function(V)
        self.p2 = Function(Q)
        self.un = Function(V)
        self.pn = Function(Q)
        self.SS = Expression("S", degree=1, S=1.0)
        a1 = inner(uu, vv)*dx
        self.A1 = assemble(a1)
        self.L1 = inner(self.u1 + self.SS*self.u2, vv)*dx
        a2 = inner(pp, qq)*dx
        self.A2 = assemble(a2)
        self.L2 = inner(self.p1 + self.SS*self.p2, qq)*dx

    def solve(self, S):
        self.SS.S = S
        b1 = assemble(self.L1)
        b2 = assemble(self.L2)
        solve(self.A1, self.un.vector(), b1, "bicgstab", "default")
        solve(self.A2, self.pn.vector(), b2, "bicgstab", "default")
        return self.un, self.pn


sav_solver_1 = SAVSolver1()
sav_solver_2 = SAVSolver2()
sav_update = SAVUpdate()
step = 0
u0 = Function(V)
t = dt
err_u = []
err_p = []
while step < 100:
    f.t = t
    u_exact.t = t
    p_exact.t = t
    (u1, p1) = sav_solver_1.solve(u0)
    (u2, p2) = sav_solver_2.solve(u0)
    sav_update.p1.assign(p1)
    sav_update.p2.assign(p2)
    sav_update.u1.assign(u1)
    sav_update.u2.assign(u2)
    S = S_solve(u0, u1, u2)
    (un,pn) = sav_update.solve(S)
    u0.assign(un)
    print(assemble(inner(un-u_exact,un-u_exact)*dx))
    print(assemble(inner(pn-p_exact,pn-p_exact)*dx))
    ufile << un
    pfile << pn
    t += dt
    step += 1
