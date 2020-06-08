from dolfin import *
import matplotlib.pyplot as plt

# Define mesh
mesh = RectangleMesh(Point(0, 0), Point(1, 1), 16, 16)

# Define function spaces W
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)

# Define trial and test functions
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

# Set parameter values
dt = 0.01
T = 100
nu = 0.01

# Define boundary conditions
noslip = DirichletBC(W.sub(0), (0, 0),"near(x[0],1) || near(x[0],0) || near(x[1],0)")
upflow = DirichletBC(W.sub(0), (1, 0), "near(x[1],1)")
pinpoint = DirichletBC(W.sub(1), 0, "near(x[0],0) && near(x[1],0)", "pointwise")
bcu = [noslip, upflow]
bcp = [pinpoint]


# Create functions
w = Function(W)
(un, pn) = w.split()

# Define variational problem
F1 = (1/dt)*inner((u-un), v)*dx + inner(grad(un)*un, v)*dx + nu*(inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Assemble matrices
A1 = assemble(a1)

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

# Use nonzero guesses - essential for CG with non-symmetric BC
parameters['krylov_solver']['nonzero_initial_guess'] = True

# Create files for storing solution
ufile = File("results/velocity.pvd")
pfile = File("results/pressure.pvd")
# Time-stepping
t = dt
w_ = Function(W)
while t < T + DOLFIN_EPS:
    # fluid velocity and pressure
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, w_.vector(), b1)
    (u_, p_) = w_.split(True)
    un.assign(u_)
    pn.assign(p_)
    # Save to file
    ufile << un
    pfile << pn
    t += dt
    print(t)

