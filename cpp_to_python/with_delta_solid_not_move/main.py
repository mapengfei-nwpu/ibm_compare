from dolfin import *
from mshr import *

from IBInterpolation import *
from IBMesh import *

import matplotlib.pyplot as plt

# Define fluid_mesh
points = [Point(0, 0, 0), Point(1, 1, 0)]
seperations = [32, 32, 0]
regular_mesh = IBMesh(points, seperations)
fluid_mesh = regular_mesh.mesh()
solid_mesh = generate_mesh(Circle(Point(0.6, 0.5), 0.2), 10)

# Define function spaces W for fluid
P2 = VectorElement("Lagrange", fluid_mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", fluid_mesh.ufl_cell(), 1)
TH = P2 * P1
W = FunctionSpace(fluid_mesh, TH)

# Define function spaces W for solid
Vs = VectorFunctionSpace(solid_mesh, "Lagrange", 2)

# Define trial and test functions for fluid
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

# Define trial and test functions for solid
us = TrialFunction(Vs)
vs = TestFunction(Vs)

# Set parameter values
dt = 0.01
T = 10
nu = 0.01
nu_s = 0.1

# Define boundary conditions
noslip = DirichletBC(W.sub(0), (0, 0),"near(x[0],1) || near(x[0],0) || near(x[1],0)")
upflow = DirichletBC(W.sub(0), (1, 0), "near(x[1],1)")
pinpoint = DirichletBC(W.sub(1), 0, "near(x[0],0) && near(x[1],0)", "pointwise")
bcu = [noslip, upflow]
bcp = [pinpoint]

# Create functions for fluid
w = Function(W)
(un, pn) = w.split()
(f, _) = Function(W).split(True)

# Create functions for solid
velocity = Function(Vs)
disp = Function(Vs)
force = Function(Vs)
disp.interpolate(Expression(("x[0]","x[1]"),degree=2))

# Define interpolation object
IB = IBInterpolation(regular_mesh, solid_mesh, disp._cpp_object)

# Define variational problem for fluid
F1 = (1/dt)*inner((u-un), v)*dx + inner(grad(un)*un, v)*dx + nu*(inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Define variational problem for solid
F2 = nu_s*inner(grad(disp), grad(vs))*dx + inner(us, vs)*dx
a2 = lhs(F2)
L2 = rhs(F2)

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

# Use nonzero guesses - essential for CG with non-symmetric BC
parameters['krylov_solver']['nonzero_initial_guess'] = True

# Create files for storing solution
ufile = File("results/velocity.pvd")
pfile = File("results/pressure.pvd")
ffile = File("results/force.pvd")

# Time-stepping
t = dt
w_ = Function(W)
while t < T + DOLFIN_EPS:
    # step 1. calculate velocity and pressure
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, w_.vector(), b1)
    (u_, p_) = w_.split(True)
    # step 2. interpolate velocity from fluid to solid
    IB.fluid_to_solid(u_._cpp_object, velocity._cpp_object)
    # step 3. calculate disp for solid and update current gauss points and dof points
    disp.vector()[:] = velocity.vector()[:]*dt + disp.vector()[:]
    IB.evaluate_current_points(disp._cpp_object)
    # step 4. calculate body force.
    b2 = assemble(L2)
    solve(A2, force.vector(), b2)
    # step 5. interpolate force from solid to fluid
    IB.solid_to_fluid(f._cpp_object, force._cpp_object)
    # step 6. update variables and save to file.
    un.assign(u_)
    pn.assign(p_)
    ufile << un
    pfile << pn
    ffile << f
    t += dt
    print(t)
