from dolfin import *
import matplotlib.pyplot as plt

mesh = RectangleMesh(Point(0,0),Point(1,1),64,64)

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Set parameter values
dt = 0.01
T = 100
nu = 0.01

# Define boundary conditions
noslip  = DirichletBC(V, (0, 0), "near(x[0],1) || near(x[0],0) || near(x[1],0)")
upflow  = DirichletBC(V, (1, 0), "x[1] > 1.0 - DOLFIN_EPS")
pinpoint = DirichletBC(Q, 0, "x[0] < DOLFIN_EPS && x[1] < 0.02")

bcu = [noslip, upflow]
bcp = [pinpoint]

# Create functions
un = Function(V)
u_ = Function(V)
pn = Function(Q)


wn = Function(V)
wn.interpolate(Expression(("x[0]","x[1]"), degree = 1))

dofmap = V.dofmap()
element = V.element()

# Define coefficients
k = Constant(dt)

# Tentative velocity step
F1 = (1/k)*inner(u - un, v)*dx + inner(grad(un)*un, v)*dx + nu*inner(grad(u), grad(v))*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u_)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u_, v)*dx - k*inner(grad(pn), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

# Use nonzero guesses - essential for CG with non-symmetric BC
parameters['krylov_solver']['nonzero_initial_guess'] = True

# Create files for storing solution
ufile = File("results/velocity.pvd")
pfile = File("results/pressure.pvd")
wfile = File("results/displacement.pvd")
yfile = File("results/det.pvd")

# Time-stepping
t = dt
i = 0
while t < T + DOLFIN_EPS:
    # Tentative velocity
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u_.vector(), b1, "bicgstab", prec)
    # Pressure correction
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    solve(A2, pn.vector(), b2, "bicgstab", prec)
    # Velocity correction
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, un.vector(), b3, "bicgstab", prec)
     # update displacement
    for cell in cells(mesh):
        dof_coordinates = element.tabulate_dof_coordinates(cell)
        dof_index = dofmap.cell_dofs(cell.index())
        for i in range(6):
            position = wn(dof_coordinates[i])
            if position[0] < 1.0 and position[0]>0.0 and position[1] < 1.0 and position[1] > 0.0:
                u0 = un(position)
            else:
                u0 = [0,0]
            wn.vector()[dof_index[i]]   += dt*u0[0]
            wn.vector()[dof_index[i+6]] += dt*u0[1]
    # Save to file
    t += dt
    ufile << un
    pfile << pn
    wfile << wn
    yn = project(det(grad(wn)),Q)
    print(yn(0.5,0.5))
    yfile << yn






