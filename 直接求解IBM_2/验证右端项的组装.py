from fenics import *
n=32
mesh=RectangleMesh(Point(0,0),Point(1,1),n,n)
V=VectorFunctionSpace(mesh,"P",2)
u=TrialFunction(V)
v=TestFunction(V)
f1=as_matrix(((1,1),(1,1)))
L=inner(f1,grad(v))*dx
force=assemble(L)
u0=Function(V)
u0.vector().set_local(force)
File("a.pvd")<<u0