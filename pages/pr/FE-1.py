import matplotlib.pyplot as plt
import numpy as np 
from fenics import *

m = 20
p = 1
c = 100

kmax = 10
kk = 1

mesh = UnitSquareMesh(m, m)
xm = mesh.coordinates()
ym = np.zeros((m+1), "float") 

V = FunctionSpace(mesh, "CG", p)
n = V.dim()-1

u = TrialFunction(V)
v = TestFunction(V)

cc = Constant([c,c])
a = dot(grad(u), grad(v))*dx + dot(cc,grad(u))*v*dx
b = u*v*dx

# Assemble stiffness form
A = PETScMatrix()
assemble(a, tensor=A)
B = PETScMatrix()
assemble(b, tensor=B)

# Create eigensolver
eigensolver = SLEPcEigenSolver(A,B)
eigensolver.parameters["spectrum"] = "smallest magnitude"

eigensolver.solve(kmax)

dataList = []
for k in range(0, kmax): 
    r, c, rx, cx = eigensolver.get_eigenpair(k)
    
    dic = {"Номер собственного значения": k+1, "Действительная часть": r, 
           "Мнимая часть": c}
    dataList.append(dic) 
    print(dic)      

    
ur = Function(V)
ur.vector()[:] = rx
ui = Function(V)
ui.vector()[:] = cx

N = 100
x = np.linspace(0,1.,N)
y = np.linspace(0,1.,N)
yy  = np.zeros((N,N)) 

for i in range(0, N): 
    for j in range(0, N): 
        pp = Point(x[i],y[j])
        yy[i,j] = ur(pp)
          
fig1 = plt.figure(1)
plt.contourf(x,y,yy) 
plt.gca().set_aspect("equal")
plt.colorbar()

for i in range(0, N): 
    for j in range(0, N): 
        pp = Point(x[i],y[j])
        yy[i,j] = ui(pp)
          
fig2 = plt.figure(2)
plt.contourf(x,y,yy) 
plt.gca().set_aspect("equal")
plt.colorbar()


plt.show()
