""" ####################################################################### 
    Title: Construction of an FE model of the spinal cord to simulate the voltage distribution caused by a dipole
    Subtitle: Functions for FEM calculations
    Author: Markus E. Oberndorfer
####################################################################### """

# ----- IMPORTS -----
from ngsolve import ElementId, VOL, TaskManager, GridFunction, CGSolver
from ngsolve import H1, CoefficientFunction, BilinearForm, LinearForm, grad, dx, Preconditioner
from math import sin, cos

# ----- FUNCTIONS -----

def FEM_Params(mesh, conductivity_vals):
    fes = H1(mesh, order = '1')     # define Finite Element Space
    u,v = fes.TnT()                 # define Trial and Test Space

    sigma_coeff = CoefficientFunction([conductivity_vals[mat] for mat in mesh.GetMaterials()])

    # create bounded X-elliptic bilinear form
    a = BilinearForm(fes)
    a += grad(v)*sigma_coeff*grad(u) * dx + 1e-8*sigma_coeff*u*v*dx
    c = Preconditioner(a, "bddc")

    # create bounded linear form
    f = LinearForm(fes)
    print('FEM ready!')

    return a, c, f, fes

def MakeDipole3D(f, x, y, z, max_mesh, alpha, beta, p):
    spc = f.space
    mp1 = spc.mesh(x,y,z)
    ei1 = ElementId(VOL, mp1.nr)    
    fel1 = spc.GetFE(ei1)
    dnums1 = spc.GetDofNrs(ei1)
    shape1 = fel1.CalcShape(*mp1.pnt)

    for d1,s1 in zip(dnums1, shape1):
        f.vec[d1] += ( p.Get()*s1 )

    x1 = x + max_mesh * cos(beta) * cos(alpha)
    y1 = y + max_mesh * cos(beta) * sin(alpha)
    z1 = z + max_mesh * sin(beta)

    if x1 > 0.01:
        x1 = 0.01
        print('x1 too large, cap at 0.01.')
    if y1 > 0.01:
        print(y1)
        y1 = 0.01
        print('y1 too large, cap at 0.01.')

    mp2 = spc.mesh(x1, y1, z1)
    ei2 = ElementId(VOL, mp2.nr)
    fel2 = spc.GetFE(ei2)
    dnums2 = spc.GetDofNrs(ei2)
    shape2 = fel2.CalcShape(*mp2.pnt)

    for d2, s2 in zip(dnums2, shape2):
        f.vec[d2] += ( - p.Get()*s2 )

def CalcGFU(dipoles,f,maxh_, a, c, fes):
    for dip in dipoles: 
        x = dip[0]
        y = dip[1]
        z = dip[2]
        alpha = dip[3]
        beta = dip[4]
        p = dip[5]
        MakeDipole3D(f, x = x, y = y, z = z, max_mesh= maxh_*1.05, alpha= alpha, beta=beta, p = p)

    # Assemble equations
    with TaskManager():
        a.Assemble()

    # Solve PDE
    gfu = GridFunction(fes)
    inv = CGSolver(a.mat, c.mat, printrates=False, precision=1e-8)
    gfu.vec.data = inv * f.vec

    return gfu, inv, f.vec
