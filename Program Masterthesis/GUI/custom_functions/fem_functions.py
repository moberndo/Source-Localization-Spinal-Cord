""" #######################################################################
    Title: Solving the source localization problem of spinal cord
    potentials for electric dipoles in the cervical spinal cord
    Subtitle: Functions for FEM calculations
    Author: Markus E. Oberndorfer
####################################################################### """

# ----- IMPORT PACKAGES -----
from ngsolve import ElementId, VOL, TaskManager, GridFunction, CGSolver
from ngsolve import H1, CoefficientFunction, BilinearForm, LinearForm
from ngsolve import grad, dx, Preconditioner
from math import sin, cos

# ----- FUNCTIONS -----
def FEM_Params(mesh, conductivity_vals):
    """
    FEM_Params calculates the linear form, the bilinear form and the
    finite element space, which are necessary for the Finite Element
    Calculations.

    Input:
        - mesh: Is an object that stores the model geometry and its
        properties. [NGSolve mesh object]
        - conductivity_vals: Dictionary that stores the conductivity
        values for all materials in the model. [dict]

    Output:
        - a: Bilinear form for the respective FE problem. [NGSolve object]
        - c: Preconditioned bilinear form. [NGSolve object]
        - f: Linear form the the respective FE problem. [NGSolve object]
        - fes: Finite element space that is used in the FE calculations.
        [NGSolve object]
    """
    fes = H1(mesh, order = '1')     # define Finite Element Space
    u,v = fes.TnT()                 # define Trial and Test Space

    sigma_coeff = CoefficientFunction([conductivity_vals[mat] for mat
                                       in mesh.GetMaterials()])

    # create bounded X-elliptic bilinear form
    a = BilinearForm(fes)
    a += grad(v)*sigma_coeff*grad(u) * dx + 1e-8*sigma_coeff*u*v*dx
    c = Preconditioner(a, "bddc")

    # create bounded linear form
    f = LinearForm(fes)
    print('FEM ready!')

    return a, c, f, fes

def MakeDipole3D(f, x, y, z, max_mesh, alpha, beta, p):
    """
    MakeDipole3D adds the dipole to the linear form f.

    Input:
        - f: linear form of the FE problem. [NGSolve object]
        - x: x location of the dipole to be added. [float]
        - y: y location of the dipole to be added. [float]
        - z: z location of the dipole to be added. [float]
        - max_mesh: Is the maximum size of one mesh element. [float]
        - alpha: alpha is the orientation angle of the dipole.
        alpha is the angle in the transverse plane [float]
        - beta: beta is the orientation angle of the dipole. beta is the
        angle in the frontal plane. [float]
        - p: Magnitude of the dipole that is added. [float]

    Output:
        - No output
    """
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

def AddDipoles(dipoles, f, maxh):
    """
    Function the iteratively adds dipoles by calling the 
    MakeDipole3D() function.

    Input:
        - dipoles: list of all possible dipoles that should be tested. Each
        list entry contains the location, direction and magnitude of one
        dipole. [list; list]
        - f: Linear form the the respective FE problem. [NGSolve object]
        - maxh: Is the maximum size of one mesh element. [float]

    Output: 
        - No output.
    """
    for dip in dipoles: 
        x = dip[0]
        y = dip[1]
        z = dip[2]
        alpha = dip[3]
        beta = dip[4]
        p = dip[5]
        MakeDipole3D(f, x = x, y = y, z = z, max_mesh= maxh*1.005,
                     alpha= alpha, beta=beta, p = p)

def CalcGFU(dipoles,f,maxh_, a, c, fes):
    """
    Calculate the GridFunction that contains the results of the inverse
    problem for the given model.
    
    Input: 
        - dipoles: List of all possible dipoles that should be tested. Each
        list entry contains the location, direction and magnitude of one
        dipole. [list; list]
        - f: Linear form the the respective FE problem. [NGSolve object]
        - maxh_: Is the maximum size of one mesh element. [float]
        - a: Bilinear form for the respective FE problem. [NGSolve object]
        - c: Preconditioned bilinear form. [NGSolve object]
        - fes: Finite element space that is used in the FE calculations.
        [NGSolve object]

    Output:
        - gfu: GridFunction that stores the results of the inverse problem.
        [NGSolve GridFunction]
        - inv: Inverse of the bilinear form. [NGSolve object]
        - f_vec: Contains the linear form of the FE problem in array form.
        [NGSolve object]
    """
    # add dipoles to linear form
    AddDipoles(dipoles, f, maxh_)

    # Assemble equations
    with TaskManager():
        a.Assemble()

    # Solve PDE
    gfu = GridFunction(fes)
    inv = CGSolver(a.mat, c.mat, printrates=False, precision=1e-8)
    gfu.vec.data = inv * f.vec

    return gfu, inv, f.vec
