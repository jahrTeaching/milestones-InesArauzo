from stl import mesh 
from numpy import deg2rad, array, size, dot, zeros, inner
from numpy.linalg import norm
from math import cos, sin, acos, pi

from matplotlib import pyplot
from mpl_toolkits import mplot3d
import vtkplotlib as vpl



def rotation_WB(alpha, beta, Data):

    alpha = deg2rad(alpha)
    beta = deg2rad(beta)

    Twb = array([[cos(alpha)*cos(beta), sin(beta), sin(alpha)*cos(beta)],
                [-cos(alpha)*sin(beta), cos(beta), -sin(alpha)*cos(beta)],
                [-sin(alpha), 0 , cos(alpha)]])   

    return dot(Twb,Data)

def rotation_BW(alpha, beta, Data):

    alpha = deg2rad(alpha)
    beta = deg2rad(beta)

    Tbw = array([[cos(alpha)*cos(beta), -cos(alpha)*sin(beta), -sin(alpha)],
                [sin(beta), cos(beta), 0],
                [sin(alpha)*cos(beta), -sin(alpha)*sin(beta), cos(alpha)]])

    return dot(Tbw,Data)

# Rotation of the STL data around an axis:


def CP_STL(stlfile: str ,alpha: float, beta: float , M: float, g: float):

    # Load STL file 

    Model = mesh.Mesh.from_file(stlfile)
    normals = Model.normals
    
    V_dir_W = array([1,0,0])
    V_dir_B = rotation_BW(alpha,beta,V_dir_W)
    #V_dir_B =  V_dir_W

    
    theta = zeros(size(normals,0))
    Cp = zeros(size(normals,0))
    Cp_mod = zeros(size(normals,0))

    Cp_max = 2/(g*M**2) * ( (((g+1)**2 * M**2)/(4*g*M**2-2*(g-1)))**(g/(g-1)) * (1-g+2*g*M**2)/(g+1) - 1 )

    for i in range(size(normals,0)):

        phi = acos(dot(V_dir_B,normals[i,:])/(norm(V_dir_B)*norm(normals[i,:])))

        if (pi/2 - phi) < 0:

            theta[i] = 0
        else: 

            theta[i] = pi/2 - phi

        Cp[i] = 2*(sin(theta[i]))**2 
        Cp_mod[i] = Cp_max*2*(sin(theta[i]))**2


    return Cp, Cp_mod

[Cp, Cp_mod] = CP_STL(stlfile='pikachu.stl', alpha = 30, beta = 0, M = 10, g = 1.4)
"""
file= open('Apollo.csv', 'w')
file.write(Cp)
file.close()
"""

mesh = mesh.Mesh.from_file('pikachu.stl')

vpl.mesh_plot(mesh, tri_scalars=Cp)

vpl.show()