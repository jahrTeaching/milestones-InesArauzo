"""
Created on Mon Nov 14 17:03:33 2022
@authors: Rafael Rivero de Nicolás, Guillermo García del Río, Inés Arauzo Andrés
"""

from numpy import deg2rad, array, size, dot, zeros, min, max, logical_and, sqrt, hstack, sin, cos, arctan, linspace, geomspace, reshape
from numpy.linalg import norm
from time import process_time

from numba import njit

import matplotlib.pyplot as plt
from matplotlib import rc # LaTeX tipography
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rc('text', usetex=True); plt.rc('font', family='serif')

import matplotlib
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

@njit(parallel=True)
def TBW_matrix(alpha: float, beta:float):

    """Creation of the rotation matrix Body Wind.

    Args:
        alpha (float): Angle of atack.          [Deg]
        beta (float): Sideslip.                 [Deg]

    Returns:
        Tbw(array): Rotation matrix.
    """    
 
    alpha = deg2rad(alpha)
    beta = deg2rad(beta)

    Tbw = array([[cos(alpha)*cos(beta), -cos(alpha)*sin(beta), -sin(alpha)],
                    [sin(beta), cos(beta), 0],
                    [sin(alpha)*cos(beta), -sin(alpha)*sin(beta), cos(alpha)]])

    return Tbw

def TWB_matrix(alpha:float, beta:float):

    """Creation of the rotation matrix Wind Body.

    Args:
        alpha (float): Angle of atack.          [Deg]
        beta (float): Sideslip.                 [Deg]

    Returns:
        Tbw(array): Rotation matrix.
    """    


    alpha = deg2rad(alpha)
    beta = deg2rad(beta)

    Twb = array([[cos(alpha)*cos(beta), sin(beta), sin(alpha)*cos(beta)],
                [-cos(alpha)*sin(beta), cos(beta), -sin(alpha)*cos(beta)],
                [-sin(alpha), 0 , cos(alpha)]])   

    return Twb
    
@njit(parallel=True)
def rotation_WB(alpha: float, beta: float, Data: array): 
    """ Coordinates change from Wind frame to Body frame.

    Args:
        alpha (float): Angle of atack.          [Deg]
        beta (float): Sideslip.                 [Deg]
        Data (array): Data in the Wind frame.

    Returns:
        Rotated array: Data in the Body frame
    """ 
    
    Twb = TWB_matrix(alpha, beta)   

    return dot(Twb,Data)
    
@njit(parallel=True)
def rotation_BW(alpha: float, beta: float, Data: array): 
    """ Coordinates change between Body and Wind frame.

    Args:
        alpha (float): Angle of atack           [Deg]
        beta (float): Sideslip                [Deg]
        Data (array): Data in the Body frame

    Returns:
        Rotated array: Data in the Wind frame
    """    
    Tbw = TBW_matrix(alpha, beta)

    return dot(Tbw,Data)

def rotation_rev(phi: float, Data: array): 
    """ Array rotation a certain angle of revolution around the X axis.

    Args:
        phi (float): Angle of revolution    [Deg]
        Data (array): Array to be rotated 

    Returns:
        Rotated Array: Rotated array
    """    
    phi = -deg2rad(phi)
    
    Rot_x = array([[1 ,   0   ,    0    ],
                [0 , cos(phi), -sin(phi)],
                [0 , sin(phi),  cos(phi)]])

    return dot(Rot_x,Data)

###########################
#### Funciones serias ####
def Analytical_Derivative(x, f):

    return ( f(x+1E-8)-f(x) ) / (1E-8)
    
@njit(parallel=True)
def Theta_calculation(x, f):

    theta = arctan( Analytical_Derivative( x, f ) )

    theta[theta<0] = 0 # betta<0 implies Cp=0 based on Newton Theory

    return theta
    
@njit(parallel=True)
def Normal(x, f):
    """Given a function f evaluated on an x partition, returns the normal vector of every partition

    Args:
        x (array): linspace in the coordinate x
        f (array): result of evaluating a function in x

    Returns:
        Normal (array): Array of normal vectors  
    """    
    
    Normals = zeros([len(x)-1, 3])
    for i in range (len(x)-1):
        dx = x[i+1] - x[i]
        dy = 0
        dz = f[i+1] - f[i]
        Normals[i, :] = (-dz, dy, dx)/norm((-dz, dy, dx))

    return Normals

# @njit(parallel=True)
# def rotation_BW(alpha: float, beta: float, Data: array): 
#     """ Coordinates change between Body and Wind frame

#     Args:
#         alpha (float): Angle of atack           [Deg]
#         beta (float): Sideslip                  [Deg]
#         Data (array): Normals array in the Body frame

#     Returns:
#         Rotated array: Data in the Wind frame
#     """    
    
#     alpha = deg2rad(alpha)
#     beta = deg2rad(beta)

#     Tbw = array([[cos(alpha)*cos(beta), -cos(alpha)*sin(beta), -sin(alpha)],
#                 [sin(beta), cos(beta), 0],
#                 [sin(alpha)*cos(beta), -sin(alpha)*sin(beta), cos(alpha)]])

#     return dot(Tbw,Data)

