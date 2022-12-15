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

# %%

#### Funciones para ploteo ####
def Plot_geometry_1D(x, y, title=r'2D Geometry', axis_labels=[r'$x$ [-]', r'$z(x)$ [-]']):

    fig, ax = plt.subplots(1,1, figsize=(9,6), constrained_layout='true')
    ax.set_xlim( x[0] , x[-1] )
    ax.set_ylim( -y[0] - x[-1] + x[0], y[0] + x[-1] - x[0] )
    ax.set_title(title, fontsize=20)
    ax.grid()
    ax.set_xlabel(axis_labels[0] ,fontsize=20)
    ax.set_ylabel(axis_labels[1] ,fontsize=20)
    ax.plot( x, y, c="blue" )
    ax.plot( x, -y, c="blue" )
    plt.plot()

def Simple_plot_1D(x, y, title=None, axis_labels=[None, None]):

    fig, ax = plt.subplots(1,1, figsize=(9,6), constrained_layout='true')
    ax.set_xlim( x[0] , x[-1] )
    ax.set_ylim( min(y)*0.95, max(y)*1.05 )
    ax.set_title(title, fontsize=20)
    ax.grid()
    ax.set_xlabel(axis_labels[0] ,fontsize=20)
    ax.set_ylabel(axis_labels[1] ,fontsize=20)
    ax.plot( x, y, c="blue" )
    plt.plot()
   
################################

#### Funciones para definir la geometría ####
def Horizontal_Parabola(x): # Libro Anderson, Cap 3, page 63

    return sqrt( (x+1)/0.769 )

def Quarter_of_a_Circle(x):

    return sqrt( x )

def Flechita(x): # Se la ha inventado rafa

    f1 = x[x<=0.5]
    f2 = 2*x[logical_and(x>0.5, x<=0.75)]-0.5
    f3 = 1-4*(x[x>0.75]-0.75)

    return hstack((f1, f2, f3))
###########################

################################################
#####          Rotation Functions          #####
################################################ 
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
    
    Normals = zeros(len(x), 3)
    for i in range (len(x)):
        dx = x[i+1] - x[i]
        dy = 0
        dz = f[i+1] - f[i]
        Normals[i, :] = (-dz, dy, dx)/norm((-dz, dy, dx))

    return Normals

def Cp_Dist_1D_by_Newton(x, f, M):
    """Calculates the pressure coefficient over a 1D shape (defined by the evaluation of the function f over the parrtition x) at Mach M 

    Args:
        x (array): linspace in the coordinate x
        f (array): result of evaluating a function in x
        M (float): Mach number

    Returns:
        _type_: _description_
    """    

    Plot_geometry_1D(x, f(x))

    Cp = 2*(sin(Theta_calculation(x, f)))**2

    Simple_plot_1D(x, Cp)

    return Cp

def Cp_Dist_1D_by_Newton_Modified(x, f, M, gamma=1.4, Mach=10):

    Plot_geometry_1D(x, f(x), title=r'2D Geometry', axis_labels=[r'$x$ [-]', r'$z(x)$ [-]'])

    Cp_max = 2/(gamma*M**2) * ( (((gamma+1)**2 * M**2)/(4*gamma*M**2-2*(gamma-1)))**(gamma/(gamma-1)) * (1-gamma+2*gamma*M**2)/(gamma+1) - 1 ) # From page 62 Andeson Hypersonic...

    Cp_modified = Cp_max*( sin(Theta_calculation(x, f)) )**2

    Simple_plot_1D(x, Cp_modified)

    return Cp_modified
########################


###############################
@njit(parallel=True)
def Revoluting_Normals(Normals_array):
    """Given the normal vector field of a 1D shape in an array, returns the  normal vector field of the 3D revoluted shape

    Args:
        Normals_array (array): 1D normal array

    Returns:
        _type_: _description_
    """    

    phi_list = linspace(0,360, 100)

    Normals_revoluted = zeros( [len(N), len(phi_list), 3] )

    for i, n in enumerate(Normals_array):

        for j, phi in enumerate(phi_list):

            Normals_revoluted[i, j, :] = rotation_rev(phi,n)

    return Normals_revoluted

# Giro de las normales para un alpha o beta

def Rotated_Normals_Loops(alpha: float, beta: float, Normals_tensor: array):
    """ Rotation of all the normals a certain angle of atack and a certain sideslip.

    Args:
        alpha (float): Angle of atack.       [Deg]
        beta (float): Sideslip.              [Deg]
        Normals_tensor (array): Tensor in which each matrix represents all the normals along a section of the revolution surface 

    Returns:
        Normals_tensor_rotated: Tensor in which each matrix represents all the normals along a section of the revolution surface
                                taking into acount the angle of atack and the sideslip. 
    """    
    t0     = process_time()
    nx     = size(Normals_tensor,0)     # Number of partitions in X = Number of matrixes in the tensor
    nphi   = size(Normals_tensor,1)     # Number of partitions in the revolution angle = Number of files in the matrix
    ncomps = size(Normals_tensor,2)     # Number of components (3 = i j k)
    
    Normals_tensor_rotated = zeros([nx, nphi, ncomps])
    
    for j in range(nx):
        for i in range(nphi):
            Normals_tensor_rotated[j,i,:] = rotation_BW(alpha,beta,Normals_tensor[j,i,:])
    
    t1 = process_time()
    print('Elapsed processtime=', t1-t0)

    return Normals_tensor_rotated




def Rotated_Normals_Giant(alpha: float, beta: float, Normals_tensor: array):
    """_summary_

    Args:
        alpha (float): Angle of atack.       [Deg]
        beta (float): Sideslip.              [Deg]
        Normals_tensor (array): Tensor in which each matrix represents all the normals along a section of the revolution surface 

    Returns:
        Output (array): Tensor in which each matrix represents all the normals 
    """    
    

    t0 = process_time()
    nx     = size(Normals_tensor,0)     # Number of partitions in X = Number of matrixes in the tensor
    nphi   = size(Normals_tensor,1)     # Number of partitions in the revolution angle = Number of files in the matrix
    ncomps = size(Normals_tensor,2)     # Number of components (3 = i j k)

    Output = zeros( [nx, nphi, ncomps] )
    
    Tbw = TBW_matrix(alpha, beta)
    # Tbw_giant = zeros([3*nphi, 3*nphi])

    # Giant matrix building
    Tbw_giant = tbw_giant_loop (zeros([3*nphi, 3*nphi]), Tbw, nphi)
    # for i in range( nphi ):

    #     Tbw_giant[3*i:3*i+3, 3*i:3*i+3] = Tbw
    
    for i in range( size(Normals_tensor,0) ):

        p_O = reshape(Output[i,:,:], [3*nphi])          # Pointer to the Output tensor. It will be a column vector
        p_T = reshape(Normals_tensor[i,:,:], [3*nphi])  # Pointer to the Input tensor. It will be a column  vector
        p_O[:] = p_O[:] + dot(Tbw_giant, p_T)           # Vector corresponding to all normals rotated for a x value

    t1 = process_time()
    print('Elapsed processtime=', t1-t0)

    return Output

###################################
@njit(parallel=True)
def tbw_giant_loop (A, Tbw, nphi):
    """_summary_
hola warra
    Args:
        A (array): tensor full of zeros where 
        Tbw (_type_): _description_
        nphi (_type_): _description_

    Returns:
        _type_: _description_
    """    
    for i in range( nphi ):

        A [3*i:3*i+3, 3*i:3*i+3] = Tbw
    return A

    






##########################################################

M_inf = 8; g = 1.4
x = linspace(-1, 1, 70)
x_1 = linspace(0, 1, 50)

# Cp = Cp_Dist_1D_by_Newton(x, Horizontal_Parabola, M_inf)

# Cp = Cp_Dist_1D_by_Newton(x_1, Quarter_of_a_Circle, M_inf)

Cpmod = Cp_Dist_1D_by_Newton_Modified(x_1, Flechita, M_inf)

plt.show()



#%%
from numpy import array, zeros, size, deg2rad, cos, sin, dot
@njit(parallel=True)
def rotation_BW(alpha: float, beta: float, Data: array): 
    """ Coordinates change between Body and Wind frame

    Args:
        alpha (float): Angle of atack           [Deg]
        beta (float): Sideslip                  [Deg]
        Data (array): Normals array in the Body frame

    Returns:
        Rotated array: Data in the Wind frame
    """    
    
    alpha = deg2rad(alpha)
    beta = deg2rad(beta)

    Tbw = array([[cos(alpha)*cos(beta), -cos(alpha)*sin(beta), -sin(alpha)],
                [sin(beta), cos(beta), 0],
                [sin(alpha)*cos(beta), -sin(alpha)*sin(beta), cos(alpha)]])

    return dot(Tbw,Data)

#t0 = process_time()
T = zeros([2,3,3])         ## El 2 son las particiones en X , el primer 3 las particiones en Phi y el último 3
                           ## es el número de componentes del vector normal (i j k)

A = array([[0.1, 0.6, -0.1],
           [0.2, 0.3, -0.3],
           [0.3, -0.2, -0.1]])

C = array([[4,0,0],
           [5,0,0],
           [6,0,0]])

T[0] = A
T[1] = C

print(T[0])

for j in range(size(T,0)):
    for i in range(size(T,1)):
        T[j,i,:] = rotation_BW(-15,0,T[j,i,:])

#t1 = process_time()
#print('Elapsed processtime=', t1-t0)

print('LA MATRIZ ROTADA ESTA BABUINO')
print(T[0])
print(T[1])

def Rotated_Normals(alpha: float, beta: float, Normals_tensor: array):
    """ Rotation of all the normals a certain angle of atack and a certain sideslip.

    Args:
        alpha (float): Angle of atack.       [Deg]
        beta (float): Sideslip.              [Deg]
        Normals_tensor (array): Tensor in which each matrix represents all the normals along a section of the revolution surface 

    Returns:
        Normals_tensor_rotated: Tensor in which each matrix represents all the normals along a section of the revolution surface
                                taking into acount the angle of atack and the sideslip. 
    """    
   
    nx     = size(Normals_tensor,0)     # Number of partitions in X = Number of matrixes in the tensor
    nphi   = size(Normals_tensor,1)     # Number of partitions in the revolution angle = Number of files in the matrix
    ncomps = size(Normals_tensor,2)     # Number of components (3 = i j k)
    
    Normals_tensor_rotated = zeros([nx, nphi, ncomps])
    
    for j in range(nx):
        for i in range(nphi):
            Normals_tensor_rotated[j,i,:] = rotation_BW(alpha,beta,Normals_tensor[j,i,:])
    

    return Normals_tensor_rotated

def Rotated_Normals_Giant(alpha: float, beta: float, Normals_tensor: array):

    t0 = process_time()
    nx     = size(Normals_tensor,0)     # Number of partitions in X = Number of matrixes in the tensor
    nphi   = size(Normals_tensor,1)     # Number of partitions in the revolution angle = Number of files in the matrix
    ncomps = size(Normals_tensor,2)     # Number of components (3 = i j k)

    Output = zeros( [nx, nphi, ncomps] )
    
    Tbw = TBW_matrix(alpha, beta)
    Twb_giant = zeros([3*nphi, 3*nphi])

    # Giant matrix building
    for i in range( nphi ):

        Twb_giant[3*i:3*i+3, 3*i:3*i+3] = TBW_matrix(alpha, beta)
    
    for i in range( size(Normals_tensor,0) ):

        p_O = reshape(Output[i,:,:], [3*nphi])         # Pointer to the Output tensor. It will be a column vector
        p_T = reshape(Normals_tensor[i,:,:], [3*nphi])  # Pointer to the Input tensor. It will be a column  vector
        p_O[:] = p_O[:] + dot(Twb_giant, p_T)          # Vector corresponding to all normals rotated for a x value

    t1 = process_time()
    print('Elapsed processtime=', t1-t0)

    return Output

A = array([[0.1, 0.6, -0.1],
           [0.2, 0.3, -0.3],
           [0.3, -0.2, -0.1]])

C = array([[4,0,0],
           [5,0,0],
           [6,0,0]])

T[0] = A
T[1] = C


X = Rotated_Normals(-15,0,T)
print('Tensor X')
print(X[0])
print(X[1])


# %% Test 1


nx = 500
nphi = 500

N = array( [linspace(0,1, nx), zeros(nx), geomspace(0.1, 1, nx) ] ).transpose() # Matriz de normales

phi_list = linspace(0,360, nphi)

Normals_revoluted = zeros( [len(N), len(phi_list), 3] )

for i, n in enumerate(N):

    for j, phi in enumerate(phi_list):

        Normals_revoluted[i, j, :] = rotation_rev(phi,n)


X = Rotated_Normals_Loops(-15,0,Normals_revoluted)

X2, time2 = Rotated_Normals_Giant(-15,0,Normals_revoluted)

print('Tensor X')
# print(X[0])
# print(X[1])


print('\n\nTensor X2')
# print(X2[0])
# print(X2[1])



if X.any() != X2.any():
    print("UYYY")


