from numpy import array, reshape, zeros, size, sqrt
from numpy.linalg import norm, eig
from scipy.optimize import fsolve, newton
from resources.Math_Operators import Jacobian


def Kepler(U, t):

    x = U[0]
    y = U[1]
    dxdt = U[2]
    dydt = U[3]
    denom = ( x**2 + y**2 )**1.5
    return array([dxdt, dydt, -x/denom, -y/denom])

def Harmonic_Oscilator(U, t):
    amplitude  = 1
    return array([U[1], -amplitude*U[0]])



def N_Body_Problem(U:array,t:array):  
    """The N-body problem is the problem of predicting the individual motions of a
       group of mass objects under the influence of their gravitational field
       The function F_NBody uses the state vector U as an input
        and returns its derivative F

    Args:
        U (array): state vector (position_i, velocity_i)
        t (array): time partition

    Returns:
        F(array): derivative of U (velocity_i, acceleration_i)
    """    

    Nb = 4                                  # Number of bodies 
    dim = 3                                 # Dimensions of the problem

    
    Us = reshape(U, (Nb, dim, 2))           # reshape U[Nb*dim*2, 1] into Us[Nb, variable(x,y,z), position(0)/velocity(1)]        
    F  = zeros(len(U))
    Fs = reshape(F, (Nb, dim, 2))           # reshape F[Nb*dim*2, 1] into Fs[Nb, variable(x,y,z), velocity(0)/acceleration(1)]      

    r  = reshape(Us[:, :, 0],(Nb, dim))     # saves the position of every body r[body_index, position]
    v  = reshape(Us[:, :, 1],(Nb, dim))     # saves the velocity of every body v[body_index, velocity]

    drdt = reshape(Fs[:, :, 0], (Nb, dim))
    dvdt = reshape(Fs[:, :, 1], (Nb, dim))

    dvdt[:, :] = 0                 

    for i in range(Nb):

        drdt[i, :] = v[i, :]

        for j in range(Nb):
            
            if j != i:                   # Only aplicable for different bodies

                d = r[j, :] - r[i, :]
                dvdt[i, :] = dvdt[i, :] + d[:]/(norm(d)**3)  

    return F

def R3BodyProblem(U:array, t:array):
    """Restrcted 3-body problem

    Args:
        U (array): state vector (position_i, velocity_i)
        t (array): time partition

    Returns:
        F(array): derivative of U (velocity_i, acceleration_i)
    """    
    mu =3.0039e-7# 0.0121505856

    x =  U[0]   
    y =  U[1]   
    vx = U[2]  
    vy = U[3]  

    d = sqrt( (x + mu)**2 + y**2 )
    r = sqrt( (x - 1 + mu)**2 + y**2 )


    ax = x + 2*vy - (1 - mu)*( x + mu )/d**3 - mu*(x - 1 + mu)/r**3
    ay = y - 2*vx - (1 - mu) * y/d**3 - mu * y/r**3

    return array( [ vx, vy, ax, ay ] )



