from Math_Operators import Newton_Raphson
from scipy.optimize import fsolve
from numpy import float64
def Euler(U, F, t, dt):
    
    return U + dt * F (U, t)

def CN( Un, F,  t, dt):

    def func(Un1):
        
        return Un1 - Un - dt/2 * ( F (Un1, t + dt) + F (Un, t) )

    U = Newton_Raphson(func, Un)

    return U

def RK4(U, F, t, dt):

    k1=F(U, t)

    k2=F(U + dt * k1/2,t + dt/2)

    k3=F(U + dt * k2/2, t + dt/2)

    k4=F(U + dt * k3,t + dt)

    return U + dt * (k1 + 2 * k2 + 2 * k3 + k4)/6

def BW_Euler(Un, F, t, dt):

    def BW_Euler_eqn(Un1):

        return Un1 - Un - dt * F(Un, t)
    U = Newton_Raphson ( BW_Euler_eqn, Un )

    return U