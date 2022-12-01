from re import U
from resources.Math_Operators import Newton_Raphson
from scipy.optimize import fsolve, newton
from numpy import array, float64


def Euler(U, F, t, dt:float64):
    Euler.__name__ = "Euler"
    return array(U + dt * F(U, t))

def CN( Un, F,  t, dt):
    CN.__name__ = "Crank Nicolson"
    def func(Un1):
        
        return Un1 - Un - dt/2 * ( F (Un1, t + dt) + F (Un, t) )

    U = newton(func, Un)

    return U

def RK4(U, F, t, dt):
    RK4.__name__ = "Runge-Kutta 4"
    k1=F(U, t)

    k2=F(U + dt * k1/2,t + dt/2)

    k3=F(U + dt * k2/2, t + dt/2)

    k4=F(U + dt * k3,t + dt)

    return array(U + dt * (k1 + 2 * k2 + 2 * k3 + k4)/6)

def BW_Euler(Un, F, t, dt):
    BW_Euler.__name__ = "Euler Inverso"
    def BW_Euler_eqn(Un1):

        return Un1 - Un - dt * F(Un, t)
    U = newton ( BW_Euler_eqn, Un )

    return array(U)

def LeapFrog (U, F, t, dt):
     LeapFrog.__name__ = "Leapfrog"
     if t == 0:
         U = U + dt*F(U, t)
     else:
        U_aux = U
        aux = F (U_aux, t)

        L_2 = int(len(U)/2)

        U_aux [L_2 :] = U_aux [L_2 :] + aux[L_2 :]*dt/2
        U_aux [: L_2] = U_aux [: L_2] + aux[: L_2]*dt

        aux = F (U_aux, t)
        U_aux [L_2 :] = U_aux [L_2 :] + aux[L_2 :]*dt/2
        U = U_aux

        
        
    
     return U


