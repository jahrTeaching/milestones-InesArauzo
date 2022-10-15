from numpy import linspace, zeros
from numpy.linalg import norm
from resources.Cauchy_Problem import Cauchy

def Error_Cauchy_Problem (Time_Scheme, Scheme_Order, F, U0, t):

    # Computes the error of the obtained solution of a Cauchy problem
    # using the Richardson extrapolation:
    # For this, the program solves the cauchy Problem on two different temporal grids
    # (t1 & t2), with different dt but the same final T. 
   
    Nt = len(t)
    NU = len(U0)
    E = zeros([Nt, NU])

    t_2 =linspace(0, t[Nt - 1], 2*Nt) #I already had a grid (t).  This one (t_2) has twice the amount of elements (Nt_2 =2*Nt), thus dt_2 = dt/2

    U_t = Cauchy(F, Time_Scheme, U0, t)
    U_t_2 =  Cauchy(F, Time_Scheme, U0, t_2)

    for i in range (Nt):
        E[i,:] = (U_t_2[2*i, :] - U_t[i, :]) / (1 - 1 / (2**Scheme_Order))
   
    return E


def Convergence_Rate(Time_Scheme, F, U0, t):
    #determines the error of the numerical solution as a function of the number of time
    #steps N. This subroutine internally integrates a sequence of refined dt_i and, by
    #means of the Richardson extrapolation, determines the error.
    
    k = 10 # Numero de ptos que cojo en la gráfica (N1, N2=2*N1, N3=2*N2...) (numero de cauchy problems que te tienes que hacer)
    Nt = len(t)
    NU = len(U0) 
    E = zeros([Nt, NU])
    log_E = zeros(k)
    log_N = zeros(k)
    lin_log_E = zeros(k)
    lin_log_N = zeros(k)
    U_t = Cauchy(F, Time_Scheme, U0, t)
    order = 0

    for i in range (k):
        t_2 =linspace(0, t[Nt - 1], 2*Nt)
        U_t_2 =  Cauchy(F, Time_Scheme, U0, t_2)
        # Tiene que haber una forma de hacerlo llamando al Cauchy problem, pero hoy ya no puc mes





    



    return [order, log_E, log_N, lin_log_E, lin_log_N ]


