from re import L
from numpy import linspace, zeros, log10, round_, array, float64
from numpy.linalg import norm
from sklearn.linear_model import LinearRegression
from resources.Cauchy_Problem import Cauchy
from resources.Time_Schemes import Euler, RK4, CN, BW_Euler, LeapFrog
from cmath import  pi, sin, cos
from mpmath import findroot


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
    T = t[Nt-1]
    NU = len(U0) 

    E = zeros([Nt, NU])
    log_diff_U21 = zeros(k)
    log_Nt = zeros(k)
    lin_log_diff_U21 = zeros(k)
    lin_log_Nt = zeros(k)
    order = 0

    U_t = Cauchy(F, Time_Scheme, U0, t)
    

    for i in range (k):

        Nt = 2*Nt
        t_2 =linspace(0, T, Nt)
        U_t_2 =  Cauchy(F, Time_Scheme, U0, t_2)
        

        log_diff_U21[i] = log10(norm(U_t_2[int(Nt - 1), :] - U_t[int(0.5*Nt - 1), :]))
        log_Nt[i] = log10(Nt)

        U_t = U_t_2
        
        for j in range(k):

         if (abs(log_diff_U21[j]) > 12):

             break
        
    l=min(j, k)


    reg= LinearRegression().fit(log_Nt[0:j+1].reshape((-1, 1)),log_diff_U21[0:j+1]) 
    order = round_(abs(reg.coef_),1)

    lin_log_Nt = log_Nt[0:j+1]
    lin_log_diff_U21 = reg.predict(log_Nt[0:j+1].reshape((-1, 1)))
        


    return [log_diff_U21, log_Nt, lin_log_diff_U21, lin_log_Nt, order]


def Stability_Region(t_scheme): 

    N = 100
    x, y = linspace(-5, 5, N), linspace(-5, 5, N)
    rho =  zeros( (N, N),  dtype = float64)

    for i in range(N): 
      for j in range(N):

          w = complex(x[i], y[j])
          r = t_scheme( 1, lambda u, t: w*u,  0, 1 )
          rho[i, j] = abs(r) 

    return rho, x, y