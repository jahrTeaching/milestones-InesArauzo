from numpy import array, linspace
import matplotlib.pyplot as  plt
from resources.Cauchy_Problem import Cauchy
from resources.Time_Schemes import Euler, CN, RK4, BW_Euler
from resources.Forcing_Sides import Kepler



def Cauchy_wrapper(T, N, t_scheme):
    t = linspace(0, T, N)
    U0 = array ([1, 0, 0, 1])
    U = Cauchy(Kepler, t_scheme, U0, t)
    return U 



t_schemes = [  CN, BW_Euler ]
t_scheme_name = [ "CN" , "BW_Euler" ]
t_scheme_fig = ["CN.png", "BW_Euler.png" ]
i = 0
for t_scheme in t_schemes:

    U_dta = Cauchy_wrapper(10, 100, t_scheme)
    U_dtb = Cauchy_wrapper(10, 10000, t_scheme)
    U_dtc = Cauchy_wrapper(10, 50, t_scheme)
    plt.figure(figsize=(5, 5), layout='constrained')
    plt.plot(U_dta[:, 0], U_dta[:, 1], label='dt=0.1') 
    plt.plot(U_dtb[:, 0], U_dtb[:, 1], label='dt=0.001') 
    plt.plot(U_dtc[:, 0], U_dtc[:, 1], label='dt=0.2')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(t_scheme_name[i])  # Add a title to the axes.
    plt.legend()
    plt.show()
    plt.savefig(t_scheme_fig[i])
    i=i+1

