from numpy import linspace, array
import matplotlib.pyplot as plt
from resources.Cauchy_Problem import Cauchy
from resources.Time_Schemes import Euler, CN, RK4, BW_Euler
from resources.Forcing_Sides import Kepler
from resources.Error_Evaluation import Error_Cauchy_Problem

N = 1000
T = 10
t = linspace(0, T, N)
U0 = array ([1, 0, 0, 1])

t_schemes = [ Euler, CN, RK4, BW_Euler ]
scheme_order = [1, 2, 4, 1]
t_scheme_name = [ "Euler", "CN" , "RK4", "BW_Euler" ]
t_scheme_fig = [ "Euler", "CN" , "RK4", "BW_Euler" ]
i = 0
for t_scheme in t_schemes:

   
    E = Error_Cauchy_Problem(t_scheme, scheme_order[i], Kepler, U0, t)
    plt.figure(figsize=(5, 5), layout='constrained')
    plt.plot(t, E[:, 0], label='X Error') 
    plt.plot(t, E[:, 1], label='Y Error')
    plt.plot(t, E[:, 2], label='U Error')
    plt.plot(t, E[:, 3], label='V Error')


    plt.xlabel('t')
    plt.ylabel('E')
    plt.title( t_scheme_name[i])  # Add a title to the axes.
    plt.legend()
    plt.show()
    plt.savefig(t_scheme_fig[i])
    i=i+1

