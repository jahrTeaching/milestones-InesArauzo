from numpy import linspace, array, round
import matplotlib.pyplot as plt
from resources.Cauchy_Problem import Cauchy
from resources.Time_Schemes import Euler, CN, RK4, BW_Euler
from resources.Forcing_Sides import Kepler
from resources.Error_Evaluation import Error_Cauchy_Problem, Convergence_Rate

from matplotlib import rc # LaTeX tipography
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rc('text', usetex=True); plt.rc('font', family='serif')

import matplotlib
matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18)

Nt = 200
T = 40
t = linspace(0, T, Nt)
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
    plt.savefig('Error_' + t_scheme_fig[i] + '.png')

    [log_diff_U21, log_Nt, lin_log_diff_U21, lin_log_Nt, order] = Convergence_Rate(t_scheme, Kepler, U0, t)
    
    
    plt.plot(log_Nt, log_diff_U21, label = t_scheme_fig[i])
    plt.plot(lin_log_Nt, lin_log_diff_U21, label = 'Linear regression')
    plt.xlabel("log(Nt)")
    plt.ylabel("log(U2-U1)")
    plt.title(t_scheme_name[i] + " order = " +str(round(order)))
    plt.legend()
    plt.plot()
    plt.grid()
    
    plt.savefig('Convergence_Rate' + t_scheme_name[i] + '.png')

    i=i+1
    