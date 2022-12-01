import matplotlib.pyplot as plt
from numpy import array, linspace, zeros, size, transpose, float64

from resources.Time_Schemes import (Euler, CN, RK4, BW_Euler, LeapFrog)
from resources.Forcing_Sides import Harmonic_Oscilator
from resources.Cauchy_Problem import Cauchy
from resources.Error_Evaluation import Stability_Region

from matplotlib import rc # LaTeX tipography
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rc('text', usetex=True); plt.rc('font', family='serif')

import matplotlib
matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18)


U0 = array([0,1])

scheme = [BW_Euler, LeapFrog]
coloring = ['r--','b--','g--']


for j in range (size(scheme)):

    T = 20                                              
    dt = array([0.1, 0.01, 0.001])                      
    
    
    for i in range(size(dt)):
        t  = linspace(0, T,int( T/dt[i]))       
        U = Cauchy(Harmonic_Oscilator, scheme[j], U0, t)
        print(U[len(t)-1,:])
        plt.title(f'Harmonic oscilator integrated with {scheme[j].__name__}')
        plt.xlabel("Tiempo (s)")
        plt.ylabel("X",rotation = 0)
        plt.grid()
        plt.plot(t,U[:,0], coloring[i], label = 'dt =' + str(dt[i]) + ' s')
        plt.legend()    
        
    # plt.savefig('Figures/'+ scheme[j].__name__ +'HarmonicOS.png')
    plt.show()
    plt.close()


    
    SR, X, Y = Stability_Region(scheme[j])

    aux = zeros(100)
    fig = plt.figure()
    ax = fig.add_subplot()
    plt.title(f'Región de estabilidad absoluta del {scheme[j].__name__} ')
    plt.plot(X, aux,'k-') 
    plt.plot(aux,Y,'k-')                         
    plt.grid()
    plt.contour(X, Y, transpose(SR), levels = [0, 1], linewidth = 2 )
    plt.contourf(X, Y, transpose(SR), levels = [0, 1], colors =['#626262'])            

    if scheme[j] != CN:
        ax.set_aspect('equal', adjustable = 'box')

    plt.xlabel("Re($\omega$)")
    plt.ylabel("Im($\omega$)",rotation = 0)
    # plt.savefig('Figures/'+ scheme[j].__name__ +'SR.png')
    plt.show()

    
 