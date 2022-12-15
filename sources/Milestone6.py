# from resources.Math_Operators import Numeric_Jacobian #, Newton_Raphson
from resources.Time_Schemes import ERK, Euler, BW_Euler, RK4, LeapFrog, CN
from resources.Cauchy_Problem import Cauchy
from resources.Forcing_Sides import R3BodyProblem, R3BodyProblem_Autonomous, Lagrange_Points_Calculation, Lagrange_Points_Stability


from numpy import zeros, array, linspace, around #, size
from numpy.random import rand
# from scipy.optimize import newton

from matplotlib import rc # LaTeX tipography
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

import matplotlib
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)


T = 1000                                # Integration duration [s]
n = int(1e6)                            # Number of Points
t = linspace(0,T,n)                     # Time array



mu = 3.0039e-7  # Earth - Sun gravitational constant


#%% Orbits around the lagrange points

NL = 5     # Number of Lagrange Points

U0 = zeros([NL,4])  # Assigning the initial values for the system resolution

U0[0,:] = array([0.8, 0.6, 0, 0])
U0[1,:] = array([0.8, -0.6, 0, 0])
U0[2,:] = array([-0.1, 0, 0, 0])
U0[3,:] = array([0.1, 0, 0, 0])
U0[4,:] = array([1.01, 0, 0, 0])

LagrangePoints = Lagrange_Points_Calculation(U0, NL)


U0_LP = zeros(4)
U0_stabLP = zeros(4)
eps = 1e-3*rand()
Lagrange_Points_List = array([1,2,3,4,5])
for k in range(5):

    selectedLP = k + 1

    if selectedLP == 5:
        label = 'L2'
    elif selectedLP == 4:
        label = 'L1'
    elif selectedLP == 3:
        label = 'L3'
    elif selectedLP == 2:
        label = 'L5'
    elif selectedLP == 1:
        label = 'L4'
    
    U0_LP[0:2] = LagrangePoints[selectedLP-1,:] + eps
    U0_LP[2:4] = eps

    U0_stabLP[0:2] = LagrangePoints[selectedLP-1,:]
    U0_stabLP[2:4] = 0

    eig = Lagrange_Points_Stability(U0_stabLP)
    print(around(eig.real,8))

    t_schemes = [Euler,RK4, CN, BW_Euler, LeapFrog, ERK]

    for j in range (len(t_schemes)):

        U_LP = Cauchy(R3BodyProblem, t_schemes[j], U0_LP, t )

        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.plot(U_LP[:,0], U_LP[:,1],'-',color = "r")
        ax1.plot(-mu, 0, 'o', color = "g")
        ax1.plot(1-mu, 0, 'o', color = "b")
        for i in range(NL):
            ax1.plot(LagrangePoints[i,0], LagrangePoints[i,1] , 'o', color = "k")

        ax2.plot(U_LP[:,0], U_LP[:,1],'-',color = "r")
        ax2.plot(LagrangePoints[selectedLP - 1,0], LagrangePoints[selectedLP - 1,1] , 'o', color = "k")

        ax1.set_xlim(-2,2)
        ax1.set_ylim(-2,2)
        ax1.set_title("Complete orbital system")
        ax2.set_title("Lagrange point")
        ax2.set_xlim(LagrangePoints[selectedLP - 1,0]-0.02,LagrangePoints[selectedLP - 1,0]+0.02)
        ax2.set_ylim(LagrangePoints[selectedLP - 1,1]-0.02,LagrangePoints[selectedLP - 1,1]+0.02)
        fig.suptitle(f"Earth-Sun Restricted 3 Body Problem ({t_schemes[j].__name__}) - orbit around {label} with t = {T}s" )
        for ax in fig.get_axes():
            ax.set(xlabel='x', ylabel='y')
            ax.grid()
            
        manager = plt.get_current_fig_manager()
        manager.full_screen_toggle()
        figure = plt.gcf()                  
        figure.set_size_inches(16, 8)       
        plt.savefig('Figures/ R3BP ' + label +' '+ t_schemes[j].__name__ +'.png', bbox_inches = 'tight')
        plt.close('all')
        plt.show()