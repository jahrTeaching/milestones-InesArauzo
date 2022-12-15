from numpy import array, linspace, zeros, size, concatenate, reshape
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import rc # LaTeX tipography
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rc('text', usetex=True); plt.rc('font', family='serif')

import matplotlib
matplotlib.rc('xtick', labelsize = 12) 
matplotlib.rc('ytick', labelsize = 12)
from mpl_toolkits.mplot3d import Axes3D

from resources.Time_Schemes import (RK4, )
from resources.Time_Schemes import (Euler, CN, RK4, BW_Euler, LeapFrog)
from resources.Forcing_Sides import N_Body_Problem
from resources.Cauchy_Problem import Cauchy
from resources.Error_Evaluation import Stability_Region



T = 80                # [s]       
dt = array([0.01])    # [s]  
Nb = 4                # Number of bodies


# Define initial conditions: 
U0 = zeros(Nb*2*3)

if Nb == 2: 
    r01 = array([1, 1, 0])
    v01 = array([-0.5, 0, 0])

    r02 = array([-1, -1, 0])
    v02 = array([0.5, 0, 0])

    r = concatenate((r01,r02),axis=0)
    v = concatenate((v01,v02),axis=0)

elif Nb == 3: 

    r01 = array([2, 2, 0])
    v01 = array([0.5, 0, -0.1])

    r02 = array([-2, -2, 0])
    v02 = array([-0.5, 0, -0.1])

    r03 = array([0, 0, 0])
    v03 = array([0 , 0, 0])

    r = concatenate((r01, r02, r03), axis=0)
    v = concatenate((v01, v02, v03), axis=0)

elif Nb == 4:

    r01 = array([2, 2, 0])
    v01 = array([-0.4, 0, 0])

    r02 = array([-2, 2, 0])
    v02 = array([0, -0.4, 0])

    r03 = array([-2, -2, 0])
    v03 = array([0.4, 0, 0])

    r04 = array([2, -2, 0])
    v04 = array([0, 0.4, 0])

    r = concatenate((r01,r02,r03,r04),axis=0)
    v = concatenate((v01,v02,v03,v04),axis=0)


for i in range(len(r)):

    U0[2*i] = r[i]
    U0[2*i + 1] = v[i]

print(U0)

t_schemes = [ Euler, CN, RK4, BW_Euler ]
#t_schemes = [RK4]

for j in range (size(t_schemes)):

    for i in range(size(dt)):

        n = int(T/dt[i])                      # Temporal steps  
        t  = linspace(0, T, n)                # Time array

        U = Cauchy(N_Body_Problem, t_schemes[j], U0, t)

        
        plt.plot(U[:, 0], U[:,2], "b")
        plt.plot(U[:, 6], U[:,8], "r")

        if Nb == 3:
            plt.plot(U[:, 12], U[:, 14], "k.")

        elif Nb == 4:
            plt.plot(U[:, 12], U[:, 14], "k")
            plt.plot(U[:, 18], U[:, 20], "purple")

        plt.xlabel("X")
        plt.ylabel("Y", rotation = 0)
        plt.title(f'2D projection of {Nb} bodies')

        plt.grid()
        plt.show()
        plt.savefig('Figures/Hito 5_' + str(Nb) + t_schemes[j].__name__ + str(dt[i])+'2D.png')

        fig = plt.figure()
        ax1 = fig.add_subplot(111,projection='3d')
        ax1.plot_wireframe(U[:,0].reshape((-1, 1)), U[:,2].reshape((-1, 1)), U[:,4].reshape((-1, 1)), color= "red", label = 'Body A')
        ax1.plot_wireframe(U[:,6].reshape((-1, 1)), U[:,8].reshape((-1, 1)), U[:,10].reshape((-1, 1)), color= "blue", label = 'Body B')
        
        if Nb == 3:
            ax1.plot_wireframe(U[:,12].reshape((-1, 1)), U[:,14].reshape((-1, 1)), U[:,16].reshape((-1, 1)), color= "black", label = 'Body C')
        elif Nb == 4:
            ax1.plot_wireframe(U[:,12].reshape((-1, 1)), U[:,14].reshape((-1, 1)), U[:,16].reshape((-1, 1)), color= "black", label = 'Body C')
            ax1.plot_wireframe(U[:,18].reshape((-1, 1)), U[:,20].reshape((-1, 1)), U[:,22].reshape((-1, 1)), color= "purple", label = 'Body D')

        plt.title(f'{Nb} bodies integrated with {t_schemes[j].__name__} and dt = {dt[i]}')
        plt.xlabel("X")
        plt.ylabel("Y",rotation = 0)
        plt.grid()

        plt.legend()
        plt.show()
        plt.savefig('Figures/Hito 5_' + str(Nb) + t_schemes[j].__name__ + str(dt[i])+'3D.png')


print(U)


Us = reshape(U,(n,4,3,2))
r = reshape( Us[:,:,:,0],(n,4,3))
        
def animate(num, data, line):
   line.set_data(data[0:2,:num])
   line.set_3d_properties(data[2,:num])
   return line
x = zeros([n+1])
y = zeros([n+1])
z = zeros([n+1])
x2 = zeros([n+1])
y2 = zeros([n+1])
z2 = zeros([n+1])
x3 = zeros([n+1])
y3 = zeros([n+1])
z3 = zeros([n+1])
x4 = zeros([n+1])
y4 = zeros([n+1])
z4 = zeros([n+1])
for i in range(n):
    x[i] = r[i,0,0]
    y[i] = r[i,0,1]
    z[i] = r[i,0,2]
    x2[i] = r[i,1,0]
    y2[i] = r[i,1,1]
    z2[i] = r[i,1,2]
    x3[i] = r[i,2,0]
    y3[i] = r[i,2,1]
    z3[i] = r[i,2,2]
    x4[i] = r[i,3,0]
    y4[i] = r[i,3,1]
    z4[i] = r[i,3,2]
data = array([x, y, z])
data2 = array([x2, y2, z2])
data3 = array([x3, y3, z3])
data4 = array([x4, y4, z4])
