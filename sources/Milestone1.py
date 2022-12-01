import numpy as np
import matplotlib.pyplot as  plt
import scipy.optimize as scp

plt.style.use ('seaborn-colorblind')


dim=4

dt=0.01
N=1000



def F(U, t):
    return np.array([U[2], U[3], -U[0]/(U[0]**2+U[1]**2)**1.5, -U[1]/(U[0]**2+U[1]**2)**1.5 ] )
#METODO DE EULER

UEuler=np.array([1,0,0,1])

xEuler=np.zeros(N)
yEuler=np.zeros(N)

xEuler[0]=UEuler[0]
yEuler[0]=UEuler[1]


#METODO RUNGE KUTTA 4
URK4=np.array([1,0,0,1])

xRK4=np.zeros(N)
yRK4=np.zeros(N)

xRK4[0]=URK4[0]
yRK4[0]=URK4[1]

def RK4(U,t,dt):
    k1=F(U,t)

    k2=F(U+dt*k1/2,t+dt/2)

    k3=F(U+dt*k2/2,t+dt/2)

    k4=F(U+dt*k3,t+dt)

    return (k1+2*k2+2*k3+k4)/6


#METODO CRANK NICOLSON
UCN=np.array([1,0,0,1])

xCN=np.zeros(N)
yCN=np.zeros(N)

xCN[0]=UCN[0]
yCN[0]=UCN[1]


def CN( Un,t,dt):

    def func(Un1):
        return Un1-Un-dt/2*(F(Un1, t+dt)+F(Un, t))

    U=scp.fsolve(func, Un)

    return U
t=0.



for i in range(1, N):
    t=t+dt

    UEuler=UEuler+dt*F(UEuler, t)
    xEuler[i]=UEuler[0]    
    yEuler[i]=UEuler[1]

    URK4=URK4+RK4(URK4, t, dt)*dt
    xRK4[i]=URK4[0]
    yRK4[i]=URK4[1]

    UCN=CN(UCN, t, dt)
    xCN[i]=UCN[0]
    yCN[i]=UCN[1]


np.savetxt('xEuler0001.csv', xEuler, delimiter=',')
np.savetxt('yEule0001r.csv', yEuler, delimiter=',')
np.savetxt('xRK40001.csv', xRK4, delimiter=',')
np.savetxt('yRK40001.csv', yRK4, delimiter=',')
np.savetxt('xCN0001.csv', xCN, delimiter=',')
np.savetxt('yCN0001.csv', yCN, delimiter=',')
plt.plot(xEuler, yEuler, 'go')
plt.title('Orbit Euler')
plt.xlabel('xEuler')
plt.ylabel('yEuler')
plt.axis('equal')
plt.show()

plt.plot(xRK4, yRK4, 'go')
plt.title('Orbit RK4')
plt.xlabel('xRK4')
plt.ylabel('yRK4')
plt.axis('equal')
plt.show()

plt.plot(xCN, yCN, 'go')
plt.title('Orbit CN')
plt.xlabel('xCN')
plt.ylabel('yCN')
plt.axis('equal')
plt.show()


