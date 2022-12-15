from re import U
from resources.Math_Operators import Newton_Raphson
from scipy.optimize import fsolve, newton
from numpy import array, float64, zeros, matmul
from numpy.linalg import norm


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

def ERK(U, F, t, dt):

    ERK.__name__ == "Embedded Runge-Kutta"
    RK_Method = "Dormand-Prince"
    tol = 1e-6

    V1 = RK_Scheme(RK_Method, "First", U, t, dt, F) 
    V2 = RK_Scheme(RK_Method, "Second", U, t, dt, F) 

    a, b, bs, c, q, Ne = Butcher_Tableau(RK_Method)

    h = min(dt, Step_Size(V1-V2, tol, dt,  min(q)))

    N_n = int(dt/h)+1
    n_dt = dt/N_n

    V1 = U
    V2 = U

    for i in range(N_n):
        time = t + i*dt/int(N_n)
        V1 = V2
        V2 = RK_Scheme(RK_Method, "First", V1, time, n_dt, F)

    U2 = V2

    ierr = 0

    return U2

def RK_Scheme(name, tag, U1, t, dt, F):
    a, b, bs, c, q, N = Butcher_Tableau(name)
    k = zeros([N, len(U1)])

    k[0,:] = F(U1, t + c[0]*dt)

    if tag=="First":
        
        for i in range(1,N):
            Up = U1
            for j in range(i):
                Up = Up + dt*a[i,j]*k[j,:]

            k[i,:] = F(Up, t + c[i]*dt)

        U2 = U1 + dt*matmul(b,k)

    elif tag == "Second":

        for i in range(1,N):
            Up = U1
            for j in range(i):
                Up = Up + dt*a[i,j]*k[j,:]

            k[i,:] = F(Up, t + c[i]*dt)

        U2 = U1 + dt*matmul(bs,k)

    return U2

def Step_Size(dU, tolerance, q, dt):

    normT = norm(dU)

    if normT > tolerance:
        step_size = dt*(tolerance/normT)**(1/(q+1))

    else:
        step_size = dt
    
    return step_size

def Butcher_Tableau(Name: str):
    """This function generates the Butcher Tableau depending on the method selected
    Args:
        Name (str): Name of the method, options are: Heun-Euler, Bogacki-Shampine, Fehlberg-RK12, Dormand-Prince, Cash-Karp and Felberg.
                    
    """
    if Name == "Heun-Euler":
        q = [2,1]        # Heun's order = 2, Eulers order = 1
        N = 2

        a = zeros([N,N-1])
        b = zeros([N])
        bst = zeros([N])
        c = zeros ([N])

        c = [0., 1.]

        a[0,:] = [  0. ]
        a[1,:] = [  1. , 0.]

        b[:]   = [0.5 , 0.5 ]
        bst[:] = [ 1.,    0.  ]

    elif  Name=="Bogacki-Shampine":
        q = [3,2]
        N = 4 

        a = zeros([N,N-1])
        b = zeros([N])
        bs = zeros([N])
        c = zeros([N])
      
        c[:] = [ 0., 1./2, 3./4, 1. ]

        a[0,:] = [  0., 0., 0.            ]
        a[1,:] = [ 1./2, 0., 0.           ]
        a[2,:] = [ 0.,	3./4, 0.    	]
        a[3,:] = [ 2./9,	1./3,	4./9 	]

        b[:]  = [ 2./9,	1./3,	4./9,	0. ]
        bs[:] = [ 7./24,	1./4,	1./3,	1./8 ]

    elif Name == 'Fehlberg-RK12':

        q = [2,1]
        N = 3 

        a = zeros([N,N-1])
        b = zeros([N])
        bst = zeros([N])
        c = zeros([N])
      
        c[:] = [ 0., 0.5, 1. ]

        a[0,:] = [  0., 0. ]
        a[1,:] = [  1./2, 0. ]
        a[2,:] = [  1./256,  255./256	]

        b[:]  = [ 1./256,	255./256,	0. ]
        bst[:] = [ 1./512,	255./256,	1./512 ]

    elif Name=="Dormand-Prince": 
        q = [5,4]
        N = 7 

        a = zeros([N,N-1])
        b = zeros([N])
        bst = zeros([N])
        c = zeros([N])
      
        c[:] = [ 0., 1./5, 3./10, 4./5, 8./9, 1., 1. ]

        a[0,:] = [          0.,           0.,           0.,         0.,           0.,     0. ]
        a[1,:] = [      1./5  ,           0.,           0.,         0.,           0.,     0. ]
        a[2,:]=  [      3./40 ,        9./40,           0.,         0.,           0.,     0. ]
        a[3,:] = [     44./45 ,      -56./15,        32./9,         0.,           0.,     0. ]
        a[4,:] = [ 19372./6561, -25360./2187,  64448./6561,  -212./729,           0.,     0. ]
        a[5,:] = [  9017./3168,    -355./33 ,  46732./5247,    49./176, -5103./18656,     0. ]
        a[6,:]=  [    35./384 ,           0.,    500./1113,   125./192, -2187./6784 , 11./84 ]

        b[:]  = [ 35./384   , 0.,   500./1113,  125./192,  -2187./6784  ,  11./84  ,     0.]
        bst[:] = [5179./57600, 0., 7571./16695,  393./640, -92097./339200, 187./2100, 1./40 ] 

    elif Name =="Cash-Karp":
        q = [5,4] 
        N = 6 

        a = zeros([N,N-1])
        b = zeros([N])
        bst = zeros([N])
        c = zeros([N])

        c[:] = [ 0., 1./5, 3./10, 3./5, 1., 7./8 ]
      
        a[1,:] = [ 0.,          0.,       0.,         0.,            0. ] 
        a[2,:] = [ 1./5,        0.,       0.,         0.,            0. ] 
        a[3,:] = [ 3./40,       9./40,    0.,         0.,            0. ] 
        a[4,:] = [ 3./10,      -9./10,    6./5,       0.,            0. ] 
        a[5,:] = [ -11./54,     5./2,    -70./27,     35./27,        0. ] 
        a[6,:] = [ 1631./55296, 175./512, 575./13824, 44275./110592, 253./4096 ] 

        b[:]  = [    37./378, 0.,     250./621,     125./594,         0., 512./1771]
        bst[:] = [2825./27648, 0., 18575./48384, 13525./55296, 277./14336,     1./4 ]

    elif Name == "Felberg":

        q = [5,4] 
        N = 6 

        a = zeros([N,N-1])
        b = zeros([N])
        bst = zeros([N])
        c = zeros([N])

        c[:] = [ 0., 1./4, 3./8, 12./13, 1., 0.5 ]
      
        a[1,:] = [ 0.,          0.,       0.,         0.,            0. ] 
        a[2,:] = [ 1./4,        0.,       0.,         0.,            0. ] 
        a[3,:] = [ 3./32,       9./32,    0.,         0.,            0. ] 
        a[4,:] = [ 1932./2197, -7200./2197 , 7296./2197 ,  0.,       0. ] 
        a[5,:] = [ 439./216,  -8 , 3680./513,     -845./4104,        0. ] 
        a[6,:] = [ -8./27, 2. , -3544./2565 , 1859./4104, -11./40 ] 

        b[:]  = [ 16./135 , 0., 6656./12825, 28561./56430 , -9./50 , 2./55]
        bst[:] = [25./216 , 0., 1408./2565, 2197./4104 , -1./5,  0 ] 



    return a, b, bst, c, q, N


