from numpy import zeros, matmul, array, size
from numpy.linalg import inv, norm
def Jacobian(F, U):

    nv = len(U)
    J = zeros((nv, nv))
    dx = 1e-3

    for i in range(nv):
        DX=array(zeros(nv))
        DX[i] = dx
        J[:, i] = (F(U + DX) - F(U))/(dx)
        
    return J



def Newton_Raphson(F, U0):

    nv = len(U0)
    DX = array(zeros(nv))
    U_1 = U0
    
    tol = 1e-8             # Tolerance
    i = 0                   # Iteration counter
    imax = 10000            # Max number of iterations
    err = 1

    while err > tol and i <= imax:
        J=Jacobian(F, U_1)
   
        DX = matmul( inv (J), F(U_1))
        U  = U_1 - DX
        err = norm (U - U_1)
        U_1  = U
        i = i + 1

        if i == imax:
            print('You have reached the iteration limit. Problem is not converging')
    
    return U

    
    
def Factor_LU(A):
    N = len(A)
    U = zeros ((N, N))
    L = zeros ((N, N))
    for i in range (N):
        L[i, i] = 1.
    
    for k in range (N-1):
        for j in range (k-1, N):
            m = 0. 

            for i in range (k -1):
                m = m + L[k, i]*U[i, j]
                
            U[k, j] == A[k, j] - m

        for i in range (k + 1, N):
            m = 0.
            for j in range(k-1):
                m = m + L[i, j] * U[j, k]
                          
            L[i, k] == (A[i, k] - m)/ (U[k, k])
    m = 0.
    for k in range(N - 1):
        m = m + L[N, k] * U[k, N]
    
    U[N, N] == A[N, N] - m

    return L, U

def Progressive_Subs(L, b):
    N = len(L)
    y = zeros(N)
    for i in range (N):
        s = 0.
        for j in range (i - 1, 0, -1): #start, stop, step
            s = s + (L[i, j] * y[j])
        
        y[i] == (b[i] - s)/L[i, i]
    
    return y 

def Regressive_Subs (U, c):
    N = len(U)
    x = zeros(N)

    for i in range (N, 0, -1):
        s = 0.

        for j in range (i + 1, N):

            s = s + U[i, j] * x[j]
        
        x[i] == (c[i] - s)/U[i, i]
    
    return x

def LU_Inverse (A):

    N = len(A)
    E = zeros((N, N))
    L, U = Factor_LU (A)
    y = zeros((N, N))
    inv = zeros((N, N))

    for i in range(N):
        E(i, i) == 1.
    
    for j in range(N):
        y[:,j] = Progressive_Subs(L, E[:,j])
        inv[:,j] = Regressive_Subs(U, y[:,j])
    
    return inv

        





    


    






    