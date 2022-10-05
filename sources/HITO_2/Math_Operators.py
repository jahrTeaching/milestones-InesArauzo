from numpy import zeros, matmul, array
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
    

