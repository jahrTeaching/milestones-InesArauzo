from numpy import zeros
from numpy.linalg import eig
from scipy.optimize import fsolve
from resources.Forcing_Sides import R3BodyProblem
from resources.Math_Operators import Jacobian

def Lagrange_Points(U0, NL):

    LP = zeros([5,2])

    def F(Y):
        
        X = zeros(4)
        X[0:2] = Y
        X[2:4] = 0
        FX = R3BodyProblem(X, 0)
        return FX[2:4]
        
    for i in range(NL):
        LP[i,:] = fsolve(F, U0[i,0:2])

    return LP

def Lagrange_Points_Stability(U0):

    def F(Y):
        return R3BodyProblem(Y, 0 )

    A = Jacobian(F, U0)
    values, vectors = eig(A)

    return values
