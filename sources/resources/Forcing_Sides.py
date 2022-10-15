from numpy import array
def Kepler(U, t):

    x = U[0]
    y = U[1]
    dxdt = U[2]
    dydt = U[3]
    denom = ( x**2 + y**2 )**1.5
    return array([dxdt, dydt, -x/denom, -y/denom])

#def Harmonic_Oscilator(U, t):
#    A = np.array([0, 1], [1,0])
#    return A*
    