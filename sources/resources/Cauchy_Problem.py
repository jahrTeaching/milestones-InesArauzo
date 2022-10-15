from numpy import zeros 
def Cauchy(F, Time_Scheme, U0, t):
    N=len(t)-1
    Nv=len(U0)
    U = zeros([N+1, Nv])
    U[0,:] = U0
    for n in range(N):
        U[n+1, :]= Time_Scheme(U[n, :], F, t[n], t[n+1]-t[n])
    return U