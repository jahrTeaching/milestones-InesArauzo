# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 17:03:33 2022

@author: Rafael Rivero de Nicolás
"""

import numpy as np
import matplotlib.pyplot as plt

def z_x(x):
    # return 0.769*(x**2)-1.0

    return np.sqrt( (x+1)/0.769 )

def Derivative(f,x):

    return ( f(x+1E-5)-f(x) ) / (1E-5)


def Betta_calculation(f,x):

    betta = np.arctan(Derivative(f,x)) # OJO CUANDO BETTA ES NEGATIVO, AHÍ LA CURVA ESTÁ DETRÁS DE UN MÁXIMO Y AHÍ NO LLEGA EL AIRE SEGÚN LA TEORÍA DE NEWTON

    return betta

M_inf = 8


x = np.linspace(0,1)
plt.plot(x,z_x(x))

betta = Betta_calculation(z_x, x)
# plt.plot( x, betta )

cp = 2*(np.sin(betta))**2
# plt.plot( x, cp )

p_pinf = 1 + cp * (0.7*M_inf**2)
# plt.plot( x, p_pinf )







