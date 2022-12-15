# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 16:08:17 2022

@author:  Rafael Rivero de Nicolás, Guillermo García del Río, Inés Arauzo Andrés
"""

from numpy import deg2rad, array, size, dot, zeros, min, max, logical_and, sqrt, hstack, sin, cos, arctan, linspace, geomspace, reshape, pi, floor, ceil, sum

import Aerodynamics as aero
import Math_Utilities as mth
import Normals_3D as nor
import Geometry as gm
import OurDecorators as pl


n_x = 10000; dx = 1/n_x; n_phi = 100; dphi = 2*pi/n_phi; # radianes

x = linspace(0, 0.95, n_x)
z = gm.Quarter_of_a_Circle(x)
dz = mth.Analytical_Derivative( x, gm.Quarter_of_a_Circle )

# z = gm.conito(x)
# dz = mth.Analytical_Derivative( x, gm.conito )

pl.Plot_geometry_1D(x, z)


N_panels = (n_x-1)*(n_phi-1) # Número de paneles totales en los que se ha mallado la superficie

# %% Función de los coeficientes

S_rev = 0
Coefs = zeros(3)

for i in range(n_x):

    S_rev = S_rev + 2*pi * z[i] * sqrt( 1 + (dz[i])**2 ) * dx


d_Surface = zeros(N_panels)

for n in range(N_panels):

    i = int(floor( n/ (n_x-1) ))
    # print(i)

    hipotenusa = sqrt( ( x[i+1] - x[i] )**2 + ( z[i+1] - z[i] )**2 )

    d_Surface[n] = dphi * ( z[i+1] + z[i] )/2 * hipotenusa

    Coefs =






S_an = 4 * pi * 1**2 / 2 # Solo para la esfera


















