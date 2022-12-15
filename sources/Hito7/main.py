"""
Created on Sun Dec 11.33.00 2022
@authors: Rafael Rivero de Nicolás, Guillermo García del Río, Inés Arauzo Andrés
"""
### IMPORTS ###

from numpy import deg2rad, array, size, dot, zeros, min, max, logical_and, sqrt, hstack, sin, cos, arctan, linspace, geomspace, reshape, pi, floor, ceil, sum

from Aerodynamics import Cp_3D_Modified, Cp_3D # Cp_Dist_1D_by_Newton, Cp_Dist_1D_by_Newton_Modified,
from Math_Utilities import Analytical_Derivative
#from Normals_3D import 
from Geometry import Quarter_of_a_Circle, conito, Flechita
#from OurDecorators import 

### GEOMETRIES AVAILABLE ###

GEOMETRY = {"Semiesfera":   Quarter_of_a_Circle, 
            "Conito":       conito, 
            "Flechita":     Flechita}

######################################
############### INPUTS ###############

Selection = 

n_x = 10       # Number of x-axis panels for a linspace partition [-]
n_phi = 10     # Number of panels for each revoluting section [-]
alpha = 0       # Angle of attack of the incident flow [DEG]
beta = 0        # Sideslip angle of the incident flow [DEG]
Mach = 5        # Mach number of the incident flow [-]
gamma = 1.4     # Adiabatic coefficent of the incident flow [-]

######################################
######################################

Geometry_selected = GEOMETRY[Selection]

x = linspace(0, 1, n_x + 1)

z = Geometry_selected(x)
dz = Analytical_Derivative( x, Geometry_selected ) #

Cps = Cp_3D(alpha,beta, n_x, n_phi, Geometry_selected)
print(Cps)

# dx = 1/n_x
# dphi = 2*pi/n_phi; # radianes



# input("Select the number of the corresponding geometry configuration:" \
#                   "\n" \
#                   "\n" \
#                   " 1 :: \t Semiesfera \n" \
#                   "2 :: conito \n" \
#                   "3 :: Flechita ")