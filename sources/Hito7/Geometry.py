# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 16:14:06 2022

@author: Rafael Rivero de Nicolás, Guillermo García del Río, Inés Arauzo Andrés
"""

from numpy import deg2rad, array, size, dot, zeros, min, max, logical_and, sqrt, hstack, sin, cos, arctan, linspace, geomspace, reshape
from numpy.linalg import norm
from time import process_time


#### Funciones para definir la geometría ####
def Horizontal_Parabola(x): # Libro Anderson, Cap 3, page 63

    return sqrt( (x+1)/0.769 )

def Quarter_of_a_Circle(x):

    return sqrt( x )

def conito(x):

    return x

def Flechita(x): # Se la ha inventado rafa

    f1 = x[x<=0.5]
    f2 = 2*x[logical_and(x>0.5, x<=0.75)]-0.5
    f3 = 1-4*(x[x>0.75]-0.75)

    return hstack((f1, f2, f3))
###########################
