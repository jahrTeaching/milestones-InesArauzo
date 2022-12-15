"""
Created on Mon Nov 14 17:03:33 2022
@authors: Rafael Rivero de Nicolás, Guillermo García del Río, Inés Arauzo Andrés
"""
from numpy import deg2rad, array, size, dot, zeros, min, max, logical_and, sqrt, hstack, sin, cos, arctan, linspace, geomspace, reshape
from numpy.linalg import norm
from time import process_time

from numba import njit

import matplotlib.pyplot as plt
from matplotlib import rc # LaTeX tipography
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rc('text', usetex=True); plt.rc('font', family='serif')

import matplotlib
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

def Plot_geometry_1D(x, y, title=r'2D Geometry', axis_labels=[r'$x$ [-]', r'$z(x)$ [-]']):

    fig, ax = plt.subplots(1,1, figsize=(9,6), constrained_layout='true')
    ax.set_xlim( x[0] , x[-1] )
    ax.set_ylim( -y[0] - x[-1] + x[0], y[0] + x[-1] - x[0] )
    ax.set_title(title, fontsize=20)
    ax.grid()
    ax.set_xlabel(axis_labels[0] ,fontsize=20)
    ax.set_ylabel(axis_labels[1] ,fontsize=20)
    ax.plot( x, y, c="blue" )
    ax.plot( x, -y, c="blue" )
    plt.plot()

def Simple_plot_1D(x, y, title=None, axis_labels=[None, None]):

    fig, ax = plt.subplots(1,1, figsize=(9,6), constrained_layout='true')
    ax.set_xlim( x[0] , x[-1] )
    ax.set_ylim( min(y)*0.95, max(y)*1.05 )
    ax.set_title(title, fontsize=20)
    ax.grid()
    ax.set_xlabel(axis_labels[0] ,fontsize=20)
    ax.set_ylabel(axis_labels[1] ,fontsize=20)
    ax.plot( x, y, c="blue" )
    plt.plot()
   
################################
