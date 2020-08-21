"""
This program plots everything
"""

from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from fenics import *
from dolfin import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import shapefile
import os
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
import math


# get mesh from file in folder mesh (if empty run mesh_generator.py)
mesh = Mesh('mesh/mesh2d.xml.gz')

# read border points from txt file
x_border = np.genfromtxt('shapefiles/switzerland_x.txt')
y_border = np.genfromtxt('shapefiles/switzerland_y.txt')

# define function space for system
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1, P1])
V = FunctionSpace(mesh, element)
W = FunctionSpace(mesh, P1)

triang = tri.Triangulation(*mesh.coordinates().reshape((-1, 2)).T,
                           triangles=mesh.cells())
bounds = np.linspace(0.0275,1.0725,40)
bounde = np.linspace(0.01,0.39,40)
boundi = np.linspace(0.01,0.39,40)
mybound1 = np.linspace(50.0, 1950.0, 40)
mybound2 = np.linspace(10000.0, 390000.0, 40)

f1 = Function(W)
f2 = Function(W)
f3 = Function(W)
f_in1 = XDMFFile("Videomaker/functions/function_S.xdmf")
f_in2 = XDMFFile("Videomaker/functions/function_E.xdmf")
f_in3 = XDMFFile("Videomaker/functions/function_I.xdmf")

for i in range(21):
    f_in1.read_checkpoint(f1, "S", i)
    f1.set_allow_extrapolation(True)
    f_in2.read_checkpoint(f2, "E", i)
    f2.set_allow_extrapolation(True)
    f_in3.read_checkpoint(f3, "I", i)
    f3.set_allow_extrapolation(True)

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.xticks([])
    plt.yticks([])
    plt.title('Percentage of initial population that is susceptible')
    plt.plot(x_border, y_border, color='r')
    Z = f1.compute_vertex_values(mesh)
    c = plot(interpolate(f1, W), mode='color', vmin=0.0, vmax=1.1)
    plt.tricontourf(triang, Z, vmin=0.0, vmax=1.1, levels=bounds, extend='both')
    plt.colorbar(c)
    plt.savefig('Videomaker/Images_S/' + str(i) + '.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.xticks([])
    plt.yticks([])
    plt.title('Percentage of initial population that is exposed')
    plt.plot(x_border, y_border, color='r')
    Z = f2.compute_vertex_values(mesh)
    c = plot(interpolate(f2, W), mode='color', vmin=0, vmax=0.4)
    plt.tricontourf(triang, Z, vmin=0.0, vmax=0.4, levels=bounde, extend='both')
    plt.colorbar(c)
    plt.savefig('Videomaker/Images_E/' + str(i) + '.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.xticks([])
    plt.yticks([])
    plt.title('Percentage of initial population that is infected')
    plt.plot(x_border, y_border, color='r')
    Z = f3.compute_vertex_values(mesh)
    c = plot(interpolate(f3, W), mode='color', vmin=0, vmax=0.4)
    plt.tricontourf(triang, Z, vmin=0.0, vmax=0.4, levels=boundi, extend='both')
    plt.colorbar(c)
    plt.savefig('Videomaker/Images_I/' + str(i) + '.jpg')
    plt.clf()

f4 = Function(W)
f5 = Function(W)
f6 = Function(W)
f_in4 = XDMFFile("Videomaker/functions/checkdiff.xdmf")
f_in5 = XDMFFile("Videomaker/functions/checkrho.xdmf")
f_in6 = XDMFFile("Videomaker/functions/checkbeta.xdmf")

f_in4.read_checkpoint(f4, "diff", 0)
f4.set_allow_extrapolation(True)

f_in5.read_checkpoint(f5, "rho", 0)
f5.set_allow_extrapolation(True)

f_in6.read_checkpoint(f6, "beta", 0)
f6.set_allow_extrapolation(True)

plt.xlabel('space [x]')
plt.ylabel('space [y]')
plt.xticks([])
plt.yticks([])
plt.title('Diffusion of neighboring municipalities')
plt.plot(x_border, y_border, color='r')
Z = f4.compute_vertex_values(mesh)
c = plot(interpolate(f4, W), vmin = 0.0, vmax = 2000.0,  mode='color')
plt.tricontourf(triang, Z, vmin = 0.0, vmax = 2000.0, levels = mybound1, extend = 'both')
plt.colorbar(c)
plt.savefig('Videomaker/Images_else/checkdiff1.jpg')
plt.clf()

plt.xlabel('space [x]')
plt.ylabel('space [y]')
plt.xticks([])
plt.yticks([])
plt.title('1 / population density')
plt.plot(x_border, y_border, color='r')
Z = f5.compute_vertex_values(mesh)
c = plot(interpolate(f5, W), vmin=0.0, vmax=400000.0, mode='color')
plt.tricontourf(triang, Z, vmin=0.0, vmax=400000.0, levels = mybound2, extend = 'both')
plt.colorbar(c)
plt.savefig('Videomaker/Images_else/checkrho1.jpg')
plt.clf()

plt.xlabel('space [x]')
plt.ylabel('space [y]')
plt.xticks([])
plt.yticks([])
plt.title('beta')
plt.plot(x_border, y_border, color='r')
Z = f6.compute_vertex_values(mesh)
c = plot(interpolate(f6, W), vmin=0.0, vmax=1.1, mode='color')
plt.tricontourf(triang, Z, vmin=0.0, vmax=1.1, levels = bounds, extend = 'both')
plt.colorbar(c)
plt.savefig('Videomaker/Images_else/checkbeta1.jpg')
plt.clf()