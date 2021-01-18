# This program calculates the diffusion functions by checking for every point in which municipality it is located
# and set it to the value from the txt files created with diffusion_reader.py


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

# turn plotting and saving off and on
plotbool = True
savebool = True

# get data
beta = 1.0
coefficients = np.genfromtxt('data/diff.txt')
densityinv = np.genfromtxt('data/rho.txt')
switzerland = gpd.read_file('data/municipalities.shp')
switzerland_new = switzerland.translate(xoff=-2000000.0, yoff=-1000000.0, zoff=0.0)
print(switzerland_new.head())

# read border points from txt file
x_border = np.genfromtxt('data/switzerland_x.txt')
y_border = np.genfromtxt('data/switzerland_y.txt')

# get mesh from file in folder mesh
mesh = Mesh('mesh/mesh2d.xml.gz')

# define function space for system
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1, P1])
V = FunctionSpace(mesh, element)
W = FunctionSpace(mesh, P1)

# checks for a Point (x,y) in which municipality it is located. If in none return -1
def containsPoint(x, y):
    ID_municipality=0
    point = Point(x, y)

    while ID_municipality<2337:         #2337
        gemeinde = switzerland_new[ID_municipality]
        containing = gemeinde.contains(point)
        if containing == True:
            return ID_municipality
        ID_municipality += 1
    return -1

histodata = []
# check for every point in which municipality it is and set values according to txt files
class MyExpression0(UserExpression):
    def eval(self, value, x):
        a = containsPoint(x[0], x[1])
        diff = coefficients[a]
        dens = densityinv[a]
        value[0] = 0.0
        value[1] = 0.0
        value[2] = 0.0
        if diff > 2000.0:
            diff = 2000.0
        if dens > 400000.0:
            dens = 400000.0
        if a > -1:
            value[0] = dens
            value[1] = diff + 500.0      # add 500 to represent non-work-related movement to avoid diffusion zero in areas with no commuters to neighbouring municipalities
            value[2] = beta
    def value_shape(self):
        return (3,)
g_ = MyExpression0()
g = interpolate(g_, V)

_g1, _g2, _g3 = g.split()

# plot functions
if plotbool == True:
    plt.title('1 / Population Density')
    plt.plot(x_border, y_border, color='r')
    c = plot(_g1)
    plt.colorbar(c)
    plt.savefig('mesh/densinv.jpg')
    plt.clf()

    plt.xticks([])
    plt.yticks([])
    plt.title('Diffusion coefficient')
    plt.plot(x_border, y_border, color='r')
    c = plot(_g2)
    bar = plt.colorbar(c)
    bar.set_ticks([])
    plt.savefig('mesh/diffusion.jpg')
    plt.clf()

    plt.title('beta')
    plt.plot(x_border, y_border, color='r')
    c = plot(_g3)
    plt.colorbar(c)
    plt.savefig('mesh/beta.jpg')
    plt.clf()

# save functions to folder difffun
if savebool == True:
    g1, g2, g3 = g.split()
    f_out = XDMFFile("difffun/rhoinv.xdmf")
    f_out.write_checkpoint(project(g1, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("difffun/diff.xdmf")
    f_out.write_checkpoint(project(g2, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("difffun/beta.xdmf")
    f_out.write_checkpoint(project(g3, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()
