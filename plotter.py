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
import datetime as dt
import matplotlib.dates as mdates


# get mesh from file in folder mesh (if empty run mesh_generator.py)
mesh = Mesh('mesh/mesh2d.xml.gz')
data_array = [1, 2, 12, 22, 32, 45, 56, 87, 120, 181, 243, 316, 365, 434, 626, 838, 1172, 1529, 1962, 2382, 2710, 3774]

# read border points from txt file
x_border = np.genfromtxt('data/switzerland_x.txt')
y_border = np.genfromtxt('data/switzerland_y.txt')

# define function space for system
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1, P1])
V = FunctionSpace(mesh, element)
W = FunctionSpace(mesh, P1)

triang = tri.Triangulation(*mesh.coordinates().reshape((-1, 2)).T,
                           triangles=mesh.cells())
smax = 1.0
emax = 0.2
imax = 0.2
res = 80
bounds = np.linspace(smax / res, smax - smax / res, res)
bounde = np.linspace(emax / res, emax - emax / res, res)
boundi = np.linspace(imax / res, imax - imax / res, res)
mybound1 = np.linspace(50.0, 1950.0, 40)
mybound2 = np.linspace(10000.0, 390000.0, 40)

f1 = Function(W)
f2 = Function(W)
f3 = Function(W)
f_in1 = XDMFFile("videomaker/functions/function_S.xdmf")
f_in2 = XDMFFile("videomaker/functions/function_E.xdmf")
f_in3 = XDMFFile("videomaker/functions/function_I.xdmf")

daily_infected_array = []
cumulative_infected_array = []
abs_array = []
cumulative_infected = 0.0
area_switzerland = 41285000000.0
for i in range(31):
    f_in1.read_checkpoint(f1, "S", i)
    f1.set_allow_extrapolation(True)
    f_in2.read_checkpoint(f2, "E", i)
    f2.set_allow_extrapolation(True)
    f_in3.read_checkpoint(f3, "I", i)
    f3.set_allow_extrapolation(True)

    daily_infected = assemble(f3*dx) / area_switzerland
    cumulative_infected += daily_infected
    daily_infected_array.append(daily_infected)
    cumulative_infected_array.append(cumulative_infected)

    #abs_array.append(daily_infected * 8570000.0)

    if i < 22:
        abs_array.append(daily_infected * 8570000.0)

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.xticks([])
    plt.yticks([])
    plt.title('Fraction of initial population that is susceptible')
    plt.plot(x_border, y_border, color='r')
    Z = f1.compute_vertex_values(mesh)
    c = plot(interpolate(f1, W), mode='color', vmin=0.0, vmax=smax)
    plt.tricontourf(triang, Z, vmin=0.0, vmax=smax, levels=bounds, extend='both')
    plt.colorbar(c)
    plt.savefig('videomaker/Images_S/' + str(i + 1000) + '.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.xticks([])
    plt.yticks([])
    plt.title('Fraction of initial population that is exposed')
    plt.plot(x_border, y_border, color='r')
    Z = f2.compute_vertex_values(mesh)
    c = plot(interpolate(f2, W), mode='color', vmin=0, vmax=emax)
    plt.tricontourf(triang, Z, vmin=0.0, vmax=emax, levels=bounde, extend='both')
    plt.colorbar(c)
    plt.savefig('videomaker/Images_E/' + str(i + 1000) + '.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.xticks([])
    plt.yticks([])
    plt.title('Fraction of initial population that is infected')
    plt.plot(x_border, y_border, color='r')
    Z = f3.compute_vertex_values(mesh)
    c = plot(interpolate(f3, W), mode='color', vmin=0, vmax=imax)
    plt.tricontourf(triang, Z, vmin=0.0, vmax=imax, levels=boundi, extend='both')
    plt.colorbar(c)
    plt.savefig('videomaker/Images_I/' + str(i + 1000) + '.jpg')
    plt.clf()

f4 = Function(W)
f5 = Function(W)
f6 = Function(W)
f_in4 = XDMFFile("videomaker/functions/checkdiff.xdmf")
f_in5 = XDMFFile("videomaker/functions/checkrho.xdmf")
f_in6 = XDMFFile("videomaker/functions/checkbeta.xdmf")

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
plt.title('Diffusion of neighbouring municipalities')
plt.plot(x_border, y_border, color='r')
Z = f4.compute_vertex_values(mesh)
c = plot(interpolate(f4, W), vmin = 0.0, vmax = 2000.0,  mode='color')
plt.tricontourf(triang, Z, vmin = 0.0, vmax = 2000.0, levels = mybound1, extend = 'both')
plt.colorbar(c)
plt.savefig('videomaker/Images_else/checkdiff.jpg')
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
plt.savefig('videomaker/Images_else/checkrho.jpg')
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
plt.savefig('videomaker/Images_else/checkbeta.jpg')
plt.clf()

plt.xlabel('time in days')
plt.ylabel('fraction of infected')
plt.title('fraction of infected individuals over time')
plt.plot(daily_infected_array)
plt.savefig('videomaker/Images_else/infected_daily.jpg')
plt.clf()

plt.xlabel('time in days')
plt.ylabel('cummulated fraction of infected')
plt.title('cummulated fraction of infected individuals over time')
plt.plot(cumulative_infected_array)
plt.savefig('videomaker/Images_else/infected_cummulated.jpg')
plt.clf()

ax = plt.gca()
plt.xlabel('time in days')
plt.ylabel('infected individuals')
plt.title('infected individuals over time')
dates = ['24/02/2020', '25/02/2020', '26/02/2020', '27/02/2020','28/02/2020', '29/02/2020', '01/03/2020', '02/03/2020', '03/03/2020', '04/03/2020', '05/03/2020', '06/03/2020', '07/03/2020', '08/03/2020', '09/03/2020', '10/03/2020', '11/03/2020', '12/03/2020', '13/03/2020', '14/03/2020', '15/03/2020', '16/03/2020']
x = [dt.datetime.strptime(d,'%d/%m/%Y').date() for d in dates]
formatter = mdates.DateFormatter("%d-%m-%Y")
ax.xaxis.set_major_formatter(formatter)
locator = mdates.DayLocator()
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
model, = plt.plot(x, abs_array, 'r')
data, = plt.plot(x, data_array, 'bo')
plt.legend([model, data], ['Model prediction', 'Data from BAG'])
plt.ylim(0.0, 8000.0)
plt.savefig('videomaker/Images_else/modelvsbag.jpg')
plt.clf()

