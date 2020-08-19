
# Approach to solve Epidemic Model with Diffusion in 2D


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
import pandas as pd
from shapely.geometry import Point
import numpy as np
import math

# define all needed parameters
T = 60.0                                # final time in days
dt = 1.0/24.0                           # time step size at beginning
theta_factor = Constant(1.1)            # factor to represent the underreporting of movement
oneoverd = Constant(1.0 / 4.0)          # one over average duration of infection
oneoverz = Constant(1.0 / 3.0)          # one over average latency period

theta = Constant(0.5)                   # theta = 0.5 means Crank-Nicolson
t = 0.0                                 # global time
timestep = [0.0]                        # array to safe the timesteps
rho = 0.1                               # safety factor
tol = 100.0                             # tolerance of the l2 norm
toldt = 0.01                            # smallest allowed timestep


# get mesh from file in folder mesh
mesh = Mesh('mesh/mesh2d.xml.gz')

# read border points from txt files
x_border = np.genfromtxt('shapefiles/switzerland_x.txt')
y_border = np.genfromtxt('shapefiles/switzerland_y.txt')

# define function space for system
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1, P1])
V = FunctionSpace(mesh, element)
W = FunctionSpace(mesh, P1)

# define test functions
v_1, v_2, v_3 = TestFunctions(V)

# define used functions
SEI_low = Function(V)
SEI_half = Function(V)
SEI_high = Function(V)
SEI_n = Function(V)
d = Function(V)

# create files and parameters for visualization output and time
name = 100000
triang = tri.Triangulation(*mesh.coordinates().reshape((-1, 2)).T,
                           triangles=mesh.cells())
bounds = np.linspace(0.0275,1.0725,40)
bounde = np.linspace(0.01,0.39,40)
boundi = np.linspace(0.01,0.39,40)

# setting Initial Conditions
SEI_0 = Expression(('1.0',
                    '0.4*exp(-0.00000001*(pow(x[0]-684000.0,2)+pow(x[1]-247000.0,2))) + '
                    '0.35*exp(-0.00000001*(pow(x[0]-600000.0,2)+pow(x[1]-165000.0,2))) + '
                    '0.4*exp(-0.00000001*(pow(x[0]-620000.0,2)+pow(x[1]-215000.0,2))) + '
                    '0.2*exp(-0.00000001*(pow(x[0]-710000.0,2)+pow(x[1]-235000.0,2)))',
                    '0.0'),
                    degree=2)
SEI_n = project(SEI_0, V)

#Diffusion, inverse population density and beta from files
f1 = Function(W)
f_in = XDMFFile("difffun/diff.xdmf")
f_in.read_checkpoint(f1, "g", 0)
f1.set_allow_extrapolation(True)
d_0 = Expression(('function',
                  'function',
                  'function'),
                 function = f1, degree = 2)
d = project(d_0, V)

f2 = Function(W)
f_in = XDMFFile("difffun/rhoinv.xdmf")
f_in.read_checkpoint(f2, "g", 0)
f2.set_allow_extrapolation(True)
rhoinv_0 = Expression(('function',
                       'function',
                       'function'),
                       function = f2, degree = 2)
rhoinv = project(rhoinv_0, V)

f3 = Function(W)
f_in = XDMFFile("difffun/beta.xdmf")
f_in.read_checkpoint(f3, "g", 0)
f3.set_allow_extrapolation(True)

beta = project(f3, W)

# check diffusion, 1/rho and beta and plot S, E, I for t=0 and safe as images
_S_0, _E_0, _I_0 = SEI_n.split()
_d_s, _d_i, _d_r = d.split()
_rhoinv1, _rhoinv2, _rhoinv3 = rhoinv.split()

mybound1 = np.linspace(50.0, 1950.0, 40)
mybound2 = np.linspace(10000.0, 390000.0, 40)

"""
plt.xlabel('space [x]')
plt.ylabel('space [y]')
plt.xticks([])
plt.yticks([])
plt.title('Diffusion of neighboring municipalities')
plt.plot(x_border, y_border, color='r')
Z = _d_s.compute_vertex_values(mesh)
c = plot(interpolate(_d_s, W), vmin = 0.0, vmax = 2000.0,  mode='color')
plt.tricontourf(triang, Z, vmin = 0.0, vmax = 2000.0, levels = mybound1, extend = 'both')
plt.colorbar(c)
plt.savefig('checkdiff1.jpg')
plt.clf()

plt.xlabel('space [x]')
plt.ylabel('space [y]')
plt.xticks([])
plt.yticks([])
plt.title('1 / population density')
plt.plot(x_border, y_border, color='r')
Z = _rhoinv1.compute_vertex_values(mesh)
c = plot(interpolate(_rhoinv1, W), vmin=0.0, vmax=400000.0, mode='color')
plt.tricontourf(triang, Z, vmin=0.0, vmax=400000.0, levels = mybound2, extend = 'both')
plt.colorbar(c)
plt.savefig('checkrho1.jpg')
plt.clf()

plt.xlabel('space [x]')
plt.ylabel('space [y]')
plt.xticks([])
plt.yticks([])
plt.title('beta')
plt.plot(x_border, y_border, color='r')
Z = beta.compute_vertex_values(mesh)
c = plot(interpolate(beta, W))
plt.tricontourf(triang, Z)
plt.colorbar(c)
plt.savefig('checkbeta1.jpg')
plt.clf()
"""

plt.xlabel('space [x]')
plt.ylabel('space [y]')
plt.xticks([])
plt.yticks([])
plt.title('Percentage of initial population that is susceptible')
plt.plot(x_border, y_border, color='r')
Z = _S_0.compute_vertex_values(mesh)
c = plot(interpolate(_S_0, W), mode='color', vmin=0.0, vmax=1.1)
plt.tricontourf(triang, Z, vmin=0.0, vmax=1.1, levels = bounds, extend = 'both')
plt.colorbar(c)
plt.savefig('Videomaker/Images_S/' + str(name) + '.jpg')
plt.clf()

plt.xlabel('space [x]')
plt.ylabel('space [y]')
plt.xticks([])
plt.yticks([])
plt.title('Percentage of initial population that is exposed')
plt.plot(x_border, y_border, color='r')
Z = _E_0.compute_vertex_values(mesh)
c = plot(interpolate(_E_0, W), mode='color', vmin=0, vmax=0.4)
plt.tricontourf(triang, Z, vmin=0.0, vmax=0.4, levels = bounde, extend = 'both')
plt.colorbar(c)
plt.savefig('Videomaker/Images_E/' + str(name) + '.jpg')
plt.clf()

plt.xlabel('space [x]')
plt.ylabel('space [y]')
plt.xticks([])
plt.yticks([])
plt.title('Percentage of initial population that is infected')
plt.plot(x_border, y_border, color='r')
Z = _I_0.compute_vertex_values(mesh)
c = plot(interpolate(_I_0, W), mode='color', vmin=0, vmax=0.4)
plt.tricontourf(triang, Z, vmin=0.0, vmax=0.4, levels = boundi, extend = 'both')
plt.colorbar(c)
plt.savefig('Videomaker/Images_I/' + str(name) + '.jpg')
plt.clf()

name += 1

# split system functions to access components
S_low, E_low, I_low = split(SEI_low)
S_half, E_half, I_half = split(SEI_half)
S_high, E_high, I_high = split(SEI_high)
S_n, E_n, I_n = split(SEI_n)
d_s, d_e, d_i = split(d)
rhoinv_s, rhoinv_e, rhoinv_i = split(rhoinv)

# time stepping
n = 0
while t < T:
    # Define variational problem for SEI_low
    F = S_low * v_1 * dx - S_n * v_1 * dx - theta * dt * (-beta * S_low * I_low * v_1 * dx + \
        theta_factor * (- rhoinv_s * d_s * dot(grad(S_low), grad(v_1)) * dx)) - \
        (1.0 - theta) * dt *(-beta * S_n * I_n * v_1 * dx + theta_factor * ( - rhoinv_s * d_s * dot(grad(S_n), grad(v_1)) * dx)) + \
        E_low * v_2 * dx - E_n * v_2 * dx - theta * dt * (beta * S_low * I_low * v_2 * dx - oneoverz * E_low * v_2 * dx + \
        theta_factor * (- rhoinv_e * d_e * dot(grad(E_low), grad(v_2)) * dx)) - \
        (1.0 - theta) * dt * (beta * S_n * I_n * v_2 * dx - oneoverz * E_n * v_2 * dx + \
        theta_factor * (- rhoinv_e * d_e * dot(grad(E_n), grad(v_2)) * dx)) + \
        I_low * v_3 * dx - I_n * v_3 * dx - theta * dt * (oneoverz * E_low * v_3 * dx - oneoverd * I_low * v_3 * dx + \
        theta_factor * (- rhoinv_i * d_i * dot(grad(I_low), grad(v_3)) * dx)) - \
        (1.0 - theta) * dt * (oneoverz * E_n * v_3 * dx - oneoverd * I_n * v_3 * dx + \
        theta_factor * (- rhoinv_i * d_i * dot(grad(I_n), grad(v_3)) * dx))

    Jac = derivative(F, SEI_low, TrialFunction(V))

    # solve variational problem for SEI_low
    solve(F == 0, SEI_low, J=Jac)

    # Define variational problem for SEI_half
    F = S_half * v_1 * dx - S_n * v_1 * dx - theta * dt / 2.0 * (-beta * S_half * I_half * v_1 * dx + \
        theta_factor * (- rhoinv_s * d_s * dot(grad(S_half), grad(v_1)) * dx)) - \
        (1.0 - theta) * dt / 2.0 * (-beta * S_n * I_n * v_1 * dx + theta_factor * (- rhoinv_s * d_s * dot(grad(S_n),grad(v_1)) * dx)) + \
        E_half * v_2 * dx - E_n * v_2 * dx - theta * dt / 2.0 * (beta * S_half * I_half * v_2 * dx - oneoverz * E_half * v_2 * dx + \
        theta_factor * (- rhoinv_e * d_e * dot(grad(E_half), grad(v_2)) * dx)) - \
        (1.0 - theta) * dt / 2.0 * (beta * S_n * I_n * v_2 * dx - oneoverz * E_n * v_2 * dx + \
        theta_factor * (- rhoinv_e * d_e * dot(grad(E_n), grad(v_2)) * dx)) + \
        I_half * v_3 * dx - I_n * v_3 * dx - theta * dt / 2.0 * (oneoverz * E_half * v_3 * dx - oneoverd * I_half * v_3 * dx + \
        theta_factor * (- rhoinv_i * d_i * dot(grad(I_half), grad(v_3)) * dx)) - \
        (1.0 - theta) * dt / 2.0 * (oneoverz * E_n * v_3 * dx - oneoverd * I_n * v_3 * dx + \
        theta_factor * (- rhoinv_i * d_i * dot(grad(I_n), grad(v_3)) * dx))

    Jac = derivative(F, SEI_half, TrialFunction(V))

    # solve variational problem for SEI_half
    solve(F == 0, SEI_half, J=Jac)

    # Define variational problem for SEI_high
    F = S_high * v_1 * dx - S_half * v_1 * dx - theta * dt / 2.0 * (-beta * S_high * I_high * v_1 * dx + \
        theta_factor * (- rhoinv_s * d_s * dot(grad(S_high), grad(v_1)) * dx)) - \
        (1.0 - theta) * dt / 2.0 * (-beta * S_half * I_half * v_1 * dx + theta_factor * (- rhoinv_s * d_s * dot(grad(S_half), grad(v_1)) * dx)) + \
        E_high * v_2 * dx - E_half * v_2 * dx - theta * dt / 2.0 * (beta * S_high * I_high * v_2 * dx - oneoverz * E_high * v_2 * dx + \
        theta_factor * (- rhoinv_e * d_e * dot(grad(E_high), grad(v_2)) * dx)) - \
        (1.0 - theta) * dt / 2.0 * (beta * S_half * I_half * v_2 * dx - oneoverz * E_half * v_2 * dx + \
        theta_factor * (- rhoinv_e * d_e * dot(grad(E_half), grad(v_2)) * dx)) + \
        I_high * v_3 * dx - I_half * v_3 * dx - theta * dt / 2.0 * (oneoverz * E_high * v_3 * dx - oneoverd * I_high * v_3 * dx + \
        theta_factor * (- rhoinv_i * d_i * dot(grad(I_high), grad(v_3)) * dx)) - \
        (1.0 - theta) * dt / 2.0 * (oneoverz * E_half * v_3 * dx - oneoverd * I_half * v_3 * dx + \
        theta_factor * (- rhoinv_i * d_i * dot(grad(I_half), grad(v_3)) * dx))

    Jac = derivative(F, SEI_high, TrialFunction(V))

    # solve variational problem for SEI_high
    solve(F == 0, SEI_high, J=Jac)

    # compute error by Richardson extrapolation
    error = errornorm(SEI_high, SEI_low, 'l2') / 3.0
    print("This is the error: ", error)

    # compute new timestep
    dt_new = pow(rho * tol / error, 0.5) * dt
    # check if dt_new gets too small
    if dt_new < toldt:
        dt_new = toldt

    # if error small enough or minimum timestep reached, update and go to the next timestep
    if (error < tol) or (dt_new == toldt):
        n += 1
        print("This is iteration number: ", n)
        print("And the time is: ", t)
        newtime = t + dt
        oldtime = t
        # interpolate between timesteps to get plots
        while (math.floor(newtime) - math.floor(t))>0:
            SEI_plot = project((SEI_high - SEI_n) * (math.floor(t) + 1.0 - oldtime) / dt + SEI_n, V)
            _S, _E, _I = SEI_plot.split()
            plt.xlabel('space [x]')
            plt.ylabel('space [y]')
            plt.xticks([])
            plt.yticks([])
            plt.title('Percentage of initial population that is susceptible')
            plt.plot(x_border, y_border, color='r')
            Z = _S.compute_vertex_values(mesh)
            c = plot(interpolate(_S, W), mode='color', vmin=0.0, vmax=1.1)
            plt.tricontourf(triang, Z, vmin=0.0, vmax=1.1, levels=bounds, extend='both')
            plt.colorbar(c)
            plt.savefig('Videomaker/Images_S/' + str(name) + '.jpg')
            plt.clf()

            plt.xlabel('space [x]')
            plt.ylabel('space [y]')
            plt.xticks([])
            plt.yticks([])
            plt.title('Percentage of initial population that is exposed')
            plt.plot(x_border, y_border, color='r')
            Z = _E.compute_vertex_values(mesh)
            c = plot(interpolate(_E, W), mode='color', vmin=0.0, vmax=0.4)
            plt.tricontourf(triang, Z, vmin=0.0, vmax=0.4, levels=bounde, extend='both')
            plt.colorbar(c)
            plt.savefig('Videomaker/Images_E/' + str(name) + '.jpg')
            plt.clf()

            plt.xlabel('space [x]')
            plt.ylabel('space [y]')
            plt.xticks([])
            plt.yticks([])
            plt.title('Percentage of initial population that is infected')
            plt.plot(x_border, y_border, color='r')
            Z = _I.compute_vertex_values(mesh)
            c = plot(interpolate(_I, W), mode='color', vmin=0.0, vmax=0.4)
            plt.tricontourf(triang, Z, vmin=0.0, vmax=0.4, levels=boundi, extend='both')
            plt.colorbar(c)
            plt.savefig('Videomaker/Images_I/' + str(name) + '.jpg')
            plt.clf()

            name += 1
            t = math.floor(t) + 1.0

        t = newtime
        SEI_n.assign(SEI_high)
        timestep.append(t)

    dt = dt_new

# print array of timesteps
print(timestep)
np.savetxt("timestep.txt", timestep, fmt="%s")
