"""
Approach to solve Epidemic Model with Diffusion in 2D

"""

from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from fenics import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from mpl_toolkits import mplot3d
import math
import dolfin as df

# define all needed parameters
T = 60.0                                # final time in days
dt = 1.0/24.0                           # time step size at beginning
r_rate = Constant(2.685*pow(10,-5))     # recruitment rate per day
d_rate_n = Constant(2.152*pow(10,-5))   # natural death rate per day
d_rate_d = Constant(0.1)                # death rate from disease 0.1
r_d = Constant(0.3)                     # recovery rate from disease 0.3
beta = Constant(0.9)                    # infection coefficient
alpha01 = Constant(0.1)                 # constants for incidence rate
alpha02 = Constant(0.02)
alpha03 = Constant(0.03)
A_const = 10.695                        # constant to model the airports
diff_deutsch= 677120000
diff_romand = 0.25*diff_deutsch         # diffusion coefficients for different areas
diff_tessin = 0.5*diff_deutsch
diff_i_deutsch= 0.8*diff_deutsch
diff_i_romand = 0.8*diff_romand
diff_i_tessin = 0.8*diff_tessin
b_austria = 8.320*pow(10,-9)*6.0        # constants to model the commuters from neighbour countries
b_germany = 2.922*pow(10,-8)*6.0
b_italy = 4.476*pow(10,-8)*6.0
b_france = 6.097*pow(10,-8)*6.0
radius = 10000                          # radius of circle in images to show sources
theta = Constant(0.5)                   # theta = 0.5 means Crank-Nicolson
t = 0.0                                 # global time
timestep = [0.0]                        # array to safe the timesteps
rho = 0.9                               # safety factor
tol = 2000.0                            # tolerance of the l2 norm

# get mesh from file in folder mesh (if empty run mesh_generator.py)
mesh = Mesh('mesh/mesh2d.xml.gz')

# define function space for system
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1, P1])
V = FunctionSpace(mesh, element)
W = FunctionSpace(mesh, P1)

# define test functions
v_1, v_2, v_3 = TestFunctions(V)

# define used functions
SIR_low = Function(V)
SIR_half = Function(V)
SIR_high = Function(V)
SIR_n = Function(V)
g_n = Function(V)
g_half = Function(V)
g = Function(V)
d = Function(V)
source = Function(V)

# create files and parameters for visualization output and time
name = 100000
triang = tri.Triangulation(*mesh.coordinates().reshape((-1, 2)).T,
                           triangles=mesh.cells())
bounds = np.linspace(0.0275,1.0725,40)
boundi = np.linspace(0.01,0.39,40)
boundr = np.linspace(0.0275,1.0725,40)

# setting Initial and Boundary Conditions
SIR_0 = Expression(('1.0',
                    '0.0',
                    '0.0'),
                    degree=2)
SIR_n = project(SIR_0, V)

# g_0 represents no movement, g_in is flow into Switzerland, g_out is flow out of Switzerland
g_0_ = Expression(('0.0',
                  '0.0',
                  '0.0'),
                  austria = b_austria, germany = b_germany, italy = b_italy, france = b_france, degree=2)
g_0 = interpolate(g_0_, V)
g_in_ = Expression(('0.0',
                   '((x[0]>750000) && (x[1]>210000)? austria : (x[0]<750000) && (x[0]>610000) && (x[1]>250000)? '
                    'germany : (x[0]<750000) && (x[0]>677000) && (x[1]<180000)? italy : (((x[0]<520000) && '
                    '(x[1]<160000)) || ((x[0]<610000) && (x[1]>160000)))? france : 0.0)',
                   '0.0'),
                   austria = b_austria, germany = b_germany, italy = b_italy, france = b_france, degree=2)
g_in = interpolate(g_in_, V)
g_out_ = Expression(('0.0',
                    '(-1.0)*((x[0]>750000) && (x[1]>210000)? austria : (x[0]<750000) && (x[0]>610000) && (x[1]>250000)? '
                    'germany : (x[0]<750000) && (x[0]>677000) && (x[1]<180000)? italy : (((x[0]<520000) && '
                    '(x[1]<160000)) || ((x[0]<610000) && (x[1]>160000)))? france : 0.0)',
                    '0.0'),
                    austria = b_austria, germany = b_germany, italy = b_italy, france = b_france, degree=2)
g_out = interpolate(g_out_, V)
g_n = interpolate(g_0, V)

#Infected individuals travel to switzerland over Zurich, Geneva and Basel Airport
source_0 = Expression(('0.0',
                       '0.54*A*exp(-b*(pow(x[0]-684202.0,2)+pow(x[1]-256909.0,2))) + '
                       '0.3*A*exp(-b*(pow(x[0]-497314.0,2)+pow(x[1]-120833.0,2))) + '
                       '0.15*A*exp(-b*(pow(x[0]-610504.0,2)+pow(x[1]-267424.0,2))) + '
                       '-1.0*(0.54*A*exp(-b*(pow(x[0]-684202.0,2)+pow(x[1]-259309.0,2))) + '
                       '0.3*A*exp(-b*(pow(x[0]-494914.0,2)+pow(x[1]-120833.0,2))) + '
                       '0.15*A*exp(-b*(pow(x[0]-608104.0,2)+pow(x[1]-267424.0,2))))',
                       '0.0'),
                      A = A_const, b = 0.0000000868, degree=2)
source = project(source_0, V)

#Diffusion in german part of switzerland is much higher than in french part and tessin
d_0 = Expression(('x[1]<206908.0? (x[1]<-2.785639*x[0]+1793034.0? d_romand: ((x[0]<743805 && x[0]>675894)? '
                  '(x[1]<0.0947945*x[0]+81841.0? d_tessin : d_deutsch) : d_deutsch)) : (x[1]<245225.0? (x[1]>1.037698*x[0]-383951.0? '
                  'd_romand : d_deutsch) : (x[1]<-0.5113657*x[0]+555276.0? d_romand : d_deutsch))',
                  'x[1]<206908.0? (x[1]<-2.785639*x[0]+1793034.0? di_romand: ((x[0]<743805 && x[0]>675894)? '
                  '(x[1]<0.0947945*x[0]+81841.0? di_tessin : di_deutsch) : di_deutsch)) : (x[1]<245225.0? (x[1]>1.037698*x[0]-383951.0? '
                  'di_romand : di_deutsch) : (x[1]<-0.5113657*x[0]+555276.0? di_romand : di_deutsch))',
                  'x[1]<206908.0? (x[1]<-2.785639*x[0]+1793034.0? d_romand: ((x[0]<743805 && x[0]>675894)? '
                  '(x[1]<0.0947945*x[0]+81841.0? d_tessin : d_deutsch) : d_deutsch)) : (x[1]<245225.0? (x[1]>1.037698*x[0]-383951.0? '
                  'd_romand : d_deutsch) : (x[1]<-0.5113657*x[0]+555276.0? d_romand : d_deutsch))'),
                  d_romand = diff_romand, d_deutsch= diff_deutsch, d_tessin = diff_tessin,
                  di_romand = diff_i_romand, di_deutsch= diff_i_deutsch, di_tessin = diff_i_tessin, degree=1)
d = project(d_0, V)

# plot bc for t=0 and safe as image
_S_0, _I_0, _R_0 = SIR_n.split()
_source_s, _source_i, _source_r = source.split()
_d_s, _d_i, _d_r = d.split()

plt.xlabel('space [x]')
plt.ylabel('space [y]')
plt.xticks([])
plt.yticks([])
plt.title('Percentage of initial population that is susceptible')
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
plt.title('Percentage of initial population that is infected')
Z = _I_0.compute_vertex_values(mesh)
c = plot(interpolate(_I_0, W), mode='color', vmin=0, vmax=0.4)
plt.tricontourf(triang, Z, vmin=0.0, vmax=0.4, levels = boundi, extend = 'both')
plt.colorbar(c)
ax = plt.gca()
circle1=plt.Circle((684202.0, 258109.0), radius, color='r', fill=False)
circle2=plt.Circle((496114.0, 120833.0), radius, color='r', fill=False)
circle3=plt.Circle((609304.0, 267424.0), radius, color='r', fill=False)
ax.add_artist(circle1)
ax.add_artist(circle2)
ax.add_artist(circle3)
plt.savefig('Videomaker/Images_I/' + str(name) + '.jpg')
plt.clf()

plt.xlabel('space [x]')
plt.ylabel('space [y]')
plt.xticks([])
plt.yticks([])
plt.title('Percentage of initial population that is recovered')
Z = _R_0.compute_vertex_values(mesh)
c = plot(interpolate(_R_0, W), mode='color', vmin=0, vmax=1.1)
plt.tricontourf(triang, Z, vmin=0.0, vmax=1.1, levels = boundr, extend = 'both')
plt.colorbar(c)
plt.savefig('Videomaker/Images_R/' + str(name) + '.jpg')
plt.clf()
name += 1

# split system functions to access components
S_low, I_low, R_low = split(SIR_low)
S_half, I_half, R_half = split(SIR_half)
S_high, I_high, R_high = split(SIR_high)
S_n, I_n, R_n = split(SIR_n)
g_n_1, g_n_2, g_n_3 = split(g_n)
g_half_1, g_half_2, g_half_3 = split(g_half)
g_1, g_2, g_3 = split(g)
d_s, d_i, d_r = split(d)
source_s, source_i, source_r = split(source)

# time stepping
n = 0
while t < T:
    # check which boundary to use at timesteps
    daytime = t - math.floor(t)
    if (daytime >= 0.25 and daytime <0.4167):
        g_n.assign(g_in)
    if (daytime < 0.25 or daytime >= 0.8333 or (daytime >= 0.4167 and daytime <0.6667)):
        g_n.assign(g_0)
    if (daytime >= 0.6667 and daytime <0.8333):
        g_n.assign(g_out)

    daytime_next = (t + dt) - math.floor((t + dt))
    if (daytime_next >= 0.25 and daytime_next < 0.4167):
        g.assign(g_in)
    if (daytime_next < 0.25 or daytime_next >= 0.8333 or (daytime_next >= 0.4167 and daytime_next < 0.6667)):
        g.assign(g_0)
    if (daytime_next >= 0.6667 and daytime_next < 0.8333):
        g.assign(g_out)

    daytime_half = (t + dt / 2.0) - math.floor((t + dt / 2.0))
    if (daytime_half >= 0.25 and daytime_half < 0.4167):
        g_half.assign(g_in)
    if (daytime_half < 0.25 or daytime_half >= 0.8333 or (daytime_half >= 0.4167 and daytime_half < 0.6667)):
        g_half.assign(g_0)
    if (daytime_half >= 0.6667 and daytime_half < 0.8333):
        g_half.assign(g_out)

    # Define variational problem for SIR_low
    F = S_low * v_1 * dx - S_n * v_1 * dx + theta * dt * (d_s * dot(grad(S_low), grad(v_1)) * dx - d_s * g_1 * v_1 * ds - \
        r_rate * v_1 * dx + d_rate_n * S_low * v_1 * dx + \
        beta * S_low * I_low / (1.0 + alpha01 * S_low + alpha02 * I_low + alpha03 * S_low * I_low) * v_1 * dx) + \
        (1.0 - theta) * dt * (d_s * dot(grad(S_n), grad(v_1)) * dx - d_s * g_n_1 * v_1 * ds - r_rate * v_1 * dx + \
        d_rate_n * S_n * v_1 * dx + \
        beta * S_n * I_n / (1.0 + alpha01 * S_n + alpha02 * I_n + alpha03 * S_n * I_n) * v_1 * dx) + \
        I_low * v_2 * dx - I_n * v_2 * dx + theta * dt * (d_i * dot(grad(I_low), grad(v_2)) * dx - d_i * g_2 * v_2 * ds - \
        beta * S_low * I_low / (1.0 + alpha01 * S_low + alpha02 * I_low + alpha03 * S_low * I_low) * v_2 * dx + \
        (d_rate_n + d_rate_d + r_d) * I_low * v_2 * dx) + (1.0 - theta) * dt * (d_i * dot(grad(I_n), grad(v_2)) * dx - \
        d_i * g_n_2 * v_2 * ds - \
        beta * S_n * I_n / (1.0 + alpha01 * S_n + alpha02 * I_n + alpha03 * S_n * I_n) * v_2 * dx + \
        (d_rate_n + d_rate_d + r_d) * I_n * v_2 * dx) + \
        R_low * v_3 * dx - R_n * v_3 * dx + theta * dt * (d_r * dot(grad(R_low), grad(v_3)) * dx - d_r * g_3 * v_3 * ds - \
        r_d * I_low * v_3 * dx + d_rate_n * R_low * v_3 * dx) + (1.0 - theta) * dt * (d_r * dot(grad(R_n), grad(v_3)) * dx - \
        d_r * g_n_3 * v_3 * ds - r_d * I_n * v_3 * dx + d_rate_n * R_n * v_3 * dx) - \
        dt * source_s * v_1 * dx - dt * source_i * v_2 * dx - dt * source_r * v_3 * dx

    Jac = derivative(F, SIR_low, TrialFunction(V))
    # solve variational problem for SIR_low
    solve(F == 0, SIR_low, J=Jac)

    # Define variational problem for SIR_half
    F = S_half * v_1 * dx - S_n * v_1 * dx + theta * dt / 2.0 * (d_s * dot(grad(S_half), grad(v_1)) * dx - d_s * g_half_1 * v_1 * ds - \
        r_rate * v_1 * dx + d_rate_n * S_half * v_1 * dx + \
        beta * S_half * I_half / (1.0 + alpha01 * S_half + alpha02 * I_half + alpha03 * S_half * I_half) * v_1 * dx) + \
        (1.0 - theta) * dt / 2.0 * (d_s * dot(grad(S_n), grad(v_1)) * dx - d_s * g_n_1 * v_1 * ds - r_rate * v_1 * dx + \
        d_rate_n * S_n * v_1 * dx + \
        beta * S_n * I_n / (1.0 + alpha01 * S_n + alpha02 * I_n + alpha03 * S_n * I_n) * v_1 * dx) + \
        I_half * v_2 * dx - I_n * v_2 * dx + theta * dt / 2.0 * (d_i * dot(grad(I_half), grad(v_2)) * dx - d_i * g_half_2 * v_2 * ds - \
        beta * S_half * I_half / (1.0 + alpha01 * S_half + alpha02 * I_half + alpha03 * S_half * I_half) * v_2 * dx + \
        (d_rate_n + d_rate_d + r_d) * I_half * v_2 * dx) + (1.0 - theta) * dt / 2.0 * (d_i * dot(grad(I_n), grad(v_2)) * dx - \
        d_i * g_n_2 * v_2 * ds - \
        beta * S_n * I_n / (1.0 + alpha01 * S_n + alpha02 * I_n + alpha03 * S_n * I_n) * v_2 * dx + \
        (d_rate_n + d_rate_d + r_d) * I_n * v_2 * dx) + \
        R_half * v_3 * dx - R_n * v_3 * dx + theta * dt / 2.0 * (d_r * dot(grad(R_half), grad(v_3)) * dx - d_r * g_half_3 * v_3 * ds - \
        r_d * I_half * v_3 * dx + d_rate_n * R_half * v_3 * dx) + (1.0 - theta) * dt / 2.0 * (d_r * dot(grad(R_n), grad(v_3)) * dx - \
        d_r * g_n_3 * v_3 * ds - r_d * I_n * v_3 * dx + d_rate_n * R_n * v_3 * dx) - \
        dt / 2.0 * source_s * v_1 * dx - dt / 2.0 * source_i * v_2 * dx - dt / 2.0 * source_r * v_3 * dx

    Jac = derivative(F, SIR_half, TrialFunction(V))
    # solve variational problem for SIR_half
    solve(F == 0, SIR_half, J=Jac)

    # Define variational problem for SIR_high
    F = S_high * v_1 * dx - S_half * v_1 * dx + theta * dt / 2.0 * (d_s * dot(grad(S_high), grad(v_1)) * dx - d_s * g_1 * v_1 * ds - \
        r_rate * v_1 * dx + d_rate_n * S_high * v_1 * dx + \
        beta * S_high * I_high / (1.0 + alpha01 * S_high + alpha02 * I_high + alpha03 * S_high * I_high) * v_1 * dx) + \
        (1.0 - theta) * dt / 2.0 * (d_s * dot(grad(S_half), grad(v_1)) * dx - d_s * g_half_1 * v_1 * ds - r_rate * v_1 * dx + \
        d_rate_n * S_half * v_1 * dx + \
        beta * S_half * I_half / (1.0 + alpha01 * S_half + alpha02 * I_half + alpha03 * S_half * I_half) * v_1 * dx) + \
        I_high * v_2 * dx - I_n * v_2 * dx + theta * dt / 2.0 * (d_i * dot(grad(I_high), grad(v_2)) * dx - d_i * g_2 * v_2 * ds - \
        beta * S_high * I_high / (1.0 + alpha01 * S_high + alpha02 * I_high + alpha03 * S_high * I_high) * v_2 * dx + \
        (d_rate_n + d_rate_d + r_d) * I_high * v_2 * dx) + (1.0 - theta) * dt / 2.0 * (d_i * dot(grad(I_half), grad(v_2)) * dx - \
        d_i * g_half_2 * v_2 * ds - \
        beta * S_half * I_half / (1.0 + alpha01 * S_half + alpha02 * I_half + alpha03 * S_half * I_half) * v_2 * dx + \
        (d_rate_n + d_rate_d + r_d) * I_half * v_2 * dx) + \
        R_high * v_3 * dx - R_half * v_3 * dx + theta * dt / 2.0 * (
        d_r * dot(grad(R_high), grad(v_3)) * dx - d_r * g_3 * v_3 * ds - \
        r_d * I_high * v_3 * dx + d_rate_n * R_high * v_3 * dx) + (1.0 - theta) * dt / 2.0 * (d_r * dot(grad(R_half), grad(v_3)) * dx - \
        d_r * g_half_3 * v_3 * ds - r_d * I_half * v_3 * dx + d_rate_n * R_half * v_3 * dx) - \
        dt / 2.0 * source_s * v_1 * dx - dt / 2.0 * source_i * v_2 * dx - dt / 2.0 * source_r * v_3 * dx

    Jac = derivative(F, SIR_high, TrialFunction(V))
    # solve variational problem for SIR_high
    solve(F == 0, SIR_high, J=Jac)

    # compute error by Richardson extrapolation
    error = errornorm(SIR_high, SIR_low, 'l2') / 3.0
    print("This is the error: ", error)

    # compute new timestep
    dt_new = pow(rho * tol / error, 0.5) * dt

    # if error small enough, update and go to the next timestep
    if error < tol:
        SIR_n.assign(SIR_high)
        n += 1
        print("This is iteration number: ", n)
        print("And the time is: ", t)
        t += dt
        timestep.append(t)

        _S, _I, _R = SIR_low.split()
        plt.xlabel('space [x]')
        plt.ylabel('space [y]')
        plt.xticks([])
        plt.yticks([])
        plt.title('Percentage of initial population that is susceptible')
        Z = _S.compute_vertex_values(mesh)
        c = plot(interpolate(_S, W), mode='color', vmin=0.0, vmax=1.1)            # , vmin=0, vmax=1.5
        plt.tricontourf(triang, Z, vmin=0.0, vmax=1.1, levels = bounds, extend = 'both')
        plt.colorbar(c)
        plt.savefig('Videomaker/Images_S/' + str(name) + '.jpg')
        plt.clf()

        plt.xlabel('space [x]')
        plt.ylabel('space [y]')
        plt.xticks([])
        plt.yticks([])
        plt.title('Percentage of initial population that is infected')
        Z = _I.compute_vertex_values(mesh)
        c = plot(interpolate(_I, W), mode='color', vmin=0.0, vmax=0.4)  # , vmin=0, vmax=1.2
        plt.tricontourf(triang, Z, vmin=0.0, vmax=0.4, levels=boundi, extend='both')
        plt.colorbar(c)
        ax = plt.gca()
        circle1 = plt.Circle((684202.0, 258109.0), radius, color='r', fill=False)
        circle2 = plt.Circle((496114.0, 120833.0), radius, color='r', fill=False)
        circle3 = plt.Circle((609304.0, 267424.0), radius, color='r', fill=False)
        ax.add_artist(circle1)
        ax.add_artist(circle2)
        ax.add_artist(circle3)
        plt.savefig('Videomaker/Images_I/' + str(name) + '.jpg')
        plt.clf()

        plt.xlabel('space [x]')
        plt.ylabel('space [y]')
        plt.xticks([])
        plt.yticks([])
        plt.title('Percentage of initial population that is recovered')
        Z = _R.compute_vertex_values(mesh)
        c = plot(interpolate(_R, W), mode='color', vmin=0.0, vmax=1.1)            # , vmin=0, vmax=2.4
        plt.tricontourf(triang, Z, vmin=0.0, vmax=1.1, levels = boundr, extend = 'both')
        plt.colorbar(c)
        plt.savefig('Videomaker/Images_R/' + str(name) + '.jpg')
        plt.clf()
        name += 1
    dt = dt_new

# print array of timesteps
print(timestep)