"""
Model for bayesian inference with korali

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
import pandas as pd
from shapely.geometry import Point
import numpy as np
import math


def model(s, X):
    beta_param = s["Parameters"][0]
    sig = s["Parameters"][1]
    s["Reference Evaluations"] = []
    s["Dispersion"] = []
    result = []

    # define all needed parameters
    T = len(X) - 1.0                            # final time in days
    dt = 4.0 / 24.0                             # time step size at beginning
    theta_factor = Constant(1.1)                # factor to represent the underreporting of movement
    beta_factor = beta_param                    # infection rate
    oneoverd = Constant(1.0 / 5.0)              # one over average duration of infection
    oneoverz = Constant(1.0 / 5.0)              # one over average latency period

    theta = Constant(0.5)                       # theta = 0.5 means Crank-Nicolson
    t = 0.0                                     # global time
    area_switzerland = 41285000000.0

    # get mesh from file in folder mesh
    mesh = Mesh('mesh/mesh2d.xml.gz')

    # get data on commuters between cantons and write to array
    array_alpha = np.genfromtxt('shapefiles/alpha.txt')
    array_beta = np.genfromtxt('shapefiles/beta.txt')

    # define function space for system
    P1 = FiniteElement('P', triangle, 1)
    element = MixedElement([P1, P1, P1])
    V = FunctionSpace(mesh, element)
    W = FunctionSpace(mesh, P1)

    # define test functions
    v_1, v_2, v_3 = TestFunctions(V)

    # define used functions
    SEI_low = Function(V)
    SEI_n = Function(V)
    d = Function(V)
    source_s_n = Function(W)
    source_e_n = Function(W)
    source_i_n = Function(W)

    # create files and parameters for visualization output and time
    name = 100000

    # setting Initial Conditions
    # 0.01 * exp(-0.00000039*(pow(x[0]-720000.0,2)+pow(x[1]-130000.0,2)))
    SEI_0 = Expression(('1.0',
                        '0.0',
                        '0.01 * exp(-0.00000039*(pow(x[0]-720000.0,2)+pow(x[1]-130000.0,2)))'),
                       degree=2)
    SEI_n = project(SEI_0, V)

    # Diffusion, inverse population density and beta from files
    f1 = Function(W)
    f_in = XDMFFile("difffun/diff.xdmf")
    f_in.read_checkpoint(f1, "g", 0)
    f1.set_allow_extrapolation(True)
    d_0 = Expression(('function',
                      'function',
                      'function'),
                     function=f1, degree=2)
    d = project(d_0, V)

    f2 = Function(W)
    f_in = XDMFFile("difffun/rhoinv.xdmf")
    f_in.read_checkpoint(f2, "g", 0)
    f2.set_allow_extrapolation(True)
    rhoinv_0 = Expression(('function',
                           'function',
                           'function'),
                          function=f2, degree=2)
    rhoinv = project(rhoinv_0, V)

    f3 = Function(W)
    f_in = XDMFFile("difffun/beta.xdmf")
    f_in.read_checkpoint(f3, "g", 0)
    f3.set_allow_extrapolation(True)
    beta = project(beta_factor * f3, W)

    # check diffusion, 1/rho and beta and plot S, E, I for t=0 and safe as images
    _S_0, _E_0, _I_0 = SEI_n.split()
    _d_s, _d_i, _d_r = d.split()
    _rhoinv1, _rhoinv2, _rhoinv3 = rhoinv.split()

    f_out = XDMFFile("Videomaker/functions/checkdiff.xdmf")
    f_out.write_checkpoint(project(_d_s, W), "diff", 0, XDMFFile.Encoding.HDF5, True)
    f_out.close()

    f_out = XDMFFile("Videomaker/functions/checkrho.xdmf")
    f_out.write_checkpoint(project(_rhoinv1, W), "rho", 0, XDMFFile.Encoding.HDF5, True)
    f_out.close()

    f_out = XDMFFile("Videomaker/functions/checkbeta.xdmf")
    f_out.write_checkpoint(project(beta, W), "beta", 0, XDMFFile.Encoding.HDF5, True)
    f_out.close()

    f_out = XDMFFile("Videomaker/functions/function_S.xdmf")
    f_out.write_checkpoint(project(_S_0, W), "S", name, XDMFFile.Encoding.HDF5, True)
    f_out.close()

    f_out = XDMFFile("Videomaker/functions/function_E.xdmf")
    f_out.write_checkpoint(project(_E_0, W), "E", name, XDMFFile.Encoding.HDF5, True)
    f_out.close()

    f_out = XDMFFile("Videomaker/functions/function_I.xdmf")
    f_out.write_checkpoint(project(_I_0, W), "I", name, XDMFFile.Encoding.HDF5, True)
    f_out.close()

    result.append(assemble(_I_0 * dx) / area_switzerland)

    name += 1

    # split system functions to access components
    S_low, E_low, I_low = split(SEI_low)
    S_n, E_n, I_n = split(SEI_n)
    d_s, d_e, d_i = split(d)
    rhoinv_s, rhoinv_e, rhoinv_i = split(rhoinv)

    # get functions and areas to represent cantons
    cantonfun = [0.0]
    if 1 == 1:
        c1 = Function(W)
        f_in = XDMFFile("cantfun/01.xdmf")
        f_in.read_checkpoint(c1, "g", 0)
        c1.set_allow_extrapolation(True)

        c2 = Function(W)
        f_in = XDMFFile("cantfun/02.xdmf")
        f_in.read_checkpoint(c2, "g", 0)
        c2.set_allow_extrapolation(True)

        c3 = Function(W)
        f_in = XDMFFile("cantfun/03.xdmf")
        f_in.read_checkpoint(c3, "g", 0)
        c3.set_allow_extrapolation(True)

        c4 = Function(W)
        f_in = XDMFFile("cantfun/04.xdmf")
        f_in.read_checkpoint(c4, "g", 0)
        c4.set_allow_extrapolation(True)

        c5 = Function(W)
        f_in = XDMFFile("cantfun/05.xdmf")
        f_in.read_checkpoint(c5, "g", 0)
        c5.set_allow_extrapolation(True)

        c6 = Function(W)
        f_in = XDMFFile("cantfun/06.xdmf")
        f_in.read_checkpoint(c6, "g", 0)
        c6.set_allow_extrapolation(True)

        c7 = Function(W)
        f_in = XDMFFile("cantfun/07.xdmf")
        f_in.read_checkpoint(c7, "g", 0)
        c7.set_allow_extrapolation(True)

        c8 = Function(W)
        f_in = XDMFFile("cantfun/08.xdmf")
        f_in.read_checkpoint(c8, "g", 0)
        c8.set_allow_extrapolation(True)

        c9 = Function(W)
        f_in = XDMFFile("cantfun/09.xdmf")
        f_in.read_checkpoint(c9, "g", 0)
        c9.set_allow_extrapolation(True)

        c10 = Function(W)
        f_in = XDMFFile("cantfun/10.xdmf")
        f_in.read_checkpoint(c10, "g", 0)
        c10.set_allow_extrapolation(True)

        c11 = Function(W)
        f_in = XDMFFile("cantfun/11.xdmf")
        f_in.read_checkpoint(c11, "g", 0)
        c11.set_allow_extrapolation(True)

        c12 = Function(W)
        f_in = XDMFFile("cantfun/12.xdmf")
        f_in.read_checkpoint(c12, "g", 0)
        c12.set_allow_extrapolation(True)

        c13 = Function(W)
        f_in = XDMFFile("cantfun/13.xdmf")
        f_in.read_checkpoint(c13, "g", 0)
        c13.set_allow_extrapolation(True)

        c14 = Function(W)
        f_in = XDMFFile("cantfun/14.xdmf")
        f_in.read_checkpoint(c14, "g", 0)
        c14.set_allow_extrapolation(True)

        c15 = Function(W)
        f_in = XDMFFile("cantfun/15.xdmf")
        f_in.read_checkpoint(c15, "g", 0)
        c15.set_allow_extrapolation(True)

        c16 = Function(W)
        f_in = XDMFFile("cantfun/16.xdmf")
        f_in.read_checkpoint(c16, "g", 0)
        c16.set_allow_extrapolation(True)

        c17 = Function(W)
        f_in = XDMFFile("cantfun/17.xdmf")
        f_in.read_checkpoint(c17, "g", 0)
        c17.set_allow_extrapolation(True)

        c18 = Function(W)
        f_in = XDMFFile("cantfun/18.xdmf")
        f_in.read_checkpoint(c18, "g", 0)
        c18.set_allow_extrapolation(True)

        c19 = Function(W)
        f_in = XDMFFile("cantfun/19.xdmf")
        f_in.read_checkpoint(c19, "g", 0)
        c19.set_allow_extrapolation(True)

        c20 = Function(W)
        f_in = XDMFFile("cantfun/20.xdmf")
        f_in.read_checkpoint(c20, "g", 0)
        c20.set_allow_extrapolation(True)

        c21 = Function(W)
        f_in = XDMFFile("cantfun/21.xdmf")
        f_in.read_checkpoint(c21, "g", 0)
        c21.set_allow_extrapolation(True)

        c22 = Function(W)
        f_in = XDMFFile("cantfun/22.xdmf")
        f_in.read_checkpoint(c22, "g", 0)
        c22.set_allow_extrapolation(True)

        c23 = Function(W)
        f_in = XDMFFile("cantfun/23.xdmf")
        f_in.read_checkpoint(c23, "g", 0)
        c23.set_allow_extrapolation(True)

        c24 = Function(W)
        f_in = XDMFFile("cantfun/24.xdmf")
        f_in.read_checkpoint(c24, "g", 0)
        c24.set_allow_extrapolation(True)

        c25 = Function(W)
        f_in = XDMFFile("cantfun/25.xdmf")
        f_in.read_checkpoint(c25, "g", 0)
        c25.set_allow_extrapolation(True)

        c26 = Function(W)
        f_in = XDMFFile("cantfun/26.xdmf")
        f_in.read_checkpoint(c26, "g", 0)
        c26.set_allow_extrapolation(True)

        cantonfun.append(c1)
        cantonfun.append(c2)
        cantonfun.append(c3)
        cantonfun.append(c4)
        cantonfun.append(c5)
        cantonfun.append(c6)
        cantonfun.append(c7)
        cantonfun.append(c8)
        cantonfun.append(c9)
        cantonfun.append(c10)
        cantonfun.append(c11)
        cantonfun.append(c12)
        cantonfun.append(c13)
        cantonfun.append(c14)
        cantonfun.append(c15)
        cantonfun.append(c16)
        cantonfun.append(c17)
        cantonfun.append(c18)
        cantonfun.append(c19)
        cantonfun.append(c20)
        cantonfun.append(c21)
        cantonfun.append(c22)
        cantonfun.append(c23)
        cantonfun.append(c24)
        cantonfun.append(c25)
        cantonfun.append(c26)

    # time stepping
    n = 0
    while t < T:
        # Define source term for long connections
        source_s_0 = Expression('0.0', degree=2)
        source_s_n = project(source_s_0, W)
        source_e_0 = Expression('0.0', degree=2)
        source_e_n = project(source_e_0, W)
        source_i_0 = Expression('0.0', degree=2)
        source_i_n = project(source_i_0, W)

        # compute number of susceptible, exposed and infected in every canton in advance and save as array
        array_S_n = [0.0]
        array_E_n = [0.0]
        array_I_n = [0.0]

        array_factor_S_n = [0.0]
        array_factor_E_n = [0.0]
        array_factor_I_n = [0.0]

        k = 1
        while k < 27:
            array_S_n.append(assemble(S_n * cantonfun[k] * dx))
            array_E_n.append(assemble(E_n * cantonfun[k] * dx))
            array_I_n.append(assemble(I_n * cantonfun[k] * dx))
            k += 1

        # calculate source terms for current time, low dt, half dt and high dt
        print("initializing sources now ...")
        ID_KT_i = 1
        while ID_KT_i < 27:
            ID_KT_j = 1
            factor_s_n = 0.0
            factor_e_n = 0.0
            factor_i_n = 0.0
            while ID_KT_j < 27:
                index = (ID_KT_i - 1) * 26 + ID_KT_j - 1
                factor_s_n += array_alpha[index] * array_S_n[ID_KT_j] - array_beta[index] * array_S_n[ID_KT_i]
                factor_e_n += array_alpha[index] * array_E_n[ID_KT_j] - array_beta[index] * array_E_n[ID_KT_i]
                factor_i_n += array_alpha[index] * array_I_n[ID_KT_j] - array_beta[index] * array_I_n[ID_KT_i]

                ID_KT_j += 1
            array_factor_S_n.append(factor_s_n)
            array_factor_E_n.append(factor_e_n)
            array_factor_I_n.append(factor_i_n)

            ID_KT_i += 1

        coeff_s = array_factor_S_n
        coeff_e = array_factor_E_n
        coeff_i = array_factor_I_n
        source_s_n = project(source_s_n + coeff_s[1] * c1 + coeff_s[2] * c2 + \
                             coeff_s[3] * c3 + coeff_s[4] * c4 + coeff_s[5] * c5 + \
                             coeff_s[6] * c6 + coeff_s[7] * c7 + coeff_s[8] * c8 + \
                             coeff_s[9] * c9 + coeff_s[10] * c10 + coeff_s[11] * c11 + \
                             coeff_s[12] * c12 + coeff_s[13] * c13 + coeff_s[14] * c14 + \
                             coeff_s[15] * c15 + coeff_s[16] * c16 + coeff_s[17] * c17 + \
                             coeff_s[18] * c18 + coeff_s[19] * c19 + coeff_s[20] * c20 + \
                             coeff_s[21] * c21 + coeff_s[22] * c22 + coeff_s[23] * c23 + \
                             coeff_s[24] * c24 + coeff_s[25] * c25 + coeff_s[26] * c26)

        source_e_n = project(source_e_n + coeff_e[1] * c1 + coeff_e[2] * c2 + \
                             coeff_e[3] * c3 + coeff_e[4] * c4 + coeff_e[5] * c5 + \
                             coeff_e[6] * c6 + coeff_e[7] * c7 + coeff_e[8] * c8 + \
                             coeff_e[9] * c9 + coeff_e[10] * c10 + coeff_e[11] * c11 + \
                             coeff_e[12] * c12 + coeff_e[13] * c13 + coeff_e[14] * c14 + \
                             coeff_e[15] * c15 + coeff_e[16] * c16 + coeff_e[17] * c17 + \
                             coeff_e[18] * c18 + coeff_e[19] * c19 + coeff_e[20] * c20 + \
                             coeff_e[21] * c21 + coeff_e[22] * c22 + coeff_e[23] * c23 + \
                             coeff_e[24] * c24 + coeff_e[25] * c25 + coeff_e[26] * c26)

        source_i_n = project(source_i_n + coeff_i[1] * c1 + coeff_i[2] * c2 + \
                             coeff_i[3] * c3 + coeff_i[4] * c4 + coeff_i[5] * c5 + \
                             coeff_i[6] * c6 + coeff_i[7] * c7 + coeff_i[8] * c8 + \
                             coeff_i[9] * c9 + coeff_i[10] * c10 + coeff_i[11] * c11 + \
                             coeff_i[12] * c12 + coeff_i[13] * c13 + coeff_i[14] * c14 + \
                             coeff_i[15] * c15 + coeff_i[16] * c16 + coeff_i[17] * c17 + \
                             coeff_i[18] * c18 + coeff_i[19] * c19 + coeff_i[20] * c20 + \
                             coeff_i[21] * c21 + coeff_i[22] * c22 + coeff_i[23] * c23 + \
                             coeff_i[24] * c24 + coeff_i[25] * c25 + coeff_i[26] * c26)

        print("source_n done ...")

        # Define variational problem for SEI_low
        F = S_low * v_1 * dx - S_n * v_1 * dx - theta * dt * (-beta * S_low * I_low * v_1 * dx + \
            theta_factor * (- rhoinv_s * d_s * dot(grad(S_low), grad(v_1)) * dx + rhoinv_s * source_s_n * v_1 * dx)) - \
            (1.0 - theta) * dt * (-beta * S_n * I_n * v_1 * dx + theta_factor * (- rhoinv_s * d_s * dot(grad(S_n), grad(v_1)) * dx + rhoinv_s * source_s_n * v_1 * dx)) + \
            E_low * v_2 * dx - E_n * v_2 * dx - theta * dt * (beta * S_low * I_low * v_2 * dx - oneoverz * E_low * v_2 * dx + \
            theta_factor * (- rhoinv_e * d_e * dot(grad(E_low), grad(v_2)) * dx + rhoinv_e * source_e_n * v_2 * dx)) - \
            (1.0 - theta) * dt * (beta * S_n * I_n * v_2 * dx - oneoverz * E_n * v_2 * dx + \
            theta_factor * (- rhoinv_e * d_e * dot(grad(E_n), grad(v_2)) * dx + rhoinv_e * source_e_n * v_2 * dx)) + \
            I_low * v_3 * dx - I_n * v_3 * dx - theta * dt * (oneoverz * E_low * v_3 * dx - oneoverd * I_low * v_3 * dx + \
            theta_factor * (- rhoinv_i * d_i * dot(grad(I_low), grad(v_3)) * dx + rhoinv_i * source_i_n * v_3 * dx)) - \
            (1.0 - theta) * dt * (oneoverz * E_n * v_3 * dx - oneoverd * I_n * v_3 * dx + \
            theta_factor * (- rhoinv_i * d_i * dot(grad(I_n), grad(v_3)) * dx + rhoinv_i * source_i_n * v_3 * dx))

        Jac = derivative(F, SEI_low, TrialFunction(V))

        # solve variational problem for SEI_low
        solve(F == 0, SEI_low, J=Jac)

        n += 1
        t += dt
        print("This is iteration number: ", n)
        print("And the time is: ", t)

        # if error small enough or minimum timestep reached, update and go to the next timestep
        if n % 6 == 0:
            _S, _E, _I = SEI_low.split()
            f_out = XDMFFile("Videomaker/functions/function_S.xdmf")
            f_out.write_checkpoint(project(_S, W), "S", name, XDMFFile.Encoding.HDF5, True)
            f_out.close()

            f_out = XDMFFile("Videomaker/functions/function_E.xdmf")
            f_out.write_checkpoint(project(_E, W), "E", name, XDMFFile.Encoding.HDF5, True)
            f_out.close()

            f_out = XDMFFile("Videomaker/functions/function_I.xdmf")
            f_out.write_checkpoint(project(_I, W), "I", name, XDMFFile.Encoding.HDF5, True)
            f_out.close()

            result.append(assemble(_I * dx) / area_switzerland)

            name += 1

        SEI_n.assign(SEI_low)

    dev = []
    for x in X:
        dev.append(sig)
    s["Reference Evaluations"] = result
    s["Dispersion"] = dev


# days until 16th of march, when the lockdown happened
def getReferencePoints():
    x = []
    time = 0.0
    while time < 22.0:
        x.append(time)
        time += 1.0
    return x

def getReferenceData():
    y = []
    bev = 8570000.0
    data_array = [1, 2, 12, 22, 32, 45, 56, 87, 120, 181, 243, 316, 365, 434, 626, 838, 1172, 1529, 1962, 2382, 2710, 3774]
    for i in data_array:
        y.append(i/bev)
    return y