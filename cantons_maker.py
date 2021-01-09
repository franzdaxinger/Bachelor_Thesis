# get txt files with commuters from canton to canton

from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from fenics import *
from dolfin import *
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.tri as tri
import shapefile
import os
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon
import numpy as np
import math
from shapely.geometry import mapping
import shapely

# get mesh from file in folder mesh
mesh = Mesh('mesh/mesh2d.xml.gz')

# define function space for system
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1, P1])
V = FunctionSpace(mesh, element)
W = FunctionSpace(mesh, P1)

# get data
data = pd.read_excel (r'data/diffusion.xlsx', sheet_name='Commune of residence perspect.',
                      usecols=['Canton number res', 'Canton number work', 'Number of employed persons'])
array_data = data.to_numpy()
number_cant = 27

# get functions and areas to represent cantons
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
    cant_area = [0.0]
    cant_area.append(assemble(c1 * dx))
    cant_area.append(assemble(c2 * dx))
    cant_area.append(assemble(c3 * dx))
    cant_area.append(assemble(c4 * dx))
    cant_area.append(assemble(c5 * dx))
    cant_area.append(assemble(c6 * dx))
    cant_area.append(assemble(c7 * dx))
    cant_area.append(assemble(c8 * dx))
    cant_area.append(assemble(c9 * dx))
    cant_area.append(assemble(c10 * dx))
    cant_area.append(assemble(c11 * dx))
    cant_area.append(assemble(c12 * dx))
    cant_area.append(assemble(c13 * dx))
    cant_area.append(assemble(c14 * dx))
    cant_area.append(assemble(c15 * dx))
    cant_area.append(assemble(c16 * dx))
    cant_area.append(assemble(c17 * dx))
    cant_area.append(assemble(c18 * dx))
    cant_area.append(assemble(c19 * dx))
    cant_area.append(assemble(c20 * dx))
    cant_area.append(assemble(c21 * dx))
    cant_area.append(assemble(c22 * dx))
    cant_area.append(assemble(c23 * dx))
    cant_area.append(assemble(c24 * dx))
    cant_area.append(assemble(c25 * dx))
    cant_area.append(assemble(c26 * dx))

# write [canton i, canton j and commuters from i to j or j to i] to array
alpha = [0.0] * (26*26)
beta = [0.0] * (26*26)
ID_home = 1
while ID_home < number_cant:
    ID_work = 1
    while ID_work < number_cant:
        number_comm = 0.0
        if ID_home < ID_work:
            for i in range(len(array_data)):
                if (array_data[i][0] == ID_home and array_data[i][1] == ID_work) or (array_data[i][1] == ID_home and array_data[i][0] == ID_work):
                    number_comm += array_data[i][2]
            index1 = (ID_home - 1) * (number_cant - 1) + ID_work - 1
            index2 = (ID_work - 1) * (number_cant - 1) + ID_home - 1
            alpha[index1] = number_comm / (cant_area[ID_home] * cant_area[ID_work])
            alpha[index2] = number_comm / (cant_area[ID_home] * cant_area[ID_work])
            beta[index1] = number_comm / (cant_area[ID_home] * cant_area[ID_home])
            beta[index2] = number_comm / (cant_area[ID_work] * cant_area[ID_work])

        ID_work += 1
    print("canton ", ID_home, " done.")
    ID_home += 1


# save to txt files in folder data
with open('data/alpha.txt', 'w') as filehandle:
    for listitem in alpha:
        filehandle.write('%s\n' % listitem)

with open('data/beta.txt', 'w') as filehandle:
    for listitem in beta:
        filehandle.write('%s\n' % listitem)