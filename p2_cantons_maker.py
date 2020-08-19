
# This program calculates the canton functions, where cantonfunction i is 1 in canton i and 0 everywhere else


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
switzerland = gpd.read_file('shapefiles/cantons.shp')
switzerland_new = switzerland.translate(xoff=-2000000.0, yoff=-1000000.0, zoff=0.0)
print(switzerland_new.head())

# get mesh from file in folder mesh
mesh = Mesh('mesh/mesh2d.xml.gz')

# define function space for system
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1])
V = FunctionSpace(mesh, element)
W = FunctionSpace(mesh, P1)

# checks for a Point (x,y) in which canton it is located. If in none return -1
def containsPoint(x, y):
    ID_KT = 0
    point = Point(x, y)

    while ID_KT < 26:         #26
        canton = switzerland_new[ID_KT]
        containing = canton.contains(point)
        if containing == True:
            return ID_KT
        ID_KT += 1
    return -1

# check for every point in which canton it is and set values 0 or 1
class MyExpression0(UserExpression):
    def eval(self, value, x):
        a = containsPoint(x[0], x[1])
        value[0] = 0.0
        value[1] = 0.0
        value[2] = 0.0
        value[3] = 0.0
        value[4] = 0.0
        value[5] = 0.0
        value[6] = 0.0
        value[7] = 0.0
        value[8] = 0.0
        value[9] = 0.0
        value[10] = 0.0
        value[11] = 0.0
        value[12] = 0.0
        value[13] = 0.0
        value[14] = 0.0
        value[15] = 0.0
        value[16] = 0.0
        value[17] = 0.0
        value[18] = 0.0
        value[19] = 0.0
        value[20] = 0.0
        value[21] = 0.0
        value[22] = 0.0
        value[23] = 0.0
        value[24] = 0.0
        value[25] = 0.0
        if a != -1:
            value[a] = 1.0
    def value_shape(self):
        return (26,)
g_ = MyExpression0()
g = interpolate(g_, V)

_g1, _g2, _g3, _g4, _g5, _g6, _g7, _g8, _g9, _g10, _g11, _g12, _g13, _g14, _g15, _g16, _g17, _g18, _g19, _g20, _g21, _g22, _g23, _g24, _g25, _g26 = g.split()

# plot functions
if plotbool == True:
    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('01')
    c = plot(_g1)
    plt.colorbar(c)
    plt.savefig('c01.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('02')
    c = plot(_g2)
    plt.colorbar(c)
    plt.savefig('c02.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('03')
    c = plot(_g3)
    plt.colorbar(c)
    plt.savefig('c03.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('04')
    c = plot(_g4)
    plt.colorbar(c)
    plt.savefig('c04.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('05')
    c = plot(_g5)
    plt.colorbar(c)
    plt.savefig('c05.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('06')
    c = plot(_g6)
    plt.colorbar(c)
    plt.savefig('c06.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('07')
    c = plot(_g7)
    plt.colorbar(c)
    plt.savefig('c07.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('08')
    c = plot(_g8)
    plt.colorbar(c)
    plt.savefig('c08.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('09')
    c = plot(_g9)
    plt.colorbar(c)
    plt.savefig('c09.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('10')
    c = plot(_g10)
    plt.colorbar(c)
    plt.savefig('c10.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('11')
    c = plot(_g11)
    plt.colorbar(c)
    plt.savefig('c11.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('12')
    c = plot(_g12)
    plt.colorbar(c)
    plt.savefig('c12.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('13')
    c = plot(_g13)
    plt.colorbar(c)
    plt.savefig('c13.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('14')
    c = plot(_g14)
    plt.colorbar(c)
    plt.savefig('c14.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('15')
    c = plot(_g15)
    plt.colorbar(c)
    plt.savefig('c15.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('16')
    c = plot(_g16)
    plt.colorbar(c)
    plt.savefig('c16.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('17')
    c = plot(_g17)
    plt.colorbar(c)
    plt.savefig('c17.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('18')
    c = plot(_g18)
    plt.colorbar(c)
    plt.savefig('c18.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('19')
    c = plot(_g19)
    plt.colorbar(c)
    plt.savefig('c19.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('20')
    c = plot(_g20)
    plt.colorbar(c)
    plt.savefig('c20.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('21')
    c = plot(_g21)
    plt.colorbar(c)
    plt.savefig('c21.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('22')
    c = plot(_g22)
    plt.colorbar(c)
    plt.savefig('c22.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('23')
    c = plot(_g23)
    plt.colorbar(c)
    plt.savefig('c23.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('24')
    c = plot(_g24)
    plt.colorbar(c)
    plt.savefig('c24.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('25')
    c = plot(_g25)
    plt.colorbar(c)
    plt.savefig('c25.jpg')
    plt.clf()

    plt.xlabel('space [x]')
    plt.ylabel('space [y]')
    plt.title('26')
    c = plot(_g26)
    plt.colorbar(c)
    plt.savefig('c26.jpg')
    plt.clf()

g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21, g22, g23, g24, g25, g26  = g.split()

# save functions in folder cantfun
if savebool == True:
    f_out = XDMFFile("cantfun/01.xdmf")
    f_out.write_checkpoint(project(g1, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/02.xdmf")
    f_out.write_checkpoint(project(g2, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/03.xdmf")
    f_out.write_checkpoint(project(g3, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/04.xdmf")
    f_out.write_checkpoint(project(g4, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/05.xdmf")
    f_out.write_checkpoint(project(g5, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/06.xdmf")
    f_out.write_checkpoint(project(g6, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/07.xdmf")
    f_out.write_checkpoint(project(g7, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/08.xdmf")
    f_out.write_checkpoint(project(g8, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/09.xdmf")
    f_out.write_checkpoint(project(g9, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/10.xdmf")
    f_out.write_checkpoint(project(g10, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/11.xdmf")
    f_out.write_checkpoint(project(g11, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/12.xdmf")
    f_out.write_checkpoint(project(g12, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/13.xdmf")
    f_out.write_checkpoint(project(g13, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/14.xdmf")
    f_out.write_checkpoint(project(g14, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/15.xdmf")
    f_out.write_checkpoint(project(g15, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/16.xdmf")
    f_out.write_checkpoint(project(g16, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/17.xdmf")
    f_out.write_checkpoint(project(g17, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/18.xdmf")
    f_out.write_checkpoint(project(g18, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/19.xdmf")
    f_out.write_checkpoint(project(g19, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/20.xdmf")
    f_out.write_checkpoint(project(g20, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/21.xdmf")
    f_out.write_checkpoint(project(g21, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/22.xdmf")
    f_out.write_checkpoint(project(g22, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/23.xdmf")
    f_out.write_checkpoint(project(g23, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/24.xdmf")
    f_out.write_checkpoint(project(g24, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/25.xdmf")
    f_out.write_checkpoint(project(g25, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()

    f_out = XDMFFile("cantfun/26.xdmf")
    f_out.write_checkpoint(project(g26, W), "g", 0, XDMFFile.Encoding.HDF5, True)  # appending to file
    f_out.close()
