
# Creates a mesh and saves it in the folder mesh

from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from fenics import *
from dolfin import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import pandas as pd
from mpl_toolkits import mplot3d
import math
from mshr import *
import xlrd

# read border points from txt files
x_border = np.genfromtxt('shapefiles/switzerland_x.txt')
y_border = np.genfromtxt('shapefiles/switzerland_y.txt')

# generate cartesian mesh with mshr
print("Generating the mesh now... please wait")
mesh2d = RectangleMesh(Point(480000.0, 70000.0), Point(840000.0, 300000.0), 20, 15, diagonal="right")
print("Mesh generation done")


# save mesh and plot in folder mesh
plt.xlabel('space [x]')
plt.ylabel('space [y]')
plt.title('mesh')
plt.plot(x_border, y_border, color='r')
plot(mesh2d)
plt.savefig('mesh2/mesh/geometry.jpg')
plt.clf()

# save the mesh as file
print("Saving file ...")
File('mesh2/mesh/mesh2d.xml.gz') << mesh2d
print("Done")