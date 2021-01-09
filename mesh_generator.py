# Creates a mesh and saves it in the folder mesh

from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from fenics import *
from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
from mshr import *

# read border points from txt files
x_border = np.genfromtxt('data/switzerland_x.txt')
y_border = np.genfromtxt('data/switzerland_y.txt')

# generate cartesian mesh with mshr
print("Generating the mesh now... please wait")
mesh2d = RectangleMesh(Point(480000.0, 70000.0), Point(840000.0, 300000.0), 40, 30, diagonal="right")
print("Mesh generation done")


# save mesh and plot in folder mesh
plt.xticks([])
plt.yticks([])
plt.title('mesh')
plt.plot(x_border, y_border, color='r')
plot(mesh2d)
plt.savefig('mesh/geometry.jpg')
plt.clf()

# save the mesh as file
print("Saving file ...")
File('mesh/mesh2d.xml.gz') << mesh2d
print("Done")