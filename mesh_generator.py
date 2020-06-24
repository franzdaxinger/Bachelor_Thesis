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

# read border points from excel file
data = pd.read_excel("test.xlsx")
print("Reading data done")
arr = data.to_numpy()
x,y = zip(*arr)

# plots used border points
plt.xlabel('space [x]')
plt.ylabel('space [y]')
plt.scatter(x,y)
plt.title('Switzerland')
plt.savefig('switzerland.jpg')
plt.clf()

# generate array containing the data as Points
arr_points = [Point(0, 0)] * (len(x))
i = 0
for i in range(len(x)):
    arr_points[i] = Point(x[i], y[i])
    i +=1

# generate mesh with mshr from polygon of border points
print("Generating the mesh now... please wait")
domain = Polygon(arr_points)
mesh2d = generate_mesh(domain,70)
print("Mesh generation done")

# save mesh and plot in folder mesh
plot(mesh2d, "2D mesh")
plt.savefig('mesh/geometry.pdf')

print("Saving file ...")
File('mesh/mesh2d.xml.gz') << mesh2d
print("Done")