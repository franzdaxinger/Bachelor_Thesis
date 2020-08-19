
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

# get data
data = pd.read_excel (r'shapefiles/diffusion.xlsx', sheet_name='Commune of residence perspect.',
                      usecols=['Canton number res', 'Canton number work', 'Number of employed persons'])
array_data = data.to_numpy()
print(array_data)

# write [canton i, canton j and commuters from i to j or j to i] to array
commuters = []
array_idhome = []
array_idwork = []
array_numcom = []
ID_home = 1
while ID_home < 27:         # 27
    ID_work = 1
    while ID_work < 27:     # 27
        number_comm = 0.0
        if ID_home < ID_work:
            for i in range(len(array_data)):
                if (array_data[i][0] == ID_home and array_data[i][1] == ID_work) or (array_data[i][1] == ID_home and array_data[i][0] == ID_work):
                    number_comm += array_data[i][2]
            commuters.append([ID_home, ID_work, number_comm])
            array_idhome.append(ID_home)
            array_idwork.append(ID_work)
            array_numcom.append(number_comm)
        ID_work += 1
    print("canton ", ID_home, " done.")
    ID_home += 1

print(commuters)

# save to txt files in folder shapefiles
with open('shapefiles/cant_commuters.txt', 'w') as filehandle:
    for listitem in commuters:
        filehandle.write('%s\n' % listitem)

with open('shapefiles/array_idhome.txt', 'w') as filehandle:
    for listitem in array_idhome:
        filehandle.write('%s\n' % listitem)

with open('shapefiles/array_idwork.txt', 'w') as filehandle:
    for listitem in array_idwork:
        filehandle.write('%s\n' % listitem)

with open('shapefiles/array_numcom.txt', 'w') as filehandle:
    for listitem in array_numcom:
        filehandle.write('%s\n' % listitem)