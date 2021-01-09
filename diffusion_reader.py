# get txt files with diffusion coefficient and 1 / population density in every municipality

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
from shapely.geometry import Point
import numpy as np
import math

# read in inhabitants per municipality
data_inhabitants = pd.read_excel (r'data/Gemeinden.xlsx.xls', sheet_name='Schweiz - Gemeinden',
                      usecols="A,C", skiprows = 8)
data_inhabitants.columns = ['ID_Gem','inhabitants']
print("processing data")
inhabitants = data_inhabitants['inhabitants'].to_numpy()
inhabitants = inhabitants[np.logical_not(np.isnan(inhabitants))]

# initialize arrays and params to save data
length = 2337
diffusion = []
densityinv = []
num_municipality = 0

data = pd.read_excel(r'data/diffusion.xlsx', sheet_name='Commune of residence perspect.',
                      usecols=['Commune number res', 'Commune number work', 'Number of employed persons'])
array_data = data.to_numpy()
municipalities = gpd.read_file('data/municipalities.shp')
id_mun = municipalities['ID_Gem'].to_numpy()

# get array with area of every municipality
area = municipalities.area.to_numpy()

# for every municipality check which are neighboring municipalities and compute diffusion to them
num_municipality = 0
while num_municipality < length:
    number = 0
    diffusion_coefficient = 0.0

    someSpecialGeometry = municipalities['geometry'].values[num_municipality]
    municipalities['isNeighbor'] = municipalities.apply(lambda row: row['geometry'].touches(someSpecialGeometry), axis=1)

    id_gem = municipalities['ID_Gem'].values[num_municipality]
    num_inhabitants = -1
    for o in range(len(data_inhabitants)):
        if data_inhabitants['ID_Gem'].values[o] == id_gem:
            num_inhabitants = o
    rhoinverse = 4830.917874        # take average pop density inverse for geometries with no municipality ID
    if num_inhabitants != -1:
        rhoinverse = area[num_municipality] / inhabitants[num_inhabitants]

    for n in range(len(array_data)):
        id_home = array_data[n][0]
        id_work = array_data[n][1]
        commuters_mun = array_data[n][2]
        num_work = 0

        if (id_home == id_gem) and (id_work != 7777):
            for k in range(len(municipalities)):
                if municipalities['ID_Gem'].values[k] == id_work:
                    num_work = k

            if municipalities['isNeighbor'].values[num_work] == True:
                mun1 = municipalities.geometry[num_municipality]
                mun2 = municipalities.geometry[num_work]
                borderline = mun1.intersection(mun2)
                blength = borderline.length
                if blength != 0:
                    number += 1
                    for i in range(len(array_data)):
                        if (array_data[i][0] == num_work) and (array_data[i][1] == num_municipality):
                            commuters_mun += array_data[i][2]
                            break
                    xdiff = municipalities.geometry.centroid.x.values[num_municipality] - \
                            municipalities.geometry.centroid.x.values[num_work]
                    ydiff = municipalities.geometry.centroid.y.values[num_municipality] - \
                            municipalities.geometry.centroid.y.values[num_work]
                    dist = math.sqrt(pow(xdiff, 2) + pow(ydiff, 2))
                    diffusion_coefficient += commuters_mun * dist / blength

    if number != 0:
        diffusion_coefficient = diffusion_coefficient / number
    diffusion.append(diffusion_coefficient)
    densityinv.append(rhoinverse)
    print("municipality ", num_municipality, " done.")
    num_municipality += 1

print("diffusion: ", diffusion)
print("densityinverse: ", densityinv)

# save as txt files in data
with open('data/diff.txt', 'w') as filehandle:
    for listitem in diffusion:
        filehandle.write('%s\n' % listitem)

with open('data/rho.txt', 'w') as filehandle:
    for listitem in densityinv:
        filehandle.write('%s\n' % listitem)