import numpy as np


def getReferenceData():
    y = []
    bev = 8570000.0
    data_array = [9.00517345645, 24.8012596045, 1.407627618, 15.5958696162, 17.5834634767, 25.4112022583, 27.3530423252, 50.4803130653, 51.3659715025, 72.5844526582]
    for i in data_array:
        y.append(i / bev)
    return y


def getReferencePoints():
    return list(np.arange(len(getReferenceData())).astype(float))
