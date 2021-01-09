import numpy as np


def getReferenceData():
    data_array = [1, 2, 12, 22, 32, 45, 56, 87, 120, 181]
    return data_array


def getReferencePoints():
    return list(np.arange(len(getReferenceData())).astype(float))
