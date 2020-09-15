import numpy as np


def getReferenceData():
    y = []
    bev = 8570000.0
    #data_array = [12, 22, 32, 45, 56, 87, 120, 181, 243, 316, 365, 434, 626, 838, 1172, 1529, 1962, 2382, 2710, 3774]
    data_array = [12, 22, 32, 45]  # XXX reduced dataset
    for i in data_array:
        y.append(i / bev)
    return y


# days until 16th of march, when the lockdown happened
def getReferencePoints():
    return list(np.arange(len(getReferenceData())).astype(float))
