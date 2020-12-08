#!/usr/bin/env python3
import os
import sys
import json
sys.path.append('./_model/')
from work import *

import numpy as np

# Opening JSON file
f = open('gen00000005.json', )
data = json.load(f)

beta_list = []
sigma_list = []
for i in data['Results']['Sample Database']:
    beta_list.append(i[0])
    sigma_list.append(i[1])

# Closing file
f.close()

# Creating new experiment
import korali
e = korali.Experiment()
k = korali.Engine()

e["Problem"]["Type"] = "Propagation"
e["Problem"]["Execution Model"] = lambda modelData: model( modelData)

e["Variables"][0]["Name"] = "Beta"
e["Variables"][0]["Precomputed Values"] = beta_list
e["Variables"][1]["Name"] = "Sigma"
e["Variables"][1]["Precomputed Values"] = sigma_list

e["Solver"]["Type"] = "Executor"
e["Solver"]["Executions Per Generation"] = 1

e["Console Output"]["Verbosity"] = "Minimal"
e["Store Sample Information"] = True
e["File Output"]["Frequency"] = 0

# Starting Korali's Engine and running experiment
k.run(e)
