#!/usr/bin/env python3
import os
import sys
sys.path.append('./_model/')
from work import *

import numpy as np

# fileName: output directory
fileName = 'samplesOut.dat'

# Prepare file
prepareOutputDir()

# Creating new experiment
import korali
e = korali.Experiment()
k = korali.Engine()

e["Problem"]["Type"] = "Propagation"
e["Problem"]["Execution Model"] = lambda modelData: model( modelData, fileName)

e["Variables"][0]["Name"] = "Beta"
e["Variables"][0]["Precomputed Values"] = list(np.arange(1.5, 3.0, 0.01))
e["Variables"][1]["Name"] = "Sigma"
e["Variables"][1]["Precomputed Values"] = list(np.arange(0, 12.0, 0.08))

e["Solver"]["Type"] = "Executor"
e["Solver"]["Executions Per Generation"] = 1

#k["Conduit"]["Type"] = "Concurrent"
#k['Conduit']['Concurrent Jobs'] = 13

e["Console Output"]["Verbosity"] = "Minimal"
e["Store Sample Information"] = True                    # new

# Starting Korali's Engine and running experiment
k.run(e)
