import sys
import korali

sys.path.append('./_model')
import model_data
from model_call import model_call

k = korali.Engine()
e = korali.Experiment()

# Setting up the reference likelihood for the Bayesian Problem
e["Problem"]["Type"] = "Bayesian/Reference"
e["Problem"]["Likelihood Model"] = "Negative Binomial"
e["Problem"]["Reference Data"] = model_data.getReferenceData()
e["Problem"]["Computational Model"] = model_call

# Configuring TMCMC parameters
e["Solver"]["Type"] = "Sampler/TMCMC"
e["Solver"]["Population Size"] = 36
#e["Solver"]["Covariance Scaling"] = 0.04
e['Solver']['Target Coefficient Of Variation'] = 0.9

# Configuring the problem's random distributions
e["Distributions"][0]["Name"] = "Uniform 0"
e["Distributions"][0]["Type"] = "Univariate/Uniform"
e["Distributions"][0]["Minimum"] = 1.5
e["Distributions"][0]["Maximum"] = +3.0

e["Distributions"][1]["Name"] = "Uniform 1"
e["Distributions"][1]["Type"] = "Univariate/Uniform"
e["Distributions"][1]["Minimum"] = 0.
e["Distributions"][1]["Maximum"] = +15.0

# Configuring the problem's variables and their prior distributions
e["Variables"][0]["Name"] = "beta_param"
e["Variables"][0]["Prior Distribution"] = "Uniform 0"

e["Variables"][1]["Name"] = "[Sigma]"
e["Variables"][1]["Prior Distribution"] = "Uniform 1"

e["Store Sample Information"] = True

# Configuring output settings
e["File Output"]["Path"] = '_korali_result_tmcmc'

# Starting Korali's Engine and running experiment
e["Console Output"]["Verbosity"] = "Detailed"

k["Conduit"]["Type"] = "Concurrent"
k['Conduit']['Concurrent Jobs'] = 36

k.run(e)
