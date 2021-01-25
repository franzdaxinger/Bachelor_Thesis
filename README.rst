Bayesian Inference for Spatial Epidemic Models
==========================================

This project is supposed to run on the ETH euler cluster. To load the required 
modules, execute first:

env2lmod \n
module load gcc/6.3.0 openmpi/3.0.1 fenics/2019.1.0_slepc \n
export OMP_NUM_THREADS=1 \n
export MPICXX=mpicxx \n


The folder /data contains the raw data used to build the epidemic model.
To generate the required data from it, the scripts must be executed in this order:

mesh_generator.py generates the mesh which is used for the computations. 

diffusion_reader.py reads in the movement of the people in switzerland from files in /data. 

diffusion_maker.py generates the non constant diffusion coefficient and inverse population density. 

cantons_reader.py reads in the shapes of the cantons of Switzerland. 

cantons_maker.py generates the used canton functions. 

Now we created all the needed data and functions to run the epidemic model. model_adaptive.py uses an 
adaptive timestep while model_const.py uses constant timesteps. The solution is saved in /videomaker/functions 
and can be visualized with plotter.py which saves images in /videomaker. To create videos showing the dynamics
of the disease, run /videomaker/videomaker.py. 

To do bayesian inference, run the script /tmcmc/run-tmcmc.py which saves the results in /tmcmc/_korali_result_tmcmc_
This result can be used to do the propagation in /propagation, by running run-execution.py and then propagate.py