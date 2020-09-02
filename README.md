The model simulates an epidemic in Switzerland and works as follows:

p2_mesh_generator.py generates a cartesian mesh

p2_diffusion_writer.py and p2_diffusion_maker.py read in files from /shapefiles and create a non constant diffusion coefficient for switzerland and functions for population density and the transmission rate beta. These functions are saved in /difffun

p2_cantons_writer.py and p2_cantons_maker.py create functions to represent every single canton and the data to take commuters over long distances into account

p2_model_adaptive.py simulates the pandemic with commuters over long distances using an adaptive time stepping

p2_model_const.py simulates the pandemic with commuters over long distances using a constant timestep

plotter.py generates plots of the solutions

Videomaker/videomaker.py creates a video from the created images
