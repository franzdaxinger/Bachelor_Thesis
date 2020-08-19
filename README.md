The model simulates an epidemic in Switzerland and works as follows:

p2_mesh_generator.py generates a cartesian mesh

p2_diffusion_writer.py and p2_diffusion_maker.py read in files from /shapefiles and create a non constant diffusion coefficient for switzerland and functions for population density and the transmission rate beta. These functions are saved in /difffun

p2_cantons_writer.py and p2_cantons_maker.py create functions to represent every single canton and the data to take commuters over long distances into account

p2_adaptive_no_lc.py is a simulation of the pandemic without the long connections between cantons

p2_model_adaptive.py simulates the pandemic with commuters over long distances

Videomaker/videomaker.py creates a video from the created images
