## skeleton for a class to read simulation data

### TODO:
###  - give the exception a good text
###  - the directory_list is a list, now it is treated as
###    a string.


import pandas as pd
import os

class MesoRDsimulation:

	def __init__(self, directory_list, cube_size):
		''' Reads a list of directories of MesoRD simulations.
	
		'''
		if (not os.path.isdir(directory_list)):
			raise NotADirectoryError('the directory is not an directory')
		self.directories = directory_list
		self.cube_size = cube_size
		pass

	def get_species(self):
		''' return a list of the different species in the simulation(s).
		'''
		
		
# for thesting purposes
if __name__ == "__main__":
	new_sim = MesoRDsimulation('simulation_example','0.01')
	print(new_sim.get_species())
	# should give an exception
	#new_sim = MesoRDsimulation('simulation_ample','0.01')