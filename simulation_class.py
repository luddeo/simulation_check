## skeleton for a class to read simulation data

### TODO:
###  - give the exception a good text
###  - add caching, try @lru_cache:
###      from functools import lru_cache

import pandas as pd
import numpy as np
import os


class MesoRDsimulation:
	def __init__(self, sim_directory, cube_size):
		''' Reads in a MesoRD simulation in sim_directory.
	
		'''
		if (not os.path.isdir(sim_directory)):
			raise NotADirectoryError('the directory is not an directory')
		self.directory = sim_directory
		self.cube_size = cube_size
		trajectories_file = os.path.join(sim_directory, 'trajectories.txt')
		reactions_file = os.path.join(sim_directory, 'reactions.txt')
		trajectories = pd.read_csv(trajectories_file,sep = ' ', header=None)
		self.reactions = pd.read_csv(reactions_file,sep = ' ', header=None).sort_values(by=[1,0])
		self.time = trajectories.iloc[:,0]
		self.IDs = trajectories.iloc[:,1::5]
		self.species = trajectories.iloc[:,2::5]
		self.species.columns =  [str(i) for i in range(len(self.species.columns))]
		
		self.x_coord = trajectories.iloc[:,3::5]
		self.x_coord.columns =  [str(i) for i in range(len(self.x_coord.columns))]
		self.y_coord = trajectories.iloc[:,4::5]
		self.y_coord.columns =  [str(i) for i in range(len(self.y_coord.columns))]
		self.z_coord = trajectories.iloc[:,5::5]
		self.z_coord.columns =  [str(i) for i in range(len(self.z_coord.columns))]
		
		self.species_order = self.species.melt().value.unique()


	def get_species(self):
		''' return a list of the different species in the simulation. Need not be sorted
		'''
		return self.species_order

	def set_species_order(self, value):
		if (np.array_equal(np.sort(self.species_order), np.sort(value))):
			self.species_order = value
		else:
			raise Exception("the elements does not match")#Need a more specific exception

	def get_pOcc_mean(self, threshold = 0):
		'''
		'''
		return self.species.iloc[self.time.values >= threshold, :].apply(pd.Series.value_counts, axis=1, normalize=True).apply(pd.Series.mean).loc[self.species_order]

	def get_pOcc_sd(self, threshold = 0):
		'''
		'''
		return self.species.iloc[self.time.values >= threshold, :].apply(pd.Series.value_counts, axis=1, normalize=True).apply(pd.Series.std).loc[self.species_order]

	def __get_D_values__(self, threshold):
		'''
		'''
		# Formula for 3D-diffusion: D = cs^2 * dx2 / (6*dt)
		# Now just using the first species in each interval, would rather use the intervals where it stays the same species
		#species.iloc[1:,:] # remove first row 
		#species.iloc[0:(-1),:] # remove last row
		# Need to compare the two and only use then both are the same.
		D_values = self.cube_size*self.cube_size* (
			self.x_coord.iloc[self.time.values >= threshold, :].diff().iloc[1:,:].pow(2) + 
			self.y_coord.iloc[self.time.values >= threshold, :].diff().iloc[1:,:].pow(2) +
			self.z_coord.iloc[self.time.values >= threshold, :].diff().iloc[1:,:].pow(2)
			).divide(6*self.time.iloc[self.time.values >= threshold].diff().iloc[1:], axis= 'rows')
		return pd.DataFrame({'species': self.species.iloc[self.time.values >= threshold, :].iloc[1:,:].melt().value, 'x':D_values.melt().value})

	def get_D_mean(self, threshold = 0):
		'''
		'''
		return self.__get_D_values__(threshold).groupby('species').mean().loc[self.species_order]

	def get_D_sd(self, threshold = 0):
		'''
		'''
		return self.__get_D_values__(threshold).groupby('species').std().loc[self.species_order]

	def __get_DT_values__(self, threshold):
		'''
		'''
		return pd.DataFrame({'species': self.reactions[self.reactions[0] >= threshold].iloc[:,2], 'diff':self.reactions[self.reactions[0] >= threshold].groupby(1)[0].diff().shift(-1)})
		#self.reactions[self.reactions[0] >= threshold].groupby(1)[0].diff().shift(-1)
		# some random things done when trying to calculate things
		#uu2 = pd.read_csv('simulation_example/reactions.txt',sep = ' ', header=None)
		#uu3 = uu2.sort_values(by=[1,0])
		#uu3[6] = uu3.groupby(1)[0].diff().shift(-1)
		#uu3.groupby(2).mean()
		#uu2.iloc[:,2].unique() # species

	def get_DT_mean(self, threshold = 0):
		'''
		'''
		return self.__get_DT_values__(threshold).groupby('species').mean().loc[self.species_order]

	def get_DT_sd(self, threshold = 0):
		'''
		'''
		return self.__get_DT_values__(threshold).groupby('species').std().loc[self.species_order]

# for testing purposes
if __name__ == "__main__":
	new_sim = MesoRDsimulation('simulation_example',0.01)
	new_sim.set_species_order(["Aa","Bb","Cc"])
	print(new_sim.get_species())
	print("pOcc")
	print(new_sim.get_pOcc_mean())
	print(new_sim.get_pOcc_mean(20))
	#print(new_sim.get_pOcc_sd())
	#print(new_sim.get_pOcc_sd(20))
	print("D")
	print(new_sim.get_D_mean())
	print(new_sim.get_D_sd())

	#print(new_sim.get_D_mean(20))
	#print(new_sim.get_D_sd(20))
	print("DT")
	print(new_sim.get_DT_mean())
	print(new_sim.get_DT_sd())

	#print(new_sim.get_pOcc_mean())
	#new_sim.set_species_order(["Aa","Bb","Cc"])
	#print(new_sim.get_species())
	#print(new_sim.get_pOcc_mean())
	# should give an exception
	#new_sim = MesoRDsimulation('simulation_ample','0.01')