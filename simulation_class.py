## skeleton for a class to read simulation data

### TODO:
###  - give the exception a good text
###  - add caching, try @lru_cache:
###      from functools import lru_cache

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


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
		self.color_list = {}
		for s in range(len(self.species_order)):
			self.color_list.update({self.species_order[s]:list(mcolors.TABLEAU_COLORS.values())[s]}) 


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
		#uu = pd.read_csv('simulation_example/trajectories.txt',sep = ' ', header=None)
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

	def plot_Occ(self, threshold = 0):
		'''
		'''
		fig, ax = plt.subplots()
		for label, data in self.species.iloc[self.time.values >= threshold, :].apply(pd.Series.value_counts, axis=1).items():
			ax.plot(self.time.iloc[self.time.values >= threshold],data,color=self.color_list[label], label=label)
		ax.set_title('Occupancy')
		ax.set_xlabel(r'time (s)')
		ax.set_ylabel(r'number of particles')
		ax.legend(loc=0)
		plt.show()
	
	def plot_trajectory(self, id, axis = "x", threshold = 0):
		'''
			axis is the axis to project down on,  practicaly just means it is the coordinate that will not be ploted.
		'''
		fig, ax = plt.subplots()
		if axis == "x":
			coord1 = self.y_coord.iloc[:,id]
			coord2 = self.z_coord.iloc[:,id]
			ax.set_xlabel(r'y-axis')
			ax.set_ylabel(r'z-axis')
		elif axis == "y":
			coord1 = self.x_coord.iloc[:,id]
			coord2 = self.z_coord.iloc[:,id]
			ax.set_xlabel(r'x-axis')
			ax.set_ylabel(r'z-axis')
		elif axis == "z":
			coord1 = self.x_coord.iloc[:,id]
			coord2 = self.y_coord.iloc[:,id]
			ax.set_xlabel(r'x-axis')
			ax.set_ylabel(r'y-axis')
		else:
			raise Exception("Not a correct axis")#Need a more specific exception
		local_species_list = self.species.iloc[:,id]
		for i in range(len(coord1)-1):
			ax.plot(coord1.iloc[i:(i+2)], coord2.iloc[i:(i+2)], color=self.color_list[local_species_list.iloc[i]])
		ax.set_aspect(1)
		#trying to make a legend, does not work like I want it to, and not that important later. So skip for now.
		#ax.legend(loc=0, labels=list(self.color_list.keys()), labelcolor=list(self.color_list.values()))
		plt.show()
		
	def plot_trajectory_radial(self, id, axis = "z", threshold = 0):
		'''
			axis is the axis that will be the non radial part
		'''
		fig, ax = plt.subplots()
		if axis == "x":
			coord1 = self.x_coord.iloc[:,id]
			coord2 = np.sqrt(self.y_coord.iloc[:,id].pow(2) + self.z_coord.iloc[:,id].pow(2))
			ax.set_xlabel(r'x-axis')
			ax.set_ylabel(r'radial')
		elif axis == "y":
			coord1 = self.y_coord.iloc[:,id]
			coord2 = np.sqrt(self.z_coord.iloc[:,id].pow(2) + self.x_coord.iloc[:,id].pow(2))
			ax.set_xlabel(r'y-axis')
			ax.set_ylabel(r'radial')
		elif axis == "z":
			coord1 = self.z_coord.iloc[:,id]
			coord2 = np.sqrt(self.y_coord.iloc[:,id].pow(2) + self.x_coord.iloc[:,id].pow(2))
			ax.set_xlabel(r'z-axis')
			ax.set_ylabel(r'radial')
		else:
			raise Exception("Not a correct axis")#Need a more specific exception
		local_species_list = self.species.iloc[:,id]
		for i in range(len(coord1)-1):
			ax.plot(coord1.iloc[i:(i+2)], coord2.iloc[i:(i+2)], color=self.color_list[local_species_list.iloc[i]])
		ax.set_aspect(1)
		#trying to make a legend, does not work like I want it to, and not that important later. So skip for now.
		#ax.legend(loc=0, labels=list(self.color_list.keys()), labelcolor=list(self.color_list.values()))
		plt.show()
		
	def plot_hist(self, id, axis = "z", threshold = 0, bins = 30):
		'''
			axis is the axis that will be the non radial part
		'''
		if axis == "x":
			coord1 = self.x_coord.iloc[:,id]
		elif axis == "y":
			coord1 = self.y_coord.iloc[:,id]
		elif axis == "z":
			coord1 = self.z_coord.iloc[:,id]
		else:
			raise Exception("Not a correct axis")#Need a more specific exception
		fig, ax = pd.DataFrame({'species': self.species[self.time.values >= threshold].iloc[:,id], 'coord':coord1[self.time.values >= threshold]}).hist(by = 'species', bins = bins)
		plt.show()

	def set_colors(self, color_set):
		'''
		'''
		if (self.color_list.keys() == color_set.keys()):
			self.color_list = color_set
		else:
			raise Exception("the keys does not match")#Need a more specific exception

