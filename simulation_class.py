### An class for reading MesoRD simulations and calculating:
###  - occupancies
###  - diffusion coefficients
###  - dwell time
### of the different species.

### Might add caching(try @lru_cache ; from functools import lru_cache)
### later efter checking some things.

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


class MesoRDsimulation:
	def __init__(self, sim_directory, cube_size):
		''' Reads in a MesoRD simulation from a directoryd.

			Parameters
			----------
			sim_directory: string
				Path to a directory containing a MesoRD simulation.

			cube_size: number
				The side size of the cubes used in the simulation.
		'''
		if (not os.path.isdir(sim_directory)):
			raise NotADirectoryError(sim_directory + ' is not an directory.')
		trajectories_file = os.path.join(sim_directory, 'trajectories.txt')
		if (not os.path.isfile(trajectories_file)):
			raise FileNotFoundError(trajectories_file + " not found.")
		reactions_file = os.path.join(sim_directory, 'reactions.txt')
		if (not os.path.isfile(reactions_file)):
			raise FileNotFoundError(reactions_file + " not found.")

		self.cube_size = cube_size
		trajectories = pd.read_csv(trajectories_file,sep = ' ', header=None)
		# Sort reaction dataframe by ID and time
		self.reactions = pd.read_csv(reactions_file,sep = ' ', header=None).sort_values(by=[1,0])
		# Order of trajectory file time, then for each particle the group of
		# ID, species, x coordinate, y coordinate and z coordinate.
		self.time = trajectories.iloc[:,0]
		# self.IDs = trajectories.iloc[:,1::5] # Not needed but keep the line for later reference
		self.species = trajectories.iloc[:,2::5]
		self.species.columns =  [str(i) for i in range(len(self.species.columns))] # Rename columns
		
		self.x_coord = trajectories.iloc[:,3::5]
		self.x_coord.columns =  [str(i) for i in range(len(self.x_coord.columns))] # Rename columns
		self.y_coord = trajectories.iloc[:,4::5]
		self.y_coord.columns =  [str(i) for i in range(len(self.y_coord.columns))] # Rename columns
		self.z_coord = trajectories.iloc[:,5::5]
		self.z_coord.columns =  [str(i) for i in range(len(self.z_coord.columns))] # Rename columns

		# Define a default order of the different species.
		self.species_order = self.species.melt().value.unique()
		
		# Define default colours for each species.
		self.color_list = {}
		for s in range(len(self.species_order)):
			self.color_list.update({self.species_order[s]:list(mcolors.TABLEAU_COLORS.values())[s]}) 


	def get_species(self):
		''' Return a list of the different species in the simulation.
			Will be in the order that species are presented.
			
			Return
			------
				list
		'''
		return self.species_order

	def set_species_order(self, value):
		''' Sets the order of the different species.

			Parameters
			----------
			value: list
				A list containing the different species.
		'''
		# Check that the species are all included.
		if (np.array_equal(np.sort(self.species_order), np.sort(value))):
			self.species_order = value
		else:
			raise ValueError("The elements in the list does not match the species.")

	def get_colors(self):
		''' Return a dictionary of the species and their colors.
			The keys are the species and the values their color.
			
			Return
			------
				dict
		'''
		return self.color_list

	def set_colors(self, color_set):
		''' Sets the colours of the different species.

			Parameters
			----------
			color_set: dict
				A dictionary having the species as keys and colors as values.
		'''
		# Check that all species are included.
		if (self.color_list.keys() == color_set.keys()):
			self.color_list = color_set
		else:
			raise ValueError("The keys does not match the species.")

	def get_pOcc_mean(self, threshold = 0):
		''' Return the mean occupany of the different species after the time
			defined by the threshold.

			Parameters
			----------
				threshold: number
					The times after the threshold are used to calculate the occupany.

			Return
			------
				panda series
		'''
		return self.species.iloc[self.time.values >= threshold, :].apply(pd.Series.value_counts, axis=1, normalize=True).apply(pd.Series.mean).loc[self.species_order]

	def get_pOcc_sd(self, threshold = 0):
		''' Return the sd of the occupany of the different species after the time
			defined by the threshold.

			Parameters
			----------
				threshold: number
					The times after the threshold are used to calculate the occupany.

			Return
			------
				panda series
		'''
		return self.species.iloc[self.time.values >= threshold, :].apply(pd.Series.value_counts, axis=1, normalize=True).apply(pd.Series.std).loc[self.species_order]

	def __get_D_values__(self, threshold):
		''' Helper function to calculate the diffusion coefficients.
			The fomula used is 
				D = <dx^2> / (6*dt)
			Where dx is the difference between timesteps as a 3D vector and dt is the timestep length.
			Here dx^2 / (6*dt) is calculated for each timestep.
			The coordinates are in units of cubes, so need to add the cube size to the formula.
			Only timesteps when the particle stays the same species are included.

			Parameters
			----------
				threshold: number
					The times after the threshold are used to calculate the diffusion.

			Return
			------
				DataFrame
		'''
		# The threshold in time.
		threshold_test = (self.time.values >= threshold)
		# Calculate the "diffusion coefficient" of each step.
		D_values = self.cube_size*self.cube_size* (
			self.x_coord.iloc[threshold_test, :].diff().iloc[1:,:].pow(2) + 
			self.y_coord.iloc[threshold_test, :].diff().iloc[1:,:].pow(2) +
			self.z_coord.iloc[threshold_test, :].diff().iloc[1:,:].pow(2)
			).divide(6*self.time.iloc[threshold_test].diff().iloc[1:], axis= 'rows')
		# Only include the ones where the the particle stay the same the whole time interval
		timestep_test = (self.species.iloc[threshold_test, :].shift(1).iloc[1:,:] == self.species.iloc[threshold_test, :].iloc[1:,:])
		return pd.DataFrame({'species': self.species.iloc[threshold_test, :].iloc[1:,:][timestep_test].melt().value, 'x':D_values[timestep_test].melt().value})

	def get_D_mean(self, threshold = 0):
		''' Return the mean diffusion of the different species after the time
			defined by the threshold.

			Parameters
			----------
				threshold: number
					The times after the threshold are used to calculate the diffusion.

			Return
			------
				panda series
		'''
		return self.__get_D_values__(threshold).groupby('species').mean().loc[self.species_order].iloc[:,0]

	def get_D_sd(self, threshold = 0):
		''' Return the sd of the diffusion of the different species after the time
			defined by the threshold.

			Parameters
			----------
				threshold: number
					The times after the threshold are used to calculate the diffusion.

			Return
			------
				panda series
		'''
		return self.__get_D_values__(threshold).groupby('species').std().loc[self.species_order].iloc[:,0]

	def __get_DT_values__(self, threshold):
		''' Helper function to calculate the dwell times. Calculated from the difference in time for
			the reactions.


			Parameters
			----------
				threshold: number
					The times after the threshold are used to calculate the dwell time.

			Return
			------
				DataFrame
		'''
		# Group by the ID, then take difference between time and shift to align with the species.
		return pd.DataFrame({'species': self.reactions[self.reactions[0] >= threshold].iloc[:,2], 'diff':self.reactions[self.reactions[0] >= threshold].groupby(1)[0].diff().shift(-1)})

	def get_DT_mean(self, threshold = 0):
		''' Return the mean dwell time of the different species after the time
			defined by the threshold.

			Parameters
			----------
				threshold: number
					The times after the threshold are used to calculate the dwell time.

			Return
			------
				panda series
		'''
		return self.__get_DT_values__(threshold).groupby('species').mean().loc[self.species_order].iloc[:,0]

	def get_DT_sd(self, threshold = 0):
		''' Return the sd of the dwell time of the different species after the time
			defined by the threshold.

			Parameters
			----------
				threshold: number
					The times after the threshold are used to calculate the dwell time.

			Return
			------
				panda series
		'''
		return self.__get_DT_values__(threshold).groupby('species').std().loc[self.species_order].iloc[:,0]

	def plot_pOcc(self):
		''' Plots the occupany of the different species over time.
		'''
		fig, ax = plt.subplots()
		# Plot the species individually for the sake of the legend.
		for label, data in self.species.apply(pd.Series.value_counts, axis=1).items():
			ax.plot(self.time,data,color=self.color_list[label], label=label)
		ax.set_title('Occupancy')
		ax.set_xlabel(r'time (s)')
		ax.set_ylabel(r'number of particles')
		ax.legend(loc=0)
		plt.show()
	
	def plot_trajectory(self, id, axis = "z", lower_time = float("-inf") , upper_time = float("inf")):
		''' Plots a trajectory projected down to a plane. The trajectory
			parts are coloured according to which species that part is.

			Parameters
			----------
				id: number
					The ID of the trajectory.
				axis: string ("x", "y" or "z")
					The axis that is normal to the plane of projection.
				lower_time: number
					The lower time threshold of trajectory plotted
				upper_time: number
					The upper time threshold of trajectory plotted
		'''
		threshold_test = (self.time.values >= lower_time) & (self.time.values < upper_time)
		fig, ax = plt.subplots()
		if axis == "x":
			coord1 = self.cube_size*self.y_coord.iloc[threshold_test,id]
			coord2 = self.cube_size*self.z_coord.iloc[threshold_test,id]
			ax.set_xlabel(r'y-axis')
			ax.set_ylabel(r'z-axis')
		elif axis == "y":
			coord1 = self.cube_size*self.x_coord.iloc[threshold_test,id]
			coord2 = self.cube_size*self.z_coord.iloc[threshold_test,id]
			ax.set_xlabel(r'x-axis')
			ax.set_ylabel(r'z-axis')
		elif axis == "z":
			coord1 = self.cube_size*self.x_coord.iloc[threshold_test,id]
			coord2 = self.cube_size*self.y_coord.iloc[threshold_test,id]
			ax.set_xlabel(r'x-axis')
			ax.set_ylabel(r'y-axis')
		else:
			raise ValueError(axis + " is not a correct axis")
		local_species_list = self.species.iloc[threshold_test,id]
		for i in range(len(coord1)-1):
			ax.plot(coord1.iloc[i:(i+2)], coord2.iloc[i:(i+2)], color=self.color_list[local_species_list.iloc[i]])
		ax.set_aspect(1)
		plt.show()
		
	def plot_trajectory_radial(self, id, axis = "x", lower_time = float("-inf") , upper_time = float("inf")):
		''' Plots the distance of the particle from the axis and the axis it self.
			Like cylindrical coordinates without the angle. The trajectory
			parts are coloured according to which species that part is.

			Parameters
			----------
				id: number
					The ID of the trajectory.
				axis: string ("x", "y" or "z")
					The axis that distance is defined from.
				lower_time: number
					The lower time threshold of trajectory plotted
				upper_time: number
					The upper time threshold of trajectory plotted
		'''
		threshold_test = (self.time.values >= lower_time) & (self.time.values < upper_time)
		fig, ax = plt.subplots()
		if axis == "x":
			coord1 = self.cube_size*self.x_coord.iloc[threshold_test,id]
			coord2 = self.cube_size*np.sqrt(self.y_coord.iloc[threshold_test,id].pow(2) + self.z_coord.iloc[threshold_test,id].pow(2))
			ax.set_xlabel(r'x-axis')
			ax.set_ylabel(r'radial')
		elif axis == "y":
			coord1 = self.cube_size*self.y_coord.iloc[threshold_test,id]
			coord2 = self.cube_size*np.sqrt(self.z_coord.iloc[threshold_test,id].pow(2) + self.x_coord.iloc[threshold_test,id].pow(2))
			ax.set_xlabel(r'y-axis')
			ax.set_ylabel(r'radial')
		elif axis == "z":
			coord1 = self.cube_size*self.z_coord.iloc[threshold_test,id]
			coord2 = self.cube_size*np.sqrt(self.y_coord.iloc[threshold_test,id].pow(2) + self.x_coord.iloc[threshold_test,id].pow(2))
			ax.set_xlabel(r'z-axis')
			ax.set_ylabel(r'radial')
		else:
			raise ValueError(axis + " is not a correct axis")
		local_species_list = self.species.iloc[threshold_test,id]
		for i in range(len(coord1)-1):
			ax.plot(coord1.iloc[i:(i+2)], coord2.iloc[i:(i+2)], color=self.color_list[local_species_list.iloc[i]])
		ax.set_aspect(1)
		plt.show()
		
	def plot_hist(self, axis = "x", bins = 30):
		''' Plots a histogram of the values of the axis coordinates. One
			subplot per species.

			Parameters
			----------
				axis: string ("x", "y" or "z")
					The axis position values are taken from.
				bins: integer
					The number of bins to use in the histogram.

		'''
		if axis == "x":
			coord1 = (self.cube_size*self.x_coord).melt()
		elif axis == "y":
			coord1 = (self.cube_size*self.y_coord).melt()
		elif axis == "z":
			coord1 = (self.cube_size*self.z_coord).melt()
		else:
			raise ValueError(axis + " is not a correct axis")
		fig, ax = pd.DataFrame({'species': self.species.melt()['value'], 'coord':coord1['value']}).hist(by = 'species', bins = bins)
		plt.show()

