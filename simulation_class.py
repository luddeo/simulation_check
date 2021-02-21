## skeleton for a class to read simulation data

### TODO:
###  - give the exception a good text
###  - Need a time threshold to funcitons giving pOcc, DT, D
###  - add caching, try lru_cache that is imported below or use somethign else.


import pandas as pd
import os
from functools import lru_cache

class MesoRDsimulation:

	def __init__(self, sim_directory, cube_size):
		''' Reads in a MesoRD simulations in sim_directory.
	
		'''
		if (not os.path.isdir(sim_directory)):
			raise NotADirectoryError('the directory is not an directory')
		self.directory = sim_directory
		self.cube_size = cube_size
		trajectories_file = os.path.join(sim_directory, 'trajectories.txt')
		reactions_file = os.path.join(sim_directory, 'reactions.txt')
		trajectories = pd.read_csv(trajectories_file,sep = ' ', header=None)
		self.time = trajectories.iloc[:,0]
		self.IDs = trajectories.iloc[:,1::5]
		self.species = trajectories.iloc[:,2::5]
		self.species.columns =  [str(i) for i in range(len(self.species.columns))]
		self.pOcc_m = self.species.apply(pd.Series.value_counts, axis=1, normalize=True).apply(pd.Series.mean)
		self.pOcc_sd = self.species.apply(pd.Series.value_counts, axis=1, normalize=True).apply(pd.Series.std)
		
		x_coord = trajectories.iloc[:,3::5]
		x_coord.columns =  [str(i) for i in range(len(x_coord.columns))]
		y_coord = trajectories.iloc[:,4::5]
		y_coord.columns =  [str(i) for i in range(len(y_coord.columns))]
		z_coord = trajectories.iloc[:,5::5]
		z_coord.columns =  [str(i) for i in range(len(z_coord.columns))]
		
		
		D_df = cube_size*cube_size*(x_coord.diff().iloc[1:,:].pow(2) + y_coord.diff().iloc[1:,:].pow(2) + z_coord.diff().iloc[1:,:].pow(2)).divide(6*self.time.diff().iloc[1:], axis= 'rows')
		D_longdf = pd.DataFrame({'species': self.species.iloc[1:,:].melt().value, 'x': D_df.melt().value})
		# Now just using the first species in each interval, would rater us ehen the interval is all the same specis
		#species.iloc[1:,:] # remove first row 
		#species.iloc[0:(-1),:] # remove last row
		# Need to compare the two and only use then both are the same.
		self.D_m = D_longdf.groupby('species').mean()
		self.D_sd = D_longdf.groupby('species').std()
		
		# some random things done when trying to calculate things
		# formula for diffusion: cs^2 * dx2 / (6*dt)
		#x_coord.divide(time,axis='rows')
		#pd.DataFrame({'species': species.melt().value, 'x': x_coord.melt().value})
		#pd.DataFrame({'species': species.melt().value, 'x': x_coord.melt().value}).groupby('species').mean()
		

		# species   get unique: .melt().iloc[:,1].unique()
		# X  b1 = uu.iloc[:,3::5].diff().iloc[1:,:].pow(2)  ; b1.columns =  [str(i) for i in range(500)]
		# .iloc[1:,:] removes the first row with NaN in it
		
		
		#uu2 = pd.read_csv('simulation_example/reactions.txt',sep = ' ', header=None)
		#uu2.iloc[:,2].unique() # species

	@lru_cache
	def get_species(self):
		''' return a list of the different species in the simulation. unsorted
		'''
		return [i for i in self.pOcc_m.index]

	@lru_cache
	def get_pOcc_mean(self):
		'''
		'''
		return self.pOcc_m
	
	def get_pOcc_sd(self):
		'''
		'''
		return self.pOcc_sd
		
	def get_D_mean(self):
		'''
		'''
		return self.D_m
	
	def get_D_sd(self):
		'''
		'''
		return self.D_sd
		
# for thesting purposes
if __name__ == "__main__":
	new_sim = MesoRDsimulation('simulation_example',0.01)
	print(new_sim.get_species())
	print("pOcc")
	print(new_sim.get_pOcc_mean())
	print(new_sim.get_pOcc_sd())
	print("D")
	print(new_sim.get_D_mean())
	print(new_sim.get_D_sd())
	# should give an exception
	#new_sim = MesoRDsimulation('simulation_ample','0.01')