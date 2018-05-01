#!/usr/bin/env python
from __future__ import division
import random
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate


class Step_characteristic():
	def __init__(self, grid, xs_data, quadrature_order): 
		self.grid = grid
		self.xs = xs_data
		self.cells = range(len(grid))
		self.directions, self.weights = np.polynomial.legendre.leggauss(quadrature_order)
		
		self.fluxes = [np.ones(len(grid)), np.ones(len(grid))]
		self.currents = [np.ones(len(grid)), np.ones(len(grid))]

		
	def update_Source(self):
		self.source = [np.zeros(len(self.grid)), np.zeros(len(self.grid))]
		for i, material in enumerate(self.grid['material']):

			sigma_s11 = self.xs.loc[material, 'sigma_s11']
			sigma_s12 = self.xs.loc[material, 'sigma_s12']
			sigma_s22 = self.xs.loc[material, 'sigma_s22']
			
			if material != 'water': #only thermal fission
				nu2 = self.xs.loc[material, 'nu2']
				sigma_f2 = self.xs.loc[material, 'sigma_f2']
				
				self.source[1][i] += self.fluxes[1][i] * nu2 * sigma_f2
				
			self.source[0][i] += self.fluxes[0][i] * sigma_s11 #within group fast scatter
			self.source[1][i] += self.fluxes[0][i] * sigma_s12 + self.fluxes[1][i] * sigma_s22 #fast down scatter thermal within scatter
			
			
	def iterate_flux(self):

		fluxes_ = [np.zeros(len(grid)), np.zeros(len(grid))]
		currents_ = [np.zeros(len(grid)), np.zeros(len(grid))]
		
		current_in = 0.005
		for group in range(len(self.fluxes)):
		
			for idx, mu in enumerate(self.directions):
				if mu > 0:
					for i in self.cells:
						material = grid['material'].iloc[i]
						sigma_t1 = xs.loc[material, 'sigma_t1']
						sigma_t2 = xs.loc[material, 'sigma_t2']
						xs_data = [sigma_t1, sigma_t2]
						
						tau = spacing * xs_data[group] / abs(mu)
						Q_cell = self.source[group][i]

						current_out = current_in * np.exp(-tau) + Q_cell / xs_data[group] * (1 - np.exp(-tau))
						current_cell = Q_cell / xs_data[group] - mu * (current_out - current_in)/(spacing * xs_data[group])
						
						fluxes_[group][i] += current_cell * self.weights[idx]
						currents_[group][i] += current_out * self.weights[idx] * self.directions[idx]

						current_in = current_out
					
	
						
				if mu < 0:
					for i in reversed(self.cells):
						material = grid['material'].iloc[i]
						sigma_t1 = xs.loc[material, 'sigma_t1']
						sigma_t2 = xs.loc[material, 'sigma_t2']
						xs_data = [sigma_t1, sigma_t2]
						
						tau = spacing * xs_data[group] / abs(mu)
						Q_cell = self.source[group][i]
						
						current_out = current_in * np.exp(-tau) + Q_cell / xs_data[group] * (1 - np.exp(-tau))
						current_cell = Q_cell / xs_data[group] - mu * (current_out - current_in)/(spacing * xs_data[group])
						
						fluxes_[group][i] += current_cell * self.weights[idx]
						currents_[group][i] += current_out * self.weights[idx] * self.directions[idx]
						
						current_in = current_out

	
		
			self.fluxes[group] = fluxes_[group]			
			self.currents[group] = currents_[group]	
		

	


	
	
grid = pd.read_csv('grid.csv')
xs = pd.read_csv('XS.csv', index_col = 'material')

global spacing
spacing = 0.15625
order = 6





flux_solver = Step_characteristic(grid, xs, 6)
for i in range(10):
	flux_solver.update_Source()
	for i in range(10):
		flux_solver.iterate_flux()


# for i in range(100):
	# for idx in range(len(directions)):
		# if directions[idx] > 0:
			# for i in cells:
				# material = grid['material'].iloc[i]
				# sigma_t1 = xs.loc[material, 'sigma_t1']
				# tau = spacing * sigma_t1 / directions[idx]
				
				# current_out = calc_current_out(current_in, Q, tau, sigma_t1)
				# current_mean = calc_mean_current(current_in, current_out, Q, directions[idx], sigma_t1)
				
				# flux[i] += current_mean * weights[idx]
				# current_cell[i] = current_mean * weights[idx] * directions[idx]
				
				# current_in = current_out
		# else:	
			# for i in reversed(cells):
				# material = grid['material'].iloc[i]
				# sigma_t1 = xs.loc[material, 'sigma_t1']
				# tau = spacing * sigma_t1 / directions[idx]
				
				# current_out = calc_current_out(current_in, Q, tau, sigma_t1)
				# current_mean = calc_mean_current(current_in, current_out, Q, directions[idx], sigma_t1)
				
				# flux[i] += current_mean * weights[idx]
				# current_cell[i] += current_mean * weights[idx] * directions[idx]
				# print('left')
				# current_in = current_out
			
	
		


#plt.plot(flux)
#plt.show()
# f = lambda x: np.cos(x)
# a = 0.0
# b = np.pi/2	