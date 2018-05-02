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
		self.fission_source = np.zeros(len(self.grid))
		
		for i, material in enumerate(self.grid['material']):

			sigma_s11 = self.xs.loc[material, 'sigma_s11']
			sigma_s12 = self.xs.loc[material, 'sigma_s12']
			sigma_s22 = self.xs.loc[material, 'sigma_s22']
			
			if material != 'water': #only thermal fission
				nu2 = self.xs.loc[material, 'nu2']
				sigma_f2 = self.xs.loc[material, 'sigma_f2']
				
				self.source[0][i] += self.fluxes[1][i] * nu2 * sigma_f2
				self.fission_source[i] += self.fluxes[1][i] * nu2 * sigma_f2
				
			self.source[0][i] += self.fluxes[0][i] * sigma_s11 #within group fast scatter
			self.source[1][i] += self.fluxes[0][i] * sigma_s12 + self.fluxes[1][i] * sigma_s22 #fast down scatter thermal within scatter
			
			
	def iterate_flux(self):
		percent_diff = 1.0
		while percent_diff > 0.01:
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
							current_cell = Q_cell / xs_data[group] - abs(mu) * (current_out - current_in)/(spacing * xs_data[group])
							
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
							current_cell = Q_cell / xs_data[group] - abs(mu) * (current_out - current_in)/(spacing * xs_data[group])
							
							fluxes_[group][i] += current_cell * self.weights[idx]
							currents_[group][i] += current_out * self.weights[idx] * self.directions[idx]
							
							current_in = current_out

		
				percent_diff = sum(((fluxes_[group] - self.fluxes[group])/(fluxes_[group]))) / len(fluxes_[group])

				self.fluxes[group] = fluxes_[group]			
				self.currents[group] = currents_[group]	
		
				print(percent_diff)

	


	
	
grid = pd.read_csv('grid.csv')
xs = pd.read_csv('XS.csv', index_col = 'material')

global spacing
spacing = 0.15625
order = 6




flux_solver = Step_characteristic(grid, xs, order)
for i in range(2):
	flux_solver.update_Source()
	flux_solver.iterate_flux()



			
	
		


plt.plot(flux_solver.fluxes[0])
plt.show()
# f = lambda x: np.cos(x)
# a = 0.0
# b = np.pi/2	