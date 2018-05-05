#!/usr/bin/env python
from __future__ import division
import random
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate


class SnTransportSolver():
	def __init__(self, grid, xs_data, quadrature_order , stopping_criteria): 
		self.grid = grid
		self.xs = xs_data
		self.cells = range(len(grid))
		self.directions, self.weights = np.polynomial.legendre.leggauss(quadrature_order)
		self.stopping_criteria = stopping_criteria
		
		self.fluxes = [np.ones(len(grid)), np.ones(len(grid))] #initial flux is uniform
		self.edge_fluxes = [np.ones(len(grid)), np.ones(len(grid))] 
		self.currents = [np.ones(len(grid)), np.ones(len(grid))]
		self.k = 1.5
		self.update_Source()
		self.converged = False

		
	def update_Source(self):
		self.scatter_source = [np.zeros(len(self.grid)), np.zeros(len(self.grid))]
		self.fission_source = [np.zeros(len(self.grid)), np.zeros(len(self.grid))]
		
		for i, material in enumerate(self.grid['material']):

			sigma_s11 = self.xs.loc[material, 'sigma_s11']
			sigma_s12 = self.xs.loc[material, 'sigma_s12']
			sigma_s22 = self.xs.loc[material, 'sigma_s22']
			
			if material != 'water': #only thermal fission
				nu2 = self.xs.loc[material, 'nu2']
				sigma_f2 = self.xs.loc[material, 'sigma_f2']
				self.fission_source[0][i] += self.fluxes[1][i] * nu2 * sigma_f2 
				
			self.scatter_source[0][i] += (self.fluxes[0][i] * sigma_s11) #within group fast scatter
			self.scatter_source[1][i] += (self.fluxes[0][i] * sigma_s12 + self.fluxes[1][i] * sigma_s22)  #fast down scatter thermal within scatter
			
		self.scatter_source = np.asarray(self.scatter_source)
		self.fission_source = np.asarray(self.fission_source) 
			
		#self.scatter_source[0] = self.scatter_source[0]/max(self.scatter_source[0])
		#self.scatter_source[1] = self.scatter_source[1]/max(self.scatter_source[1])
		
		#self.fission_source[0] = self.fission_source[0]/max(self.fission_source[0])
		
		return self.scatter_source, self.fission_source
			
			
	def iterate_flux(self):
		percent_diff = 100
		
		boundary_currents = [np.zeros(len(self.directions)), np.zeros(len(self.directions))]
		
		while percent_diff > self.stopping_criteria: #stopping criteria
			fluxes_ = [np.zeros(len(grid)), np.zeros(len(grid))]
			edge_fluxes_ = [np.zeros(len(grid)), np.zeros(len(grid))]
			currents_ = [np.zeros(len(grid)), np.zeros(len(grid))]
			

			for group in range(len(self.fluxes)):
				for idx, mu in enumerate(self.directions):
					if mu < 0:
						for i in reversed(self.cells):
							material = grid['material'].iloc[i]
							sigma_t1 = xs.loc[material, 'sigma_t1']
							sigma_t2 = xs.loc[material, 'sigma_t2']
							xs_data = [sigma_t1, sigma_t2]
							
							tau = spacing * xs_data[group] / abs(mu)
							Q_cell = (self.scatter_source[group][i] + self.fission_source[group][i]/self.k) / 2
							
							current_in = boundary_currents[group][idx]
							
							current_out = current_in * np.exp(-tau) + Q_cell / xs_data[group] * (1 - np.exp(-tau))
							current_cell = Q_cell / xs_data[group] - abs(mu) * (current_out - current_in)/(spacing * xs_data[group])
							
							fluxes_[group][i] += current_cell * self.weights[idx]
							edge_fluxes_[group][i] += current_out * self.weights[idx]
							currents_[group][i] += current_out * self.weights[idx] * self.directions[idx]
							
							boundary_currents[group][idx] = current_out
							
					if mu > 0:
						for i in self.cells:
							material = grid['material'].iloc[i]
							sigma_t1 = xs.loc[material, 'sigma_t1']
							sigma_t2 = xs.loc[material, 'sigma_t2']
							xs_data = [sigma_t1, sigma_t2]
							
							tau = spacing * xs_data[group] / abs(mu)
							Q_cell =  (self.scatter_source[group][i] + self.fission_source[group][i]/self.k) / 2
							
							current_in = boundary_currents[group][idx]

							current_out = current_in * np.exp(-tau) + Q_cell / xs_data[group] * (1 - np.exp(-tau))
							current_cell = Q_cell / xs_data[group] - abs(mu) * (current_out - current_in)/(spacing * xs_data[group])
							
							fluxes_[group][i] += current_cell * self.weights[idx]
							edge_fluxes_[group][i] += current_out * self.weights[idx]
							currents_[group][i] += current_out * self.weights[idx] * self.directions[idx]
							
							boundary_currents[group][idx] = current_out
						
				#stopping criteria

				percent_diff = np.amax(abs(fluxes_[group] - self.fluxes[group]))
				self.fluxes[group] = fluxes_[group]		
				self.edge_fluxes[group] = edge_fluxes_[group]	
				self.currents[group] = currents_[group]	
				
	def power_iteration(self):
		fission_old = self.fission_source[0]
		source_new, fission_new = self.update_Source()
		#print(sum(source_new[0]), sum(source_new[1]))
		
		k_new = self.k * sum(fission_new[0]) / sum(fission_old)
		
		
		k_diff = abs(k_new - self.k)/k_new
		F_diff = abs(np.nanmax((fission_new[0] - fission_old) / fission_new[0]))
		print(k_diff, F_diff)
		self.k = k_new
		if k_diff < self.stopping_criteria and F_diff < self.stopping_criteria:
			self.converged = True
			print('converged')



		
		


	
grid = pd.read_csv('grid.csv')
xs = pd.read_csv('XS.csv', index_col = 'material')

global spacing
spacing = 0.15625
order = 20
percent_diff = .0001



solver = SnTransportSolver(grid, xs, order, percent_diff)
power_iterations = 0
while True:
	solver.iterate_flux()
	solver.power_iteration()
	power_iterations += 1
	if solver.converged == True:
		print('k:', solver.k, power_iterations, 'iterations')
		break


plt.plot(range(len(solver.fluxes[0])), solver.edge_fluxes[0])
plt.plot(range(len(solver.fluxes[0])), solver.edge_fluxes[1])
plt.show()
# f = lambda x: np.cos(x)
# a = 0.0
# b = np.pi/2	