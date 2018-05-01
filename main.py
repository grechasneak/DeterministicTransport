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

	def gen_fissionSource(self): 
		fission_sourceF = []
		for material in self.grid['material']:
			nu2 = self.xs.loc[material, 'nu2']
			sigma_f2 = self.xs.loc[material, 'sigma_f2']
			fission_sourceF.append(nu2 * sigma_f2)
		return np.asarray(fission_sourceF) * self.fluxes[1]	
		
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
		self.update_Source()
		fluxes_ = [np.zeros(len(grid)), np.zeros(len(grid))]
		currents_ = [np.zeros(len(grid)), np.zeros(len(grid))]
		
		for group in range(len(fluxes)):
		
			for idx, mu in enumerate(self.directions)):
				if mu > 0:
					for i in self.cells:
						material = grid['material'].iloc[i]
						sigma_t1 = xs.loc[material, 'sigma_t1']
						sigma_t2 = xs.loc[material, 'sigma_t2']
						
						tau_1 = spacing * sigma_t1 / mu
						tau_2 = spacing * sigma_t2 / mu
					
					# current_out = current_in * np.exp(-tau) + Q_cell / sigma_t * (1 - np.exp(-tau))
					# current_mean = calc_mean_current(current_in, current_out, Q, directions[idx], sigma_t1)
					
					# flux[i] += current_mean * weights[idx]
					# current_cell[i] = current_mean * weights[idx] * directions[idx]
					
					# current_in = current_out
				else:	
					for i in reversed(cells):
						material = grid['material'].iloc[i]
					sigma_t1 = xs.loc[material, 'sigma_t1']
					# tau = spacing * sigma_t1 / directions[idx]
					
					# current_out = calc_current_out(current_in, Q, tau, sigma_t1)
					# current_mean = calc_mean_current(current_in, current_out, Q, directions[idx], sigma_t1)
					
					# flux[i] += current_mean * weights[idx]
					# current_cell[i] += current_mean * weights[idx] * directions[idx]
					# print('left')
					# current_in = current_out
		

	
	def calc_current_out(current_in, Q_cell, tau, sigma_t):
		return current_in * np.exp(-tau) + Q_cell / sigma_t * (1 - np.exp(-tau))
	
	def calc_mean_current(current_in, current_out, Q_cell, mu, sigma_t):
		return Q_cell / sigma_t - mu * (current_out - current_in)/(spacing * sigma_t)
	

	


def gen_fissionSource(grid, xs, flux): 
	fission_sourceF = []
	for material in grid['material']:
		nu2 = xs.loc[material, 'nu2']
		sigma_f2 = xs.loc[material, 'sigma_f2']
		fission_sourceF.append(nu2 * sigma_f2)
	return np.asarray(fission_sourceF) * flux
	
grid = pd.read_csv('grid.csv')
xs = pd.read_csv('XS.csv', index_col = 'material')

global spacing
spacing = 0.15625
order = 20
directions, weights = np.polynomial.legendre.leggauss(order)



Q = 1/100
cells = range(128)

current_in = 0.001

flux = np.zeros(128)
current_cell = np.zeros(128)

flux_solver = Step_characteristic(grid, xs, 20)

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