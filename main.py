#!/usr/bin/env python
from __future__ import division
import random
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate


def gaussian_integration(f, a, b, degree):
	# Gauss-Legendre (default interval is [-1, 1])
	degree = 10
	x, w = np.polynomial.legendre.leggauss(degree)
	# Translate x values from the interval [-1, 1] to [a, b]
	t = 0.5 * (x + 1) * (b - a) + a
	gauss = sum(w * f(t)) * 0.5 * (b - a)
	return gauss
	
	
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


for i in range(100):
	for idx in range(len(directions)):
		if directions[idx] > 0:
			for i in cells:
				material = grid['material'].iloc[i]
				sigma_t1 = xs.loc[material, 'sigma_t1']
				tau = spacing * sigma_t1 / directions[idx]
				
				current_out = calc_current_out(current_in, Q, tau, sigma_t1)
				current_mean = calc_mean_current(current_in, current_out, Q, directions[idx], sigma_t1)
				
				flux[i] += current_mean * weights[idx]
				current_cell[i] = current_mean * weights[idx] * directions[idx]
				
				current_in = current_out
		else:	
			for i in reversed(cells):
				material = grid['material'].iloc[i]
				sigma_t1 = xs.loc[material, 'sigma_t1']
				tau = spacing * sigma_t1 / directions[idx]
				
				current_out = calc_current_out(current_in, Q, tau, sigma_t1)
				current_mean = calc_mean_current(current_in, current_out, Q, directions[idx], sigma_t1)
				
				flux[i] += current_mean * weights[idx]
				current_cell[i] += current_mean * weights[idx] * directions[idx]
				print('left')
				current_in = current_out
			
	
		


plt.plot(flux)
plt.show()
# f = lambda x: np.cos(x)
# a = 0.0
# b = np.pi/2	