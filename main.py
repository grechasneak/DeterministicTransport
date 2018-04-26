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
	
def calc_flux(currents, weights):
	return sum(currents * weights.reshape(len(weights), -1))
	
def calc_average_current(currents, weights, directions):
	return (currents * (weights.reshape(len(weights), -1) * directions.reshape(len(directions), -1)))
	
def generate_currentArray(directions):
	currents = []
	for i in directions:
		currents.append(np.zeros(128))
	return np.asarray(currents)
	
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
order = 5
directions, weights = np.polynomial.legendre.leggauss(order)



Q = 1/10000
cells = range(128)

current_in = 0.1


current_out = np.zeros(128)
current_mean = np.zeros(128)

currents = generate_currentArray(directions)
currents_mean = generate_currentArray(directions)

for idx in range(len(directions)):
	for i in cells:
		material = grid['material'].iloc[i]
		sigma_t1 = xs.loc[material, 'sigma_t1']
		tau = spacing * sigma_t1 / abs(directions[idx])
		currents[idx][i] = calc_current_out(current_in, Q, directions[idx], sigma_t1) 
		currents_mean[idx][i] = calc_mean_current(current_in, currents[idx][i], Q, tau, sigma_t1)
		current_in = currents[idx][i]



flux = calc_flux(currents, weights)
average_current = calc_average_current(currents, weights, directions)
# f = lambda x: np.cos(x)
# a = 0.0
# b = np.pi/2	