#!/usr/bin/env python
from __future__ import division
import random
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate


grid = pd.read_csv('grid.csv')
xs = pd.read_csv('XS.csv', index_col = 'material')



material = grid.iloc[0, 1]
sigma_t1 = xs.loc[material, 'sigma_t1']
sigma_s11 = xs.loc[material, 'sigma_s11']
sigma_t2 = xs.loc[material, 'sigma_t2']
sigma_s22 = xs.loc[material, 'sigma_s22']

f = lambda x: np.cos(x)
a = 0.0
b = np.pi/2	
	
def gaussian_integration(f, a, b, degree):
	# Gauss-Legendre (default interval is [-1, 1])
	degree = 10
	x, w = np.polynomial.legendre.leggauss(degree)
	# Translate x values from the interval [-1, 1] to [a, b]
	t = 0.5*(x + 1)*(b - a) + a
	gauss = sum(w * f(t)) * 0.5*(b - a)
	return gauss

	
	
def step_characteristic():
	pass

	
	
# For comparison
quad, quad_err = integrate.quad(f, a, b)
gauss = gaussian_integration(f, a, b, 10)
print(quad, quad_err)
print(gauss)
print(gauss - quad)

