#!/usr/bin/env python
from __future__ import division
import random
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

grid = pd.read_csv('grid.csv')
xs = pd.read_csv('XS.csv', index_col = 'material')



material = grid.iloc[1, 1]
sigma_t1 = xs.loc[material, 'sigma_t1']
sigma_s11 = xs.loc[material, 'sigma_s11']
sigma_t2 = xs.loc[material, 'sigma_t2']
sigma_s22 = xs.loc[material, 'sigma_s22']