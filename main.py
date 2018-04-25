#!/usr/bin/env python
from __future__ import division
import random
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

grid = pd.read_csv('grid.csv')
xs = pd.read_csv('XS.csv', index_col = 'material')