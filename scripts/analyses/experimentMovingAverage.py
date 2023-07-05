# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 09:11:43 2023

@author: danie
"""


import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))
from general import *
from parseTable import *

#chose a species
species = 'ri'

#chose replicates to summarize
experiments = ['bhri', 'btri', 'bhbtri']


#location of the experiment data folder
strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', species)


#chose the state ('live', 'dead', 'glucose', 'acetate', 'pH', etc.)
measuredState = 'live'

#chose the state type for the plot
stateType = 'cells'

#chose the regular interval
intervals = 4

#get the data
stFile = parseTable(os.path.join(strainSummaryFolder, measuredState +  '.tsv'))
df_state = getDFdict(stFile, measuredState, False)

labels = ['ri1', 'ri2', 'ri3']
colors = ['#00ff26', '#003eff', '#ff0000']

s = summarizeExperiments(df_state, 'ri_live', experiments, interval = intervals)

makeExperimentPlot(species, measuredState, stateType, stFile, labels, colors)
plt.plot(s, lw = 5, ls = '--', color = 'k')
plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'movingAverage.png'), dpi = 100)


#get spline
sp = get_spline('live', strainSummaryFolder, labels, df_state = s)

#define an arbitrary time interval
t = np.linspace(0, 120, 10000)

#get spline points
live = sp(t)

#add to plot
plt.plot(t, live, lw = 8, color='#FF6700', alpha = 0.75)

plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'spline.png'), dpi = 100)