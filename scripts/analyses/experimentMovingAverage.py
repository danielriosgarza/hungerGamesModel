# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 09:11:43 2023

@author: danie
"""


import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))
from general import *

 
    
ri_live = parseTable(os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'ri', 'dead.tsv'))

live_df = getDFdict(ri_live, 'dead', False)

experiments = ['bhri', 'btri', 'bhbtri']

labels = ['ri1', 'ri2', 'ri3']
colors = ['#00ff26', '#003eff', '#ff0000']

s = summarizeExperiments(live_df, 'ri_live', experiments, interval = 8)
sp = get_spline('dead', 'something', 'dead' , df_state=s)
t = np.linspace(0,120, 1000)

makeExperimentPlot('ri', 'dead', 'cells', experiments, labels, colors)
#plt.plot(s)
plt.plot(t, sp(t))