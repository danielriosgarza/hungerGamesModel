# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 13:34:26 2023

@author: danie
"""

import os
import sys
from pathlib import Path

from aquarel import load_theme

theme = load_theme("boxy_light")
theme.apply()


import matplotlib.pyplot as plt

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))
from general import *


########### bh ###############

species = 'bh'
experiments = ['bhwc1', 'bhwctreh']
labels = ['wc', 'wc+treh']
colors = ['#00ff26', '#003eff']
figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

states = ['glucose',
          'pyruvate',
          'trehalose',
          'acetate',
          'lactate']


stTypes = ['metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite']

for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, v + 'wc+treh.png'), dpi = 150)
    plt.show()
    
for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, 'logos', v + 'wc+treh.png'), dpi = 50)
    plt.show()

####################################################


