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
experiments = ['bhbt', 'bhri', 'bhbtri']
labels = ['bh1', 'bh2', 'bh3']
colors = ['#00ff26', '#003eff', '#ff0000']
figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

states = ['live',
          'dead',
          'pH',
          'glucose',
          'pyruvate',
          'trehalose',
          'acetate',
          'lactate']


stTypes = ['cells',
           'cells',
           'pH',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite']

for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, v + '.png'), dpi = 150)
    plt.show()
    
for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, 'logos', v + '.png'), dpi = 50)
    plt.show()

####################################################

################### bt #############################

species = 'bt'
experiments = ['bhbt', 'btri', 'bhbtri']
labels = ['bt1', 'bt2', 'bt3']
colors = ['#00ff26', '#003eff', '#ff0000']
figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

states = ['live',
          'dead',
          'pH',
          'glucose',
          'pyruvate',
          'succinate',
          'formate',
          'acetate',
          'lactate']


stTypes = ['cells',
           'cells',
           'pH',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite']

for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, v + '.png'), dpi = 150)
    plt.show()
    
for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, 'logos', v + '.png'), dpi = 50)
    plt.show()






############## ri #####################

species = 'ri'
experiments = ['bhri', 'btri', 'bhbtri']
labels = ['ri1', 'ri2', 'ri3']
colors = ['#00ff26', '#003eff', '#ff0000']
figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

states = ['live',
          'dead',
          'pH',
          'glucose',
          'pyruvate',
          'butyrate',
          'acetate',
          'lactate']


stTypes = ['cells',
           'cells',
           'pH',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite']

for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, v + '.png'), dpi = 150)
    plt.show()
    
for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, 'logos', v + '.png'), dpi = 50)
    plt.show()

#####################################

############## bhbt #####################

species = 'bhbt'
experiments = ['bhbt']
labels = ['bhbt']
colors = ['#003eff']
figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

states = ['live_bh',
          'live_bt',
          'dead',
          'pH',
          'trehalose',
          'glucose',
          'pyruvate',
          'succinate',
          'acetate',
          'lactate',
          'formate']


stTypes = ['cells',
           'cells',
           'cells',
           'pH',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite']

for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, v + '.png'), dpi = 150)
    plt.show()
    
for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, 'logos', v + '.png'), dpi = 50)
    plt.show()

############## bhri #####################

species = 'bhri'
experiments = ['bhri']
labels = ['bhri']
colors = ['#003eff']
figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

states = ['live_bh',
          'live_ri',
          'dead',
          'pH',
          'trehalose',
          'glucose',
          'pyruvate',
          'butyrate',
          'acetate',
          'lactate']


stTypes = ['cells',
           'cells',
           'cells',
           'pH',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite']

for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, v + '.png'), dpi = 150)
    plt.show()
    
for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, 'logos', v + '.png'), dpi = 50)
    plt.show()

############## btri #####################

species = 'btri'
experiments = ['btri']
labels = ['btri']
colors = ['#003eff']
figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

states = ['live_bt',
          'live_ri',
          'dead',
          'pH',
          'butyrate',
          'glucose',
          'pyruvate',
          'succinate',
          'acetate',
          'lactate',
          'formate']


stTypes = ['cells',
           'cells',
           'cells',
           'pH',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite']

for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, v + '.png'), dpi = 150)
    plt.show()
    
for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, 'logos', v + '.png'), dpi = 50)
    plt.show()

############## bhbtri #####################

species = 'bhbtri'
experiments = ['bhbtri']
labels = ['bhbtri']
colors = ['#003eff']
figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

states = ['live_bh',
          'live_bt',
          'live_ri',
          'dead',
          'pH',
          'trehalose',
          'glucose',
          'pyruvate',
          'succinate',
          'butyrate',
          'acetate',
          'lactate',
          'formate']


stTypes = ['cells',
           'cells',
           'cells',
           'cells',
           'pH',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite']

for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, v + '.png'), dpi = 150)
    plt.show()
    
for i,v in enumerate(states):
    pH = makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors)
    plt.savefig(os.path.join(figPath, 'logos', v + '.png'), dpi = 50)
    plt.show()

#####################################




