# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 14:56:14 2023

@author: danie
"""

import os
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))

from parseTable import *
from general import *

species = 'bh'
experiments = ['bhbt', 'bhri', 'bhbtri']
labels = ['bh1', 'bh2', 'bh3']
colors = ['#00ff26', '#003eff', '#ff0000']

params = getPramsFromFile('bh', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bh.tsv'))

databaseName = 'modelDB_bhbtri_bh.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

database = os.path.join(databaseFolder, databaseName)

conn = create_connection(os.path.join(databaseFolder, databaseName))
assignBhParams(params, conn)

measuredStates = ['live', 
                  'dead',
                  'pH',
                  
                  'trehalose',
                  'pyruvate',
                  'glucose',
                  'lactate',
                  'acetate']

intervals = [4,
             12,
             4,
             
             4,
             4,
             4,
             16,
             4]


bh1 = simulateExperiment('bh', 
                       'bhbt', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=intervals)

bh2 = simulateExperiment('bh', 
                       'bhri', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=intervals)

bh3 = simulateExperiment('bh', 
                       'bhbtri', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=intervals)







states = ['live',
          'dead',
          'pH',
          'trehalose',
          'pyruvate',
          'glucose',
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

figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

for i,v in enumerate(states):
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bh1, bh2, bh3], alpha=0.5)
    #plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 150)
    #plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
    #plt.show()
    

measuredStates = [
                  'trehalose',
                  'pyruvate',
                  'glucose',
                  'lactate',
                  'acetate']

bhwctreh = simulateExperiment('bh', 
                       'bhwctreh', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=intervals,
                       starttime=17,
                       endtime = 90)

bhwc1 = simulateExperiment('bh', 
                       'bhwc1', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=intervals)
