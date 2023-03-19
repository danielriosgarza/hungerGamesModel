# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 15:48:54 2023

@author: danie
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))

from general import *



species = 'bh'
experiments = ['bhbt', 'bhri', 'bhbtri']
labels = ['bh1', 'bh2', 'bh3']
colors = ['#00ff26', '#003eff', '#ff0000']
figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

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



params = getPramsFromFile('bh', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bh0.tsv'))

databaseName = 'modelDB_bhA.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')


database = os.path.join(databaseFolder, databaseName)

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


r1 = simulateExperiment(species, experiments, params, database, ['live',
                                                               'trehalose',
                                                               'pyruvate',
                                                               'glucose',
                                                               'lactate',
                                                               'acetate'], combined = True, intervals = intervals )






for i,v in enumerate(states):
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [r1, r1, r1], alpha=0.05)
    
    
    stFile = parseTable(os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh', v + '.tsv'))
    
    st_df = getDFdict(stFile, v, False)
    
    s = summarizeExperiments(st_df, v, experiments, interval = intervals[i])
    sp = get_spline('dead', 'something', 'dead' , df_state=s)
    t = np.linspace(0,120, 1000)

    plt.plot(s)
    plt.plot(t, sp(t))
