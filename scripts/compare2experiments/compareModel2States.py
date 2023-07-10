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




bh1 = simulateExperiment('bh', 
                       'bhbt', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=None)

bh2 = simulateExperiment('bh', 
                       'bhri', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=None)

bh3 = simulateExperiment('bh', 
                       'bhbtri', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=None)







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
    plt.show()
    


####################################################################################

species = 'bt'
experiments = ['bhbt', 'btri', 'bhbtri']
labels = ['bt1', 'bt2', 'bt3']
colors = ['#00ff26', '#003eff', '#ff0000']

params = getPramsFromFile('bt', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bt.tsv'))

databaseName = 'modelDB_bhbtri_bh.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

database = os.path.join(databaseFolder, databaseName)

conn = create_connection(os.path.join(databaseFolder, databaseName))
assignBtParams(params, conn)






measuredStates = ['live',
          
          
          'pyruvate',
          'glucose',
          'acetate',
          'lactate',
          'formate',
          'succinate']

bt1 = simulateExperiment('bt', 
                       'bhbt', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=None)

bt2 = simulateExperiment('bt', 
                       'btri', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=None)

bt3 = simulateExperiment('bt', 
                       'bhbtri', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=None)

figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

states = ['live',
          'dead',
          'pH',
          
          'pyruvate',
          'glucose',
          'acetate',
          'lactate',
          'formate',
          'succinate']


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
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bt1, bt2, bt3], alpha=0.5)
    #plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 150)
    #plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
    plt.show()
    
####################################################################################

species = 'ri'
experiments = ['bhri', 'btri', 'bhbtri']
labels = ['ri1', 'ri2', 'ri3']
colors = ['#00ff26', '#003eff', '#ff0000']

params = getPramsFromFile('ri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'ri.tsv'))

databaseName = 'modelDB_bhbtri_bh.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

database = os.path.join(databaseFolder, databaseName)

conn = create_connection(os.path.join(databaseFolder, databaseName))
assignRiParams(params, conn)






measuredStates = ['live',
          
          'pyruvate',
          'glucose',
          'acetate',
          'lactate',
          'butyrate']

ri1 = simulateExperiment('ri', 
                       'bhri', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=None)

ri2 = simulateExperiment('ri', 
                       'btri', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=None)

ri3 = simulateExperiment('ri', 
                       'bhbtri', 
                       params, 
                       database, 
                       measuredStates, 
                       combined=False, 
                       intervals=None)

figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

states = ['live',
          'dead',
          'pH',
          
          'pyruvate',
          'glucose',
          'acetate',
          'lactate',
          'butyrate']


stTypes = ['cells',
            'cells',
            'pH',
            'metabolite',
            'metabolite',
            'metabolite',
            'metabolite',
            'metabolite']


for i,v in enumerate(states):
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [ri1, ri2, ri3], alpha=0.5)
    plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 150)
    plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
    plt.show()
