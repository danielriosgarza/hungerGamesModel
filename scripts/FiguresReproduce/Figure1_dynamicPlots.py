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

params = getPramsFromFile('bh', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'allParamsFitted.tsv'))

databaseName = 'modelDB_bhbtri_test.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

database = os.path.join(databaseFolder, databaseName)

conn = create_connection(os.path.join(databaseFolder, databaseName))
assignBhParams(params, conn)

measuredStates = ['live', 
                  
                  'trehalose',
                  'pyruvate',
                  'glucose',
                  'lactate',
                  'acetate']




bh1 = simulateExperiment(group = 'bh', 
                       experimentLabel = 'bhbt', 
                       dbPath = database, 
                       measuredStates = measuredStates, 
                       )

bh2 = simulateExperiment(group = 'bh', 
                       experimentLabel = 'bhri', 
                       dbPath = database, 
                       measuredStates = measuredStates
                       )

bh3 = simulateExperiment(group = 'bh', 
                       experimentLabel = 'bhbtri', 
                       dbPath = database, 
                       measuredStates = measuredStates
                       )







states = ['live',
          'pH',
          'trehalose',
          'pyruvate',
          'glucose',
          'acetate',
          'lactate']

stTypes = ['cells',
           'pH',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite',
           'metabolite']

figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

for i,v in enumerate(states):
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bh1, bh2, bh3], alpha=0.5, legend=False)
    plt.savefig(os.path.join(figPath, v + 'noLegend_model.png'), dpi = 300)
    plt.show()
    


####################################################################################

species = 'bt'
experiments = ['bhbt', 'btri', 'bhbtri']
labels = ['bt1', 'bt2', 'bt3']
colors = ['#00ff26', '#003eff', '#ff0000']

params = getPramsFromFile('bt', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'allParamsFitted.tsv'))

databaseName = 'modelDB_bhbtri.sqlite3'

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

bt1 = simulateExperiment(group = 'bt', 
                       experimentLabel = 'bhbt', 
                       dbPath = database, 
                       measuredStates = measuredStates
                       )

bt2 = simulateExperiment(group = 'bt', 
                       experimentLabel = 'btri', 
                       dbPath = database, 
                       measuredStates = measuredStates
                       )

bt3 = simulateExperiment(group = 'bt', 
                       experimentLabel = 'bhbtri', 
                       dbPath = database, 
                       measuredStates = measuredStates
                       )


figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')

states = ['live',
          'pH',
          
          'pyruvate',
          'glucose',
          'acetate',
          'lactate',
          'formate',
          'succinate']


stTypes = ['cells',
            'pH',
            'metabolite',
            'metabolite',
            'metabolite',
            'metabolite',
            'metabolite',
            'metabolite']



for i,v in enumerate(states):
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bt1, bt2, bt3], alpha=0.5, legend=False)
    plt.savefig(os.path.join(figPath, v + 'noLegend_model.png'), dpi = 300)
    plt.show()
    


    
####################################################################################

species = 'ri'
experiments = ['bhri', 'btri', 'bhbtri']
labels = ['ri1', 'ri2', 'ri3']
colors = ['#00ff26', '#003eff', '#ff0000']

params = getPramsFromFile('ri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'allParamsFitted.tsv'))

databaseName = 'modelDB_bhbtri.sqlite3'

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

ri1 = simulateExperiment(group = 'ri', 
                       experimentLabel = 'bhri', 
                       dbPath = database, 
                       measuredStates = measuredStates
                       )

ri2 = simulateExperiment(group = 'ri', 
                       experimentLabel = 'btri', 
                       dbPath = database, 
                       measuredStates = measuredStates
                       )

ri3 = simulateExperiment(group = 'ri', 
                       experimentLabel = 'bhbtri', 
                       dbPath = database, 
                       measuredStates = measuredStates
                       )


states = ['live',
          'pH',
          
          'pyruvate',
          'glucose',
          'acetate',
          'lactate',
          'butyrate']


stTypes = ['cells',
            'pH',
            'metabolite',
            'metabolite',
            'metabolite',
            'metabolite',
            'metabolite']

figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')



for i,v in enumerate(states):
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [ri1, ri2, ri3], alpha=0.5, legend=False)
    #plt.savefig(os.path.join(figPath, v + 'noLegend_model.png'), dpi = 300)
    plt.show()
    
####################################################################################

