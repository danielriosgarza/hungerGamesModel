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

databaseName = 'modelDB_bhbtri.sqlite3'

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
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bh1, bh2, bh3], alpha=0.5, legend=False)
    plt.savefig(os.path.join(figPath, v + 'noLegend_model.png'), dpi = 300)
    plt.show()
    


for i,v in enumerate(states):
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bh1, bh2, bh3], alpha=0.5, legend=True)
    #plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
    #plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
    plt.show()



####################################################################################

species = 'bt'
experiments = ['bhbt', 'btri', 'bhbtri']
labels = ['bt1', 'bt2', 'bt3']
colors = ['#00ff26', '#003eff', '#ff0000']

params = getPramsFromFile('bhbtri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'btri.tsv'))

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
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bt1, bt2, bt3], alpha=0.5, legend=False)
    plt.savefig(os.path.join(figPath, v + 'noLegend_model.png'), dpi = 300)
    plt.show()
    


for i,v in enumerate(states):
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bt1, bt2, bt3], alpha=0.5, legend=True)
    plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
    #plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
    plt.show()

    
####################################################################################

species = 'ri'
experiments = ['bhri', 'btri', 'bhbtri']
labels = ['ri1', 'ri2', 'ri3']
colors = ['#00ff26', '#003eff', '#ff0000']

params = getPramsFromFile('bhbtri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'btri.tsv'))

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

figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')



for i,v in enumerate(states):
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [ri1, ri2, ri3], alpha=0.5, legend=False)
    plt.savefig(os.path.join(figPath, v + 'noLegend_model.png'), dpi = 300)
    plt.show()
    


for i,v in enumerate(states):
    makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [ri1, ri2, ri3], alpha=0.5, legend=True)
    plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
    plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
    plt.show()

####################################################################################

species = 'bhbt'
experiments = ['bhbt']
labels = ['bhbt']
colors = ['#003eff']

#params = getPramsFromFile('ri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'ri.tsv'))

databaseName = 'modelDB_bhbtri.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

database = os.path.join(databaseFolder, databaseName)

#conn = create_connection(os.path.join(databaseFolder, databaseName))
#assignRiParams(params, conn)






measuredStates = ['live_bh',
                  'live_bt',
          
            'trehalose',
          'pyruvate',
          'glucose',
          'acetate',
          'lactate',
          'succinate',
          'formate']

bhbt = simulateExperiment(group = 'bhbt', 
                       experimentLabel = 'bhbt', 
                       dbPath = database, 
                       measuredStates = measuredStates
                       )


states = ['live_bh',
          'live_bt',
          'dead',
          
          'pH',
          
          'trehalose',
        'pyruvate',
        'glucose',
        'acetate',
        'lactate',
        'succinate',
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

figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')





for i,v in enumerate(states):
    if (v=='live_bh') or (v=='live_bt'):
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bhbt, None, None], alpha=0.5, legend=False, ylim = [0,13])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()
    
    elif v=='pH':
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bhbt, None, None], alpha=0.5, legend=False, ylim = [0,8])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()
    
    elif (v=='glucose') or (v=='pyruvate') or (v=='trehalose'):
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bhbt, None, None], alpha=0.5, legend=False, ylim = [0,9])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()
    
    else:
 
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bhbt, None, None], alpha=0.5, legend=False, ylim = [0,26])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()


    
    
####################################################################################

species = 'bhri'
experiments = ['bhri']
labels = ['bhri']
colors = ['#003eff']

#params = getPramsFromFile('ri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'ri.tsv'))

databaseName = 'modelDB_bhbtri.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

database = os.path.join(databaseFolder, databaseName)

#conn = create_connection(os.path.join(databaseFolder, databaseName))
#assignRiParams(params, conn)






measuredStates = ['live_bh',
                  'live_ri',
          
            'trehalose',
          'pyruvate',
          'glucose',
          'acetate',
          'lactate',
          'butyrate'
          ]

bhri = simulateExperiment(group = 'bhri', 
                       experimentLabel = 'bhri', 
                       dbPath = database, 
                       measuredStates = measuredStates
                       )


states = ['live_bh',
          'live_ri',
          'dead',
          
          'pH',
          
          'trehalose',
        'pyruvate',
        'glucose',
        'acetate',
        'lactate',
        'butyrate']


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
            ]

figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')


for i,v in enumerate(states):
    if (v=='live_bh') or (v=='live_ri'):
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bhri, None, None], alpha=0.5, legend=False, ylim = [0,13])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()
    
    elif v=='pH':
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bhri, None, None], alpha=0.5, legend=False, ylim = [0,8])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()
    
    elif (v=='glucose') or (v=='pyruvate') or (v=='trehalose'):
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bhri, None, None], alpha=0.5, legend=False, ylim = [0,9])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()
    
    else:
 
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bhri, None, None], alpha=0.5, legend=False, ylim = [0,26])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()    
    
####################################################################################

species = 'btri'
experiments = ['btri']
labels = ['btri']
colors = ['#003eff']

#params = getPramsFromFile('ri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'ri.tsv'))

databaseName = 'modelDB_bhbtri.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

database = os.path.join(databaseFolder, databaseName)

#conn = create_connection(os.path.join(databaseFolder, databaseName))
#assignRiParams(params, conn)






measuredStates = ['live_bt',
                  'live_ri',
          
            
          'pyruvate',
          'glucose',
          'acetate',
          'lactate',
          'butyrate',
          'succinate',
          'formate'
          ]

btri = simulateExperiment(group = 'btri', 
                       experimentLabel = 'btri', 
                       dbPath = database, 
                       measuredStates = measuredStates
                       )


states = ['live_bt',
          'live_ri',
          'dead',
          
          'pH',
          
          
        'pyruvate',
        'glucose',
        'acetate',
        'lactate',
        'butyrate',
        'succinate',
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
            'metabolite'
            ]

figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')
for i,v in enumerate(states):
    if (v=='live_bt') or (v=='live_ri'):
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [btri, None, None], alpha=0.5, legend=False, ylim = [0,13])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()
    
    elif v=='pH':
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [btri, None, None], alpha=0.5, legend=False, ylim = [0,8])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()
    
    elif (v=='glucose') or (v=='pyruvate') or (v=='trehalose'):
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [btri, None, None], alpha=0.5, legend=False, ylim = [0,9])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()
    
    else:
 
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [btri, None, None], alpha=0.5, legend=False, ylim = [0,26])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()   

####################################################################################

species = 'bhbtri'
experiments = ['bhbtri']
labels = ['bhbtri']
colors = ['#003eff']

#params = getPramsFromFile('ri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'ri.tsv'))

databaseName = 'modelDB_bhbtri.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

database = os.path.join(databaseFolder, databaseName)

#conn = create_connection(os.path.join(databaseFolder, databaseName))
#assignRiParams(params, conn)






measuredStates = ['live_bh',
                  'live_bt',
                  'live_ri',
                  
                  'trehalose',
          
            
          'pyruvate',
          'glucose',
          'acetate',
          'lactate',
          'butyrate',
          'succinate',
          'formate'
          ]

bhbtri = simulateExperiment(group = 'bhbtri', 
                       experimentLabel = 'bhbtri', 
                       dbPath = database, 
                       measuredStates = measuredStates
                       )


states = ['live_bh',
          'live_bt',
          'live_ri',
          'dead',
          
          'pH',
          
          'trehalose',
          
          
        'pyruvate',
        'glucose',
        'acetate',
        'lactate',
        'butyrate',
        'succinate',
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
            'metabolite'
            ]

figPath = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', species+'Experiments')
for i,v in enumerate(states):
    if (v=='live_bh') or (v=='live_bt') or (v=='live_ri'):
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bhbtri, None, None], alpha=0.5, legend=False, ylim = [0,13])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()
    
    elif v=='pH':
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bhbtri, None, None], alpha=0.5, legend=False, ylim = [0,8])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()
    
    elif (v=='glucose') or (v=='pyruvate') or (v=='trehalose'):
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bhbtri, None, None], alpha=0.5, legend=False, ylim = [0,9])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()
    
    else:
 
        makeExperimentPlot(species, v, stTypes[i], experiments, labels, colors, simulObj = [bhbtri, None, None], alpha=0.5, legend=False, ylim = [0,26])
        plt.savefig(os.path.join(figPath, v + '_model.png'), dpi = 300)
        plt.savefig(os.path.join(figPath, 'logos', v + '_model.png'), dpi = 50)
        plt.show()   