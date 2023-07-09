# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 18:16:15 2023

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




coculture = 'bhbtri'

states = ['live_bh', 
          'live_bt', 
          'live_ri',
          'trehalose',
          'pyruvate',
          'glucose',
          'butyrate',
          'succinate',
          'acetate',
          'lactate',
          'formate']



intervals = 4

strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', coculture)


initialStates = {}


for i,v in enumerate(states):
    stFile  = parseTable(os.path.join(strainSummaryFolder, v + '.tsv'))
    df_state  = getDFdict(stFile, v, True)[coculture]
    initialStates[v] = get_initialState(v, strainSummaryFolder, coculture, combined=False, interval = intervals) 







databaseName = 'modelDB_bhbtri_bh.sqlite3'
databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

#Load the database
db = get_database(os.path.join(databaseFolder, databaseName))

r_species = genericSimulation(db)




for state in states:
    makeExperimentPlot(coculture, 
              state, 
              stateType = 'metabolite',
              experiments = [coculture],
              lables = [coculture],
              colors = ['#003eff'],
              simulObj = [r_species],
              alpha=1)
    plt.show()
#plt.ylim(0,1)

states = ['live_bh', 
          'live_bt', 
          'live_ri',
          ]



for state in states:
    makeExperimentPlot(coculture, 
              state, 
              stateType = 'metabolite',
              experiments = [coculture],
              lables = [coculture],
              colors = ['#003eff'],
              simulObj = [r_species],
              alpha=1)
plt.show()