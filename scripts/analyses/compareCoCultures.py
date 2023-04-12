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
plt.style.use('seaborn-bright')
import seaborn as sns


sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))

from parseTable import *
from general import *



coculture = 'bhbtri'

state = 'butyrate'


strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', coculture)

stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
        
stateDF = getDFdict(stateTable, state, True)[coculture]
#sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF) 



makeExperimentPlot(coculture, 
             state, 
             stateType = 'metabolite',
             experiments = [coculture],
             lables = [coculture],
             colors = ['#003eff'],
             simulObj = [None],
             alpha=1)
#plt.ylim(0,1)
