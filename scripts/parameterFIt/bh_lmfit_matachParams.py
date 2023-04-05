# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 08:54:49 2022

@author: u0139894
"""

from pathlib import Path
import os
import sys

from scipy.interpolate import PchipInterpolator as CubicSpline
from lmfit import minimize, Parameters, fit_report


sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'core'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'db'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))


from mainClasses import *
from parseTable import *
from updateParameters import *
from readModelDB import *
from loadParameters import *
from general import *


def pseudoHuberLoss(y_true, y_pred, delta = 0.50):
    """
    Compute the Pseudo-Huber loss between y_true and y_pred with a given delta value.
    """
    
    choice = np.random.choice(np.arange(len(y_true)), size = int(1.0*len(y_true)), replace=False)
    
    y_t = y_true[choice]
    y_p = y_pred[choice]
    
    # Compute the squared error
    error = y_t - y_p
    squared_error = np.square(error)

    # Compute the loss for small errors
    small_error = delta**2 * (np.sqrt(1 + squared_error / delta**2) - 1)

    # For values of error greater than delta, use the linear part of the Pseudo-Huber loss
    large_error = delta * (np.sqrt(squared_error) - 0.5 * delta)

    # Return the mean loss
    return np.mean(np.where(np.abs(error) <= delta, small_error, large_error))

def sResidual(y_true, y_pred):
    choice = np.random.choice(np.arange(len(y_true)), size = int(0.3*len(y_true)), replace=False)
    y_t = y_true[choice]
    y_p = y_pred[choice]
    
    return np.sqrt(np.mean((y_t-y_p)**2))
    
    

####################################################################
def writeOutput(lmfit_params, outputFile):
    with open(outputFile, 'w') as f:
        f.write('species\tparameter\tvalue\tmin\tmax\n')
        
        for i in lmfit_params:
            f.write(species + '\t' + i + '\t' + str(lmfit_params[i].value) + '\t' + str(lmfit_params[i].min) + '\t' + str(lmfit_params[i].max) + '\n')




def distance(lmfit_params, database, initialStates, measuredStates, splines, experimentLabel, species, intervals, combined = True):
    
    conn = create_connection(database)
    
    assignBhParams(lmfit_params, conn)
    
    db = get_database(database)
    
    r = simulateExperiment(species, 
                           experimentLabel, 
                           lmfit_params, 
                           database, 
                           initialStates,
                           combined=combined, 
                           intervals=intervals)
    
    
    distances = []
    
    for i in measuredStates:
        if i=='live':
            distances.append(10*pseudoHuberLoss(splines['live'](r.time_simul), r.cellActive_dyn[0]))
        
        elif i=='dead':
            
            distances.append(pseudoHuberLoss(splines['dead'](r.time_simul), r.cellInactive_dyn[0]))
        
        elif i=='pH':
            distances.append(pseudoHuberLoss(splines['pH'](r.time_simul), r.pH_simul))
        
        elif i=='glucose':
            distances.append(10*pseudoHuberLoss(splines[i](r.time_simul), r.met_simul[r.metabolome.metabolites.index(i)]))
        
        elif i=='acetate':
            distances.append(pseudoHuberLoss(splines[i](r.time_simul), r.met_simul[r.metabolome.metabolites.index(i)]))
        
        else:
            distances.append(pseudoHuberLoss(splines[i](r.time_simul), r.met_simul[r.metabolome.metabolites.index(i)]))
    
    
    
    objV = sum(distances)
   
    ssrSum = np.round(objV,3)
    
    print(ssrSum)
    if len(evals)>0:    
        if ssrSum<min(evals):
            writeOutput(lmfit_params, outputFile)
            plt.plot(evals)
            plt.title(str(ssrSum))
            plt.show()
    evals.append(ssrSum)
    return objV



##############SETUP###########################################
ipH_path = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bhbtri_ipH4.tsv') 

species = 'bh'

experimentLabel = ['bhbt', 'bhri', 'bhbtri']

strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh')

inputParams = getPramsFromFile('bh', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bh.tsv'))


outputFile = os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bh.tsv')


databaseName = 'modelDB_bhbtri_bh.sqlite3'

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

initialStates = ['live',
                 'trehalose',
                 'pyruvate',
                 'glucose',
                 'lactate',
                 'acetate'
                 ]

#####################################################################







####get splines
splines = {}

for i,v in enumerate(measuredStates):
    
    stFile = parseTable(os.path.join(strainSummaryFolder, v +  '.tsv'))
    df_state = getDFdict(stFile, v, False)
    summ_state = summarizeExperiments(df_state, v, experimentLabel, interval = intervals[i])
    splines[v] = get_spline(v, 'nothing', experimentLabel, df_state = summ_state)



evals = []






    
out = minimize(distance, params=inputParams, method='powell', kws = {'database' : database,
                                                                  'initialStates' : initialStates,
                                                                  'measuredStates' : measuredStates,
                                                                  'splines': splines,
                                                                  'experimentLabel':experimentLabel,
                                                                  'species' : species,
                                                                  'intervals': intervals,
                                                                  'combined': True,
                                                                  })
    
