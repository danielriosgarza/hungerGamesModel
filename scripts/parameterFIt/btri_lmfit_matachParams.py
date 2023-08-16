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
from scipy.stats import pearsonr

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
    choice = np.random.choice(np.arange(len(y_true)), size = int(1.0*len(y_true)), replace=False)
    y_t = y_true[choice]
    y_p = y_pred[choice]
    
    return np.sqrt(np.mean((y_t-y_p)**2))
    
def custom_loss(y_true, y_pred):
    # Compute MSE
    mse_loss = sResidual(y_true, y_pred)
    
    # Compute Pearson correlation coefficient
    corr, _ = pearsonr(y_true, y_pred)
    
    # Return a combination of MSE and negative correlation
    return (5*mse_loss) - corr     

####################################################################
def writeOutput(lmfit_params, outputFile):
    with open(outputFile, 'w') as f:
        f.write('species\tparameter\tvalue\tmin\tmax\n')
        
        for i in lmfit_params:
            f.write(species + '\t' + i + '\t' + str(lmfit_params[i].value) + '\t' + str(lmfit_params[i].min) + '\t' + str(lmfit_params[i].max) + '\n')




def distance(lmfit_params, database, initialStates, measuredStates, splines, experimentLabel, species, intervals, combined = True):
    
    conn = create_connection(database)
    
    assignBtParams(lmfit_params, conn)
    assignRiParams(lmfit_params, conn)
    
    
    
    
    db = get_database(database)
    
    r = simulateExperiment(group = species, 
                           experimentLabel = experimentLabel, 
                           dbPath = database, 
                           measuredStates = initialStates,
                           combined=combined, 
                           intervals=intervals)
    
    r2 = simulateExperiment(group = 'bt', 
                           experimentLabel = experimentLabel_bt, 
                           dbPath = database, 
                           measuredStates = initialStates_bt,
                           combined=combined, 
                           intervals=intervals_bt)
    
    r3 = simulateExperiment(group = 'ri', 
                           experimentLabel = experimentLabel_ri, 
                           dbPath = database, 
                           measuredStates = initialStates_ri,
                           combined=combined, 
                           intervals=intervals_ri)
    
    
    distances = []
    
    
    
    
    for i in measuredStates_bt:
        if i=='live':
            distances.append(custom_loss(splines_bt['live'](r2.time_simul), np.sum(r2.cellActive_dyn,axis=0)))
        
        elif i=='dead':
            
            distances.append(pseudoHuberLoss(splines_bt['dead'](r2.time_simul), np.sum(r2.cellInactive_dyn,axis=0)))
        
        elif i=='pH':
            distances.append(pseudoHuberLoss(splines_bt['pH'](r2.time_simul), r2.pH_simul))
        
        else:
            distances.append(pseudoHuberLoss(splines_bt[i](r2.time_simul), r2.met_simul[r2.metabolome.metabolites.index(i)]))
    
    for i in measuredStates_ri:
        if i=='live':
            distances.append(custom_loss(splines_ri['live'](r3.time_simul), np.sum(r3.cellActive_dyn,axis=0)))
            
        
        elif i=='dead':
            
            distances.append(pseudoHuberLoss(splines_ri['dead'](r3.time_simul), np.sum(r3.cellInactive_dyn,axis=0)))
        
        elif i=='pH':
            distances.append(pseudoHuberLoss(splines_ri['pH'](r3.time_simul), r3.pH_simul))
        
        else:
            distances.append(pseudoHuberLoss(splines_ri[i](r3.time_simul), r3.met_simul[r3.metabolome.metabolites.index(i)]))
    
    
    
    for i in measuredStates:
        if i=='live_bh':
            distances.append(custom_loss(splines['live_bh'](r.time_simul), r.cellActive_dyn[0]))
        
        elif i=='live_bt':
            distances.append(custom_loss(splines['live_bt'](r.time_simul), r.cellActive_dyn[1]))
            
        elif i=='live_ri':
            distances.append(custom_loss(splines['live_ri'](r.time_simul), r.cellActive_dyn[2]))
            
        
        elif i=='dead':
            
            distances.append(pseudoHuberLoss(splines['dead'](r.time_simul), np.sum(r.cellInactive_dyn,axis=0)))
        
        elif i=='pH':
            distances.append(pseudoHuberLoss(splines['pH'](r.time_simul), r.pH_simul))
        
        
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

species = 'bhbtri'

experimentLabel = ['bhbtri']

strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bhbtri')

inputParams = getPramsFromFile('btri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'btri.tsv'))



outputFile = os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'btri.tsv')



databaseName = 'modelDB_bhbtri_bh.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')


database = os.path.join(databaseFolder, databaseName)

measuredStates = ['live_bh',
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


intervals = [4,
             4,
             4,
             4,
             4,
             
             4,
             4,
             4,
             4,
             4,
             4,
             4,
             4]

initialStates = ['live_bh',
                 'live_bt',
                 'live_ri',
                 'trehalose',
                 'pyruvate',
                 'glucose',
                 'lactate',
                 'acetate',
                 'butyrate',
                 'succinate',
                 'formate'
                 ]


#####################################################################


####get splines
splines = {}

for i,v in enumerate(measuredStates):
    
    stFile = parseTable(os.path.join(strainSummaryFolder, v +  '.tsv'))
    df_state = getDFdict(stFile, v, False)
    summ_state = summarizeExperiments(df_state, v, experimentLabel, interval = intervals[i])
    splines[v] = get_spline(v, 'nothing', experimentLabel, df_state = summ_state)




#################bt###################################
strainSummaryFolder_bt = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bt')
experimentLabel_bt = ['bhbt', 'btri', 'bhbtri']


measuredStates_bt = ['live', 
                  'dead',
                  'pH',
                  
                  
                  'pyruvate',
                  'glucose',
                  'lactate',
                  'acetate',
                  'formate',
                  'succinate']


intervals_bt = [4,
             12,
             4,
             
             4,
             4,
             16,
             4,
             4,
             4]

initialStates_bt = ['live',
                 'pyruvate',
                 'glucose',
                 'lactate',
                 'acetate',
                 'formate',
                 'succinate'
                 ]

####get splines
splines_bt = {}

for i,v in enumerate(measuredStates_bt):
    
    stFile = parseTable(os.path.join(strainSummaryFolder_bt, v +  '.tsv'))
    df_state = getDFdict(stFile, v, False)
    summ_state = summarizeExperiments(df_state, v, experimentLabel_bt, interval = intervals_bt[i])
    splines_bt[v] = get_spline(v, 'nothing', experimentLabel_bt, df_state = summ_state)
    
 ##########################################################################################


##############ri################################################
experimentLabel_ri = ['bhri', 'btri', 'bhbtri']

strainSummaryFolder_ri = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'ri')

measuredStates_ri = ['live', 
                  'dead',
                  'pH',
                  
                  
                  'pyruvate',
                  'glucose',
                  'lactate',
                  'acetate',
                  'butyrate']


intervals_ri = [4,
             12,
             4,
             
             4,
             4,
             16,
             4,
             4]

initialStates_ri = ['live',
                 'pyruvate',
                 'glucose',
                 'lactate',
                 'acetate',
                 'butyrate'
                 ]

####get splines
splines_ri = {}

for i,v in enumerate(measuredStates_ri):
    
    stFile = parseTable(os.path.join(strainSummaryFolder_ri, v +  '.tsv'))
    df_state = getDFdict(stFile, v, False)
    summ_state = summarizeExperiments(df_state, v, experimentLabel_ri, interval = intervals_ri[i])
    splines_ri[v] = get_spline(v, 'nothing', experimentLabel_ri, df_state = summ_state)




################################################################   







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
    
