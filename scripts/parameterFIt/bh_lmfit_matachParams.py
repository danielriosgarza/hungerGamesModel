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
    
def rmsle(y_true, y_pred):
    # Ensure the predicted and true values have no negative values
    assert (y_true >= 0).all() 
    assert (y_pred >= 0).all()
    
    # Compute RMSLE
    log_error = np.log1p(y_pred) - np.log1p(y_true)
    return np.sqrt(np.mean(log_error**2))    



def custom_loss(y_true, y_pred):
    # Compute MSE
    mse_loss = sResidual(y_true, y_pred)
    
    # Compute Pearson correlation coefficient
    corr, _ = pearsonr(y_true, y_pred)
    
    # Return a combination of MSE and negative correlation
    return mse_loss - corr    

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
    
    r = simulateExperiment(group = species, 
                           experimentLabel = experimentLabel, 
                           dbPath = database, 
                           measuredStates = initialStates,
                           combined=combined, 
                           intervals=intervals)
    
    
    r2 = genericSimulation(db)
    
    measuredStates_treh = [
                      'trehalose',
                      'pyruvate',
                      'glucose',
                      'lactate',
                      'acetate']

   
    
    
    distances = []
    distances.append(pseudoHuberLoss(coSpline(r2.time_simul), r2.cellActive_dyn[0]))
    distances.append(pseudoHuberLoss(coSpline_treh(r2.time_simul), r2.met_simul[r.metabolome.metabolites.index('trehalose')]))
    
    #distances.append(sResidual(coSpline_treh_treh(bhwctreh.time_simul), bhwctreh.met_simul[r.metabolome.metabolites.index('trehalose')]))
    #distances.append(sResidual(5*coSpline_treh_glc(bhwctreh.time_simul), bhwctreh.met_simul[r.metabolome.metabolites.index('glucose')]))
    #distances.append(pseudoHuberLoss(coSpline_treh_ace(bhwctreh.time_simul), bhwctreh.met_simul[r.metabolome.metabolites.index('acetate')]))
    #distances.append(pseudoHuberLoss(coSpline_treh_lac(bhwctreh.time_simul), bhwctreh.met_simul[r.metabolome.metabolites.index('lactate')]))
    
    for i in measuredStates:
        if i=='live':
            distances.append(5*custom_loss(splines['live'](r.time_simul), np.sum(r.cellActive_dyn,axis=0)))
            
            
        
        elif i=='dead':
            
            distances.append(pseudoHuberLoss(splines['dead'](r.time_simul), np.sum(r.cellInactive_dyn,axis=0)))
        
        elif i=='pH':
            distances.append(pseudoHuberLoss(splines['pH'](r.time_simul), r.pH_simul))
        
        elif i=='glucose':
            distances.append(5*custom_loss(splines[i](r.time_simul), r.met_simul[r.metabolome.metabolites.index(i)]))
        
        elif i=='lactate':
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

inputParams = getPramsFromFile('bh', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bh_2.txt'))


outputFile = os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bh_2.txt')


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



strainSummaryFolder_cc = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bhbtri')
stFile  = parseTable(os.path.join(strainSummaryFolder_cc, 'live_bh' + '.tsv'))
df_state  = getDFdict(stFile, 'live_bh', False)['bhbtri']

timeV = np.array(df_state.index)

stateMean = np.array(df_state.mean(axis=1))


coSpline = CubicSpline(timeV, stateMean, extrapolate=False)

strainSummaryFolder_cc = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bhbtri')
stFile  = parseTable(os.path.join(strainSummaryFolder_cc, 'trehalose' + '.tsv'))
df_state  = getDFdict(stFile, 'trehalose', False)['bhbtri']

timeV = np.array(df_state.index)

stateMean = np.array(df_state.mean(axis=1))


coSpline_treh = CubicSpline(timeV, stateMean, extrapolate=False)



strainSummaryFolder_cc = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh')
stFile  = parseTable(os.path.join(strainSummaryFolder_cc, 'trehalose' + '.tsv'))
df_state  = getDFdict(stFile, 'trehalose', False)['bhwctreh']

timeV = np.array(df_state.index)

stateMean = np.array(df_state.mean(axis=1))


coSpline_treh_treh = CubicSpline(timeV, stateMean, extrapolate=False)


strainSummaryFolder_cc = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh')
stFile  = parseTable(os.path.join(strainSummaryFolder_cc, 'glucose' + '.tsv'))
df_state  = getDFdict(stFile, 'glucose', False)['bhwctreh']

timeV = np.array(df_state.index)

stateMean = np.array(df_state.mean(axis=1))


coSpline_treh_glc = CubicSpline(timeV, stateMean, extrapolate=False)



strainSummaryFolder_cc = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh')
stFile  = parseTable(os.path.join(strainSummaryFolder_cc, 'acetate' + '.tsv'))
df_state  = getDFdict(stFile, 'acetate', False)['bhwctreh']

timeV = np.array(df_state.index)

stateMean = np.array(df_state.mean(axis=1))


coSpline_treh_ace = CubicSpline(timeV, stateMean, extrapolate=False)


strainSummaryFolder_cc = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh')
stFile  = parseTable(os.path.join(strainSummaryFolder_cc, 'lactate' + '.tsv'))
df_state  = getDFdict(stFile, 'lactate', False)['bhwctreh']

timeV = np.array(df_state.index)

stateMean = np.array(df_state.mean(axis=1))


coSpline_treh_lac = CubicSpline(timeV, stateMean, extrapolate=False)

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
    
