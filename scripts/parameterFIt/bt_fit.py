# -*- coding: utf-8 -*-
"""
Created on Sun Jul 20 17:56:31 2025

@author: drgarza
"""

from pathlib import Path
import os
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator as CubicSpline
from lmfit import minimize, Parameters, fit_report
from scipy.stats import pearsonr

import matplotlib.pyplot as plt


sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'core'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'db'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))


from mainClasses import *
from parseTable import *
from updateParameters import *
from readModelDB import *
from loadParameters import *
from general import *


def writeOutput(lmfit_params, outputFile):
    with open(outputFile, 'w') as f:
        f.write('species\tparameter\tvalue\tmin\tmax\n')

        for i in lmfit_params:
            f.write(species + '\t' + i + '\t' + str(lmfit_params[i].value) + '\t' + str(lmfit_params[i].min) + '\t' + str(lmfit_params[i].max) + '\n')

     


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


def distance(lmfit_params, database, measuredStates, splines, experimentLabel, species, intervals, combined = True):
  '''
  lmfit_params: lmfit_params object containing the parameter id and upper and lower bound
  database: database that will be populated with the parameters
  initialStates: initial state used to fit
  measuredStates: measured states used to fit
  splines: dictionary of splines representing the experimental data
  species: which species we are fitting
  intervals: intervals to use in the moving average
  combined: wether to use the combined experimental data or not
  '''
  #connect to the database
  conn = create_connection(database)

  #assign the parameters to from the lmfit obj to the database
  assignBtParams(lmfit_params, conn)

  #load the database
  db = get_database(database)


  #simulate the experiments with the model parameters and initial values from the experiments
  r = simulateExperiment(group = species,
                         experimentLabel = experimentLabel,
                         dbPath = database,
                         measuredStates = measuredStates,
                         combined=combined,
                         intervals=intervals)



  #get the distances from the experimental group splines and the models simulation for each measured state
  distances = []

  for i in measuredStates:
    if i == 'live':
      distances.append(5*pseudoHuberLoss(splines['live'](r.time_simul), np.sum(r.cellActive_dyn,axis=0)))
    elif i == 'dead':
      distances.append(pseudoHuberLoss(splines['dead'](r.time_simul), np.sum(r.cellInactive_dyn,axis=0)))
    elif i == 'pH':
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

#file with acids concentration and pH values to setup the elastic net
ipH_path = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bhbtri_ipH4.tsv')

#species to fit the parameters
species = 'bt'

#experiments to fit the parameters on
experimentLabel = ['bhbt', 'btri', 'bhbtri']

#folder from where to get the experimental data
strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bt')

# load the parameters into an lmfit obj
inputParams = lmfit_params = getPramsFromFile('bt', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'allParams2Fit.tsv'))

#database to write the parameters
databaseName = 'modelDB_fitting_bt.sqlite3'
databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')
database = os.path.join(databaseFolder, databaseName)


#the states we measured that relevant for Bacteroides thetaiotaomicron
measuredStates = ['live',
                  'dead',
                  'pH',
                  'pyruvate',
                  'glucose',
                  'lactate',
                  'acetate',
                  'formate',
                  'succinate']






# used to compute the moving average
intervals = [4,
             12,
             4,
             4,
             4,
             16,
             4,
             4,
             4]


#file to store the results
outputFile = os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bt_powell_optimization.txt')



#####################################################################


####get splines
splines = {}

for i,v in enumerate(measuredStates):

    stFile = parseTable(os.path.join(strainSummaryFolder, v +  '.tsv'))
    df_state = getDFdict(stFile, v, False)
    summ_state = summarizeExperiments(df_state, v, experimentLabel, interval = intervals[i])
    splines[v] = get_spline(v, 'nothing', experimentLabel, df_state = summ_state)


evals = []



out = minimize(distance,
               params=inputParams,
               method='powell',
               kws = {'database' : database,
                      'measuredStates' : measuredStates,
                      'splines': splines,
                      'experimentLabel':experimentLabel,
                      'species' : species,
                      'intervals': intervals,
                      'combined': True,
                      }
               )



