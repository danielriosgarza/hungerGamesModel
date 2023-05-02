# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 06:04:24 2023

@author: danie
"""

from pathlib import Path
import os
import sys

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'db'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))

from readModelDB import *
from mainClasses import *
from loadParameters import *
from general import *
import matplotlib.pyplot as plt
plt.style.use('seaborn-bright')
import seaborn as sns

from scipy.spatial.distance import cosine
import numpy as np


import plotly.io as pio
pio.renderers.default='browser'





def getDistance(d1 = 0.615, d2 = 0, d3 = 0.615, pt = 24, plot = False):
    
    #pH profile
    ipH_path = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bhbtri_ipH4.tsv') 
    
    databaseName = 'modelDB_bhbtri.sqlite3'
    
    databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')
    
    #update database with parameters from a file
    ##########################################
    
    #create a database connection
    conn = create_connection(os.path.join(databaseFolder, databaseName))
    
    #load the parameter file (parameter files are located at "/files/params" )
    bh_params = getPramsFromFile('bh', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bh.tsv'))
    bt_params = getPramsFromFile('bt', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bt.tsv'))
    ri_params = getPramsFromFile('ri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'ri.tsv'))
    
    
    #assign these parameters (depending on the strain, use the specific function)
    assignBhParams(bh_params, conn)
    assignBtParams(bt_params, conn)
    assignRiParams(ri_params, conn)
    
    #Load database
    db = get_database(os.path.join(databaseFolder, databaseName))
    
    
    #getStarting pH
    wc = createMetabolome(db, 'wc')
    predictpH = getpH(wc.metabolites, ipH_path)
    pH =  predictpH(wc.get_concentration())
    
    #get the feed media and the reactor media
    wc_feed = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    wc_reactor = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    
    
    
    #get the feed obj. Make it sterile
    feed_microbiome = Microbiome({'bh':createBacteria(db, 'bh', 'wc'),
                                  'bt':createBacteria(db, 'bt', 'wc'),
                                  'ri':createBacteria(db, 'ri', 'wc')})
    feed_microbiome.subpopD['xa'].count = 0
    feed_microbiome.subpopD['xe'].count = 0
    feed_microbiome.subpopD['xi'].count = 0
    
    #create the reactor obj, with starting populations
    reactor_microbiome = Microbiome({'bh':createBacteria(db, 'bh', 'wc'),
                                     'bt':createBacteria(db, 'bt', 'wc'),
                                     'ri':createBacteria(db, 'ri', 'wc')})
    reactor_microbiome.subpopD['xa'].count = 0.01
    reactor_microbiome.subpopD['xe'].count = 0.01
    reactor_microbiome.subpopD['xi'].count = 0.01

    
   
     
    batchA = Pulse(wc_feed, feed_microbiome, 0, 2000, 10000, 0, 0, d1,d1)
    
    batchB = Pulse(wc_feed, feed_microbiome, 2000, 2000 + pt, 10000, 0, 0, d2,d2)
    
    batchC = Pulse(wc_feed, feed_microbiome, 2000 + pt, 4000, 10000, 0, 0, d3,d3)
    
    
    #simulate
    reactor = Reactor(reactor_microbiome, wc_reactor,[
                                                      batchA,
                                                      batchB,
                                                      batchC
                                                      
                                                       ], 15)
    reactor.simulate()
    
    if plot:
        reactor.makePlots()
    
    
    #dist = cosine(reactor.subpop_simul.T[10000], reactor.subpop_simul.T[-1])
    
    dist = reactor.met_simul[reactor.metabolome.metabolites.index('trehalose')][20000]
    
    print('\n', dist, '\n')
    return dist



times = np.arange(3, 48, 1)

departure = np.linspace(0,0.5, 100)

results = {i:[] for i in times}


for time in times:
    for dilution in departure:
        print(time, '\t', dilution)
        
        results[time].append(getDistance(d1 = 0.61+dilution, d2 = 0, d3 = 0.61+dilution, pt = time))
        

from pylab import *


for t in times:
    scatter([t]*100, departure, c=results[t], cmap=cm.coolwarm, vmin=0, vmax=0.6,s=20)
colorbar(label='trehalose concentration after perturbation')
xlabel('perturbation duration (h)')
ylabel('distance from threshold ($h^{-1}$)')

xlim(10,48)
ylim(-0.01,.32)
savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'perturbationTIme.png'), dpi = 600)