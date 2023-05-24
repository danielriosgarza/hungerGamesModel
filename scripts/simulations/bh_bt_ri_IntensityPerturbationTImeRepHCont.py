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
    wc_feed = createMetabolome(db, 'wc', pH, pHFunc=None)
    wc_reactor = createMetabolome(db, 'wc', pH, pHFunc=None)
    
    
    
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

    
   
     
    batchA = Pulse(wc_feed, feed_microbiome, 0, 2000, 100, 0, 0, d1,d1)
    
    batchB = Pulse(wc_feed, feed_microbiome, 2000, 2000 + pt, 100, 0, 0, d2,d2)
    
    batchC = Pulse(wc_feed, feed_microbiome, 2000 + pt, 3000, 100, 0, 0, d3,d3)
    
    
    
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
    
    
    d = {'time' : reactor.time_simul,
         'xa' : reactor.subpop_simul[reactor.microbiome.subpops.index('xa')],
         'xb' : reactor.subpop_simul[reactor.microbiome.subpops.index('xb')],
         'xe' : reactor.subpop_simul[reactor.microbiome.subpops.index('xe')],
         'xf' : reactor.subpop_simul[reactor.microbiome.subpops.index('xf')],
         'xi' : reactor.subpop_simul[reactor.microbiome.subpops.index('xi')],
         'xj' : reactor.subpop_simul[reactor.microbiome.subpops.index('xj')],
         'pH' : reactor.pH_simul,
         'trehalose' : reactor.met_simul[reactor.metabolome.metabolites.index('trehalose')],
         'glucose' : reactor.met_simul[reactor.metabolome.metabolites.index('glucose')],
         'pyruvate' : reactor.met_simul[reactor.metabolome.metabolites.index('pyruvate')],
         'butyrate' : reactor.met_simul[reactor.metabolome.metabolites.index('butyrate')],
         'acetate' : reactor.met_simul[reactor.metabolome.metabolites.index('acetate')],
         'lactate' : reactor.met_simul[reactor.metabolome.metabolites.index('lactate')],
         'succinate' : reactor.met_simul[reactor.metabolome.metabolites.index('succinate')],
         'formate' : reactor.met_simul[reactor.metabolome.metabolites.index('formate')]
         }
    
    #log2(xb/xa)
    #dist = np.log2(reactor.subpop_simul[1][20000]/reactor.subpop_simul[0][20000])
    #print('\n', dist, '\n')
    return d



times = np.arange(3, 48, 1)

departure = np.linspace(0,0.5, 300)

results = {i:{} for i in times}


for time in times:
    
    counter = 0
    for dilution in departure:
        print(time, '\t', counter)
        
        results[time][counter] = getDistance(d1 = 0.89 + dilution, d2 = 0, d3 = 0.89 + dilution, pt = time)
        counter+=1
        



states = list(results[3][0].keys())[1::]

cosdict = {i:{} for i in times}

for time in times:
    
    for i in range(len(departure)):
        cosdict[time][i] = cosine(np.array([results[time][i][st][99] for st in states]), np.array([results[time][i][st][299] for st in states]) )
    
befdict = {i:{} for i in times}
for time in times:
    for i in range(len(departure)):
        befdict[time][i] = np.array([results[time][i][st][99] for st in states])

aftdict = {i:{} for i in times}
for time in times:
    for i in range(len(departure)):
        aftdict[time][i] = np.array([results[time][i][st][299] for st in states])



from pylab import *


cosMat = []

for time in times:
    v = np.array([cosdict[time][i] for i in range(len(departure))])
    cosMat.append(v)

cosMat = np.array(cosMat)

pcolormesh(cosMat, cmap=cm.coolwarm)
tp = [0,20,40,60,80,99]

xticks(tp, np.round((0.89+departure[tp])/15,3))
colorbar(label = 'cosine distance to return community')
ylabel('feed perturbation duration (h)')
xlabel('dilution rate ($h^{-1}$)')
#savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'perturbDuration.png'), dpi = 600)
show()



treMatBef = []

for time in times:
    v = np.array([befdict[time][i][states.index('trehalose')] for i in range(len(departure))])
    treMatBef.append(v)

treMatBef = np.array(treMatBef)

pcolormesh(treMatBef.T, cmap=cm.coolwarm)
colorbar(label = 'trehalose concentration')
xlabel('feed perturbation duration (h)')
ylabel('distance from threshold ($h^{-1}$)')
show()

treMatAft = []

for time in times:
    v = np.array([aftdict[time][i][states.index('trehalose')] for i in range(len(departure))])
    treMatAft.append(v)

treMatAft = np.array(treMatAft)

pcolormesh(treMatAft.T, cmap=cm.coolwarm)
colorbar(label = 'trehalose concentration')
xlabel('feed perturbation duration (h)')
ylabel('distance from threshold ($h^{-1}$)')
show()

xbMatBef = []

for time in times:
    v = np.array([befdict[time][i][states.index('xb')] for i in range(len(departure))])
    xbMatBef.append(v)

xbMatBef = np.array(xbMatBef)

pcolormesh(xbMatBef.T, cmap=cm.coolwarm)
colorbar(label = 'xb concentration')
xlabel('feed perturbation duration (h)')
ylabel('distance from threshold ($h^{-1}$)')
show()

xbMatAft = []

for time in times:
    v = np.array([aftdict[time][i][states.index('xb')] for i in range(len(departure))])
    xbMatAft.append(v)

xbMatAft = np.array(xbMatAft)

pcolormesh(xbMatAft.T, cmap=cm.coolwarm)
colorbar(label = 'xb concentration')
xlabel('feed perturbation duration (h)')
ylabel('distance from threshold ($h^{-1}$)')
show()

def heatMaps(state, vmin = None, vmax =None):
    bef = []

    for time in times:
        v = np.array([befdict[time][i][states.index(state)] for i in range(len(departure))])
        bef.append(v)

    bef = np.array(bef)
    
    if vmin is not None:
        pcolormesh(bef.T, cmap=cm.coolwarm, vmin = vmin, vmax = vmax)
    
    else:
        pcolormesh(bef.T, cmap=cm.coolwarm)
    colorbar(label = state + ' concentration')
    xlabel('feed perturbation duration (h)')
    ylabel('distance from threshold ($h^{-1}$)')
    show()

    aft = []

    for time in times:
        v = np.array([aftdict[time][i][states.index(state)] for i in range(len(departure))])
        aft.append(v)

    aft = np.array(aft)
    
    if vmin is not None:
        pcolormesh(aft.T, cmap=cm.coolwarm, vmin = vmin, vmax = vmax)
    else:
        pcolormesh(aft.T, cmap=cm.coolwarm)
    colorbar(label = state + ' concentration')
    xlabel('feed perturbation duration (h)')
    ylabel('distance from threshold ($h^{-1}$)')
    show()



#######twoStates############
twoStates = []
for i in range(len(departure)):
    q = np.array([results[10][i][z][99] for z in states])
    twoStates.append(q)
twoStates = np.array(twoStates).T


bh = twoStates[0] + twoStates[1]
bt = twoStates[2] + twoStates[3]
ri = twoStates[4] + twoStates[5]

mStates = ['bh', 'bt', 'ri'] + states[6::]

l = [bh,bt,ri]

for i in twoStates[6::]:
    l.append(i)

twoStates = np.array(l)
for i,v in enumerate(twoStates):
    twoStates[i] = twoStates[i]/max(twoStates[i])
    
pcolormesh(twoStates, cmap=cm.coolwarm, vmin = 0, vmax=1, antialiased=True)
yticks(np.arange(len(mStates))+0.5, mStates)
tp = [0,20,40,60,80,99]

xticks(tp, np.round((0.4+departure[tp])/15,3))
colorbar()
xlabel('dilution rate ($h^{-1}$)')
tight_layout()
#savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'twoStates.png'), dpi = 600)
show()
##########################################

#####hysteresis###################

stA = []
for i in range(len(departure)):
    q = np.array([results[30][i][z][100] for z in states])
    stA.append(q)
stA = np.array(stA).T

stB = []
for i in range(len(departure)):
    q = np.array([results[30][i][z][299] for z in states])
    stB.append(q)
stB = np.array(stB).T

bhA = stA[0] + stA[1]
btA = stA[2] + stA[3]
riA = stA[4] + stA[5]

bhB = stB[0] + stB[1]
btB = stB[2] + stB[3]
riB = stB[4] + stB[5]

lA = [bhA,btA,riA]
for i in stA[6::]:
    lA.append(i)
stA = np.array(lA)
for i,v in enumerate(stA):
    stA[i] = stA[i]/max(stA[i])


lB = [bhB,btB,riB]
for i in stB[6::]:
    lB.append(i)
stB = np.array(lB)
for i,v in enumerate(stB):
    stB[i] = stB[i]/max(stB[i])

hysteric = np.concatenate([stA, stB])
pcolormesh(hysteric, cmap=cm.coolwarm, vmin = 0, vmax=1, antialiased=True)
yticks(np.arange(len(mStates + mStates))+0.5, mStates+mStates)
tp = [0,20,40,60,80,99]

xticks(tp, np.round((0.4+departure[tp])/15,3))
colorbar()
xlabel('dilution rate ($h^{-1}$)')
tight_layout()
#savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'hysteresis.png'), dpi = 600)
show()

########DrJekyllMrHyde
twoStates = []
for i in range(len(departure)):
    q = np.array([results[10][i][z][99] for z in states])
    twoStates.append(q)
twoStates = np.array(twoStates).T


bh = twoStates[0] + twoStates[1]
bt = twoStates[2] + twoStates[3]
ri = twoStates[4] + twoStates[5]

mStates = ['bh', 'bt', 'ri'] + states[6::]

l = [bh,bt,ri]

for i in twoStates[6::]:
    l.append(i)

twoStates = np.array(l)

twoStates = twoStates[[0,1,2]]

mst = ['bh', 'bt', 'ri']
    
pcolormesh(twoStates, cmap=cm.coolwarm, antialiased=True)
yticks(np.arange(len(mst))+0.5, mst)
tp = [0,20,40,60,80,99]

xticks(tp, np.round((0.4+departure[tp])/15,3))
colorbar()
xlabel('dilution rate ($h^{-1}$)')
tight_layout()
savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'DrJekyllMrHyde.png'), dpi = 600)
show()





# from pylab import *
# from aquarel import load_theme

# theme = load_theme("boxy_light")
# theme.apply()


# for t in times:
#     scatter([t]*100, departure, c=results[t], cmap=cm.coolwarm, vmin=-30, vmax=-6,s=20, alpha = .5)
# colorbar(label='log2(xb/xa)')
# xlabel('perturbation duration (h)')
# ylabel('distance from threshold ($h^{-1}$)')

# #xlim(10,48)
# #ylim(-0.01,.32)
# savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'perturbationTimeSubpops.png'), dpi = 600)