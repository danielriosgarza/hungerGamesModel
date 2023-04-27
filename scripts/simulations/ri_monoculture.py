# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 16:54:49 2023

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




import plotly.io as pio
pio.renderers.default='browser'




def getStates(reactor, state):
    
    states = {}
    
    
    for st in state:
        if st == 'xi':
            states['xi'] = reactor.subpop_simul[0][-1]
        
        if st == 'xj':
            states['xj'] = reactor.subpop_simul[1][-1]
        
        if st == 'xk':
            states['xk'] = reactor.subpop_simul[2][-1]
        
        if st == 'xl':
            states['xl'] = reactor.subpop_simul[3][-1]
        
        
        if st == 'pH':
            states['pH'] = reactor.pH_simul[-1]
            
        
        if st == 'glucose':
            states['glucose'] = reactor.met_simul[reactor.metabolome.metabolites.index('glucose')][-1]
            
        if st == 'pyruvate':
            states['pyruvate'] = reactor.met_simul[reactor.metabolome.metabolites.index('pyruvate')][-1]
            
        if st == 'lactate':
            states['lactate'] = reactor.met_simul[reactor.metabolome.metabolites.index('lactate')][-1]
        
        
        if st == 'acetate':
            states['acetate'] = reactor.met_simul[reactor.metabolome.metabolites.index('acetate')][-1]
            
        
        if st == 'butyrate':
            states['butyrate'] = reactor.met_simul[reactor.metabolome.metabolites.index('butyrate')][-1]
    return states


#pH profile
ipH_path = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bhbtri_ipH4.tsv') 

databaseName = 'modelDB_bhbtri.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

#update database with parameters from a file
##########################################

#create a database connection



states = ['xi', 'xj', 'xk', 'xl', 'pH', 'glucose', 'pyruvate', 'lactate', 'acetate', 'butyrate']

sts = {i:[] for i in states}

glucose  = np.linspace(0, 20, 300)


for i, v in enumerate(glucose):
    
    print(i)
    
    conn = create_connection(os.path.join(databaseFolder, databaseName))
    
    #load the parameter file (parameter files are located at "/files/params" )
    ri_params = getPramsFromFile('ri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'ri.tsv'))
    
    
    #assign these parameters (depending on the strain, use the specific function)
    assignRiParams(ri_params, conn)
    
    
    #Load database
    db = get_database(os.path.join(databaseFolder, databaseName))
    
    
    #getStarting pH
    wc = createMetabolome(db, 'wc')
    predictpH = getpH(wc.metabolites, ipH_path)
    pH =  predictpH(wc.get_concentration())
    
    #get the feed media and the reactor media
    wc_feed = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    wc_feed.metD['glucose'].update(v)
    wc_reactor = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    # wc_reactor.metD['lactate'].update(10.0)
    # wc_reactor.metD['glucose'].update(0.0)
    
    
    #get the feed obj. Make it sterile
    feed_microbiome = Microbiome({'ri':createBacteria(db, 'ri', 'wc')})
    feed_microbiome.subpopD['xi'].count = 0
    
    #create the reactor obj, with starting populations
    reactor_microbiome = Microbiome({'ri':createBacteria(db, 'ri', 'wc')})
    #reactor_microbiome.subpopD['xa'].count = 0.01
    #reactor_microbiome.subpopD['xb'].count = 0.01
    batchA = Pulse(wc_feed, feed_microbiome, 0, 1000, 100000, 0, 0, 5.5,5.5)
    
    #simulate
    reactor = Reactor(reactor_microbiome, wc_reactor,[batchA], 60)
    reactor.simulate()
    #reactor.makePlots()
    st = getStates(reactor, states)
    
    for si in st:
        sts[si].append(st[si])





# conn = create_connection(os.path.join(databaseFolder, databaseName))

# #load the parameter file (parameter files are located at "/files/params" )
# ri_params = getPramsFromFile('ri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'ri.tsv'))


# #assign these parameters (depending on the strain, use the specific function)
# assignRiParams(ri_params, conn)


# #Load database
# db = get_database(os.path.join(databaseFolder, databaseName))


# #getStarting pH
# wc = createMetabolome(db, 'wc')
# predictpH = getpH(wc.metabolites, ipH_path)
# pH =  predictpH(wc.get_concentration())

# #get the feed media and the reactor media
# wc_feed = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
# wc_reactor = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
# # wc_reactor.metD['lactate'].update(10.0)
# # wc_reactor.metD['glucose'].update(0.0)


# #get the feed obj. Make it sterile
# feed_microbiome = Microbiome({'ri':createBacteria(db, 'ri', 'wc')})
# feed_microbiome.subpopD['xi'].count = 0

# #create the reactor obj, with starting populations
# reactor_microbiome = Microbiome({'ri':createBacteria(db, 'ri', 'wc')})
# #reactor_microbiome.subpopD['xa'].count = 0.01
# #reactor_microbiome.subpopD['xb'].count = 0.01
# batchA = Pulse(wc_feed, feed_microbiome, 0, 500, 100000, 0, 0, 1,1)

# batchB = Pulse(wc_feed, feed_microbiome, 500, 1000, 100000, 0, 0, 10,10)

# batchC = Pulse(wc_feed, feed_microbiome, 1500, 2000, 100000, 0, 0, 20,20)

# batchD = Pulse(wc_feed, feed_microbiome, 2500, 3000, 100000, 0, 0, 10,10)

# batchE = Pulse(wc_feed, feed_microbiome, 3500, 4000, 100000, 0, 0, 1,1)

# batchF = Pulse(wc_feed, feed_microbiome, 4500, 5000, 100000, 0, 0, 10,10)

# batchG = Pulse(wc_feed, feed_microbiome, 5500, 6000, 100000, 0, 0, 20,20)

# batchH = Pulse(wc_feed, feed_microbiome, 6500, 7000, 100000, 0, 0, 1,1)

# #simulate
# reactor = Reactor(reactor_microbiome, wc_reactor,[batchA, batchB, batchC, batchD, batchE, batchF, batchG, batchH], 60)
# reactor.simulate()





# #########fermentation acids ###################
# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.met_simul[reactor.metabolome.metabolites.index('acetate')],
#                 color = '#003eff',
#                 legend = 'acetate (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = 'concentration (mM)',
#                 title = 'Ri acids',
#                 linestyle = '-')


# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.met_simul[reactor.metabolome.metabolites.index('lactate')],
#                 color = '#00ffaf',
#                 legend = 'lactate (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = 'concentration (mM)',
#                 title = None,
#                 linestyle = '-')


# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.met_simul[reactor.metabolome.metabolites.index('butyrate')],
#                 color = '#ff00a1',
#                 legend = 'butyrate (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = 'concentration (mM)',
#                 title = None,
#                 linestyle = '-')



# plt.show()


# ########Carbon Consumption
# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.met_simul[reactor.metabolome.metabolites.index('pyruvate')],
#                 color = '#ff8900',
#                 legend = 'pyr (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = 'concentration (mM)',
#                 title = 'carbon consumption',
#                 linestyle = '-')


# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.met_simul[reactor.metabolome.metabolites.index('lactate')],
#                 color = '#00ffaf',
#                 legend = 'lactate (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = 'concentration (mM)',
#                 title = None,
#                 linestyle = '-')


# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.met_simul[reactor.metabolome.metabolites.index('glucose')],
#                 color = '#ff0000',
#                 legend = 'glucose (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = 'concentration (mM)',
#                 title = None,
#                 linestyle = '-')

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.met_simul[reactor.metabolome.metabolites.index('acetate')],
#                 color = '#003eff',
#                 legend = 'acetate (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = 'concentration (mM)',
#                 title = None,
#                 linestyle = '-')


# plt.show()
# #############Subpopulations


# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.subpop_simul[0],
#                 color = '#3998FF',
#                 legend = 'live_xi (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '--')

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.subpop_simul[1],
#                 color = '#45A4FFA6',
#                 legend = 'live_xj (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '--')


# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.cellActive_dyn[0],
#                 color = '#1F51FF',
#                 legend = 'live_cells (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '-')



# plt.show()

# ################pH

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.pH_simul,
#                 color = '#39FF14',
#                 legend = 'pH (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = 'pH',
#                 title = None,
#                 linestyle = '-')


# plt.ylim(5,7)