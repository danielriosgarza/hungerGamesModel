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




import plotly.io as pio
pio.renderers.default='browser'



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



bh, bt, ri = [],[],[]


for i in np.linspace(0, 25, 25):
    
    print(i)
        
    
    #getStarting pH
    wc = createMetabolome(db, 'wc')
    predictpH = getpH(wc.metabolites, ipH_path)
    pH =  predictpH(wc.get_concentration())
    
    #get the feed media and the reactor media
    wc_feed = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    wc_reactor = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    #wc_reactor.metD['mannose'].update(0)
    # wc_reactor.metD['glucose'].update(0.0)
    
    
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
    #reactor_microbiome.subpopD['xa'].count = 0.00
    #reactor_microbiome.subpopD['xb'].count = 0.01
    
    batchA = Pulse(wc_feed, feed_microbiome, 0, 12, 10000, 0, 0, 0,0)
    
    chemostat = Pulse(wc_feed, feed_microbiome, 12, 1500, 10000, 0, 0, i,i)
    
    #simulate
    reactor = Reactor(reactor_microbiome, wc_reactor,[batchA, chemostat], 100)
    reactor.simulate()
    reactor.makePlots()
    bh.append(reactor.cellActive_dyn[0][-1])
    bt.append(reactor.cellActive_dyn[1][-1])
    ri.append(reactor.cellActive_dyn[2][-1])
    
    print(bh[-1], bt[-1], ri[-1])
    


# #########fermentation acids ###################
# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.met_simul[reactor.metabolome.metabolites.index('acetate')],
#                 color = '#003eff',
#                 legend = 'acetate (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = 'concentration (mM)',
#                 title = 'acids',
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

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.met_simul[reactor.metabolome.metabolites.index('formate')],
#                 color = '#00c6ff',
#                 legend = 'formate (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = 'concentration (mM)',
#                 title = None,
#                 linestyle = '-')

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.met_simul[reactor.metabolome.metabolites.index('succinate')],
#                 color = '#00ff26',
#                 legend = 'succinate (simul)',
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


# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.met_simul[reactor.metabolome.metabolites.index('trehalose')],
#                 color = '#8900ff',
#                 legend = 'trehalose (simul)',
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


# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.met_simul[reactor.metabolome.metabolites.index('trehalose')],
#                 color = '#8900ff',
#                 legend = 'trehalose (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = 'concentration (mM)',
#                 title = None,
#                 linestyle = '-')


# plt.show()

# ####################Subpopulations

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.subpop_simul[0],
#                 color = '#FF10F0',
#                 legend = 'live_xa (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '--')

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.subpop_simul[1],
#                 color = '#FF2E2EC9',
#                 legend = 'live_xb (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '--')


# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.cellActive_dyn[0],
#                 color = '#FF6EC7',
#                 legend = 'live_bh_cells (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '-')

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.subpop_simul[4],
#                 color = '#FF7727',
#                 legend = 'live_xe (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '--')

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.subpop_simul[5],
#                 color = '#FFDB9E',
#                 legend = 'live_xf (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '--')


# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.cellActive_dyn[1],
#                 color = '#FF5F1F',
#                 legend = 'live_bt_cells (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '-')

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.subpop_simul[8],
#                 color = '#3998FF',
#                 legend = 'live_xi (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '--')

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.subpop_simul[9],
#                 color = '#45A4FFA6',
#                 legend = 'live_xj (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '--')


# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.cellActive_dyn[2],
#                 color = '#1F51FF',
#                 legend = 'live_ri_cells (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '-')



# plt.show()


