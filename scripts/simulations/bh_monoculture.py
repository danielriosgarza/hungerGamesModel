# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 11:39:38 2023

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

from aquarel import load_theme

theme = load_theme("boxy_light")
theme.apply()

import matplotlib.pyplot as plt




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


#assign these parameters (depending on the strain, use the specific function)
assignBhParams(bh_params, conn)


#Load database
db = get_database(os.path.join(databaseFolder, databaseName))


#getStarting pH
wc = createMetabolome(db, 'wc')
predictpH = getpH(wc.metabolites, ipH_path)
pH =  predictpH(wc.get_concentration())

#get the feed media and the reactor media
wc_feed = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
wc_reactor = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
#wc_reactor.metD['pyruvate'].update(0.0)
#wc_reactor.metD['trehalose'].update(8.0)
#wc_reactor.metD['glutamate'].update(9.0)

#get the feed obj. Make it sterile
feed_microbiome = Microbiome({'bh':createBacteria(db, 'bh', 'wc')})
feed_microbiome.subpopD['xa'].count = 0

#create the reactor obj, with starting populations
reactor_microbiome = Microbiome({'bh':createBacteria(db, 'bh', 'wc')})
#reactor_microbiome.subpopD['xa'].count = 0.01
#reactor_microbiome.subpopD['xb'].count = 0.01
batchA = Pulse(wc_feed, feed_microbiome, 0, 120, 100000, 0, 0, 0,0)

#simulate
reactor = Reactor(reactor_microbiome, wc_reactor,[batchA], 60)
reactor.simulate()
#reactor.makePlots()



#########makes a lot of acetate ###################
makeKineticPlot(x = reactor.time_simul,
                y = reactor.met_simul[reactor.metabolome.metabolites.index('acetate')],
                color = '#003eff',
                legend = 'acetate (simul)',
                xlabel = 'time (h)',
                ylabel = 'concentration (mM)',
                title = 'Bh acids',
                linestyle = '-')


makeKineticPlot(x = reactor.time_simul,
                y = reactor.met_simul[reactor.metabolome.metabolites.index('lactate')],
                color = '#00ffaf',
                legend = 'lactate (simul)',
                xlabel = 'time (h)',
                ylabel = 'concentration (mM)',
                title = None,
                linestyle = '-')



plt.show()

########Carbon Consumption
makeKineticPlot(x = reactor.time_simul,
                y = reactor.met_simul[reactor.metabolome.metabolites.index('pyruvate')],
                color = '#ff8900',
                legend = 'pyr (simul)',
                xlabel = 'time (h)',
                ylabel = 'concentration (mM)',
                title = 'carbon consumption',
                linestyle = '-')


makeKineticPlot(x = reactor.time_simul,
                y = reactor.met_simul[reactor.metabolome.metabolites.index('trehalose')],
                color = '#8900ff',
                legend = 'treh (simul)',
                xlabel = 'time (h)',
                ylabel = 'concentration (mM)',
                title = None,
                linestyle = '-')


makeKineticPlot(x = reactor.time_simul,
                y = reactor.met_simul[reactor.metabolome.metabolites.index('glucose')],
                color = '#ff0000',
                legend = 'glucose (simul)',
                xlabel = 'time (h)',
                ylabel = 'concentration (mM)',
                title = None,
                linestyle = '-')

makeKineticPlot(x = reactor.time_simul,
                y = reactor.met_simul[reactor.metabolome.metabolites.index('glutamate')],
                color = '#00B8FF',
                legend = 'glutamate (simul)',
                xlabel = 'time (h)',
                ylabel = 'concentration (mM)',
                title = None,
                linestyle = '-')


plt.show()







#############Subpopulations


makeKineticPlot(x = reactor.time_simul,
                y = reactor.subpop_simul[0],
                color = '#FF10F0',
                legend = 'live_xa (simul)',
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '--')

makeKineticPlot(x = reactor.time_simul,
                y = reactor.subpop_simul[1],
                color = '#FF2E2EC9',
                legend = 'live_xb (simul)',
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '--')


makeKineticPlot(x = reactor.time_simul,
                y = reactor.cellActive_dyn[0],
                color = '#1F51FF',
                legend = 'live_cells (simul)',
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-')



plt.show()

################pH

makeKineticPlot(x = reactor.time_simul,
                y = reactor.pH_simul,
                color = '#39FF14',
                legend = 'pH (simul)',
                xlabel = 'time (h)',
                ylabel = 'pH',
                title = None,
                linestyle = '-')

plt.ylim(5,7)
plt.show()