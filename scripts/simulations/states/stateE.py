# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 06:04:24 2023

@author: danie
"""

from pathlib import Path
import os
import sys


sys.path.append(os.path.join(Path(os.getcwd()).parents[1], 'core'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[1], 'db'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[1], 'compare2experiments'))

from readModelDB import *
from mainClasses import *
from loadParameters import *
from general import *




import plotly.io as pio
pio.renderers.default='browser'



#pH profile
ipH_path = os.path.join(Path(os.getcwd()).parents[2], 'files', 'strainSummaries', 'bhbtri_ipH4.tsv') 

databaseName = 'modelDB_bhbtri.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[2], 'files', 'dbs')

#update database with parameters from a file
##########################################

#create a database connection
conn = create_connection(os.path.join(databaseFolder, databaseName))

#load the parameter file (parameter files are located at "/files/params" )
bh_params = getPramsFromFile('bh', os.path.join(Path(os.getcwd()).parents[2], 'files', 'params', 'bh.tsv'))
bt_params = getPramsFromFile('bhbtri', os.path.join(Path(os.getcwd()).parents[2], 'files', 'params', 'btri.tsv'))
ri_params = getPramsFromFile('bhbtri', os.path.join(Path(os.getcwd()).parents[2], 'files', 'params', 'btri.tsv'))


#assign these parameters (depending on the strain, use the specific function)
assignBhParams(bh_params, conn)
assignBtParams(bt_params, conn)
assignRiParams(ri_params, conn)

#Load database
db = get_database(os.path.join(databaseFolder, databaseName))

def mockpHfunc(metObj, pH=7.0):
    def pHfunc(metObj):
        return pH
    return pHfunc





#getStarting pH
wc = createMetabolome(db, 'wc')
predictpH = getpH(wc.metabolites, ipH_path)
fixedpHA = mockpHfunc(wc.metabolites,pH=5.53)
fixedpHB = mockpHfunc(wc.metabolites,pH=3.5)
fixedpHC = mockpHfunc(wc.metabolites,pH=3.5)
pH =  predictpH(wc.get_concentration())

#get the feed media and the reactor media
wc_feed = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
wc_feedA = createMetabolome(db, 'wc', pH, pHFunc=fixedpHA)
wc_feedB = createMetabolome(db, 'wc', pH, pHFunc=fixedpHB)

wc_reactor = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
wc_reactorA = createMetabolome(db, 'wc', pH, pHFunc=fixedpHA)
wc_reactorB = createMetabolome(db, 'wc', pH, pHFunc=fixedpHB)


#wc_feed_pH = createMetabolome(db, 'wc', 2.5, pHFunc=predictpH)
#wc_reactor.metD['trehalose'].update(0.1)
#wc_feed.metD['trehalose'].update(0.000)
#wc_reactor.metD['pyruvate'].update(0.1)
#wc_feed.metD['pyruvate'].update(0.000)
#wc_reactor.metD['mannose'].update(0.1)
#wc_feed.metD['mannose'].update(0.000)


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
reactor_microbiome.subpopD['xa'].count = 0.003
reactor_microbiome.subpopD['xb'].count = 0.00
reactor_microbiome.subpopD['xe'].count = 0.003
reactor_microbiome.subpopD['xf'].count = 0.0
reactor_microbiome.subpopD['xi'].count = 0.003
reactor_microbiome.subpopD['xj'].count = 0.


d = 1.0
d2 = 1.0
d3 = 0.85




batchA = Pulse(wc_feed, feed_microbiome, 0, 600, 100, 0, 0, d,d)

batchB = Pulse(wc_feed, feed_microbiome, 600, 700, 100, 0, 0, d2,d2)

batchC = Pulse(wc_feed, feed_microbiome, 700, 3000, 100, 0, 0, d,d)


#simulate
reactorA = Reactor(reactor_microbiome, wc_reactorA,[
                                                  batchA,
                                                  
                                                  
                                                  
                                                   ], 15)



reactorA.simulate()
#reactorA.makePlots()

reactorB = Reactor(reactorA.microbiome, wc_reactorA,[
                                                  
                                                  batchB,
                                                  
                                                   ], 15)

reactorB.simulate()
#reactorB.makePlots()

reactorC = Reactor(reactorB.microbiome, wc_reactorA,[
                                                  batchC
                                                  
                                                  
                                                   ], 15)



reactorC.simulate()
#reactorC.makePlots()












# ####################Subpopulations

makeKineticPlot(x = reactorA.time_simul*0.1,
                y = reactorA.cellActive_dyn[0],
                color = '#FF10F0',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-')


makeKineticPlot(x = reactorB.time_simul*0.1,
                y = reactorB.cellActive_dyn[0],
                color = '#FF10F0',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-')

makeKineticPlot(x = reactorC.time_simul*0.1,
                y = reactorC.cellActive_dyn[0],
                color = '#FF10F0',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-')




makeKineticPlot(x = reactorA.time_simul*0.1,
                y = reactorA.cellActive_dyn[1],
                color = '#ff8300',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-')


makeKineticPlot(x = reactorB.time_simul*0.1,
                y = reactorB.cellActive_dyn[1],
                color = '#ff8300',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-')


makeKineticPlot(x = reactorC.time_simul*0.1,
                y = reactorC.cellActive_dyn[1],
                color = '#ff8300',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-')


makeKineticPlot(x = reactorA.time_simul*0.1,
                y = reactorA.cellActive_dyn[2],
                color = '#00B8FF',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-')

makeKineticPlot(x = reactorB.time_simul*0.1,
                y = reactorB.cellActive_dyn[2],
                color = '#00B8FF',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-')


makeKineticPlot(x = reactorC.time_simul*0.1,
                y = reactorC.cellActive_dyn[2],
                color = '#00B8FF',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-')


title = 'stateE'
plt.savefig(os.path.join(Path(os.getcwd()).parents[2], 'files', 'Figures', 'multistability', title + '.png'), transparent=True, dpi=600)


plt.show()


