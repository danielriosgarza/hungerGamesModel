# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 06:04:24 2023

@author: danie
"""

from pathlib import Path
import os
import sys


sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'core'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'db'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))

from readModelDB import *
from mainClasses import *
from loadParameters import *
from general import *




import plotly.io as pio
pio.renderers.default='browser'



def mockpHfunc(metObj, pH=7.0):
    def pHfunc(metObj):
        return pH
    return pHfunc


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
bt_params = getPramsFromFile('bhbtri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'btri.tsv'))
ri_params = getPramsFromFile('bhbtri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'btri.tsv'))


#assign these parameters (depending on the strain, use the specific function)
assignBhParams(bh_params, conn)
assignBtParams(bt_params, conn)
assignRiParams(ri_params, conn)

#Load database
db = get_database(os.path.join(databaseFolder, databaseName))

#getStarting pH
wc = createMetabolome(db, 'wc')
predictpH = getpH(wc.metabolites, ipH_path)
fixedpHA = mockpHfunc(wc.metabolites,pH=5.60)
fixedpHC = mockpHfunc(wc.metabolites,pH=5.47)
pH =  predictpH(wc.get_concentration())

#get the feed media and the reactor media
wc_feedA = createMetabolome(db, 'wc', pH, pHFunc=fixedpHA)
wc_feedC = createMetabolome(db, 'wc', pH, pHFunc=fixedpHC)

wc_reactorA = createMetabolome(db, 'wc', pH, pHFunc=fixedpHA)
wc_reactorC = createMetabolome(db, 'wc', pH, pHFunc=fixedpHC)


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
reactor_microbiome.subpopD['xe'].count = 0.003
reactor_microbiome.subpopD['xi'].count = 0.003


d = 1.0




batchA = Pulse(wc_feedA, feed_microbiome, 0, 600, 100, 0, 0, d,d)

batchB = Pulse(wc_feedC, feed_microbiome, 600, 3000, 100, 0, 0, d,d)


#simulate
reactorA = Reactor(reactor_microbiome, wc_reactorA,[batchA], 15)



reactorA.simulate()
#reactorA.makePlots()



for i in wc_reactorA.metD:
    wc_reactorC.metD[i].update(wc_reactorA.metD[i].concentration)
    

reactorB = Reactor(reactorA.microbiome, wc_reactorC,[batchB], 15)

reactorB.simulate()





makeKineticPlot(x = reactorA.time_simul*0.1,
                y = reactorA.cellActive_dyn[0],
                color = '#FF10F0',
                legend = 'Blautia hydrogenotrophica',
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-',
                legendSize = 10)

makeKineticPlot(x = reactorB.time_simul*0.1,
                y = reactorB.cellActive_dyn[0],
                color = '#FF10F0',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-',
                legendSize = 10)




makeKineticPlot(x = reactorA.time_simul*0.1,
                y = reactorA.cellActive_dyn[1],
                color = '#ff8300',
                legend = 'Bacteroides thetaiotaomicron',
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-',
                legendSize = 10)

makeKineticPlot(x = reactorB.time_simul*0.1,
                y = reactorB.cellActive_dyn[1],
                color = '#ff8300',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-',
                legendSize = 10)



makeKineticPlot(x = reactorA.time_simul*0.1,
                y = reactorA.cellActive_dyn[2],
                color = '#00B8FF',
                legend = 'Roseburia instestinalis',
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-',
                legendSize = 10)

makeKineticPlot(x = reactorB.time_simul*0.1,
                y = reactorB.cellActive_dyn[2],
                color = '#00B8FF',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-',
                legendSize = 10)


title = 'stateA'
#plt.savefig(os.path.join(Path(os.getcwd()).parents[2], 'files', 'Figures', 'multistability', title + '.png'), transparent=True, dpi=600)
plt.show()

