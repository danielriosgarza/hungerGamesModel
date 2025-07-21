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



#pH profile
ipH_path = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bhbtri_ipH4.tsv') 

databaseName = 'modelDB_bhbtri.sqlite3'

databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

#update database with parameters from a file
##########################################

#create a database connection
conn = create_connection(os.path.join(databaseFolder, databaseName))

#load the parameter file (parameter files are located at "/files/params" )
bh_params = getPramsFromFile('bh', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'allParamsFitted.tsv'))
bt_params = getPramsFromFile('bt', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'allParamsFitted.tsv'))
ri_params = getPramsFromFile('ri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'allParamsFitted.tsv'))


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
reactor_microbiome.subpopD['xa'].count = 0.003
reactor_microbiome.subpopD['xe'].count = 0.003
reactor_microbiome.subpopD['xi'].count = 0.003


d = 0.615
perturb = 0.0




batchA = Pulse(wc_feed, feed_microbiome, 0, 2400, 100, 0, 0, d,d)

batchB = Pulse(wc_feed, feed_microbiome, 2400, 2640, 100, 0, 0, perturb,perturb)

batchC = Pulse(wc_feed, feed_microbiome, 2640, 5040, 100, 0, 0, d,d)


#simulate
reactorA = Reactor(reactor_microbiome, wc_reactor,[
                                                  batchA,
                                                  batchB,
                                                  batchC
                                                  
                                                  
                                                  
                                                   ], 15)



reactorA.simulate()
#reactorA.makePlots()





# ####################Subpopulations

makeKineticPlot(x = reactorA.time_simul*0.1,
                y = reactorA.cellActive_dyn[0],
                color = '#FF10F0',
                legend = 'Blautia hydrogenotrophica',
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


makeKineticPlot(x = reactorA.time_simul*0.1,
                y = reactorA.cellActive_dyn[2],
                color = '#00B8FF',
                legend = 'Roseburia instestinalis',
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-',
                legendSize = 10)


title = 'stateA'
#plt.savefig(os.path.join(Path(os.getcwd()).parents[2], 'files', 'Figures', 'multistability', title + '.png'), transparent=True, dpi=600)
plt.show()


