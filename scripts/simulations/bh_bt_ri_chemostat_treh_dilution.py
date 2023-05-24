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




import plotly.io as pio
pio.renderers.default='browser'





def getState(reactor):
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
    
    return d

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



xa_1, xb_1, xi_1, xj_1, xe_1, xf_1, treh_1, glc_1, pyr_1,phvar_1 = [],[],[],[],[],[],[],[],[],[]
xa_2, xb_2, xi_2, xj_2, xe_2, xf_2, treh_2, glc_2, pyr_2,phvar_2 = [],[],[],[],[],[],[],[],[],[]





counter = 0

for i in np.linspace(0.1,2.3, 50):
    
    print(counter)
    counter+=1
        
    
    #getStarting pH
    wc = createMetabolome(db, 'wc')
    predictpH = getpH(wc.metabolites, ipH_path)
    pH =  predictpH(wc.get_concentration())
    
    #get the feed media and the reactor media
    wc_feed = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    wc_reactor = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    #wc_reactor.metD['trehalose'].update(i)
    #wc_feed.metD['trehalose'].update(i)
    
    
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
    reactor_microbiome.subpopD['xb'].count = 0.00
    
    
    
    chemostatA = Pulse(wc_feed, feed_microbiome, 0, 700, 100, 0, 0, 0.406,0.406)
    chemostatB = Pulse(wc_feed, feed_microbiome, 700, 2000, 100, 0, 0, i,i)
    
    #simulate
    reactor = Reactor(reactor_microbiome, wc_reactor,[chemostatA, chemostatB], 15)
    reactor.simulate()
    #reactor.makePlots()
    xa_1.append(reactor.subpop_simul[reactor.microbiome.subpops.index('xa')][-1])
    xb_1.append(reactor.subpop_simul[reactor.microbiome.subpops.index('xb')][-1])
    xe_1.append(reactor.subpop_simul[reactor.microbiome.subpops.index('xe')][-1])
    xf_1.append(reactor.subpop_simul[reactor.microbiome.subpops.index('xf')][-1])
    xi_1.append(reactor.subpop_simul[reactor.microbiome.subpops.index('xi')][-1])
    xj_1.append(reactor.subpop_simul[reactor.microbiome.subpops.index('xj')][-1])
    treh_1.append(reactor.met_simul[reactor.metabolome.metabolites.index('trehalose')][-1])
    glc_1.append(reactor.met_simul[reactor.metabolome.metabolites.index('glucose')][-1])
    pyr_1.append(reactor.met_simul[reactor.metabolome.metabolites.index('pyruvate')][-1])
    phvar_1.append(reactor.pH)




counter = 0

for i in np.linspace(0.1,2.3, 50):
    
    print(counter)
    counter+=1
        
    
    #getStarting pH
    wc = createMetabolome(db, 'wc')
    predictpH = getpH(wc.metabolites, ipH_path)
    pH =  predictpH(wc.get_concentration())
    
    #get the feed media and the reactor media
    wc_feed = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    wc_reactor = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    #wc_reactor.metD['trehalose'].update(i)
    #wc_feed.metD['trehalose'].update(i)
    
    
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
    reactor_microbiome.subpopD['xb'].count = 0.00
    
    
    
    chemostatA = Pulse(wc_feed, feed_microbiome, 0, 700, 100, 0, 0, 1.29387755, 1.29387755)
    chemostatB = Pulse(wc_feed, feed_microbiome, 700, 2000, 100, 0, 0, i,i)
    
    #simulate
    reactor = Reactor(reactor_microbiome, wc_reactor,[chemostatA, chemostatB], 15)
    reactor.simulate()
    #reactor.makePlots()
    xa_2.append(reactor.subpop_simul[reactor.microbiome.subpops.index('xa')][-1])
    xb_2.append(reactor.subpop_simul[reactor.microbiome.subpops.index('xb')][-1])
    xe_2.append(reactor.subpop_simul[reactor.microbiome.subpops.index('xe')][-1])
    xf_2.append(reactor.subpop_simul[reactor.microbiome.subpops.index('xf')][-1])
    xi_2.append(reactor.subpop_simul[reactor.microbiome.subpops.index('xi')][-1])
    xj_2.append(reactor.subpop_simul[reactor.microbiome.subpops.index('xj')][-1])
    treh_2.append(reactor.met_simul[reactor.metabolome.metabolites.index('trehalose')][-1])
    glc_2.append(reactor.met_simul[reactor.metabolome.metabolites.index('glucose')][-1])
    pyr_2.append(reactor.met_simul[reactor.metabolome.metabolites.index('pyruvate')][-1])
    phvar_2.append(reactor.pH)


    
    
    
from pylab import *
s1 = np.array([xa_1, xb_1, xi_1, xj_1, xe_1, xf_1, treh_1, glc_1, pyr_1,phvar_1]).T
s2 = np.array([xa_2, xb_2, xi_2, xj_2, xe_2, xf_2, treh_2, glc_2, pyr_2,phvar_2]).T

cos_1 = np.array([cosine(s1[0], i) for i in s1])
cos_2 = np.array([cosine(s2[0], i) for i in s2])

x =  np.linspace(0.1,2.3, 50)

plot(x, cos_2, 'o-', color='r');plot(x,cos_1, 'o-', color='b')
show()

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


