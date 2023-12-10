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
bt_params = getPramsFromFile('bhbtri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'btri.tsv'))
ri_params = getPramsFromFile('bhbtri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'btri.tsv'))


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
fixedpHA = mockpHfunc(wc.metabolites,pH=5.6)
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
reactor_microbiome.subpopD['xa'].count = 0.01
reactor_microbiome.subpopD['xb'].count = 0.00
reactor_microbiome.subpopD['xe'].count = 0.01
reactor_microbiome.subpopD['xf'].count = 0.0
reactor_microbiome.subpopD['xi'].count = 0.01
reactor_microbiome.subpopD['xj'].count = 0.


d = 1.0
d2 = 1.0
d3 = 0.85




batchA = Pulse(wc_feedA, feed_microbiome, 0, 600, 100, 0, 0, d,d)

batchB = Pulse(wc_feedA, feed_microbiome, 600, 700, 100, 0, 0, d2,d2)

batchC = Pulse(wc_feedA, feed_microbiome, 700, 1200, 100, 0, 0, d,d)


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


title = 'stateD'
plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', title + '.png'), transparent=True, dpi=600)


plt.show()

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.subpop_simul[1],
#                 color = '#FF2E2EC9',
#                 legend = 'Bh glucose suppop',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '-')








# b = reactor.cellActive_dyn.T[-1]

# bac_composition = b/sum(b)

# bac_labels = ['Bh', 'Bt', 'Ri']

# bac_colors = ['#FF10F0', '#ff8300', '#00B8FF']



# pyru = reactor.met_simul[reactor.metabolome.metabolites.index('pyruvate')]
# gluc = reactor.met_simul[reactor.metabolome.metabolites.index('glucose')]
# treh = reactor.met_simul[reactor.metabolome.metabolites.index('trehalose')]
# mann = reactor.met_simul[reactor.metabolome.metabolites.index('mannose')]
# acet = reactor.met_simul[reactor.metabolome.metabolites.index('acetate')]
# lact = reactor.met_simul[reactor.metabolome.metabolites.index('lactate')]
# succ = reactor.met_simul[reactor.metabolome.metabolites.index('succinate')]
# buty = reactor.met_simul[reactor.metabolome.metabolites.index('butyrate')]

# metsA = np.array([pyru[0]*wc_feed.metD['pyruvate'].carbons,
#         gluc[0]*wc_feed.metD['glucose'].carbons,
#         treh[0]*wc_feed.metD['trehalose'].carbons,
#         mann[0]*wc_feed.metD['mannose'].carbons,
#         acet[0]*wc_feed.metD['acetate'].carbons,
#         lact[0]*wc_feed.metD['lactate'].carbons,
#         succ[0]*wc_feed.metD['succinate'].carbons,
#         buty[0]*wc_feed.metD['butyrate'].carbons
#         ])



# metsB = np.array([pyru[-1]*wc_feed.metD['pyruvate'].carbons,
#         gluc[-1]*wc_feed.metD['glucose'].carbons,
#         treh[-1]*wc_feed.metD['trehalose'].carbons,
#         mann[-1]*wc_feed.metD['mannose'].carbons,
#         acet[-1]*wc_feed.metD['acetate'].carbons,
#         lact[-1]*wc_feed.metD['lactate'].carbons,
#         succ[-1]*wc_feed.metD['succinate'].carbons,
#         buty[-1]*wc_feed.metD['butyrate'].carbons
#         ])

# met_composition = metsB/sum(metsB) - metsA/sum(metsA)

# met_labels = ['Pyru', 'Gluc', 'Treh', 'Mann','Acet', 'Lact', 'Succ', 'Buty']

# met_colors = ['#ff8900', '#ff0000', '#8900ff', '#024059', '#003eff', '#00ffaf', '#00ff26', '#ff00a1']


# plot_stacked_bar_charts(bac_composition, bac_labels, bac_colors, 'Rel_Ab', met_composition, met_labels, met_colors, 'Rel_C_mM')



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
#                 legend = 'Bh trehalose suppop',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '-')

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.subpop_simul[1],
#                 color = '#FF2E2EC9',
#                 legend = 'Bh glucose suppop',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '-')

# plt.title('dilution rate {0:.3f}'.format(d/15))
# #savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'BhSubpopsDilutionRate.png'), dpi = 600)
# plt.show()


# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.cellActive_dyn[0],
#                 color = '#FF6EC7',
#                 legend = 'live_bh_cells (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '-')

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.cellActive_dyn[1],
#                 color = '#FF5F1F',
#                 legend = 'live_bt_cells (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '-')


# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.cellActive_dyn[2],
#                 color = '#1F51FF',
#                 legend = 'live_ri_cells (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$10^5$ cells/uL',
#                 title = None,
#                 linestyle = '-')

# plt.title('dilution rate {0:.3f}'.format(d/15))
# plt.tight_layout()
# #plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', '30hperturbCells.png'), dpi = 600)

# plt.show()

# makeKineticPlot(x = reactor.time_simul,
#                 y = reactor.met_simul[reactor.metabolome.metabolites.index('trehalose')],
#                 color = '#8900ff',
#                 legend = 'trehalose (simul)',
#                 xlabel = 'time (h)',
#                 ylabel = '$mM$',
#                 title = None,
#                 linestyle = '-')

# plt.title('dilution rate {0:.3f}'.format(d/15))
# plt.tight_layout()
# #plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', '30hperturbTreh.png'), dpi = 600)
# plt.show()

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


