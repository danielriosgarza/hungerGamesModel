# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:11:52 2023

@author: drgarza
"""
from pathlib import Path
import os
import sys

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'db'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))
import cmasher as cmr

from readModelDB import *
from mainClasses import *
from loadParameters import *
from general import *
from tqdm import tqdm



import plotly.io as pio
pio.renderers.default='browser'

from matplotlib.colors import ListedColormap


def create_heatmap(data, row_labels, x_positions, x_values, xlabel, ylabel, title, fileName = None):
    fig, ax = plt.subplots()
    heatmap = ax.imshow(data, aspect='auto')
    ax.grid(False)

    # Set row labels (y-axis)
    ax.set_yticks(np.arange(data.shape[0]))
    
    yticklabels = ax.set_yticklabels(row_labels)
    yticklabels[-1].set_fontstyle('italic')
    yticklabels[-2].set_fontstyle('italic')
    yticklabels[-3].set_fontstyle('italic')

    # Set specific x-ticks and labels (x-axis)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_values)

    # Set labels and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    # Add grid
    ax.set_xticks(np.arange(-.5, data.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, data.shape[0], 1), minor=True)
    ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", size=0)

    # Add colorbar
    cbar = ax.figure.colorbar(heatmap, ax=ax)
    plt.tight_layout()
    
    if fileName is not None:
        plt.savefig(fileName, transparent=True, dpi=600)

    plt.show()


def mockpHfunc(metObj, pH=7.0):
    def pHfunc(metObj):
        return pH
    return pHfunc

def makeSimulation(pHControl = None, 
                   dilutionRate = 0.5, 
                   bhA = 0.0, 
                   bhB = 0.01,
                   bt = 0.01,
                   ri = 0.01,
                   pyruvate = None,
                   glucose = None,
                   trehalose = None,
                   mannose = None):
    
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
    
    if pHControl is None:
        predictpH = getpH(wc.metabolites, ipH_path)
    else:
        predictpH = mockpHfunc(wc.metabolites,pH=pHControl)
        
        
    pH =  predictpH(wc.get_concentration())

    #get the feed media and the reactor media
    wc_feed = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    wc_reactor = createMetabolome(db, 'wc', pH, pHFunc=predictpH)

    #wc_feed_pH = createMetabolome(db, 'wc', 2.5, pHFunc=predictpH)
    
    if pyruvate is not None:
        wc_reactor.metD['pyruvate'].update(pyruvate)
        wc_feed.metD['pyruvate'].update(pyruvate)
        
    if glucose is not None:
        wc_reactor.metD['glucose'].update(glucose)
        wc_feed.metD['glucose'].update(glucose)

    if trehalose is not None:
        wc_reactor.metD['trehalose'].update(trehalose)
        wc_feed.metD['trehalose'].update(trehalose)

    if mannose is not None:
        wc_reactor.metD['mannose'].update(trehalose)
        wc_feed.metD['mannose'].update(trehalose)        
  

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
    reactor_microbiome.subpopD['xa'].count = bhA
    reactor_microbiome.subpopD['xb'].count = bhB
    reactor_microbiome.subpopD['xe'].count = bt
    reactor_microbiome.subpopD['xi'].count = ri




    batchA = Pulse(wc_feed, feed_microbiome, 0, 2000, 100, 0, 0, dilutionRate,dilutionRate)

    #simulate
    reactor = Reactor(reactor_microbiome, wc_reactor,[
                                                  batchA
                                                   ], 15)
    reactor.simulate()
    #reactor.makePlots()



    b = reactor.cellActive_dyn.T[-1]

    bac_composition = b#/sum(b)

    bac_labels = ['Bh', 'Bt', 'Ri']

    bac_colors = ['#FF10F0', '#ff8300', '#00B8FF']



    pyru = reactor.met_simul[reactor.metabolome.metabolites.index('pyruvate')]
    gluc = reactor.met_simul[reactor.metabolome.metabolites.index('glucose')]
    treh = reactor.met_simul[reactor.metabolome.metabolites.index('trehalose')]
    mann = reactor.met_simul[reactor.metabolome.metabolites.index('mannose')]
    acet = reactor.met_simul[reactor.metabolome.metabolites.index('acetate')]
    lact = reactor.met_simul[reactor.metabolome.metabolites.index('lactate')]
    form = reactor.met_simul[reactor.metabolome.metabolites.index('formate')]
    succ = reactor.met_simul[reactor.metabolome.metabolites.index('succinate')]
    buty = reactor.met_simul[reactor.metabolome.metabolites.index('butyrate')]

    metsA = np.array([pyru[0],
        gluc[0],
        treh[0],
        mann[0],
        acet[0],
        lact[0],
        form[0],
        succ[0],
        buty[0]
        ])



    metsB = np.array([pyru[-1],
        gluc[-1],
        treh[-1],
        mann[-1],
        acet[-1],
        lact[-1],
        form[-1],
        succ[-1],
        buty[-1]
        ])

    met_composition = metsB/sum(metsB) - metsA/sum(metsA)

    met_labels = ['Pyru', 'Gluc', 'Treh', 'Mann','Acet', 'Lact','Form', 'Succ', 'Buty']

    met_colors = ['#ff8900', '#ff0000', '#8900ff', '#024059', '#003eff', '#00ffaf', '#00c6ff', '#00ff26', '#ff00a1']


    plot_stacked_bar_charts(bac_composition, bac_labels, bac_colors, 'Rel_Ab', met_composition, met_labels, met_colors, 'Rel_C_mM')
    return bac_composition,metsA,metsB


#pH_points = np.linspace(5,6.5,30)
dilution_rate_points = np.linspace(0,3,150)


bhbtri2_mA = []
bhbtri2_mB = []
bhbtri2_b = []

for d in tqdm(dilution_rate_points):
    
    
    
    
    bac, metA, metB = makeSimulation(pHControl=None,
                                      dilutionRate=d,
                                     
                                      bhA = 0.003,
                                      bhB =0.0,
                                      bt = 0.003,
                                      ri = 0.003)
    bhbtri2_b.append(bac)
    bhbtri2_mA.append(metA)
    bhbtri2_mB.append(metB)



bh_bhbtri = np.array([i[0] for i in bhbtri2_b])
#bh_bhbtri = bh_bhbtri/max(bh_bhbtri)

bt_bhbtri = np.array([i[1] for i in bhbtri2_b])
#bt_bhbtri = bt_bhbtri/max(bt_bhbtri)

ri_bhbtri = np.array([i[2] for i in bhbtri2_b])
#ri_bhbtri = ri_bhbtri/max(ri_bhbtri)


pyru_bhbtri = np.array([i[0] for i in bhbtri2_mB])
#pyru_bhbtri = pyru_bhbtri/max(pyru_bhbtri)

gluc_bhbtri = np.array([i[1] for i in bhbtri2_mB])
#gluc_bhbtri = gluc_bhbtri/max(gluc_bhbtri)


treh_bhbtri = np.array([i[2] for i in bhbtri2_mB])
#treh_bhbtri = treh_bhbtri/max(treh_bhbtri)

mann_bhbtri = np.array([i[3] for i in bhbtri2_mB])
#mann_bhbtri = mann_bhbtri/max(mann_bhbtri)

acet_bhbtri = np.array([i[4] for i in bhbtri2_mB])
#acet_bhbtri = acet_bhbtri/max(acet_bhbtri)


lact_bhbtri = np.array([i[5] for i in bhbtri2_mB])
#lact_bhbtri = lact_bhbtri/max(lact_bhbtri)

form_bhbtri = np.array([i[6] for i in bhbtri2_mB])
#succ_bhbtri = succ_bhbtri/max(succ_bhbtri)


succ_bhbtri = np.array([i[7] for i in bhbtri2_mB])
#succ_bhbtri = succ_bhbtri/max(succ_bhbtri)


buty_bhbtri = np.array([i[8] for i in bhbtri2_mB])
#buty_bhbtri = buty_bhbtri/max(buty_bhbtri)


dataM = np.array([pyru_bhbtri, gluc_bhbtri, treh_bhbtri, mann_bhbtri, acet_bhbtri, lact_bhbtri, form_bhbtri, succ_bhbtri, buty_bhbtri, bh_bhbtri, bt_bhbtri, ri_bhbtri])

rows = ['pyruvate', 
        'glucose', 
        'trehalose',
        'mannose',
        'acetate',
        'lactate',
        'succinate',
        'formate',
        'butyrate',
        'Blautia hydrogenotrophica',
        'Bacteroides thetaiotaomicron',
        'Roseburia intinalis'
        ]

create_heatmap(dataM, rows, np.linspace(0,149, 6), np.round(np.linspace(0, 3/15, 6),3), 'dilution rate($h^{-1})$', None, None, fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'dilution.png'))