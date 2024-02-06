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

from sklearn.decomposition import PCA
pca = PCA(n_components=1)


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
                   dilutionFactor = 0.5, 
                   bh = 0.0, 
                   bt = 0.01,
                   ri = 0.01,
                   perturbation = False):
    
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
    
    if perturbation:
        predictpHB = mockpHfunc(wc.metabolites,pH=5.6)
        pHB =  predictpHB(wc.get_concentration())

        #get the feed media and the reactor media
        wc_feedB = createMetabolome(db, 'wc', pHB, pHFunc=predictpHB)
        wc_reactorB = createMetabolome(db, 'wc', pHB, pHFunc=predictpHB)
     
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
    reactor_microbiome.subpopD['xa'].count = bh
    reactor_microbiome.subpopD['xe'].count = bt
    reactor_microbiome.subpopD['xi'].count = ri
    
    if perturbation:
        batchA = Pulse(wc_feed, feed_microbiome, 0, 600, 100, 0, 0, dilutionFactor,dilutionFactor)

        batchB = Pulse(wc_feed, feed_microbiome, 600, 5040, 100, 0, 0, dilutionFactor,dilutionFactor)

        
        
        #simulate
        reactorA = Reactor(reactor_microbiome, wc_reactorB,[
                                                          batchA,
                                                          
                                                          
                                                          
                                                           ], 15)



        reactorA.simulate()
        #reactorA.makePlots()



        for i in wc_reactorB.metD:
            wc_reactor.metD[i].update(wc_reactorB.metD[i].concentration)
            

        reactor = Reactor(reactorA.microbiome, wc_reactor,[
                                                          
                                                          batchB,
                                                          
                                                           ], 15)

        reactor.simulate()

        
        
        
        
    else:
        
        batchA = Pulse(wc_feed, feed_microbiome, 0, 5040, 100, 0, 0, dilutionFactor,dilutionFactor)
    
        #simulate
        reactor = Reactor(reactor_microbiome, wc_reactor,[
                                                      batchA
                                                       ], 15)
        reactor.simulate()
    
    b = reactor.cellActive_dyn.T[-1]
    bac_composition = b
    pyru = reactor.met_simul[reactor.metabolome.metabolites.index('pyruvate')]
    gluc = reactor.met_simul[reactor.metabolome.metabolites.index('glucose')]
    treh = reactor.met_simul[reactor.metabolome.metabolites.index('trehalose')]
    mann = reactor.met_simul[reactor.metabolome.metabolites.index('mannose')]
    acet = reactor.met_simul[reactor.metabolome.metabolites.index('acetate')]
    form = reactor.met_simul[reactor.metabolome.metabolites.index('formate')]
    lact = reactor.met_simul[reactor.metabolome.metabolites.index('lactate')]
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
    pH = reactor.pH_simul[-1]
    return bac_composition,metsA,metsB, pH




    
    
simulation_points = 150 

pH_points = np.linspace(5,6.5,simulation_points)


bhbtri2_mA = []
bhbtri2_mB = []
bhbtri2_b = []

for pH in tqdm(pH_points):
    
    
    bac, metA, metB,_ = makeSimulation(pHControl=pH,
                                      dilutionFactor=1.0,
                                     
                                      bh = 0.003,
                                      bt = 0.003,
                                      ri = 0.003)
    bhbtri2_b.append(bac)
    bhbtri2_mA.append(metA)
    bhbtri2_mB.append(metB)



bh_bhbtri = np.array([i[0] for i in bhbtri2_b])
bt_bhbtri = np.array([i[1] for i in bhbtri2_b])
ri_bhbtri = np.array([i[2] for i in bhbtri2_b])

pyru_bhbtri = np.array([i[0] for i in bhbtri2_mB])
gluc_bhbtri = np.array([i[1] for i in bhbtri2_mB])
treh_bhbtri = np.array([i[2] for i in bhbtri2_mB])
mann_bhbtri = np.array([i[3] for i in bhbtri2_mB])
acet_bhbtri = np.array([i[4] for i in bhbtri2_mB])
lact_bhbtri = np.array([i[5] for i in bhbtri2_mB])
form_bhbtri = np.array([i[6] for i in bhbtri2_mB])
succ_bhbtri = np.array([i[7] for i in bhbtri2_mB])
buty_bhbtri = np.array([i[8] for i in bhbtri2_mB])


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
        'Roseburia intestinalis'
        ]

create_heatmap(dataM, 
               rows, 
               np.linspace(0,simulation_points-1, 6), 
               np.round(np.linspace(5, 6.5, 6),3), 
               'pH', 
               None, 
               None, 
               fileName = None) #os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'dilution.png'))



pc = pca.fit_transform(dataM.T).flatten()


bhbtri1_mA = []
bhbtri1_mB = []
bhbtri1_b = []

for pH in tqdm(pH_points):
    
    
    
    
    bac, metA, metB,_ = makeSimulation(pHControl=pH,
                                      dilutionFactor=1.0,
                                      bh = 0.003,
                                      bt = 0.003,
                                      ri = 0.003,
                                      perturbation=True)
    bhbtri1_b.append(bac)
    bhbtri1_mA.append(metA)
    bhbtri1_mB.append(metB)



bh_bhbtri_p = np.array([i[0] for i in bhbtri1_b])
bt_bhbtri_p = np.array([i[1] for i in bhbtri1_b])
ri_bhbtri_p = np.array([i[2] for i in bhbtri1_b])

pyru_bhbtri_p = np.array([i[0] for i in bhbtri1_mB])
gluc_bhbtri_p = np.array([i[1] for i in bhbtri1_mB])
treh_bhbtri_p = np.array([i[2] for i in bhbtri1_mB])
mann_bhbtri_p = np.array([i[3] for i in bhbtri1_mB])
acet_bhbtri_p = np.array([i[4] for i in bhbtri1_mB])
lact_bhbtri_p = np.array([i[5] for i in bhbtri1_mB])
form_bhbtri_p = np.array([i[6] for i in bhbtri1_mB])
succ_bhbtri_p = np.array([i[7] for i in bhbtri1_mB])
buty_bhbtri_p = np.array([i[8] for i in bhbtri1_mB])


dataM_p = np.array([pyru_bhbtri_p, 
                    gluc_bhbtri_p, 
                    treh_bhbtri_p, 
                    mann_bhbtri_p, 
                    acet_bhbtri_p, 
                    lact_bhbtri_p, 
                    form_bhbtri_p, 
                    succ_bhbtri_p, 
                    buty_bhbtri_p, 
                    bh_bhbtri_p, 
                    bt_bhbtri_p, 
                    ri_bhbtri_p])

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
        'Roseburia intestinalis'
        ]

create_heatmap(dataM_p, 
               rows, 
               np.linspace(0,149, 6), 
               np.round(np.linspace(5, 6.5, 6),3), 
               'pH', 
               None, 
               None, 
               fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'pH_perturbation.png'))

pc_p = pca.fit_transform(dataM_p.T).flatten()



pH_values = np.linspace(5,6.5,simulation_points)
plt.plot(pH_values[pc>5], pc[pc>5], 'o', color='red')

plt.plot(pH_values[pc_p<5], pc_p[pc_p<5], 'o', color='red')
plt.vlines(x = max(pH_values[pc_p>2]), 
            ymin = min(pc[pc<2])-1,
            ymax = max(pc_p[pc_p>2])+1,
            linestyles = 'dashed',
            color='k')

plt.vlines(x = min(pH_values[pc<2]), 
            ymin = min(pc[pc<2])-1,
            ymax = max(pc_p[pc_p>2])+1,
            linestyles = 'dashed',
            color='k')
plt.xlabel('pH')
plt.ylabel('first principle component score')
plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'pH_hysterysis.png'), transparent=True, dpi=600)
plt.show()