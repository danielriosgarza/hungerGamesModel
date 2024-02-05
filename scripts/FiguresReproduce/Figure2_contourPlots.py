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
import cmasher as cmr

from readModelDB import *
from mainClasses import *
from loadParameters import *
from general import *
from tqdm import tqdm



import plotly.io as pio
pio.renderers.default='browser'

from matplotlib.colors import ListedColormap






def getComposition(compVec):
    rounded_comp = np.round(compVec, 10)
    binary_comp = np.where(rounded_comp > 0.000000001, 1, 0)
    
    # Convert array to a tuple of its elements
    binary_comp_tuple = tuple(binary_comp.tolist())

    composition_dict = {
        (1, 0, 0): 1,
        (0, 1, 0): 2,
        (0, 0, 1): 3,
        (1, 0, 1): 4,
        (1, 1, 0): 5,
        (0, 1, 1): 6,
        (1, 1, 1): 7,
    }
    
    return composition_dict.get(binary_comp_tuple, 0)

    
    
def makeContour(x, y, z, xlabel, ylabel, title, vmin = 0, vmax = 1, cbar=False, cmap = 'cmr.lavender'):
    cmap = cmap  
    cmap = plt.get_cmap(cmap)   # MPL
    xd, yd = np.meshgrid(x, y)
    zd = np.array(z).reshape(len(x), len(y))

    # Create a contour plot
    plt.figure(figsize=(5, 5))
    contour = plt.contourf(xd, yd, zd, cmap=cmap, vmin = vmin, vmax = vmax)

    # Add colorbar and set fontsize
    if cbar:
        colorbar = plt.colorbar(contour)
        colorbar.ax.tick_params(labelsize=14)  # Set the fontsize for colorbar labels

    # Add labels and title
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    # plt.title(title)

    # Show the plot
    plt.tight_layout()
    if title is not None:
        plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'contours', title),transparent=True, dpi=600)
    plt.show()










def makeDiscreteContour(x, y, z, xlabel, ylabel, title):
    # Define colors for 0, a placeholder for 1, and then 2-7
    custom_colors = ['#121313',  # Color for 0
                     '#121313',  # Placeholder color for 1 (not in data)
                     '#ff8300',  # Color for 2
                     '#00B8FF',  # Color for 3
                     '#984ea3',  # Color for 4
                     '#e41a1c',  # Color for 5
                     '#377eb8',  # Color for 6
                     '#4daf4a']  # Color for 7
    custom_cmap = ListedColormap(custom_colors)

    xd, yd = np.meshgrid(x, y)
    zd = np.array(z).reshape(len(x), len(y))

    plt.figure(figsize=(5, 5))
    # Levels covering 0 to 7 (including placeholder for 1)
    levels = np.arange(-0.5, 8.5, 1)
    contour = plt.contourf(xd, yd, zd, levels=levels, cmap=custom_cmap, extend='both')
    cbar = plt.colorbar(contour, ticks=np.arange(0, 8))  # Ticks for 0 to 7

    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    
    plt.tight_layout()
    if title is not None:
        plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'contours', title + '.png'), transparent=True, dpi=600)
    plt.show()








def mockpHfunc(metObj, pH=7.0):
    def pHfunc(metObj):
        return pH
    return pHfunc


def makeSimulation(pHControl = None, 
                   dilutionFactor = 0.5, 
                   bh = 0.0, 
                   bt = 0.01,
                   ri = 0.01):
    
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
    batchA = Pulse(wc_feed, feed_microbiome, 0, 6000, 100, 0, 0, dilutionFactor,dilutionFactor)

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
    lact = reactor.met_simul[reactor.metabolome.metabolites.index('lactate')]
    succ = reactor.met_simul[reactor.metabolome.metabolites.index('succinate')]
    buty = reactor.met_simul[reactor.metabolome.metabolites.index('butyrate')]

    metsA = np.array([pyru[0],
        gluc[0],
        treh[0],
        mann[0],
        acet[0],
        lact[0],
        succ[0],
        buty[0]
        ])
    metsB = np.array([pyru[-1],
        gluc[-1],
        treh[-1],
        mann[-1],
        acet[-1],
        lact[-1],
        succ[-1],
        buty[-1]
        ])
    
    return bac_composition,metsA,metsB





pH_points = np.linspace(5,7.0,100)
dilution_rate_points = np.linspace(0,3,100)





bhbtri2_mA = []
bhbtri2_mB = []
bhbtri2_b = []

for ph in tqdm(pH_points):
    for d in tqdm(dilution_rate_points):
        
        
        
        
        bac, metA, metB = makeSimulation(pHControl=ph,
                                          dilutionRate=d,
                                         
                                          bhA = 0.003,
                                          bhB =0.0,
                                          bt = 0.003,
                                          ri = 0.003)
        bhbtri2_b.append(bac)
        bhbtri2_mA.append(metA)
        bhbtri2_mB.append(metB)




pH_points = np.linspace(5,7.0,100)
dilution_rate_points = np.linspace(0,3,100)/15



bh_bhbtri = bh_bhbtri/np.max(bh_bhbtri)
makeContour(dilution_rate_points, 
            pH_points, 
            bh_bhbtri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'bh_bhbtri.png'
            )


bt_bhbtri = np.array([i[1] for i in bhbtri2_b])



bt_bhbtri = bt_bhbtri/np.max(bt_bhbtri)
makeContour(dilution_rate_points, 
            pH_points, 
            bt_bhbtri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'bt_bhbtri.png'
            )


ri_bhbtri = np.array([i[2] for i in bhbtri2_b])



ri_bhbtri = ri_bhbtri/np.max(ri_bhbtri)
makeContour(dilution_rate_points, 
            pH_points, 
            ri_bhbtri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'ri_bhbtri.png'
            )


pyru_bhbtri = np.array([i[0] for i in bhbtri2_mB])
gluc_bhbtri = np.array([i[1] for i in bhbtri2_mB])
treh_bhbtri = np.array([i[2] for i in bhbtri2_mB])
mann_bhbtri = np.array([i[3] for i in bhbtri2_mB])


acet_bhbtri = np.array([i[4] for i in bhbtri2_mB])
lact_bhbtri = np.array([i[5] for i in bhbtri2_mB])
succ_bhbtri = np.array([i[6] for i in bhbtri2_mB])
buty_bhbtri = np.array([i[7] for i in bhbtri2_mB])

comp = np.array([getComposition(i) for i in bhbtri2_b])

makeDiscreteContour(dilution_rate_points, 
            pH_points, 
            comp,   
            'dilution rate ($h^{-1}$)', 
            'pH',
            'coexistence_bhbtri2.png')



pyru_bhbtri = pyru_bhbtri/np.max(pyru_bhbtri)
makeContour(dilution_rate_points, 
            pH_points, 
            pyru_bhbtri,   
            'dilution rate ($h^{-1}$)', 
            'pH',
            'pyruvate_bhbtri.png'
            )



gluc_bhbtri = gluc_bhbtri/np.max(gluc_bhbtri)
makeContour(dilution_rate_points, 
            pH_points, 
            gluc_bhbtri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'glucose_bhbtri.png'
            )



treh_bhbtri = treh_bhbtri/np.max(treh_bhbtri)
makeContour(dilution_rate_points, 
            pH_points, 
            treh_bhbtri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'trehalose_bhbtri.png'
            )



mann_bhbtri = mann_bhbtri/np.max(mann_bhbtri)
makeContour(dilution_rate_points, 
            pH_points, 
            mann_bhbtri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'mannose_bhbtri.png'
            )

acet_bhbtri = acet_bhbtri/np.max(acet_bhbtri)
makeContour(dilution_rate_points, 
            pH_points, 
            acet_bhbtri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'acetate_bhbtri.png'
            )


lact_bhbtri = lact_bhbtri/np.max(lact_bhbtri)
makeContour(dilution_rate_points, 
            pH_points, 
            lact_bhbtri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'lactate_bhbtri.png'
            )




succ_bhbtri = succ_bhbtri/np.max(succ_bhbtri)
makeContour(dilution_rate_points, 
            pH_points, 
            succ_bhbtri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'succinate_bhbtri.png'
            )


buty_bhbtri = buty_bhbtri/np.max(buty_bhbtri)
makeContour(dilution_rate_points, 
            pH_points, 
            buty_bhbtri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'butyrate_bhbtri.png'
            )

makeContour(dilution_rate_points, 
            pH_points, 
            buty_bhbtri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'butyrate_bhbtri_cbar.png',
            cbar=True
            )
