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
    cmap = cmr.rainforest  
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
    
    if title is not None:
        plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'contours', title + '.png'), transparent=True, dpi=600)
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

    met_composition = metsB/sum(metsB) - metsA/sum(metsA)

    met_labels = ['Pyru', 'Gluc', 'Treh', 'Mann','Acet', 'Lact', 'Succ', 'Buty']

    met_colors = ['#ff8900', '#ff0000', '#8900ff', '#024059', '#003eff', '#00ffaf', '#00ff26', '#ff00a1']


    plot_stacked_bar_charts(bac_composition, bac_labels, bac_colors, 'Rel_Ab', met_composition, met_labels, met_colors, 'Rel_C_mM')
    return bac_composition,metsA,metsB



pH_points = np.linspace(5,6.5,30)
dilution_rate_points = np.linspace(0,3,30)



bh_mA = []
bh_mB = []
bh_b = []

for ph in tqdm(pH_points):
    for d in tqdm(dilution_rate_points):
        bac, metA, metB = makeSimulation(pHControl=ph,
                                          dilutionRate=d,
                                          bhA = 0,
                                          bhB =0.0,
                                          bt = 0.0,
                                          ri = 0.0)
        bh_b.append(bac)
        bh_mA.append(metA)
        bh_mB.append(metB)


bt_mA = []
bt_mB = []
bt_b = []

for ph in tqdm(pH_points):
    for d in tqdm(dilution_rate_points):
        bac, metA, metB = makeSimulation(pHControl=ph,
                                          dilutionRate=d,
                                          bhA = 0,
                                          bhB =0.0,
                                          bt = 0.01,
                                          ri = 0.0)
        bt_b.append(bac)
        bt_mA.append(metA)
        bt_mB.append(metB)




ri_mA = []
ri_mB = []
ri_b = []

for ph in tqdm(pH_points):
    for d in tqdm(dilution_rate_points):
        bac, metA, metB = makeSimulation(pHControl=ph,
                                          dilutionRate=d,
                                          bhA = 0,
                                          bhB =0.0,
                                          bt = 0.0,
                                          ri = 0.01)
        ri_b.append(bac)
        ri_mA.append(metA)
        ri_mB.append(metB)







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




pH_points = np.linspace(5,6.5,30)
dilution_rate_points = np.linspace(0,3,30)/15


pyru_bh = np.array([i[0] for i in bh_mB])
gluc_bh = np.array([i[1] for i in bh_mB])
treh_bh = np.array([i[2] for i in bh_mB])


acet_bh = np.array([i[4] for i in bh_mB])
lact_bh = np.array([i[5] for i in bh_mB])

bh_bh = np.array([i[0] for i in bh_b])
bh_bhbtri = np.array([i[0] for i in bhbtri2_b])


pyru_bh = pyru_bh/np.max(pyru_bh)
makeContour(dilution_rate_points, 
            pH_points, 
            pyru_bh,   
            'dilution rate ($h^{-1}$)', 
            'pH',
            'pyruvate_bh.png'
            )



gluc_bh = gluc_bh/np.max(gluc_bh)
makeContour(dilution_rate_points, 
            pH_points, 
            gluc_bh, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'glucose_bh.png'
            )


treh_bh = treh_bh/np.max(treh_bh)
makeContour(dilution_rate_points, 
            pH_points, 
            treh_bh, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'trehalose_bh.png'
            )



acet_bh = acet_bh/np.max(acet_bh)
makeContour(dilution_rate_points, 
            pH_points, 
            acet_bh, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'acetate_bh.png'
            )


lact_bh = lact_bh/np.max(lact_bh)
makeContour(dilution_rate_points, 
            pH_points, 
            lact_bh, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'lactate_bh.png'
            )


bh_bh = bh_bh/np.max(bh_bh)
makeContour(dilution_rate_points, 
            pH_points, 
            bh_bh, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'bh_bh.png'
            )

bh_bhbtri = bh_bhbtri/np.max(bh_bhbtri)
makeContour(dilution_rate_points, 
            pH_points, 
            bh_bhbtri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'bh_bhbtri.png'
            )


pyru_bt = np.array([i[0] for i in bt_mB])
gluc_bt = np.array([i[1] for i in bt_mB])
mann_bt = np.array([i[3] for i in bt_mB])


acet_bt = np.array([i[4] for i in bt_mB])
lact_bt = np.array([i[5] for i in bt_mB])
succ_bt = np.array([i[6] for i in bt_mB])

bt_bt = np.array([i[1] for i in bt_b])
bt_bhbtri = np.array([i[1] for i in bhbtri2_b])


pyru_bt = pyru_bt/np.max(pyru_bt)
makeContour(dilution_rate_points, 
            pH_points, 
            pyru_bt,   
            'dilution rate ($h^{-1}$)', 
            'pH',
            'pyruvate_bt.png'
            )



gluc_bt = gluc_bt/np.max(gluc_bt)
makeContour(dilution_rate_points, 
            pH_points, 
            gluc_bt, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'glucose_bt.png'
            )


mann_bt = mann_bt/np.max(mann_bt)
makeContour(dilution_rate_points, 
            pH_points, 
            mann_bt, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'mannose_bt.png'
            )



acet_bt = acet_bt/np.max(acet_bt)
makeContour(dilution_rate_points, 
            pH_points, 
            acet_bt, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'acetate_bt.png'
            )


lact_bt = lact_bt/np.max(lact_bt)
makeContour(dilution_rate_points, 
            pH_points, 
            lact_bt, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'lactate_bt.png'
            )


succ_bt = succ_bt/np.max(succ_bt)
makeContour(dilution_rate_points, 
            pH_points, 
            succ_bt, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'succinate_bt.png'
            )


bt_bt = bt_bt/np.max(bt_bt)
makeContour(dilution_rate_points, 
            pH_points, 
            bt_bt, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'bt_bt.png'
            )

bt_bhbtri = bt_bhbtri/np.max(bt_bhbtri)
makeContour(dilution_rate_points, 
            pH_points, 
            bt_bhbtri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'bt_bhbtri.png'
            )


pyru_ri = np.array([i[0] for i in ri_mB])
gluc_ri = np.array([i[1] for i in ri_mB])


acet_ri = np.array([i[4] for i in ri_mB])
lact_ri = np.array([i[5] for i in ri_mB])
buty_ri = np.array([i[7] for i in ri_mB])

ri_ri = np.array([i[2] for i in ri_b])
ri_bhbtri = np.array([i[2] for i in bhbtri2_b])


pyru_ri = pyru_ri/np.max(pyru_ri)
makeContour(dilution_rate_points, 
            pH_points, 
            pyru_ri,   
            'dilution rate ($h^{-1}$)', 
            'pH',
            'pyruvate_ri.png'
            )



gluc_ri = gluc_ri/np.max(gluc_ri)
makeContour(dilution_rate_points, 
            pH_points, 
            gluc_ri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'glucose_ri.png'
            )




acet_ri = acet_ri/np.max(acet_ri)
makeContour(dilution_rate_points, 
            pH_points, 
            acet_ri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'acetate_ri.png'
            )


lact_ri = lact_ri/np.max(lact_ri)
makeContour(dilution_rate_points, 
            pH_points, 
            lact_ri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'lactate_ri.png'
            )


buty_ri = buty_ri/np.max(buty_ri)
makeContour(dilution_rate_points, 
            pH_points, 
            buty_ri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'butyrate_ri.png'
            )


ri_ri = ri_ri/np.max(ri_ri)
makeContour(dilution_rate_points, 
            pH_points, 
            ri_ri, 
            'dilution rate ($h^{-1}$)', 
            'pH',
            'ri_ri.png'
            )

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
