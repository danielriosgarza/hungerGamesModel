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
from tqdm import tqdm



import plotly.io as pio
pio.renderers.default='browser'




def getComposition(compVec):
    rounded_comp = np.round(compVec, 2)
    binary_comp = np.where(rounded_comp > 0, 1, 0)
    
    composition_dict = {
        (1, 0, 0): 1,
        (0, 1, 0): 2,
        (0, 0, 1): 3,
        (1, 0, 1): 4,
        (1, 1, 0): 5,
        (0, 1, 1): 6,
        (1, 1, 1): 7,
    }
    
    return composition_dict.get(tuple(binary_comp), 0)
    
    
def makeContour(x, y, z, label, xlabel, ylabel, title):
    
    xd, yd = np.meshgrid(x, y)

    zd= np.array(z).reshape(len(x),len(y))


    # Create a contour plot
    plt.figure(figsize=(10, 8))
    contour = plt.contourf(xd, yd, zd, cmap='viridis')

    # Add colorbar
    plt.colorbar(contour, label=label)

    # Add labels and title
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    # Show the plot
    plt.show()

def mockpHfunc(metObj, pH=7.0):
    def pHfunc(metObj):
        return pH
    return pHfunc

def makeSimulation(pHControl = None, dilutionRate = 0.5, filePath = None):
    
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
    reactor_microbiome.subpopD['xi'].count = 0.01




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

    metsA = np.array([pyru[0]*wc_feed.metD['pyruvate'].carbons,
        gluc[0]*wc_feed.metD['glucose'].carbons,
        treh[0]*wc_feed.metD['trehalose'].carbons,
        mann[0]*wc_feed.metD['mannose'].carbons,
        acet[0]*wc_feed.metD['acetate'].carbons,
        lact[0]*wc_feed.metD['lactate'].carbons,
        succ[0]*wc_feed.metD['succinate'].carbons,
        buty[0]*wc_feed.metD['butyrate'].carbons
        ])



    metsB = np.array([pyru[-1]*wc_feed.metD['pyruvate'].carbons,
        gluc[-1]*wc_feed.metD['glucose'].carbons,
        treh[-1]*wc_feed.metD['trehalose'].carbons,
        mann[-1]*wc_feed.metD['mannose'].carbons,
        acet[-1]*wc_feed.metD['acetate'].carbons,
        lact[-1]*wc_feed.metD['lactate'].carbons,
        succ[-1]*wc_feed.metD['succinate'].carbons,
        buty[-1]*wc_feed.metD['butyrate'].carbons
        ])

    met_composition = metsB/sum(metsB) - metsA/sum(metsA)

    met_labels = ['Pyru', 'Gluc', 'Treh', 'Mann','Acet', 'Lact', 'Succ', 'Buty']

    met_colors = ['#ff8900', '#ff0000', '#8900ff', '#024059', '#003eff', '#00ffaf', '#00ff26', '#ff00a1']


    plot_stacked_bar_charts(bac_composition, bac_labels, bac_colors, 'Rel_Ab', met_composition, met_labels, met_colors, 'Rel_C_mM')
    return bac_composition,met_composition


m = []
b = []

for ph in tqdm(np.linspace(5,8,25)):
    for d in tqdm(np.linspace(0,1,25)):
        bac, met = makeSimulation(pHControl=ph,dilutionRate=d)
        b.append(bac)
        m.append(met)


pH_points = np.linspace(5,8,25)
dilution_rate_points = np.linspace(0,1,25)

bh= [i[0] for i in b]

makeContour(pH_points, 
            dilution_rate_points, 
            bh, 
            'Bh', 
            'pH', 
            'dilution factor',
            'Contour Plot of Bh vs. pH and Dilution Rate'
            )

bt= [i[1] for i in b]
makeContour(pH_points, 
            dilution_rate_points, 
            bt, 
            'Bt', 
            'pH', 
            'dilution factor',
            'Contour Plot of Bt vs. pH and Dilution Rate'
            )

ri= [i[2] for i in b]
makeContour(pH_points, 
            dilution_rate_points, 
            ri, 
            'Ri', 
            'pH', 
            'dilution factor',
            'Contour Plot of Ri vs. pH and Dilution Rate'
            )

composition = [getComposition(i) for i in b]

makeContour(pH_points, 
            dilution_rate_points, 
            composition, 
            'composition', 
            'pH', 
            'dilution factor',
            'Contour Plot of composition vs. pH and Dilution Rate'
            )



treh = [i[2] for i in m]

makeContour(pH_points, 
            dilution_rate_points, 
            treh, 
            'trehalose', 
            'pH', 
            'dilution factor',
            'Contour Plot of trehalose vs. pH and Dilution Rate'
            )
