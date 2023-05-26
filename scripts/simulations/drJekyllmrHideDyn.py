# -*- coding: utf-8 -*-
"""
Created on Sat May 20 12:11:08 2023

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
bt_params = getPramsFromFile('bt', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'bt.tsv'))
ri_params = getPramsFromFile('ri', os.path.join(Path(os.getcwd()).parents[1], 'files', 'params', 'ri.tsv'))


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
reactor_microbiome.subpopD['xa'].count = 0.01
reactor_microbiome.subpopD['xe'].count = 0.01
reactor_microbiome.subpopD['xi'].count = 0.01
reactor_microbiome.subpopD['xb'].count = 0.00

d = 0.9
d2 = 0
d3 = 1.4

b = []


t = 0
for i in range(50):
    
    if i%2==0:
        t0=t
        t=t0+24
        b.append(Pulse(wc_feed, feed_microbiome, t0, t, 100, 0, 0, d2,d2))
    else:
        t0 = t
        t=t0+24
        b.append(Pulse(wc_feed, feed_microbiome, t0, t, 100, 0, 0, d,d))



#simulate
reactor = Reactor(reactor_microbiome, wc_reactor,b, 15)
reactor.simulate()
reactor.makePlots()


