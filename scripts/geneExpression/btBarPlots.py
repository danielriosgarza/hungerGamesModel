# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:10:05 2023

@author: danie
"""

import os
from pathlib import Path
import numpy as np
import scipy.stats as sts
import cobra
from aquarel import load_theme
theme = load_theme("minimal_light")
theme.apply()

import matplotlib.pyplot as plt

from parseGenExpData import *

geneFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bt', 'genes')
modelFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bt', 'gsmms') 

figuresFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'btExperiments')


model = cobra.io.read_sbml_model(os.path.join(modelFolder, 'bt_final.xml'))

btGE = GeneExpr(geneFolder = geneFolder, 
                tpmFile = 'bt_tpm.txt', 
                groupLabels = ['t04', 't12', 't36', 'btri_bt_t12', 'btri_bt_t36', 'mono_bt', 'btri_bt'], 
                groupIDX = [[4,5,6],
                            [7,8,9],
                            [10,11,12],
                            [13,14,15],
                            [19,20,21],
                            [7,8,9,10,11,12],
                            [13,14,15,16,17,18,19,20,21]
                            ], 
                groupComparison = {('t04', 't12'):'t04vst12_deseq.txt',
                                   ('t04', 't36'): 't04vst36_deseq.txt',
                                   ('t12', 't36'): 't12vst36_deseq.txt',
                                   ('t12', 'btri_bt_t12'):'t12vs_btri_bt_t12_deseq.txt',
                                   ('t36', 'btri_bt_t36'):'t36vs_btri_bt_t36_deseq.txt',
                                   ('mono_bt', 'btri_bt'):'t12t36vs_btri_bt_t12t24t36_deseq.txt'},
                            
                
                featuresFile =  'bt_BVBRC_genome_feature.txt', 
                sbmlModel = model,
                rootID = None
                )   


fcRange = 8
# #########################t04 vs 12 #########################
gA = 't04'
gB = 't12'

reactionList = []

with open(os.path.join(geneFolder, 'wcReactions.txt')) as f:

    for line in f:
        reactionList.append(line.strip())


reactionList = np.array(reactionList)

group = (gA, gB)

x,p,g = extracReactions(btGE, reactionList, group)

sorter = np.argsort(x)

x = x[sorter]
g = g[sorter]
reactionList = reactionList[sorter]
p = p[sorter]

labels = reactionList.copy()
labels[labels=='FRD2'] = 'FRD'
labels[labels=='ACKr'] = 'ACK'
labels[labels=='LDH_L'] = 'LDH'
labels[labels=='POR4i'] = 'POR4'
labels[labels=='PPCKr'] = 'PPCKr'
labels[labels=='PTAr'] = 'PTA'


c= np.array(['#6B6ACF']*len(x))

c[(p < 0.01) & (x > 0)] = '#001eff'

c[(p < 0.01) & (x < 0)] = '#44A043'




plt.rcParams["figure.figsize"] = (3 ,9)
fig, ax = plt.subplots()

ax.vlines(0,0, len(x), color = 'k', lw=1, ls='--')

ax.barh(y = np.arange(len(x)), width = x, color=c, edgecolor='k', zorder=5)

ax.hlines(np.arange(len(x)),-fcRange, fcRange, color = 'k', lw=.5, ls='--')    

ax.set_yticks(np.arange(len(x)), labels=labels, fontsize=10)



ax.set_xlim(-fcRange, fcRange)
ax.set_xlabel('gene expression FC')

ax.text(-fcRange/2, len(x) + .2, gA)
ax.text(fcRange/2, len(x) + .2, gB)

plt.rcParams["figure.figsize"] = (3 ,9)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)





#########################t04 vs 36 #########################
gA = 't04'
gB = 't36'


group = (gA, gB)

x,p,g = extracReactions(btGE, reactionList, group)



c= np.array(['#6B6ACF']*len(x))

c[(p < 0.01) & (x > 0)] = '#bd00ff'

c[(p < 0.01) & (x < 0)] = '#44A043'





fig, ax = plt.subplots()

ax.vlines(0,0, len(x), color = 'k', lw=1, ls='--')

ax.barh(y = np.arange(len(x)), width = x, color=c, edgecolor='k', zorder=5)

ax.hlines(np.arange(len(x)),-fcRange, fcRange, color = 'k', lw=.5, ls='--')    

ax.set_yticks(np.arange(len(x)), labels=labels, fontsize=10)



ax.set_xlim(-fcRange, fcRange)
ax.set_xlabel('gene expression FC')

ax.text(-fcRange/2, len(x) + .2, gA)
ax.text(fcRange/2, len(x) + .2, gB)

plt.rcParams["figure.figsize"] = (3 ,9)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)

#########################t12 vs 36 #########################
gA = 't12'
gB = 't36'


group = (gA, gB)

x,p,g = extracReactions(btGE, reactionList, group)



c= np.array(['#6B6ACF']*len(x))


c[(p < 0.01) & (x > 0)] = '#bd00ff'

c[(p < 0.01) & (x < 0)] = '#001eff'







fig, ax = plt.subplots()

ax.vlines(0,0, len(x), color = 'k', lw=1, ls='--')

ax.barh(y = np.arange(len(x)), width = x, color=c, edgecolor='k', zorder=5)

ax.hlines(np.arange(len(x)),-fcRange, fcRange, color = 'k', lw=.5, ls='--')    


ax.set_yticks(np.arange(len(x)), labels=labels, fontsize=10)



ax.set_xlim(-fcRange, fcRange)
ax.set_xlabel('gene expression FC')

ax.text(-fcRange/2, len(x) + .2, gA)
ax.text(fcRange/2, len(x) + .2, gB)

plt.rcParams["figure.figsize"] = (3 ,9)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)


#########################'t12' vs 'btri_bt_t12'#########################
gA = 't12'
gB = 'btri_bt_t12'


group = (gA, gB)

x,p,g = extracReactions(btGE, reactionList, group)



c= np.array(['#6B6ACF']*len(x))



c[(p < 0.05) & (x > 0)] = '#001eff'

c[(p < 0.05) & (x < 0)] = '#44A043'





fig, ax = plt.subplots()

ax.vlines(0,0, len(x), color = 'k', lw=1, ls='--')

ax.barh(y = np.arange(len(x)), width = x, color=c, edgecolor='k', zorder=5)

ax.hlines(np.arange(len(x)),-fcRange, fcRange, color = 'k', lw=.5, ls='--')    


ax.set_yticks(np.arange(len(x)), labels=labels, fontsize=10)



ax.set_xlim(-fcRange, fcRange)
ax.set_xlabel('gene expression FC')

ax.text(-7.5, len(x) + .2, gA, rotation=15)
ax.text(0.5, len(x) + .2, gB, rotation=15)

plt.rcParams["figure.figsize"] = (3 ,9)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)





#########################'t36' vs 'btri_bt_t36'#########################
gA = 't36'
gB = 'btri_bt_t36'


group = (gA, gB)

x,p,g = extracReactions(btGE, reactionList, group)



c= np.array(['#6B6ACF']*len(x))



c[(p < 0.05) & (x > 0)] = '#001eff'

c[(p < 0.05) & (x < 0)] = '#44A043'





fig, ax = plt.subplots()

ax.vlines(0,0, len(x), color = 'k', lw=1, ls='--')

ax.barh(y = np.arange(len(x)), width = x, color=c, edgecolor='k', zorder=5)

ax.hlines(np.arange(len(x)),-fcRange, fcRange, color = 'k', lw=.5, ls='--')    


ax.set_yticks(np.arange(len(x)), labels=labels, fontsize=10)



ax.set_xlim(-fcRange, fcRange)
ax.set_xlabel('gene expression FC')

ax.text(-7.5, len(x) + .2, gA, rotation=15)
ax.text(0.5, len(x) + .2, gB, rotation=15)

plt.rcParams["figure.figsize"] = (3 ,9)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)





#########################'mono_ri' vs 'btri_ri'#########################
gA = 'mono_bt'
gB = 'btri_bt'


group = (gA, gB)

x,p,g = extracReactions(btGE, reactionList, group)



c= np.array(['#6B6ACF']*len(x))



c[(p < 0.05) & (x > 0)] = '#001eff'

c[(p < 0.05) & (x < 0)] = '#44A043'





fig, ax = plt.subplots()

ax.vlines(0,0, len(x), color = 'k', lw=1, ls='--')

ax.barh(y = np.arange(len(x)), width = x, color=c, edgecolor='k', zorder=5)

ax.hlines(np.arange(len(x)),-fcRange, fcRange, color = 'k', lw=.5, ls='--')    


ax.set_yticks(np.arange(len(x)), labels=labels, fontsize=10)



ax.set_xlim(-fcRange, fcRange)
ax.set_xlabel('gene expression FC')

ax.text(-7.5, len(x) + .2, gA, rotation=15)
ax.text(0.5, len(x) + .2, gB, rotation=15)

plt.rcParams["figure.figsize"] = (3 ,9)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)


