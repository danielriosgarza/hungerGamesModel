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
import matplotlib.pyplot as plt
plt.style.use('seaborn-bright')

from parseGenExpData import *

geneFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'ri', 'genes')
modelFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'ri', 'gsmms') 

figuresFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'riExperiments')


model = cobra.io.read_sbml_model(os.path.join(modelFolder, 'ri_final.xml'))

riGE = GeneExpr(geneFolder = geneFolder, 
                tpmFile = 'ri_tpm.txt', 
                groupLabels = ['t04', 't12', 't48', 'btri_ri_t12', 'btri_ri_t36', 'mono_ri', 'btri_ri'], 
                groupIDX = [[4,5,6],
                            [7,8,9],
                            [10,11,12],
                            [13,14,15],
                            [19,20,21],
                            [7,8,9,10,11,12],
                            [13,14,15,16,17,18,19,20,21]
                            ], 
                groupComparison = {('t04', 't12'):'t04vst12_deseq.txt',
                                   ('t04', 't48'): 't04vst48_deseq.txt',
                                   ('t12', 't48'): 't12vst48_deseq.txt',
                                   ('btri_ri_t12', 'btri_ri_t36'):'btri_ri_t12vst36_deseq.txt',
                                   ('mono_ri', 'btri_ri'):'btri__ri_t12t48vst12t24t36_deseq.txt',
                                   ('t12', 'btri_ri_t12'):'ri_t12vsbtri_ri_t12_deseq.txt'},
                            
                
                featuresFile =  'ri_BVBRC_genome_feature.txt', 
                sbmlModel = model,
                rootID = 'patric',
                species='ri'
                )   


fcRange = 6

#########################t04 vs 12 #########################
gA = 't04'
gB = 't12'

reactionList = []

with open(os.path.join(geneFolder, 'wcReactions.txt')) as f:
    for line in f:
        reactionList.append(line.strip())


reactionList = np.array(reactionList)

group = (gA, gB)

x,p,g = extracReactions(riGE, reactionList, group)

sorter = np.argsort(x)

x = x[sorter]
g = g[sorter]
reactionList = reactionList[sorter]
p = p[sorter]


labels = reactionList.copy()


labels[labels=='PFK(ppi)'] = 'PFK'
labels[labels=='FDNADOX_H'] = 'OXrd'
labels[labels=='LDH_L'] = 'LDH'
labels[labels=='ACACT1r'] = 'ACACT'
labels[labels=='ACKr'] = 'ACK'
labels[labels=='ACOAD1i'] = 'ACOAD'
labels[labels=='ECOAH1'] = 'ECOAH'
labels[labels=='HACD1'] = 'HACD'


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

ax.text(-fcRange/2, len(x) + .2, gA)
ax.text(fcRange/2, len(x) + .2, gB)

plt.rcParams["figure.figsize"] = (3 ,8)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)




#########################t04 vs 48 #########################
gA = 't04'
gB = 't48'


group = (gA, gB)

x,p,g = extracReactions(riGE, reactionList, group)



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

ax.text(-fcRange/2, len(x) + .2, gA)
ax.text(fcRange/2, len(x) + .2, gB)

plt.rcParams["figure.figsize"] = (3 ,8)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)

#########################t12 vs 48 #########################
gA = 't12'
gB = 't48'


group = (gA, gB)

x,p,g = extracReactions(riGE, reactionList, group)



c= np.array(['#6B6ACF']*len(x))



c[(p < 0.05) & (x > 0)] = '#BD00FF'

c[(p < 0.05) & (x < 0)] = '#001EFF'





fig, ax = plt.subplots()

ax.vlines(0,0, len(x), color = 'k', lw=1, ls='--')

ax.barh(y = np.arange(len(x)), width = x, color=c, edgecolor='k', zorder=5)

ax.hlines(np.arange(len(x)),-fcRange, fcRange, color = 'k', lw=.5, ls='--')    


ax.set_yticks(np.arange(len(x)), labels=labels, fontsize=10)



ax.set_xlim(-fcRange, fcRange)
ax.set_xlabel('gene expression FC')

ax.text(-fcRange/2, len(x) + .2, gA)
ax.text(fcRange/2, len(x) + .2, gB)

plt.rcParams["figure.figsize"] = (3 ,8)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)


#########################'btri_ri_t12' vs 'btri_ri_t36'#########################
gA = 'btri_ri_t12'
gB = 'btri_ri_t36'


group = (gA, gB)

x,p,g = extracReactions(riGE, reactionList, group)



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

plt.rcParams["figure.figsize"] = (3 ,8)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)

#########################'mono_ri' vs 'btri_ri'#########################
gA = 'mono_ri'
gB = 'btri_ri'


group = (gA, gB)

x,p,g = extracReactions(riGE, reactionList, group)



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

plt.rcParams["figure.figsize"] = (3 ,8)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)



#########################'t12' vs 'btri_ri_t12'#########################
gA = 't12'
gB = 'btri_ri_t12'


group = (gA, gB)

x,p,g = extracReactions(riGE, reactionList, group)



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

plt.rcParams["figure.figsize"] = (3 ,8)
plt.tight_layout()

plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)
