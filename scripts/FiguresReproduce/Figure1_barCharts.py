# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:10:05 2023

@author: danie
"""
import sys
import os
from pathlib import Path
import numpy as np
import scipy.stats as sts
import cobra
from aquarel import load_theme
theme = load_theme("minimal_light")
theme.apply()


import matplotlib.pyplot as plt


project_root = Path(os.getcwd()).parents[0]
sys.path.append(os.path.join(project_root, 'geneExpression'))
from parseGenExpData import *


#################################################################################################################
###############################Blautia hydrogenotryphica#########################################################
#################################################################################################################



geneFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh', 'genes')
modelFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bh', 'gsmms') 
figuresFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'bhExperiments')


#load the GSMM
model = cobra.io.read_sbml_model(os.path.join(modelFolder, 'bh_final.xml'))

#load the gene expression class
bhGE = GeneExpr(geneFolder = geneFolder, 
                tpmFile = 'bh_tpm.txt', 
                groupLabels = ['t14', 't32', 't72', 't32 and t72'], 
                groupIDX = [[4,5,6],
                            [7,8,9],
                            [10,11,12],
                            [7,8,9,10,11,12]
                            ], 
                groupComparison = {('t14', 't32'):'t14vst32_deseq.txt',
                                   ('t14', 't72'): 't14vst72_deseq.txt',
                                   ('t32', 't72'): 't32vst72_deseq.txt'}, 
                featuresFile =  'bh_BVBRC_features.txt', 
                sbmlModel = model,
                species='bh'
                )   

fcRange = 16


#########################t14 vs 32 #########################
gA = 't14'
gB = 't32'

reactionList = []

with open(os.path.join(geneFolder, 'wcReactions.txt')) as f:
    for line in f:
        reactionList.append(line.strip())


reactionList = np.array(reactionList)

group = (gA, gB)

x,p,g = extracReactions(bhGE, reactionList, group)

sorter = np.argsort(x)

x = x[sorter]
g = g[sorter]
reactionList = reactionList[sorter]
p = p[sorter]


labels = reactionList.copy()

#Change labels for better readability
labels[labels=='CODH_ACS'] = 'CODH'
labels[labels=='OOR2r'] = 'OOR2'
labels[labels=='OOR2r'] = 'OOR2'
labels[labels=='PFK(ppi)'] = 'PFK'
labels[labels=='HYDFDNrfdx'] = 'HydABC'
labels[labels=='GLUt2r'] = 'GLUt'
labels[labels=='Rnf'] = 'RNF'
labels[labels=='LDH_L'] = 'LDH'
labels[labels=='LDH_L'] = 'LDH'


#########################t14 vs 32 #########################
gA = 't14'
gB = 't32'


group = (gA, gB)

x,p,g = extracReactions(bhGE, reactionList, group)

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

ax.text(-fcRange/2, len(x) + .1, gA)
ax.text(fcRange/2, len(x) + .1, gB)

plt.rcParams["figure.figsize"] = (3 ,9)
plt.tight_layout()

#plt.savefig(os.path.join(figuresFolder,  'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)

plt.show()

#########################t14 vs 72 #########################
gA = 't14'
gB = 't72'


group = (gA, gB)

x,p,g = extracReactions(bhGE, reactionList, group)

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

ax.text(-fcRange/2, len(x) + .1, gA)
ax.text(fcRange/2, len(x) + .1, gB)

plt.rcParams["figure.figsize"] = (3 ,9)
plt.tight_layout()

#plt.savefig(os.path.join(figuresFolder,  'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)

plt.show()
#########################t32 vs 72 #########################
gA = 't32'
gB = 't72'


group = (gA, gB)

x,p,g = extracReactions(bhGE, reactionList, group)

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

ax.text(-fcRange/2, len(x) + .1, gA)
ax.text(fcRange/2, len(x) + .1, gB)

plt.rcParams["figure.figsize"] = (3 ,9)
plt.tight_layout()

#plt.savefig(os.path.join(figuresFolder,  'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)
plt.show()

#################################################################################################################
#################################################################################################################
#################################################################################################################



#################################################################################################################
##############################Bacteroides thetaiotaomicron#######################################################
#################################################################################################################

geneFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bt', 'genes')
modelFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bt', 'gsmms') 

figuresFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'btExperiments')


model = cobra.io.read_sbml_model(os.path.join(modelFolder, 'bt_final.xml'))

btGE = GeneExpr(geneFolder = geneFolder, 
                tpmFile = 'bt_tpm.txt', 
                groupLabels = ['t04', 't12', 't36'], 
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
                                   },
                            
                
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

#plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)
plt.show()

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

#plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)
plt.show()
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

#plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)
plt.show()

#################################################################################################################
#################################################################################################################
#################################################################################################################


#################################################################################################################
##############################Roseburia intestinalis#############################################################
#################################################################################################################
geneFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'ri', 'genes')
modelFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'ri', 'gsmms') 

figuresFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'riExperiments')


model = cobra.io.read_sbml_model(os.path.join(modelFolder, 'ri_final.xml'))

riGE = GeneExpr(geneFolder = geneFolder, 
                tpmFile = 'ri_tpm.txt', 
                groupLabels = ['t04', 't12', 't48'], 
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
                                   },
                            
                
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

#plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)
plt.show()



#########################t04 vs 48 #########################
gA = 't04'
gB = 't48'


group = (gA, gB)

x,p,g = extracReactions(riGE, reactionList, group)



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

#plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)
plt.show()

#########################t12 vs 48 #########################
gA = 't12'
gB = 't48'


group = (gA, gB)

x,p,g = extracReactions(riGE, reactionList, group)



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

#plt.savefig(os.path.join(figuresFolder, 'geneExp_' + gA + 'vs' + gB + '.png'), dpi = 300)
plt.show()



