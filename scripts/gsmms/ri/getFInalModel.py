# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 15:19:41 2023

@author: danie
"""

import os
from pathlib import Path

import cobra
from cobra import Reaction
import numpy as np



in_modelPath = os.path.join(Path(os.getcwd()).parents[2], 'files', 'strainSummaries', 'ri', 'gsmms')
out_modelPath = os.path.join(in_modelPath, 'efms')
model = cobra.io.read_sbml_model(os.path.join(in_modelPath, 'Roseburia_intestinalis_L1_82.xml'))
geneFolder = os.path.join(Path(os.getcwd()).parents[2], 'files', 'strainSummaries', 'ri', 'genes')

agora2patric = {}

agora2patric['536231.5.peg.2404'] = 'uknown'

with open(os.path.join(in_modelPath, 'agora2patric.tsv')) as f:
    for line in f:
        a = line.strip().split('\t')
        agora2patric[a[0]] = a[1]



for reac in model.reactions:
    if reac.genes != frozenset():
        
        genes = []
        
        for gene in reac.genes:
            genes.append(gene.id)
        genes.sort(key=len, reverse=True)
        rule = reac.gene_reaction_rule[:]
        
        
        
        for gene in genes:
            rule = rule.replace(gene, agora2patric[gene])[:]
        
        reac.gene_reaction_rule = rule[:]


# for i in model.reactions:
#     if i.id in toPatric:
#         model.reactions.get_by_id(i.id).gene_reaction_rule = toPatric[i.id]

cobra.manipulation.delete.remove_genes(model,['uknown'])

reactionList = []

with open(os.path.join(geneFolder, 'wcReactions.txt')) as f:
    for line in f:
        reactionList.append(line.strip())
        
#
for i in reactionList:
    if model.reactions.get_by_id(i).genes == frozenset():
        print(i)

#BTCOAACCOAT

model.reactions.get_by_id('BTCOAACCOAT').gene_reaction_rule = '536231.75.peg.425'


#ECOAH1
model.reactions.get_by_id('ECOAH1').gene_reaction_rule = '536231.75.peg.3357'

#FDNADOX_H
model.reactions.get_by_id('FDNADOX_H').gene_reaction_rule = '( 536231.75.peg.2207 and 536231.75.peg.2208 and 536231.75.peg.2209 and 536231.75.peg.2210 and 536231.75.peg.2211 and 536231.75.peg.2212)'

modelOut = 'ri_final.xml'

cobra.io.write_sbml_model(model, os.path.join(in_modelPath, modelOut))