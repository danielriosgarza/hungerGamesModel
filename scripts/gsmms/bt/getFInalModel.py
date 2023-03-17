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


patricPath = os.path.join(Path(os.getcwd()).parents[2], 'files', 'strainSummaries', 'bt', 'genes', 'bt_BVBRC_genome_feature.txt')
in_modelPath = os.path.join(Path(os.getcwd()).parents[2], 'files', 'strainSummaries', 'bt', 'gsmms')
out_modelPath = os.path.join(in_modelPath, 'efms')
model = cobra.io.read_sbml_model(os.path.join(in_modelPath, 'iAH991.xml'))

#add RNF
'2.0 fdxrd[c] + 2.0 h[c] + nad[c] <=> 2.0 fdxox[c] + h[e] + nadh[c]'
fdxrd_c = model.metabolites.get_by_id('fdxrd[c]')
fdxox_c = model.metabolites.get_by_id('fdxox[c]')
h_c = model.metabolites.get_by_id('h[c]')
h_e = model.metabolites.get_by_id('h[e]')
nad_c =  model.metabolites.get_by_id('nad[c]')
nadh_c =  model.metabolites.get_by_id('nadh[c]')

reaction = Reaction('RNF')

reaction.name = 'ferredoxin---NAD+ oxidoreductase (H+-transporting)'


reaction.add_metabolites({fdxrd_c:-2, h_c:-2, nad_c:-1, fdxox_c:2, h_e:2, nadh_c:1})
reaction.lower_bound=-1000
reaction.upper_bound=1000
reaction.gene_reaction_rule = '( BT_0617 and BT_0619 and BT_0620 and BT_0622 and BT_0618 and BT_0621 )'

model.add_reactions([reaction])

model.repair()
model.optimize()


modelOut = 'bt_final.xml'

cobra.io.write_sbml_model(model, os.path.join(in_modelPath, modelOut))