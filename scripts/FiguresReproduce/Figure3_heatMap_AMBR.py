# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 14:05:22 2025

@author: drgarza
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Ellipse
import numpy as np


from aquarel import load_theme

theme = load_theme("boxy_light")
theme.apply()

#from sklearn.manifold import MDS as PCA
from sklearn.covariance import MinCovDet
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
pca = PCA(n_components=2, svd_solver ='randomized', power_iteration_normalizer='QR', iterated_power=10000)


def get_data(fileName):
    
    strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'ambr')
    
    
    
    with open(os.path.join(strainSummaryFolder, fileName)) as f:
        labels = f.readline().strip().split('\t')[6::]
        experiments = []
        conditions = []
        description = []
        time = []
        timeCondition = []
        replicate = []
        exclude= [ ]
        data = {i:[] for i in labels if i not in exclude}
        
        for line in f:
            a = line.strip().split('\t')
            experiments.append(a[0])
            conditions.append(a[3])
            description.append(a[2])
            
            for i in range(len(labels)):
                if labels[i] not in exclude:
                    data[labels[i]].append(float(a[i+6]))
            time.append(float(a[4]))
            timeCondition.append(float(a[5]))
            
            ambrExp = a[0].split('_')
            replicate.append(ambrExp[0] + '_' + ambrExp[1] + '_' + a[4])
            
    
    dataM = np.array([np.array(data[i]) for i in labels if i not in exclude]).T
    return labels, np.array(experiments), np.array(description), conditions, time, np.array(timeCondition), replicate, dataM





def makeHeatMap(data, labels, figPath = None):
    data = data.T
    size_x = data.shape[0]
    fig, ax = plt.subplots()
    heatmap = ax.imshow(data, aspect='auto', vmin = 0, vmax =12)
    ax.grid(False)
    ax.set_yticks(np.arange(data.shape[0]))
    yticklabels = ax.set_yticklabels(labels)
    
    
    # Set specific x-ticks and labels (x-axis)
    
    
    #ax.set_xticks([0.0, 13.0, 22.0, 28.0])
    #ax.set_xticklabels([0.0, 13.0, 22.0, 28.0])
    
        # # Set labels and title
        # ax.set_xlabel(xlabel)
        # ax.set_ylabel(ylabel)
        # ax.set_title(title)
        
    # Add grid
    ax.set_xticks(np.arange(-.5, data.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, data.shape[0], 1), minor=True)
    ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", size=0)
    
    # Add colorbar
    cbar = ax.figure.colorbar(heatmap, ax=ax)
    plt.tight_layout()
    if figPath is not None:
        plt.savefig(figPath, dpi = 600)
    
    plt.show()
    
    
    
#############################HeatMap###############


labels, experiments, description, conditions, time, timeCond, replicate, dataM = get_data('ambrSS_alt.txt')


ambr_8 = np.array([True if i[0:3]=='A8_' else False for i in experiments])
not_ambr_8 = ~ambr_8

idx_control = (description == 'control conditions') & not_ambr_8
idx_control_ambr_8 = (description == 'control conditions') & ambr_8

feed = description == 'post feed perturbation'

pHfeed = (description == 'post pH/feed perturbation') & not_ambr_8
pHfeed_ambr_8 = (description == 'post pH/feed perturbation') & ambr_8

control = dataM[idx_control]
pH_feed = dataM[pHfeed]
feed = dataM[feed]

control_ambr_8 = dataM[idx_control_ambr_8]
pH_feed_ambr_8 = dataM[pHfeed_ambr_8]


# Save the figure
fileName1 = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'composition_heatMap_controlb.png')

fileName2 = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'composition_heatMap_pH_feedb.png')

fileName3 = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'composition_heatMap_feedb.png')

fileName4 = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'composition_heatMap_controlb_ambr8.png')

fileName5 = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'composition_heatMap_pH_feedb_ambr8.png')

makeHeatMap(control, labels, fileName1)
makeHeatMap(pH_feed, labels, fileName2)
makeHeatMap(feed, labels, fileName3)

makeHeatMap(control_ambr_8, labels, fileName4)
makeHeatMap(pH_feed_ambr_8, labels, fileName5)
