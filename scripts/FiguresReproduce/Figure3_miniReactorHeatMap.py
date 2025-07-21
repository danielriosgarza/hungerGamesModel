# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 07:31:19 2023

@author: drgarza
"""

import os
import sys
from pathlib import Path

from aquarel import load_theme

theme = load_theme("boxy_light")
theme.apply()

import matplotlib.pyplot as plt

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))
from general import *


def create_heatmap(data, row_labels, x_positions, x_values, xlabel, ylabel, title, fileName=None):
    fig, ax = plt.subplots()
    
    # Use interpolation='none' to avoid color bleeding
    heatmap = ax.imshow(data, aspect='auto', interpolation='none')
    ax.grid(False)

    # Set row labels (y-axis)
    ax.set_yticks(np.arange(data.shape[0]))
    
    yticklabels = ax.set_yticklabels(row_labels)
    yticklabels[-1].set_fontstyle('italic')
    yticklabels[-2].set_fontstyle('italic')
    yticklabels[-3].set_fontstyle('italic')

    # Set specific x-ticks and labels (x-axis)
    #ax.set_xticks(x_positions)
    #ax.set_xticklabels(x_values, rotation=90, fontsize=4)

    # Set labels and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    # Add grid
    ax.set_xticks(np.arange(-.5, data.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, data.shape[0], 1), minor=True)
    ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", size=0)

    # Add colorbar
    cbar = ax.figure.colorbar(heatmap, ax=ax)
    plt.tight_layout()
    
    if fileName is not None:
        plt.savefig(fileName, transparent=True, dpi=600)

    plt.show()







strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'ambr')



with open(os.path.join(strainSummaryFolder, 'ambrSS.txt')) as f:
    labels = f.readline().strip().split('\t')[5::]
    experiments = []
    condition = []
    data = {i:[] for i in labels}
    
    for line in f:
        a = line.strip().split('\t')
        experiments.append(a[0])
        condition.append(a[1])
        for i in range(len(labels)):
            data[labels[i]].append(float(a[i+5]))



cols = []

for i in condition:
    if 'E1_c' in i:
        cols.append('red')
    
    if 'E1_p' in i:
        cols.append('blue')
    
    if 'E2_c' in i:
        cols.append('green')
    if 'E2_p' in i:
        cols.append('purple')
        
    if i=='E2_Co':
        cols.append('black')
dataM = np.array([np.array(data[i]) for i in labels]).T

create_heatmap(dataM.T, labels, np.arange(len(experiments)), condition, None, None, None, fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambrHM.png'))



