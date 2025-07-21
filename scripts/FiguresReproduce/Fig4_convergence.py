# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 11:10:47 2025

@author: drgarza
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Ellipse
import numpy as np

from scipy.spatial.distance import cosine
from scipy.spatial.distance import braycurtis
from scipy.spatial.distance import correlation


from aquarel import load_theme
import ternary

theme = load_theme("boxy_light")
theme.apply()


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


# Add jitter to the data
def add_jitter(arr, scale=0.5):
    np.random.seed(666)
    return arr + np.random.uniform(-scale, scale, size=arr.shape)

def get_cosine_d(chunk_array, w = 2):
    distances = []
    
    #chunk_array = [(i.T/i.max(axis=1)).T for i in control_data_chunks]
    for i in range(len(chunk_array)):
        
        nd = chunk_array[i] / np.max(chunk_array[i], axis=1, keepdims=True)
        cosine_distances = np.array([
            braycurtis(nd[z], nd[z+w])
            for z in range(len(nd) - w)
        ])
        distances.append(cosine_distances[:])
    
    return distances

def plot_cosine_chunks(time_chunks, cosine_chunks, highlight=100, fileName=None, ball_color="tab:blue"):
    plt.figure(figsize=(8, 5))

    for i in range(len(cosine_chunks)):
        y = cosine_chunks[i]
        t = time_chunks[i][0:len(y)]  # Match X and Y dynamically
        plt.plot(t, y, color='k', linewidth=0.5, zorder=1)
        plt.scatter(t, y, color=ball_color, edgecolors='black', s=50, zorder=2)

    t_max = max(np.concatenate(time_chunks).flatten())
    plt.axvspan(highlight, t_max, color='lightgray', alpha=0.5, zorder=0)

    plt.ylim(0, 1)
    plt.axvline(highlight, linestyle='--', color='gray', linewidth=1.5)
    plt.xlabel('Time (h)', fontsize=20)
    plt.ylabel('Bray-Curtis dissimilarity', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()

    if fileName is not None:
        plt.savefig(fileName, dpi=600, transparent=True)
    plt.show()



labels, experiments, description, conditions, time, timeCond, replicate, dataM = get_data('ambrAll.txt')

ambr_8 = np.array([True if i[0:3]=='A8_' else False for i in experiments])
ambr_6 = np.array([True if i[0:3]=='A6_' else False for i in experiments])
not_ambr_8 = ~ambr_8
not_ambr_6 = ~ambr_6

idx_control = (description == 'control conditions') & not_ambr_8
idx_control_ambr_8 = (description == 'control conditions') & ambr_8

idx_feed = (description == 'post feed perturbation')

idx_pHfeed = (description == 'post pH/feed perturbation') & not_ambr_8
idx_pHfeed_ambr_8 = (description == 'post pH/feed perturbation') & ambr_8

#Split the control samples by their time points

control_data = dataM[idx_control]

control_time = np.array(time)[idx_control]

control_zero_indices = np.where(control_time == 0)[0]
control_zero_indices = np.append(control_zero_indices, len(control_time))

control_time_chunks = [control_time[control_zero_indices[i]:control_zero_indices[i+1]] for i in range(len(control_zero_indices)-1) if control_zero_indices[i+1] > control_zero_indices[i]]

control_data_chunks = [control_data[control_zero_indices[i]:control_zero_indices[i+1]] for i in range(len(control_zero_indices)-1) if control_zero_indices[i+1] > control_zero_indices[i]]

control_cosine = get_cosine_d(control_data_chunks)


fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'convergence_control.png')

plot_cosine_chunks(control_time_chunks, control_cosine, 84, fileName, '#00FF00')



#Split the feed samples by their time points

feed_data = dataM[idx_feed]

feed_time = np.array(time)[idx_feed]

feed_time -= 132
#feed_time *=feed_time>=0

feed_zero_indices = np.where(feed_time == 0)[0]
feed_zero_indices = np.append(feed_zero_indices, len(feed_time))

feed_time_chunks = [feed_time[feed_zero_indices[i]:feed_zero_indices[i+1]] for i in range(len(feed_zero_indices)-1) if feed_zero_indices[i+1] > feed_zero_indices[i]]

feed_data_chunks = [feed_data[feed_zero_indices[i]:feed_zero_indices[i+1]] for i in range(len(feed_zero_indices)-1) if feed_zero_indices[i+1] > feed_zero_indices[i]]

feed_cosine = get_cosine_d(feed_data_chunks)

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'convergence_feed.png')

plot_cosine_chunks(feed_time_chunks, feed_cosine, 78, fileName, '#FF007F')


#Split the feed samples by their time points

pHfeed_data = dataM[idx_pHfeed]

pHfeed_time = np.array(time)[idx_pHfeed]

pHfeed_time -= 176
#pHfeed_time *=feed_time>=0

pHfeed_zero_indices = np.where(pHfeed_time == 0)[0]
pHfeed_zero_indices = np.append(pHfeed_zero_indices, len(pHfeed_time))

pHfeed_time_chunks = [pHfeed_time[pHfeed_zero_indices[i]:pHfeed_zero_indices[i+1]] for i in range(len(pHfeed_zero_indices)-1) if pHfeed_zero_indices[i+1] > pHfeed_zero_indices[i]]

pHfeed_data_chunks = [pHfeed_data[pHfeed_zero_indices[i]:pHfeed_zero_indices[i+1]] for i in range(len(pHfeed_zero_indices)-1) if pHfeed_zero_indices[i+1] > pHfeed_zero_indices[i]]

pHfeed_cosine = get_cosine_d(pHfeed_data_chunks)

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'convergence_pHfeed.png')

plot_cosine_chunks(pHfeed_time_chunks, pHfeed_cosine, 88, fileName, '#FFD700')
    
plt.show()

# =========================
# Control conditions (ambr_8)
# =========================

control_ambr8_data = dataM[idx_control_ambr_8]
control_ambr8_time = np.array(time)[idx_control_ambr_8]

# Identify starting indices
control_ambr8_zero_indices = np.where(control_ambr8_time == 0)[0]
control_ambr8_zero_indices = np.append(control_ambr8_zero_indices, len(control_ambr8_time))

# Split into chunks
control_ambr8_time_chunks = [
    control_ambr8_time[control_ambr8_zero_indices[i]:control_ambr8_zero_indices[i+1]]
    for i in range(len(control_ambr8_zero_indices)-1) if control_ambr8_zero_indices[i+1] > control_ambr8_zero_indices[i]
]
control_ambr8_data_chunks = [
    control_ambr8_data[control_ambr8_zero_indices[i]:control_ambr8_zero_indices[i+1]]
    for i in range(len(control_ambr8_zero_indices)-1) if control_ambr8_zero_indices[i+1] > control_ambr8_zero_indices[i]
]

control_ambr8_cosine = get_cosine_d(control_ambr8_data_chunks)

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'convergence_control_ambr8.png')

plot_cosine_chunks(control_ambr8_time_chunks, control_ambr8_cosine, 94, fileName, '#00FF00')  # Blue color


# =========================
# pH/feed perturbation (ambr_8)
# =========================

pHfeed_ambr8_data = dataM[idx_pHfeed_ambr_8]
pHfeed_ambr8_time = np.array(time)[idx_pHfeed_ambr_8]

# Align zero point
pHfeed_ambr8_time -= 268  # Same shift as before

pHfeed_ambr8_zero_indices = np.where(pHfeed_ambr8_time == 0)[0]
pHfeed_ambr8_zero_indices = np.append(pHfeed_ambr8_zero_indices, len(pHfeed_ambr8_time))

# Split into chunks
pHfeed_ambr8_time_chunks = [
    pHfeed_ambr8_time[pHfeed_ambr8_zero_indices[i]:pHfeed_ambr8_zero_indices[i+1]]
    for i in range(len(pHfeed_ambr8_zero_indices)-1) if pHfeed_ambr8_zero_indices[i+1] > pHfeed_ambr8_zero_indices[i]
]
pHfeed_ambr8_data_chunks = [
    pHfeed_ambr8_data[pHfeed_ambr8_zero_indices[i]:pHfeed_ambr8_zero_indices[i+1]]
    for i in range(len(pHfeed_ambr8_zero_indices)-1) if pHfeed_ambr8_zero_indices[i+1] > pHfeed_ambr8_zero_indices[i]
]

pHfeed_ambr8_cosine = get_cosine_d(pHfeed_ambr8_data_chunks, w=1)

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'convergence_pHfeed_ambr8.png')

plot_cosine_chunks(pHfeed_ambr8_time_chunks, pHfeed_ambr8_cosine, 60, fileName, '#FFD700')  #