# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 11:10:47 2025
@author: drgarza
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import braycurtis
from aquarel import load_theme

theme = load_theme("boxy_light")
theme.apply()


# ===========================================
# Load Data
# ===========================================
def get_data(fileName):
    strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'ambr')

    with open(os.path.join(strainSummaryFolder, fileName)) as f:
        labels = f.readline().strip().split('\t')[6::]
        experiments, conditions, description, time, timeCondition, replicate = [], [], [], [], [], []
        exclude = []
        data = {i: [] for i in labels if i not in exclude}

        for line in f:
            a = line.strip().split('\t')
            experiments.append(a[0])
            conditions.append(a[3])
            description.append(a[2])

            for i in range(len(labels)):
                if labels[i] not in exclude:
                    data[labels[i]].append(float(a[i + 6]))
            time.append(float(a[4]))
            timeCondition.append(float(a[5]))

            ambrExp = a[0].split('_')
            replicate.append(ambrExp[0] + '_' + ambrExp[1] + '_' + a[4])

    dataM = np.array([np.array(data[i]) for i in labels if i not in exclude]).T
    return labels, np.array(experiments), np.array(description), conditions, time, np.array(timeCondition), replicate, dataM


# ===========================================
# Compute Bray-Curtis Distances
# ===========================================
def get_cosine_d(chunk_array, w=2):
    distances = []
    for i in range(len(chunk_array)):
        nd = chunk_array[i] / np.max(chunk_array[i], axis=1, keepdims=True)
        cosine_distances = np.array([
            braycurtis(nd[z], nd[z + w])
            for z in range(len(nd) - w)
        ])
        distances.append(cosine_distances[:])
    return distances


# ===========================================
# Plot function
# ===========================================
def plot_cosine_chunks(time_chunks, cosine_chunks, highlight=100, fileName=None, ball_color="tab:blue", y_label = 'Bray-Curtis dissimilarity (t vs t+2)'):
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
    plt.ylabel(y_label, fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()

    if fileName is not None:
        plt.savefig(fileName, dpi=600, transparent=True)
    plt.show()


# ===========================================
# Improved chunking logic
# ===========================================
def split_chunks(time_array, data_array):
    time_array = np.array(time_array)
    data_array = np.array(data_array)

    # Detect replicate breaks when time decreases
    breaks = [0]
    for i in range(1, len(time_array)):
        if time_array[i] < time_array[i - 1]:
            breaks.append(i)
    breaks.append(len(time_array))

    time_chunks = []
    data_chunks = []

    for i in range(len(breaks) - 1):
        start, end = breaks[i], breaks[i + 1]
        t_chunk = time_array[start:end]
        d_chunk = data_array[start:end]

        # Sort and remove duplicates
        sort_idx = np.argsort(t_chunk)
        t_chunk = t_chunk[sort_idx]
        d_chunk = d_chunk[sort_idx]
        unique_idx = np.unique(t_chunk, return_index=True)[1]
        t_chunk = t_chunk[unique_idx]
        d_chunk = d_chunk[unique_idx]

        time_chunks.append(t_chunk)
        data_chunks.append(d_chunk)

    return time_chunks, data_chunks


# ===========================================
# Main Workflow
# ===========================================
labels, experiments, description, conditions, time, timeCond, replicate, dataM = get_data('ambrAll.txt')

ambr_8 = np.array([True if i[0:3] == 'A8_' else False for i in experiments])
ambr_6 = np.array([True if i[0:3] == 'A6_' else False for i in experiments])
not_ambr_8 = ~ambr_8
not_ambr_6 = ~ambr_6

idx_control = (description == 'control conditions') & not_ambr_8
idx_control_ambr_8 = (description == 'control conditions') & ambr_8
idx_feed = (description == 'post feed perturbation')
idx_pHfeed = (description == 'post pH/feed perturbation') & not_ambr_8
idx_pHfeed_ambr_8 = (description == 'post pH/feed perturbation') & ambr_8


# =========================
# Control (not ambr_8)
# =========================
control_data = dataM[idx_control]
control_time = np.array(time)[idx_control]
control_time_chunks, control_data_chunks = split_chunks(control_time, control_data)

control_cosine = get_cosine_d(control_data_chunks, w=2)
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'convergence_control.png')
plot_cosine_chunks(control_time_chunks, control_cosine, 84, fileName, '#00FF00')


# =========================
# Feed perturbation
# =========================
feed_data = dataM[idx_feed]
feed_time = np.array(time)[idx_feed] - 132
feed_time_chunks, feed_data_chunks = split_chunks(feed_time, feed_data)

feed_cosine = get_cosine_d(feed_data_chunks, w=2)
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'convergence_feed.png')
plot_cosine_chunks(feed_time_chunks, feed_cosine, 78, fileName, '#FF007F')


# =========================
# pH/feed perturbation (not ambr_8)
# =========================
pHfeed_data = dataM[idx_pHfeed]
pHfeed_time = np.array(time)[idx_pHfeed] - 176
pHfeed_time_chunks, pHfeed_data_chunks = split_chunks(pHfeed_time, pHfeed_data)

pHfeed_cosine = get_cosine_d(pHfeed_data_chunks, w=2)
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'convergence_pHfeed.png')
plot_cosine_chunks(pHfeed_time_chunks, pHfeed_cosine, 88, fileName, '#FFD700')


# =========================
# Control (ambr_8)
# =========================
control_ambr8_data = dataM[idx_control_ambr_8]
control_ambr8_time = np.array(time)[idx_control_ambr_8]
control_ambr8_time_chunks, control_ambr8_data_chunks = split_chunks(control_ambr8_time, control_ambr8_data)

control_ambr8_cosine = get_cosine_d(control_ambr8_data_chunks, w=1)  # ✅ w=1 for sparse
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'convergence_control_ambr8.png')
plot_cosine_chunks(control_ambr8_time_chunks, control_ambr8_cosine, 94, fileName, '#00FF00')


# =========================
# pH/feed perturbation (ambr_8)
# =========================
pHfeed_ambr8_data = dataM[idx_pHfeed_ambr_8]
pHfeed_ambr8_time = np.array(time)[idx_pHfeed_ambr_8] - 268
pHfeed_ambr8_time_chunks, pHfeed_ambr8_data_chunks = split_chunks(pHfeed_ambr8_time, pHfeed_ambr8_data)

pHfeed_ambr8_cosine = get_cosine_d(pHfeed_ambr8_data_chunks, w=1)
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'convergence_pHfeed_ambr8.png')
plot_cosine_chunks(pHfeed_ambr8_time_chunks, pHfeed_ambr8_cosine, 60, fileName, '#FFD700', y_label = 'Bray-Curtis dissimilarity (t vs t+1)')

plt.show()
