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

#from sklearn.manifold import MDS as PCA
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
PC = PCA(n_components=2)


import matplotlib.pyplot as plt

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))
from general import *


def create_heatmap(data, row_labels, x_positions, x_values, xlabel, ylabel, title, fileName = None):
    fig, ax = plt.subplots()
    heatmap = ax.imshow(data, aspect='auto')
    ax.grid(False)

    # Set row labels (y-axis)
    ax.set_yticks(np.arange(data.shape[0]))
    
    yticklabels = ax.set_yticklabels(row_labels)
    yticklabels[-1].set_fontstyle('italic')
    yticklabels[-2].set_fontstyle('italic')
    yticklabels[-3].set_fontstyle('italic')

    # Set specific x-ticks and labels (x-axis)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_values, rotation =90)

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




def plot_pca_scatter(data, point_labels=None, labels=None, colors=None, title=None, xlabel='PC1', ylabel='PC2', jitter_amount=0.1, label_offset=0.02, point_size=100, fileName=None):
    """
    Plots a scatter plot of the first two principal components of the given data with jitter and point labels.

    :param data: 2D numpy array or list of lists with the data.
    :param point_labels: Optional; List of short labels for each data point.
    :param labels: Optional; List of labels for each data point.
    :param colors: Optional; Color for each data point or a single color for all.
    :param title: Optional; Title of the plot.
    :param xlabel: Optional; Label for the x-axis.
    :param ylabel: Optional; Label for the y-axis.
    :param jitter_amount: Amount of jitter to apply (default 0.1).
    :param label_offset: Offset for the point labels in the x-direction.
    :param point_size: Size of the points in the scatter plot.
    :param fileName: Optional; If provided, the plot will be saved to this file.
    """
    # Scale the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data)

    # Perform PCA
    pca = PCA(n_components=2)
    transformed_data = pca.fit_transform(scaled_data)

    # Add jitter
    jittered_data = transformed_data + np.random.normal(0, jitter_amount, transformed_data.shape)

    # Explained variance
    explained_variance = pca.explained_variance_ratio_
    print(f"Explained Variance: PC1: {explained_variance[0]:.2f}, PC2: {explained_variance[1]:.2f}")

    # Plotting
    plt.figure(figsize=(8, 6))
    if labels is not None:
        for i, label in enumerate(set(labels)):
            subset = jittered_data[labels == label]
            plt.scatter(subset[:, 0], subset[:, 1], s=point_size, 
                        label=label, c=(np.array(colors)[labels == label] if colors is not None else None))
        plt.legend()
    else:
        plt.scatter(jittered_data[:, 0], jittered_data[:, 1], color=colors, s=point_size)

    # Overlay point labels with an offset
    if point_labels is not None:
        for (x, y), label in zip(jittered_data, point_labels):
            plt.text(x + label_offset, y, label, fontsize=14)

    # Labels and title
    plt.xlabel(f'{xlabel} ({explained_variance[0]*100:.2f} %)', fontsize=16)
    plt.ylabel(f'{ylabel} ({explained_variance[1]*100:.2f} %)', fontsize=16)
    plt.title(title)

    # Adjust plot margins
    plt.xlim(min(jittered_data[:, 0])-1, max(jittered_data[:, 0])+1)
    plt.ylim(min(jittered_data[:, 1])-1, max(jittered_data[:, 1])+1)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

    plt.grid(True)
    plt.tight_layout()
    if fileName is not None:
        plt.savefig(fileName, transparent=True, dpi=600)
    plt.show()












strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'ambr')



with open(os.path.join(strainSummaryFolder, 'ambrAll.txt')) as f:
    labels = f.readline().strip().split('\t')[1:]
    experiments = []
    data = {i:[] for i in labels}
    
    for line in f:
        a = line.strip().split('\t')
        experiments.append(a[0])
        for i in range(len(labels)):
            data[labels[i]].append(float(a[i+1]))


dataM = np.array([np.array(data[i]) for i in labels]).T

pca = PC.fit_transform(dataM)

create_heatmap(dataM.T, labels, np.arange(len(experiments)), experiments, None, None, None, fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambrHM.png'))

plot_pca_scatter(dataM, point_labels=experiments, jitter_amount=0.4, fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambrPC.png'))


