# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 04:36:03 2025

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
import ternary

theme = load_theme("boxy_light")
theme.apply()

#from sklearn.manifold import MDS as PCA
from sklearn.covariance import MinCovDet
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
pca = PCA(n_components=2, svd_solver ='randomized', power_iteration_normalizer='QR', iterated_power=10000)
#%matplotlib qt

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




labels, experiments, description, conditions, time, timeCond, replicate, dataM = get_data('ambrAll.txt')

ambr_8 = np.array([True if i[0:3]=='A8_' else False for i in experiments])
not_ambr_8 = ~ambr_8

scaler = StandardScaler()
scaled_data = scaler.fit_transform(dataM)
transformed_data = pca.fit_transform(scaled_data)
explained_variance = pca.explained_variance_ratio_

idx_control = (description == 'control conditions') & not_ambr_8
idx_control_ambr_8 = (description == 'control conditions') & ambr_8

feed = description == 'post feed perturbation'

pHfeed = (description == 'post pH/feed perturbation') & not_ambr_8
pHfeed_ambr_8 = (description == 'post pH/feed perturbation') & ambr_8

Ri = dataM.T[-1]
Bt = dataM.T[-2]
Bh = dataM.T[-3]

x_pc = transformed_data[:, 0]
y_pc = transformed_data[:, 1]


steadstate_control = (description == 'control conditions') & (not_ambr_8) & (timeCond>99) 

control_composition = np.array([Ri[steadstate_control],
                                Bt[steadstate_control],
                                Bh[steadstate_control]]).T

control_composition = np.array([i/sum(i) for i in control_composition])

steadstate_control_ambr_8 = (description == 'control conditions') & (ambr_8) & (timeCond>99) 

control_composition_ambr_8 = np.array([Ri[steadstate_control_ambr_8],
                                Bt[steadstate_control_ambr_8],
                                Bh[steadstate_control_ambr_8]]).T

control_composition_ambr_8 = np.array([i/sum(i) for i in control_composition_ambr_8])


steadstate_feed = (description == 'post feed perturbation') & (not_ambr_8) & (timeCond>85) 

feed_composition = np.array([Ri[steadstate_feed],
                                Bt[steadstate_feed],
                                Bh[steadstate_feed]]).T

feed_composition = np.array([i/sum(i) for i in feed_composition])

steadstate_pHfeed = (description == 'post pH/feed perturbation') & (not_ambr_8) & (timeCond>99) 

pHfeed_composition = np.array([Ri[steadstate_pHfeed],
                                Bt[steadstate_pHfeed],
                                Bh[steadstate_pHfeed]]).T

pHfeed_composition = np.array([i/sum(i) for i in pHfeed_composition])

steadstate_pHfeed_ambr_8 = (description == 'post pH/feed perturbation') & (ambr_8) & (timeCond>64) 

pHfeed_composition_ambr_8 = np.array([Ri[steadstate_pHfeed_ambr_8],
                                Bt[steadstate_pHfeed_ambr_8],
                                Bh[steadstate_pHfeed_ambr_8]]).T

pHfeed_composition_ambr_8 = np.array([i/sum(i) for i in pHfeed_composition_ambr_8])


#Ternary plot
#####################################################################

scale = 1.0  # since data is normalized to 1
fig, tax = ternary.figure(scale=scale)


# Remove the outer rectangular frame
ax = tax.get_axes()
ax.set_frame_on(False)
for spine in ax.spines.values():
    spine.set_visible(False)
tax.clear_matplotlib_ticks()  # Removes ticks outside the ternary plot


tax.boundary(linewidth=2.0)
tax.gridlines(multiple=0.25, color="k")

# Scatter plot the data
tax.scatter(control_composition_ambr_8,
            s=250,
            marker='o', 
            color='#00FF00', 
            label="control",
            edgecolors='k',
            linewidths=0.5)


# # Label axes
# tax.bottom_axis_label("R. intestinalis", fontsize=20)
# tax.left_axis_label("B. thetaiotaomicron", fontsize=20)
# tax.right_axis_label("B. hydrogenotrophica", fontsize=20)

# Draw and show the plot
#tax.legend(fontsize=25)
tax.ticks(axis='lbr', multiple=0.25, linewidth=1, tick_formats="%.1f", fontsize=0)
plt.tight_layout()

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'ternary_control_ambr_8.png')

plt.savefig(fileName, dpi=600, transparent = True)
plt.show()
############################################

#Ternary plot
#####################################################################

scale = 1.0  # since data is normalized to 1
fig, tax = ternary.figure(scale=scale)


# Remove the outer rectangular frame
ax = tax.get_axes()
ax.set_frame_on(False)
for spine in ax.spines.values():
    spine.set_visible(False)
tax.clear_matplotlib_ticks()  # Removes ticks outside the ternary plot


tax.boundary(linewidth=2.0)
tax.gridlines(multiple=0.25, color="k")

# Scatter plot the data
tax.scatter(pHfeed_composition_ambr_8,
            s=250,
            marker = 's', 
            color='#FFD700', 
            label="pH + feed pert.",
            edgecolors='k',
            linewidths=0.5)

# # Label axes
# tax.bottom_axis_label("R. intestinalis", fontsize=20)
# tax.left_axis_label("B. thetaiotaomicron", fontsize=20)
# tax.right_axis_label("B. hydrogenotrophica", fontsize=20)

# Draw and show the plot
#tax.legend(fontsize=25)
tax.ticks(axis='lbr', multiple=0.25, linewidth=1, tick_formats="%.1f", fontsize=0)
plt.tight_layout()

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'ternary_pHfeed_ambr_8.png')

plt.savefig(fileName, dpi=600, transparent = True)
plt.show()
############################################



###############pca############################

fig, ax = plt.subplots(figsize=(8, 6), facecolor='none')
ax.set_facecolor('#dddddd')
background = ax.scatter(x_pc,
                        y_pc,
                        c='white', 
                        s =500,
                        alpha=0.01)
control_time_ambr_8 = ax.scatter(x_pc[idx_control_ambr_8], 
                          y_pc[idx_control_ambr_8],
                          marker = 'o',
                          c = timeCond[idx_control_ambr_8],
                          cmap = 'Purples',
                          s = 500,
                          vmin = 0,
                          vmax = 120)

control_data_ambr_8 = ax.scatter(x_pc[idx_control_ambr_8], 
                          y_pc[idx_control_ambr_8],
                          edgecolors='k',
                          linewidths=0.5,
                          c='#00FF00', 
                          marker = 'o',
                          s = 150, zorder =100)


# Labels
ax.set_xlabel(f"PC1 ({explained_variance[0]*100:.2f}% variance explained)", fontsize=20)
ax.set_ylabel(f"PC2 ({explained_variance[1]*100:.2f}% variance explained)", fontsize=20)

plt.tight_layout()
# Ensure the figure background remains transparent
fig.patch.set_alpha(0)

# Save the figure
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'pca_control_ambr_8.png')

plt.savefig(fileName, dpi=600, transparent = False)
plt.show()


fig, ax = plt.subplots(figsize=(8, 6))
ax.set_facecolor('#dddddd')

background = ax.scatter(x_pc,
                        y_pc,
                        c='white', 
                        s =100,
                        alpha=0.01)


control_time = ax.scatter(x_pc[idx_control_ambr_8], 
                          y_pc[idx_control_ambr_8],
                          marker = 'o',
                          c = timeCond[idx_control_ambr_8],
                          cmap = 'Purples',
                          s = 500,
                          vmin = 0,
                          vmax = 120)

control_data = ax.scatter(x_pc[idx_control_ambr_8], 
                          y_pc[idx_control_ambr_8],
                          edgecolors='k',
                          linewidths=0.5,
                          marker = 'o',
                          c='#00FF00', 
                          s = 150, zorder =10)

pHfeed_time = ax.scatter(x_pc[pHfeed_ambr_8], 
                          y_pc[pHfeed_ambr_8],
                          c = timeCond[pHfeed_ambr_8],
                          marker = 'o',
                          cmap = 'Purples',
                          s = 500,
                          vmin = 0,
                          vmax = 120,zorder=20)

pHfeed_data = ax.scatter(x_pc[pHfeed_ambr_8], 
                          y_pc[pHfeed_ambr_8],
                          edgecolors='k',
                          linewidths=0.5,
                          marker = 's',
                          c='#FFD700', 
                          s = 150, zorder =30)

# Labels
ax.set_xlabel(f"PC1 ({explained_variance[0]*100:.2f}% variance explained)", fontsize=20)
ax.set_ylabel(f"PC2 ({explained_variance[1]*100:.2f}% variance explained)", fontsize=20)

plt.tight_layout()
fig.patch.set_alpha(0)
# Save the figure
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'pca_control_pHfeed_ambr_8.png')

plt.savefig(fileName, dpi=600, transparent = False)
plt.show()
