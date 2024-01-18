# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 11:09:27 2023

@author: drgarza
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Ellipse
import numpy as np


import scipy.spatial as sps

from aquarel import load_theme

theme = load_theme("boxy_light")
theme.apply()

#from sklearn.manifold import MDS as PCA
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import umap

PC = PCA(n_components=2)



def get_data(fileName):
    
    strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'ambr')
    
    
    
    with open(os.path.join(strainSummaryFolder, fileName)) as f:
        labels = f.readline().strip().split('\t')[4::]
        experiments = []
        conditions = []
        data = {i:[] for i in labels}
        
        for line in f:
            a = line.strip().split('\t')
            experiments.append(a[0])
            conditions.append(a[2])
            for i in range(len(labels)):
                data[labels[i]].append(float(a[i+4]))
    
    dataM = np.array([np.array(data[i]) for i in labels]).T
    return labels, experiments, conditions, dataM













labels, experiments, conditions, dataM = get_data('ambrAll.txt')




#Background PCA

scaler = StandardScaler()
scaled_data = scaler.fit_transform(dataM)
pca = PCA(n_components=2)

transformed_data = pca.fit_transform(scaled_data)


transformed_data = umap.UMAP(n_neighbors=50, n_components=2, metric = 'canberra', min_dist=0.5).fit(scaled_data)

transformed_data = transformed_data.embedding_[:]

explained_variance = [0,0]#pca.explained_variance_ratio_

#pinit

initInd_c1 = [i for i in range(len(conditions)) if conditions[i] == 'E1_C1']# or conditions[i] == 'E2_c1' ]
initInd_c2 = [i for i in range(len(conditions)) if conditions[i] == 'E1_C2']# or conditions[i] == 'E2_c2' ]
initInd_c3 = [i for i in range(len(conditions)) if conditions[i] == 'E1_C3']# or conditions[i] == 'E2_c2' ]
initInd_c4 = [i for i in range(len(conditions)) if conditions[i] == 'E1_C4']# or conditions[i] == 'E2_c2' ]
initInd_c5 = [i for i in range(len(conditions)) if conditions[i] == 'E1_C5']# or conditions[i] == 'E2_c2' ]
initInd_c6 = [i for i in range(len(conditions)) if conditions[i] == 'E1_C6']# or conditions[i] == 'E2_c2' ]
initInd_c7 = [i for i in range(len(conditions)) if conditions[i] == 'E1_C7']# or conditions[i] == 'E2_c2' ]
initInd_c8 = [i for i in range(len(conditions)) if conditions[i] == 'E1_C8']# or conditions[i] == 'E2_c2' ]
initInd_c9 = [i for i in range(len(conditions)) if conditions[i] == 'E1_C9']# or conditions[i] == 'E2_c2' ]
initInd_c10 = [i for i in range(len(conditions)) if conditions[i] == 'E1_C10']# or conditions[i] == 'E2_c2' ]
initInd_c11 = [i for i in range(len(conditions)) if conditions[i] == 'E1_C11']# or conditions[i] == 'E2_c2' ]
initInd_c12 = [i for i in range(len(conditions)) if conditions[i] == 'E1_C12']# or conditions[i] == 'E2_c2' ]

initInd_p1 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P1']# or conditions[i] == 'E2_c2' ]
initInd_p2 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P2']# or conditions[i] == 'E2_c2' ]
initInd_p3 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P3']# or conditions[i] == 'E2_c2' ]
initInd_p4 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P4']# or conditions[i] == 'E2_c2' ]
initInd_p5 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P5']# or conditions[i] == 'E2_c2' ]
initInd_p6 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P6']# or conditions[i] == 'E2_c2' ]
initInd_p7 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P7']# or conditions[i] == 'E2_c2' ]
initInd_p8 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P8']# or conditions[i] == 'E2_c2' ]
initInd_p9 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P9']# or conditions[i] == 'E2_c2' ]
initInd_p10 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P10']# or conditions[i] == 'E2_c2' ]
initInd_p11 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P11']# or conditions[i] == 'E2_c2' ]
initInd_p12 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P12']# or conditions[i] == 'E2_c2' ]
initInd_p13 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P13']# or conditions[i] == 'E2_c2' ]
initInd_p14 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P14']# or conditions[i] == 'E2_c2' ]
initInd_p15 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P15']# or conditions[i] == 'E2_c2' ]
initInd_p16 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P16']# or conditions[i] == 'E2_c2' ]
initInd_p17 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P17']# or conditions[i] == 'E2_c2' ]
initInd_p18 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P18']# or conditions[i] == 'E2_c2' ]
initInd_p19 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P19']# or conditions[i] == 'E2_c2' ]
initInd_p20 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P20']# or conditions[i] == 'E2_c2' ]
initInd_p21 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P21']# or conditions[i] == 'E2_c2' ]
initInd_p22 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P22']# or conditions[i] == 'E2_c2' ]
initInd_p23 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P23']# or conditions[i] == 'E2_c2' ]
initInd_p24 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P24']# or conditions[i] == 'E2_c2' ]
initInd_p25 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P25']# or conditions[i] == 'E2_c2' ]
initInd_p26 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P26']# or conditions[i] == 'E2_c2' ]
initInd_p27 = [i for i in range(len(conditions)) if conditions[i] == 'E1_P27']# or conditions[i] == 'E2_c2' ]


initInd_c1b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C1' ]
initInd_c2b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C2' ]
initInd_c3b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C3' ]
initInd_c4b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C4' ]
initInd_c5b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C5' ]
initInd_c6b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C6' ]
initInd_c7b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C7' ]
initInd_c8b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C8' ]
initInd_c9b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C9' ]
initInd_c10b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C10' ]
initInd_c11b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C11' ]
initInd_c12b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C12' ]
initInd_c13b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C13' ]
initInd_c14b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C14' ]
initInd_c15b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C15' ]
initInd_c16b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C16' ]
initInd_c17b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C17' ]
initInd_c18b = [i for i in range(len(conditions)) if conditions[i] == 'E2_c18' ]
initInd_c19b = [i for i in range(len(conditions)) if conditions[i] == 'E2_c19' ]
initInd_c21b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C21' ]
initInd_c23b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C23' ]
initInd_c25b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C25' ]
initInd_c27b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C27' ]
initInd_c29b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C29' ]
initInd_c31b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C31' ]



initInd_p1b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P1' ]
initInd_p2b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P2' ]
initInd_p3b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P3' ]
initInd_p4b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P4' ]
initInd_p5b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P5' ]
initInd_p6b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P6' ]
initInd_p7b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P7' ]
initInd_p8b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P8' ]
initInd_p9b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P9' ]
initInd_p10b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P10' ]
initInd_p11b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P11' ]
initInd_p12b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P12' ]
initInd_p13b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P13' ]
initInd_p14b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P14' ]
initInd_p15b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P15' ]
initInd_p16b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P16' ]
initInd_p17b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P17' ]
initInd_p18b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P18' ]
initInd_p19b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P19' ]
initInd_p20b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P20' ]
initInd_p21b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P21' ]
initInd_p22b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P22' ]
initInd_p23b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P23' ]
initInd_p24b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P24' ]
initInd_p25b = [i for i in range(len(conditions)) if conditions[i] == 'E2_P25' ]




##############################


def update(frame):
    # Extract the x and y coordinates of the current group
    x_coords = transformed_data[arbitrary_groups[frame], 0]
    y_coords = transformed_data[arbitrary_groups[frame], 1]

    # Update the data of the highlight scatter plot
    highlight_scatter.set_offsets(np.c_[x_coords, y_coords])
    
    ax.set_title(frame_titles[frame])
    
    plt.tight_layout()

    return scatter, highlight_scatter



#Path to non-perturbed state
#####################################################################################################################
plt.figure(figsize=(8, 6))
fig, ax = plt.subplots()
scatter = ax.scatter(transformed_data.T[0], transformed_data.T[1], c='#39FF14', s=100, alpha = 0.1, linewidths=0.1, zorder=1)
highlight_scatter = ax.scatter([], [], color='#1F51FF', alpha=1.0, zorder=2)

plt.xlabel(f'UMAP1', fontsize=16)
plt.ylabel(f'UMAP2', fontsize=16)

# ell1 = Ellipse(xy=(0.4, 0.1),
#               width=2.1, height=3.6,
#               angle=-20, edgecolor='black', fc='None', lw=0.5)



# ell2 = Ellipse(xy=(-1.8, -0.8),
#               width=3.5, height=1.5,
#               angle= 80, edgecolor='black', fc='None', lw=0.5)


# ell3 = Ellipse(xy=(3.5, 2.0),
#               width=3.5, height=1.5,
#               angle=np.degrees(30), edgecolor='black', fc='None', lw=0.5)


# ax.add_patch(ell1)
# ax.add_patch(ell2)
# ax.add_patch(ell3)



# Example arbitrary groups of indices
arbitrary_groups = [
    initInd_c1,
    initInd_c2,
    initInd_c3,
    initInd_c4,
    initInd_c5,
    initInd_c6,
    initInd_c7,
    initInd_c8,
    initInd_c9,
    initInd_c10,
    initInd_c11,
    initInd_c12
    
]

frame_titles = ['0 h',
                '8 h',
                '16 h',
                '24 h',
                '32 h',
                '40 h',
                '48 h',
                '56 h',
                '64 h',
                '72 h',
                '80 h',
                '88 h',]


ani = FuncAnimation(fig, update, frames=len(arbitrary_groups), interval=1000, blit=True)

plt.tight_layout()

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'ambrControlExperiment1_bac.gif')
ani.save(fileName, writer='imagemagick')
plt.show()
# ######################################################################################################################

# ###################################################################################################################
plt.figure(figsize=(8, 6))
fig, ax = plt.subplots()
scatter = ax.scatter(transformed_data.T[0], transformed_data.T[1], c='#39FF14', s=100, alpha = 0.1, linewidths=0.1)
highlight_scatter = ax.scatter([], [], color='#bc13fe', alpha=1.0, zorder=2)
plt.xlabel(f'UMAP1', fontsize=16)
plt.ylabel(f'UMAP2', fontsize=16)

# ell1 = Ellipse(xy=(0.4, 0.1),
#               width=2.1, height=3.6,
#               angle=-20, edgecolor='black', fc='None', lw=0.5)



# ell2 = Ellipse(xy=(-1.8, -0.8),
#               width=3.5, height=1.5,
#               angle= 80, edgecolor='black', fc='None', lw=0.5)


# ell3 = Ellipse(xy=(3.5, 2.0),
#               width=3.5, height=1.5,
#               angle=np.degrees(30), edgecolor='black', fc='None', lw=0.5)


# ax.add_patch(ell1)
# ax.add_patch(ell2)
# ax.add_patch(ell3)



# Example arbitrary groups of indices
arbitrary_groups = [
    initInd_c1b,
    initInd_c2b,
    initInd_c3b,
    initInd_c4b,
    initInd_c5b,
    initInd_c6b,
    initInd_c7b,
    initInd_c8b,
    initInd_c9b,
    initInd_c10b,
    initInd_c11b,
    initInd_c12b,
    initInd_c13b,
    initInd_c14b,
    initInd_c15b,
    initInd_c16b,
    initInd_c17b,
    initInd_c18b,
    initInd_c19b,
    initInd_c21b,
    initInd_c23b,
    initInd_c25b,
    initInd_c27b,
    initInd_c29b,
    initInd_c31b,
    

    
]



frame_titles = ['0 h',
                '6 h',
                '12 h',
                '18 h',
                '24 h',
                '30 h',
                '36 h',
                '42 h',
                '48 h',
                '54 h',
                '60 h',
                '66 h',
                '72 h',
                '78 h',
                '84 h',
                '90 h',
                '96 h',
                '102 h',
                '108 h',
                '120 h',
                '132 h',
                '144 h',
                '156 h',
                '168 h',
                '180 h',
                '192 h',
                '204 h',
                '216 h',
                '228 h',
                '240 h',
                '252 h',
                
                ]


ani = FuncAnimation(fig, update, frames=len(arbitrary_groups), interval=1000, blit=True)

plt.tight_layout()
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'ambrControlExperiment2_bac.gif')
ani.save(fileName, writer='imagemagick')
plt.show()
# ###################################################################################################################



plt.figure(figsize=(8, 6))
fig, ax = plt.subplots()
scatter = ax.scatter(transformed_data.T[0], transformed_data.T[1], c='#39FF14', s=100, alpha = 0.1, linewidths=0.1)
highlight_scatter = ax.scatter([], [], color='#1F51FF', alpha=1.0, zorder=2)
plt.xlabel(f'UMAP1', fontsize=16)
plt.ylabel(f'UMAP2', fontsize=16)

# ell1 = Ellipse(xy=(0.4, 0.1),
#               width=2.1, height=3.6,
#               angle=-20, edgecolor='black', fc='None', lw=0.5)



# ell2 = Ellipse(xy=(-1.8, -0.8),
#               width=3.5, height=1.5,
#               angle= 80, edgecolor='black', fc='None', lw=0.5)


# ell3 = Ellipse(xy=(3.5, 2.0),
#               width=3.5, height=1.5,
#               angle=np.degrees(30), edgecolor='black', fc='None', lw=0.5)


# ax.add_patch(ell1)
# ax.add_patch(ell2)
# ax.add_patch(ell3)



# Example arbitrary groups of indices
arbitrary_groups = [
    initInd_p1,
    initInd_p2,
    initInd_p3,
    initInd_p4,
    initInd_p5,
    initInd_p6,
    initInd_p7,
    initInd_p8,
    initInd_p9,
    initInd_p10,
    initInd_p11,
    initInd_p12,
    initInd_p13,
    initInd_p14,
    initInd_p15,
    initInd_p16,
    initInd_p17,
    initInd_p18,
    initInd_p19,
    initInd_p20,
    initInd_p21,
    initInd_p22,
    initInd_p23,
    initInd_p24,
    initInd_p25,
    initInd_p26,
    initInd_p27

    
]

frame_titles = ['96 h',
                '104 h',
                '112 h',
                '120 h',
                '128 h',
                '136 h',
                '144 h',
                '152 h',
                '160 h',
                '168 h',
                '176 h',
                '184 h',
                '192 h',
                '200 h',
                '208 h',
                '216 h',
                '224 h',
                '232 h',
                '240 h',
                '248 h',
                '256 h',
                '264 h',
                '272 h',
                '280 h',
                '288 h',
                '296 h',
                '308 h',
                
                ]


ani = FuncAnimation(fig, update, frames=len(arbitrary_groups), interval=1000, blit=True)
plt.tight_layout()

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'ambrPerturbExperiment1_bac.gif')
ani.save(fileName, writer='imagemagick')
plt.show()


# ###################################################################################################################



plt.figure(figsize=(8, 6))
fig, ax = plt.subplots()
scatter = ax.scatter(transformed_data.T[0], transformed_data.T[1], c='#39FF14', s=100, alpha = 0.1, linewidths=0.1)
highlight_scatter = ax.scatter([], [], color='#bc13fe', alpha=1.0, zorder=2)
plt.xlabel(f'UMAP1', fontsize=16)
plt.ylabel(f'UMAP2', fontsize=16)

# ell1 = Ellipse(xy=(0.4, 0.1),
#               width=2.1, height=3.6,
#               angle=-20, edgecolor='black', fc='None', lw=0.5)



# ell2 = Ellipse(xy=(-1.8, -0.8),
#               width=3.5, height=1.5,
#               angle= 80, edgecolor='black', fc='None', lw=0.5)


# ell3 = Ellipse(xy=(3.5, 2.0),
#               width=3.5, height=1.5,
#               angle=np.degrees(30), edgecolor='black', fc='None', lw=0.5)


# ax.add_patch(ell1)
# ax.add_patch(ell2)
# ax.add_patch(ell3)



# Example arbitrary groups of indices
arbitrary_groups = [
    initInd_p1b,
    initInd_p2b,
    initInd_p3b,
    initInd_p4b,
    initInd_p5b,
    initInd_p6b,
    initInd_p7b,
    initInd_p8b,
    initInd_p9b,
    initInd_p10b,
    initInd_p11b,
    initInd_p12b,
    initInd_p13b,
    initInd_p14b,
    initInd_p15b,
    initInd_p16b,
    initInd_p17b,
    initInd_p18b,
    initInd_p19b,
    initInd_p20b,
    initInd_p21b,
    initInd_p22b,
    initInd_p23b,
    initInd_p24b,
    initInd_p25b

    
]


frame_titles = ['114 h',
                '120 h',
                '126 h',
                '132 h',
                '138 h',
                '144 h',
                '150 h',
                '156 h',
                '162 h',
                '168 h',
                '174 h',
                '180 h',
                '186 h',
                '192 h',
                '198 h',
                '204 h',
                '210 h',
                '216 h',
                '222 h',
                '228 h',
                '234 h',
                '240 h',
                '246 h',
                '252 h',
                '258 h'
                ]

ani = FuncAnimation(fig, update, frames=len(arbitrary_groups), interval=1000, blit=True)

plt.tight_layout()

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'ambrPerturbExperiment2_bac.gif')
ani.save(fileName, writer='imagemagick')
plt.show()
