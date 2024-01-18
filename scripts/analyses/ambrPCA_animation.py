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
        labels = f.readline().strip().split('\t')[4::]
        experiments = []
        conditions = []
        exclude= [ ]
        data = {i:[] for i in labels if i not in exclude}
        
        for line in f:
            a = line.strip().split('\t')
            experiments.append(a[0])
            conditions.append(a[2])
            for i in range(len(labels)):
                if labels[i] not in exclude:
                    data[labels[i]].append(float(a[i+4]))
    
    dataM = np.array([np.array(data[i]) for i in labels if i not in exclude]).T
    return labels, experiments, conditions, dataM






def getComposition(compVec):
    
    binary_comp = np.argsort(compVec)
    
    # Convert array to a tuple of its elements
    binary_comp_tuple = tuple(binary_comp.tolist())

    composition_dict = {
        (0, 1, 2): 1,
        (0, 2, 1): 2,
        (1, 0, 2): 3,
        (1, 2, 0): 4,
        (2, 0, 1): 5,
        (2, 1, 0): 6,
        
    }
    
    return composition_dict.get(binary_comp_tuple, 0)




def getEllipse(points, scale_factor=2.5):
    """
    Plot the smallest enclosing ellipse around the given points.

    Arguments:
    points -- np.array of points
    """
    # Robustly estimate covariance
    robust_cov = MinCovDet().fit(points)
    center = robust_cov.location_
    covariance = robust_cov.covariance_

    # Eigenvalues and eigenvectors of the covariance matrix
    vals, vecs = np.linalg.eigh(covariance)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:,order]

    # Calculate the angle and dimensions of the ellipse
    angle = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * np.sqrt(vals) * scale_factor
    
    return center, width, height, angle


def minimum_bounding_ellipse(data):
    """
    Finds the minimum bounding ellipse for a set of points without using OpenCV.

    Parameters:
        data (numpy.ndarray): A 2D array of points.

    Returns:
        numpy.ndarray: The center of the ellipse.
        float: The semi-major axis of the ellipse.
        float: The semi-minor axis of the ellipse.
        double: The angle of rotation of the ellipse.
    """

    # Check if data is a 2D array
    if data.ndim != 2:
        raise ValueError("Data must be a 2D array")

    # Find the mean of the data
    mean = np.mean(data, axis=0)

    # Create a matrix of pairwise distances
    distances = np.zeros((data.shape[0], data.shape[0]))
    for i, row in enumerate(data):
        for j, other_row in enumerate(data):
            distances[i, j] = np.linalg.norm(row - other_row)

    # Form the covariance matrix
    covariance = np.zeros((2, 2))
    for i, row in enumerate(data):
        covariance += np.outer((row - mean), (row - mean))

    # Calculate the eigenvalues and eigenvectors of the covariance matrix
    eigenvalues, eigenvectors = np.linalg.eig(covariance)

    # Find the semi-major and semi-minor axes
    semi_major_axis = np.sqrt(eigenvalues[0])
    semi_minor_axis = np.sqrt(eigenvalues[1])

    # Find the angle of rotation
    angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))

    # Calculate the center of the ellipse
    ellipse_center = mean

    return ellipse_center, semi_major_axis, semi_minor_axis, angle


def min_volume_enclosing_ellipse(points, tolerance=0.01):
    """
    Find the minimum volume enclosing ellipse of a set of points.
    
    Arguments:
    points -- 2D numpy array of points.
    tolerance -- Tolerance for the stopping criterion.

    Returns:
    A tuple containing the center and the axes lengths of the ellipse.
    """
    N, d = points.shape
    Q = np.vstack([points.T, np.ones(N)])
    err = 1.0 + tolerance
    u = np.ones(N) / N

    while err > tolerance:
        X = np.dot(np.dot(Q, np.diag(u)), Q.T)
        M = np.diag(np.dot(np.dot(Q.T, np.linalg.inv(X)), Q))
        j = np.argmax(M)
        maximum = M[j]
        step_size = (maximum - d - 1) / ((d + 1) * (maximum - 1))
        new_u = (1 - step_size) * u
        new_u[j] += step_size
        err = np.linalg.norm(new_u - u)
        u = new_u

    center = np.dot(points.T, u)
    A = np.linalg.inv(np.dot(np.dot(points.T, np.diag(u)), points) - np.outer(center, center)) / d

    U, s, _ = np.linalg.svd(A)
    radii = 1 / np.sqrt(s)
    
    angle = np.arctan2(U[0, 1], U[0, 0]) * 180 / np.pi

    return center,  2*radii[0], 2*radii[1],angle


labels, experiments, conditions, dataM = get_data('ambrAll.txt')





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
initInd_c18b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C18' ]
initInd_c19b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C19' ]
initInd_c20b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C20' ]
initInd_c21b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C21' ]
initInd_c22b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C22' ]
initInd_c23b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C23' ]
initInd_c24b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C24' ]
initInd_c25b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C25' ]
initInd_c26b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C26' ]
initInd_c27b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C27' ]
initInd_c28b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C28' ]
initInd_c29b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C29' ]
initInd_c30b = [i for i in range(len(conditions)) if conditions[i] == 'E2_C30' ]
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


#Background PCA




scaler = StandardScaler()
scaled_data = scaler.fit_transform(dataM)


transformed_data = pca.fit_transform(scaled_data)

explained_variance = pca.explained_variance_ratio_

stateA = initInd_c10 + initInd_c11 + initInd_c12  +[i for i in range(len(conditions)) if ((conditions[i] == 'E1_P27') or (conditions[i] == 'E1_P26') or (conditions[i] == 'E1_P25')) and (experiments[i]=='A6_C2_V12') ]


stateB =  initInd_c17b + initInd_c18b + initInd_c19b #initInd_c29b + initInd_c30b + initInd_c31b 

stateC = initInd_p23b + initInd_p24b + initInd_p25b

samplesD = ['A6_C2_V7', 'A6_C2_V8', 'A6_C2_V9', 'A6_C2_V11']
stateD = [i for i in range(len(conditions)) if ((conditions[i] == 'E1_P25') or (conditions[i] == 'E1_P26') or (conditions[i] == 'E1_P27')) and (experiments[i] in samplesD) ]

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
scatter = ax.scatter(transformed_data.T[0], transformed_data.T[1], c='#39FF14', s=5, alpha = 0.5, linewidths=0.1, zorder=1)
highlight_scatter = ax.scatter([], [], color='#1F51FF', alpha=1.0, zorder=2, s=10)

plt.xlabel(f'PC1 ({explained_variance[0]*100:.2f} %)', fontsize=16)
plt.ylabel(f'PC2 ({explained_variance[1]*100:.2f} %)', fontsize=16)

c1, w1, h1, a1 = min_volume_enclosing_ellipse(transformed_data[stateA])
c2, w2, h2, a2 = min_volume_enclosing_ellipse(transformed_data[stateB])
c3, w3, h3, a3 = min_volume_enclosing_ellipse(transformed_data[stateC])
c4, w4, h4, a4 = min_volume_enclosing_ellipse(transformed_data[stateD])
ell1 = Ellipse(xy=c1,
              width=w1, height=h1,
              angle=a1, edgecolor='black', fc='None', lw=0.5)



ell2 = Ellipse(xy=c2,
              width=w2, height=h2,
              angle= a2, edgecolor='black', fc='None', lw=0.5)


ell3 = Ellipse(xy=c3,
              width=w3, height=h3,
              angle=a3, edgecolor='black', fc='None', lw=0.5)

ell4 = Ellipse(xy=c4,
              width=w4, height=h4,
              angle=a4, edgecolor='black', fc='None', lw=0.5)


ax.add_patch(ell1)
ax.add_patch(ell2)
ax.add_patch(ell3)
ax.add_patch(ell4)



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

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'ambrControlExperiment1.gif')
ani.save(fileName, writer='imagemagick')
plt.show()
# ######################################################################################################################

# ###################################################################################################################
plt.figure(figsize=(8, 6))
fig, ax = plt.subplots()
scatter = ax.scatter(transformed_data.T[0], transformed_data.T[1], c='#39FF14', s=5, alpha = 0.5, linewidths=0.1, zorder=1)
highlight_scatter = ax.scatter([], [], color='#bc13fe', alpha=1.0, zorder=2, s=10)
plt.xlabel(f'PC1 ({explained_variance[0]*100:.2f} %)', fontsize=16)
plt.ylabel(f'PC2 ({explained_variance[1]*100:.2f} %)', fontsize=16)

c1, w1, h1, a1 = min_volume_enclosing_ellipse(transformed_data[stateA])
c2, w2, h2, a2 = min_volume_enclosing_ellipse(transformed_data[stateB])
c3, w3, h3, a3 = min_volume_enclosing_ellipse(transformed_data[stateC])
c4, w4, h4, a4 = min_volume_enclosing_ellipse(transformed_data[stateD])
ell1 = Ellipse(xy=c1,
              width=w1, height=h1,
              angle=a1, edgecolor='black', fc='None', lw=0.5)



ell2 = Ellipse(xy=c2,
              width=w2, height=h2,
              angle= a2, edgecolor='black', fc='None', lw=0.5)


ell3 = Ellipse(xy=c3,
              width=w3, height=h3,
              angle=a3, edgecolor='black', fc='None', lw=0.5)

ell4 = Ellipse(xy=c4,
              width=w4, height=h4,
              angle=a4, edgecolor='black', fc='None', lw=0.5)


ax.add_patch(ell1)
ax.add_patch(ell2)
ax.add_patch(ell3)
ax.add_patch(ell4)



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
    initInd_c20b,
    initInd_c21b,
    initInd_c22b,
    initInd_c23b,
    initInd_c24b,
    initInd_c25b,
    initInd_c26b,
    initInd_c27b,
    initInd_c28b,
    initInd_c29b,
    initInd_c30b,
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
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'ambrControlExperiment2.gif')
ani.save(fileName, writer='imagemagick')
plt.show()
# ###################################################################################################################



plt.figure(figsize=(8, 6))
fig, ax = plt.subplots()
scatter = ax.scatter(transformed_data.T[0], transformed_data.T[1], c='#39FF14', s=5, alpha = 0.5, linewidths=0.1, zorder=1)
highlight_scatter = ax.scatter([], [], color='#1F51FF', alpha=1.0, zorder=2, s=10)
plt.xlabel(f'PC1 ({explained_variance[0]*100:.2f} %)', fontsize=16)
plt.ylabel(f'PC2 ({explained_variance[1]*100:.2f} %)', fontsize=16)

c1, w1, h1, a1 = min_volume_enclosing_ellipse(transformed_data[stateA])
c2, w2, h2, a2 = min_volume_enclosing_ellipse(transformed_data[stateB])
c3, w3, h3, a3 = min_volume_enclosing_ellipse(transformed_data[stateC])
c4, w4, h4, a4 = min_volume_enclosing_ellipse(transformed_data[stateD])
ell1 = Ellipse(xy=c1,
              width=w1, height=h1,
              angle=a1, edgecolor='black', fc='None', lw=0.5)



ell2 = Ellipse(xy=c2,
              width=w2, height=h2,
              angle= a2, edgecolor='black', fc='None', lw=0.5)


ell3 = Ellipse(xy=c3,
              width=w3, height=h3,
              angle=a3, edgecolor='black', fc='None', lw=0.5)

ell4 = Ellipse(xy=c4,
              width=w4, height=h4,
              angle=a4, edgecolor='black', fc='None', lw=0.5)


ax.add_patch(ell1)
ax.add_patch(ell2)
ax.add_patch(ell3)
ax.add_patch(ell4)



# Example arbitrary groups of indices
arbitrary_groups = [
    initInd_c12,
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

frame_titles = ['88 h',
                '96 h',
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

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'ambrPerturbExperiment1.gif')
ani.save(fileName, writer='imagemagick')
plt.show()


# ###################################################################################################################



plt.figure(figsize=(8, 6))
fig, ax = plt.subplots()
scatter = ax.scatter(transformed_data.T[0], transformed_data.T[1], c='#39FF14', s=5, alpha = 0.5, linewidths=0.1, zorder=1)
highlight_scatter = ax.scatter([], [], color='#bc13fe', alpha=1.0, zorder=2, s=10)
plt.xlabel(f'PC1 ({explained_variance[0]*100:.2f} %)', fontsize=16)
plt.ylabel(f'PC2 ({explained_variance[1]*100:.2f} %)', fontsize=16)

c1, w1, h1, a1 = min_volume_enclosing_ellipse(transformed_data[stateA])
c2, w2, h2, a2 = min_volume_enclosing_ellipse(transformed_data[stateB])
c3, w3, h3, a3 = min_volume_enclosing_ellipse(transformed_data[stateC])
c4, w4, h4, a4 = min_volume_enclosing_ellipse(transformed_data[stateD])
ell1 = Ellipse(xy=c1,
              width=w1, height=h1,
              angle=a1, edgecolor='black', fc='None', lw=0.5)



ell2 = Ellipse(xy=c2,
              width=w2, height=h2,
              angle= a2, edgecolor='black', fc='None', lw=0.5)


ell3 = Ellipse(xy=c3,
              width=w3, height=h3,
              angle=a3, edgecolor='black', fc='None', lw=0.5)

ell4 = Ellipse(xy=c4,
              width=w4, height=h4,
              angle=a4, edgecolor='black', fc='None', lw=0.5)


ax.add_patch(ell1)
ax.add_patch(ell2)
ax.add_patch(ell3)
ax.add_patch(ell4)



# Example arbitrary groups of indices
arbitrary_groups = [
    initInd_c19b,
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


frame_titles = ['108 h',
                '114 h',
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

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'ambrPerturbExperiment2.gif')
ani.save(fileName, writer='imagemagick')
plt.show()
