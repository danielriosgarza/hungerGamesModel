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

scaler = StandardScaler()
scaled_data = scaler.fit_transform(dataM)


transformed_data = pca.fit_transform(scaled_data)

explained_variance = pca.explained_variance_ratio_


stateA = initInd_c10 + initInd_c11 + initInd_c12  +[i for i in range(len(conditions)) if ((conditions[i] == 'E1_P27') or (conditions[i] == 'E1_P26') or (conditions[i] == 'E1_P25')) and (experiments[i]=='A6_C2_V12') ]


stateB =  initInd_c17b + initInd_c18b + initInd_c19b #initInd_c29b + initInd_c30b + initInd_c31b 

samplesC = ['A7_C2_V7', 'A7_C2_V8', 'A7_C2_V9', 'A7_C2_V10', 'A7_C2_V11', 'A7_C2_V12']
stateC = [i for i in range(len(conditions)) if ((conditions[i] == 'E2_P23') or (conditions[i] == 'E2_P24') or (conditions[i] == 'E2_P25')) and (experiments[i] in samplesC) ]

samplesD = ['A6_C2_V7', 'A6_C2_V8', 'A6_C2_V9', 'A6_C2_V11']
stateD = [i for i in range(len(conditions)) if ((conditions[i] == 'E1_P25') or (conditions[i] == 'E1_P26') or (conditions[i] == 'E1_P27')) and (experiments[i] in samplesD) ]




#Path to non-perturbed state
#####################################################################################################################
plt.figure(figsize=(8, 6))
fig, ax = plt.subplots()
scatter = ax.scatter(transformed_data.T[0], transformed_data.T[1], c='#39FF14', s=5, alpha = 0.5, linewidths=0.1, zorder=1)
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
plt.tight_layout()

ax.scatter(transformed_data.T[0][initInd_c1], transformed_data.T[1][initInd_c1], c='#1F51FF', s=20, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c2], transformed_data.T[1][initInd_c2], c='#1F51FF', s=20, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c3], transformed_data.T[1][initInd_c3], c='#1F51FF', s=20, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c4], transformed_data.T[1][initInd_c4], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c5], transformed_data.T[1][initInd_c5], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c6], transformed_data.T[1][initInd_c6], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c7], transformed_data.T[1][initInd_c7], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c8], transformed_data.T[1][initInd_c8], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c9], transformed_data.T[1][initInd_c9], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c10], transformed_data.T[1][initInd_c10], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c11], transformed_data.T[1][initInd_c11], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c12], transformed_data.T[1][initInd_c12], c='#1F51FF', s=10, alpha = 1.)


fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'ambrE1_C0.png')
plt.savefig(fileName, transparent=True, dpi=600)
plt.show()
######################################################################################################################


#Path after perturbation
###################################################################################################################
plt.figure(figsize=(8, 6))
fig, ax = plt.subplots()
scatter = ax.scatter(transformed_data.T[0], transformed_data.T[1], c='#39FF14', s=5, alpha = 0.5, linewidths=0.1, zorder=1)
plt.xlabel(f'PC1 ({explained_variance[0]*100:.2f} %)', fontsize=16)
plt.ylabel(f'PC2 ({explained_variance[1]*100:.2f} %)', fontsize=16)

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
plt.tight_layout()

ax.scatter(transformed_data.T[0][initInd_c12], transformed_data.T[1][initInd_c12], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p1], transformed_data.T[1][initInd_p1], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p2], transformed_data.T[1][initInd_p2], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p3], transformed_data.T[1][initInd_p3], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p4], transformed_data.T[1][initInd_p4], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p5], transformed_data.T[1][initInd_p5], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p6], transformed_data.T[1][initInd_p6], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p7], transformed_data.T[1][initInd_p7], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p8], transformed_data.T[1][initInd_p8], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p9], transformed_data.T[1][initInd_p9], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p10], transformed_data.T[1][initInd_p10], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p10], transformed_data.T[1][initInd_p11], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p10], transformed_data.T[1][initInd_p12], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p13], transformed_data.T[1][initInd_p13], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p14], transformed_data.T[1][initInd_p14], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p15], transformed_data.T[1][initInd_p15], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p16], transformed_data.T[1][initInd_p16], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p17], transformed_data.T[1][initInd_p17], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p18], transformed_data.T[1][initInd_p18], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p19], transformed_data.T[1][initInd_p19], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p20], transformed_data.T[1][initInd_p20], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p21], transformed_data.T[1][initInd_p21], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p22], transformed_data.T[1][initInd_p22], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p23], transformed_data.T[1][initInd_p23], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p24], transformed_data.T[1][initInd_p24], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p25], transformed_data.T[1][initInd_p25], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p26], transformed_data.T[1][initInd_p26], c='#1F51FF', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p27], transformed_data.T[1][initInd_p27], c='#1F51FF', s=10, alpha = 1.)

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'ambrE1_P1.png')
plt.savefig(fileName, transparent=True, dpi=600)
plt.show()

###################################################################################################################


plt.figure(figsize=(8, 6))
fig, ax = plt.subplots()
scatter = ax.scatter(transformed_data.T[0], transformed_data.T[1], c='#39FF14', s=5, alpha = 0.5, linewidths=0.1, zorder=1)
plt.xlabel(f'PC1 ({explained_variance[0]*100:.2f} %)', fontsize=16)
plt.ylabel(f'PC2 ({explained_variance[1]*100:.2f} %)', fontsize=16)

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
plt.tight_layout()



ax.scatter(transformed_data.T[0][initInd_c1b], transformed_data.T[1][initInd_c1b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c2b], transformed_data.T[1][initInd_c2b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c3b], transformed_data.T[1][initInd_c3b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c4b], transformed_data.T[1][initInd_c4b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c5b], transformed_data.T[1][initInd_c5b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c6b], transformed_data.T[1][initInd_c6b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c7b], transformed_data.T[1][initInd_c7b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c8b], transformed_data.T[1][initInd_c8b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c9b], transformed_data.T[1][initInd_c9b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c10b], transformed_data.T[1][initInd_c10b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c11b], transformed_data.T[1][initInd_c11b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c12b], transformed_data.T[1][initInd_c12b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c13b], transformed_data.T[1][initInd_c13b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c14b], transformed_data.T[1][initInd_c14b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c15b], transformed_data.T[1][initInd_c15b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c16b], transformed_data.T[1][initInd_c16b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c17b], transformed_data.T[1][initInd_c17b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c18b], transformed_data.T[1][initInd_c18b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_c19b], transformed_data.T[1][initInd_c19b], c='#bc13fe', s=10, alpha = 1.)
# ax.scatter(transformed_data.T[0][initInd_c21b], transformed_data.T[1][initInd_c21b], c='#bc13fe', s=10, alpha = 1.)
# ax.scatter(transformed_data.T[0][initInd_c23b], transformed_data.T[1][initInd_c23b], c='#bc13fe', s=10, alpha = 1.)
# ax.scatter(transformed_data.T[0][initInd_c25b], transformed_data.T[1][initInd_c25b], c='#bc13fe', s=10, alpha = 1.)
# ax.scatter(transformed_data.T[0][initInd_c27b], transformed_data.T[1][initInd_c27b], c='#bc13fe', s=10, alpha = 1.)
# ax.scatter(transformed_data.T[0][initInd_c29b], transformed_data.T[1][initInd_c29b], c='#bc13fe', s=10, alpha = 1.)
# ax.scatter(transformed_data.T[0][initInd_c31b], transformed_data.T[1][initInd_c31b], c='#bc13fe', s=10, alpha = 1.)

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'ambrE2_C2.png')
plt.savefig(fileName, transparent=True, dpi=600)
plt.show()
##############################################################################################################################



plt.figure(figsize=(8, 6))
fig, ax = plt.subplots()
scatter = ax.scatter(transformed_data.T[0], transformed_data.T[1], c='#39FF14', s=5, alpha = 0.5, linewidths=0.1, zorder=1)
plt.xlabel(f'PC1 ({explained_variance[0]*100:.2f} %)', fontsize=16)
plt.ylabel(f'PC2 ({explained_variance[1]*100:.2f} %)', fontsize=16)

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
plt.tight_layout()



ax.scatter(transformed_data.T[0][initInd_c19b], transformed_data.T[1][initInd_c19b], c='#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p1b], transformed_data.T[1][initInd_p1b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p2b], transformed_data.T[1][initInd_p2b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p3b], transformed_data.T[1][initInd_p3b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p4b], transformed_data.T[1][initInd_p4b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p5b], transformed_data.T[1][initInd_p5b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p6b], transformed_data.T[1][initInd_p6b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p7b], transformed_data.T[1][initInd_p7b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p8b], transformed_data.T[1][initInd_p8b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p9b], transformed_data.T[1][initInd_p9b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p10b], transformed_data.T[1][initInd_p10b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p11b], transformed_data.T[1][initInd_p11b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p12b], transformed_data.T[1][initInd_p12b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p13b], transformed_data.T[1][initInd_p13b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p14b], transformed_data.T[1][initInd_p14b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p15b], transformed_data.T[1][initInd_p15b], c = '#bc13fe', s=10, alpha = 1.)


ax.scatter(transformed_data.T[0][initInd_p16b], transformed_data.T[1][initInd_p16b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p17b], transformed_data.T[1][initInd_p17b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p18b], transformed_data.T[1][initInd_p18b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p19b], transformed_data.T[1][initInd_p19b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p20b], transformed_data.T[1][initInd_p20b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p21b], transformed_data.T[1][initInd_p21b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p22b], transformed_data.T[1][initInd_p22b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p23b], transformed_data.T[1][initInd_p23b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p24b], transformed_data.T[1][initInd_p24b], c = '#bc13fe', s=10, alpha = 1.)
ax.scatter(transformed_data.T[0][initInd_p25b], transformed_data.T[1][initInd_p25b], c = '#bc13fe', s=10, alpha = 1.)


fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'ambrE2_P2.png')
plt.savefig(fileName, transparent=True, dpi=600)
plt.show()


# Labels and title
