# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 16:29:22 2023

@author: danie
"""

from pathlib import Path
import os
import sys
import sdeint
from tqdm import tqdm
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import beta
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from aquarel import load_theme

# Adjusting the import path
project_root = Path(os.getcwd()).parents[0]
sys.path.append(os.path.join(project_root, 'core'))
from phenomenologicalModelGeneric import *

# PCA setup
pca = PCA(n_components=2, whiten=1)

# Theme setup for plots
theme = load_theme("boxy_light")
theme.apply()

def get_interaction_v(N, a, b):
    """
    Generates interaction values using a beta distribution.
    """
    x = beta.rvs(a, b, size=N)
    choice = np.random.choice(np.arange(N), size = int(N/2), replace=False)
    x[choice] = 0.0
    return 2 * x - 1
    

def random_deriv_dict(N, a=50, b=51):
    """
    Creates a list of dictionaries representing derivative functions.
    """
    dicts = []
    for i in range(N):
        interactions = get_interaction_v(N, a, b)
        interactions[i] = -1
        dicts.append({
            'index': i,
            'growthRate': np.random.uniform(),
            'source': [],
            'sink': [],
            'interactions': interactions
        })
    return dicts

def ode(t,N, epsilon, derivs):
    """
    Represents the Ordinary Differential Equation for the system.
    """
    derivatives = np.zeros(len(N))
    for i in range(len(N)):
        derivatives[i] = max(-N[i], derivs[i](N, epsilon))
    return derivatives

def diffusion(N, t):
    """
    Represents the diffusion component of the system.
    """
    return np.diag(N * 0.1)



def get_ss(solObj):
    y = solObj.T
    sp1 = y[0] + y[1]
    sp2 = y[2] + y[3]
    sp3 = y[4] + y[5]
    
    v = [sp1[-1], sp2[-1], sp3[-1]]
    
    
    for i in y[6::]:
        v.append(i[-1])
    
    return np.array(v)

def close(func, *args):
    def newfunc(x, t):
        return func(x, t, *args)
    return newfunc


N = 50
Nc = 500
negInt1 = 90
negInt2 = 65
negInt3 = 55




derivDicts = random_deriv_dict(N)

fxa_xb = get_hill(0.1, 9, 0.1, HillType.ACTIVATION)
fxb_xa = get_hill(0.1, 10, 1.5, HillType.INHIBITION)

fxe_xf = get_hill(1, 1, 0.99, HillType.ACTIVATION)
fxf_xe = get_hill(0.1, 10, 0.0001, HillType.INHIBITION)

fxi_xj = get_hill(0.1, 9, 0.025, HillType.ACTIVATION)
fxj_xi = get_hill(0.1, 10, 0.01, HillType.INHIBITION)

d1 = derivDicts[0]
d2 = derivDicts[1]
d3 = derivDicts[2]
d4 = derivDicts[3]
d5 = derivDicts[4]
d6 = derivDicts[5]

#xa

d1['growthRate'] = 0.19
d1['source'] = [(1, fxb_xa)]
d1['sink'] = [(0, fxa_xb)]
d1['interactions'][1] = 0

ssp = get_interaction_v(N, 50, negInt1)

ssp[0] = 0 #interaction with its other phenotype
ssp[1] = -1 #self interaction

#xb
d2['growthRate'] = 0.97
d2['source'] = [(0, fxa_xb)]
d2['sink'] = [(1, fxb_xa)]

derivDicts[0] = d1.copy()
derivDicts[1] = d2.copy()



#xe
d3['growthRate'] = 0.92
d3['source'] = [(3, fxf_xe)]
d3['sink'] = [(2, fxe_xf)]
d3['interactions'][3] = 0

ssp2 = get_interaction_v(N, 50, negInt2)

ssp2[2] = 0 #interaction with its other phenotype
ssp2[3] = -1 #self interaction

#xf
d4['growthRate'] = 0.3
d4['source'] = [(2, fxe_xf)]
d4['sink'] = [(3, fxf_xe)]

derivDicts[2] = d3.copy()
derivDicts[3] = d4.copy()

#xi
d4['growthRate'] = 0.70
d4['source'] = [(5, fxj_xi)]
d4['sink'] = [(4, fxi_xj)]
d4['interactions'][3] = 0

ssp3 = get_interaction_v(N, 50, negInt3)

ssp3[4] = 0 #interaction with its other phenotype
ssp3[5] = -1 #self interaction

#xj
d5['growthRate'] = 0.02
d5['source'] = [(4, fxi_xj)]
d5['sink'] = [(5, fxj_xi)]

derivDicts[4] = d4.copy()
derivDicts[5] = d5.copy()



for i in range(len(ssp)):
    derivDicts[i]['interactions'][1] = ssp[i]

for i in range(len(ssp2)):
    derivDicts[i]['interactions'][2] = ssp2[i]

for i in range(len(ssp3)):
    derivDicts[i]['interactions'][4] = ssp3[i]


derivs = [build_derivative(derivDicts[i]) for i in range(N)]









envP = []
ss = []

for i in tqdm(range(Nc)):

    
    N0 = np.random.rand(N) * 0.1
    epsilon = np.random.uniform(low=0, high=0.25)
    envP.append(epsilon)

    args = (epsilon, derivs)
    #sol = sdeint.stratHeun(close(ode, *args), diffusion, N0, np.linspace(0, 100, 1000))
    sol = solve_ivp(ode, (0,1000), y0 = N0, t_eval=np.linspace(0,100,10000), args=args)
    
    #fsol = get_ss(sol)
    
    #if sum(np.isnan(fsol))==0:
    #    ss.append(fsol)
    ss.append(sol['y'].T[-1])

ss = np.array(ss)
plt.pcolormesh(ss.T, cmap=cm.coolwarm)
plt.colorbar(label='abundance')
plt.show()

pc = pca.fit_transform(ss)
explained_variance = pca.explained_variance_ratio_

plt.scatter(pc.T[0], pc.T[1], s=10, c=envP, cmap=cm.coolwarm)
plt.xlabel(f'PC1 ({explained_variance[0] * 100:.2f} %)', fontsize=16)
plt.ylabel(f'PC2 ({explained_variance[1] * 100:.2f} %)', fontsize=16)
plt.colorbar()
plt.show()

sspA = get_interaction_v(10000, 50, 51)
sspB = get_interaction_v(10000, 50, negInt1)

plt.hist(sspA, bins=100, density=True, histtype='step', color='k')
plt.hist(sspA, bins=100, density=True, histtype='bar', color='gray')
plt.hist(sspB, bins=100, density=True, histtype='step', color='red')
plt.hist(sspB, bins=100, density=True, histtype='bar', color='red', alpha=0.5)
plt.show()