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
    #choice = np.random.choice(np.arange(N), size = int(N*.8), replace=False)
    #x[choice] = 0.0
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
            'growthRate': np.random.uniform(low = 0.1, high =1.3),
            'source': [],
            'sink': [],
            'interactions': interactions
        })
    return dicts

def ode(N, t, epsilon, derivs):
    """
    Represents the Ordinary Differential Equation for the system.
    """
    derivatives = np.zeros(len(N))
    for i in range(len(N)):
        derivatives[i] = max(-N[i], derivs[i](N, epsilon))
    return derivatives

def diffusion(N, t, coeff):
    """
    Represents the diffusion component of the system with a given diffusion coefficient.
    """
    return np.diag(N * coeff)

def make_diffusion_function(coeff):
    """
    Creates a diffusion function with a fixed coefficient.
    """
    def diffusion_function(N, t):
        return diffusion(N, t, coeff)
    return diffusion_function


def get_ss(solObj):
    y = solObj.T
    sp1 = y[0] + y[1]
    #sp2 = y[2] + y[3]
    #sp3 = y[4] + y[5]
    
    v = [sp1[-1]]#, sp2[-1], sp3[-1]]
    
    
    for i in y[2::]:
        v.append(i[-1])
    
    return np.array(v)

def close(func, *args):
    def newfunc(x, t):
        return func(x, t, *args)
    return newfunc

def main(N = 50, 
         Nc = 1000, 
         negInt1 = 75, 
         gr = 0.3, 
         diffusion_coeff = 0.2,
         root = None):
    
    
    
    derivDicts = random_deriv_dict(N)
    
    fxa_xb = get_hill(0.1, 10, 0.1, HillType.ACTIVATION)
    fxb_xa = get_hill(0.1, 10, 0.5, HillType.INHIBITION)
    
    
    d1 = derivDicts[0]
    d2 = derivDicts[1]
    
    
    #xa
    
    d1['growthRate'] = 0.1
    d1['source'] = [(1, fxb_xa)]
    d1['sink'] = [(0, fxa_xb)]
    d1['interactions'][1] = 0
    
    ssp = get_interaction_v(N, 50, negInt1)
    
    ssp[0] = 0 #interaction with its other phenotype
    ssp[1] = -1 #self interaction
    
    #xb
    d2['growthRate'] = gr
    d2['source'] = [(0, fxa_xb)]
    d2['sink'] = [(1, fxb_xa)]
    
    derivDicts[0] = d1.copy()
    derivDicts[1] = d2.copy()
    
    
    
    
    
    
    
    for i in range(len(ssp)):
        derivDicts[i]['interactions'][1] = ssp[i]
    
    
    
    
    derivs = [build_derivative(derivDicts[i]) for i in range(N)]
    
    
    
    
    
    
    
    
    
    envP = []
    ss = []
    
    
    
    while len(ss)<Nc:
        
        N0 = np.random.rand(N) * 0.1
        epsilon = np.random.uniform(low=0, high=0.25)
        
    
        args = (epsilon, derivs)
        current_diffusion = make_diffusion_function(diffusion_coeff)
        sol = sdeint.stratHeun(close(ode, *args), current_diffusion, N0, np.linspace(0, 100, 1000))
        #sol = solve_ivp(ode, (0,1000), y0 = N0, t_eval=np.linspace(0,100,10000), args=args)
        
        fsol = get_ss(sol)
        
        if sum(np.isnan(fsol))==0:
            envP.append(epsilon)
            ss.append(fsol)
        #ss.append(sol['y'].T[-1])
    
    ss = np.array(ss)
    sorter = np.argsort(envP)
    plt.pcolormesh((ss[sorter]).T, cmap=cm.coolwarm)
    plt.colorbar(label='abundance')
    
    if root is not None:
        plt.tight_layout()
        plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'phenomenological', 'heatMap_' + root + '.png'), transparent=True, dpi=600)
    plt.show()

    pc = pca.fit_transform(ss)
    explained_variance = pca.explained_variance_ratio_
    
    plt.scatter(pc.T[0], pc.T[1], s=10, c=envP, cmap=cm.coolwarm)
    plt.xlabel(f'PC1 ({explained_variance[0] * 100:.2f} %)', fontsize=16)
    plt.ylabel(f'PC2 ({explained_variance[1] * 100:.2f} %)', fontsize=16)
    plt.colorbar()
    if root is not None:
        plt.tight_layout()
        plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'phenomenological', 'pca_' + root + '.png'), transparent=True, dpi=600)
    
    plt.show()
    
    sspA = get_interaction_v(10000, 50, 51)
    sspB = get_interaction_v(10000, 50, negInt1)
    
    plt.hist(sspA, bins=100, density=True, histtype='step', color='k')
    plt.hist(sspA, bins=100, density=True, histtype='bar', color='gray')
    plt.hist(sspB, bins=100, density=True, histtype='step', color='red')
    plt.hist(sspB, bins=100, density=True, histtype='bar', color='red', alpha=0.5)
    if root is not None:
        plt.tight_layout()
        plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'phenomenological', 'dist_' + root + '.png'), transparent=True, dpi=600)
    plt.show()
    
    plt.scatter(envP, ss.T[0])
    plt.show()
    
    plt.hist(ss.T[0], density=True, histtype='step', color='red')
    plt.hist(ss.T[0], density=True, histtype='bar', color='red', alpha=0.5)
    

    
    if root is not None:
        plt.tight_layout()
        plt.savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'phenomenological', 'bidist_' + root + '.png'), transparent=True, dpi=600)
    plt.show()



for i in tqdm(np.linspace(51,100, 5)):
    main(Nc=100, negInt1=i, gr=0.1, diffusion_coeff=0.40, root = 'NC_' + str(np.round(i)) + '_gr_0.5_DC_0.2')
    