# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 09:18:24 2023

@author: drgarza
"""

from pathlib import Path
import os
import sys
from scipy.integrate import solve_ivp
from tqdm import tqdm

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from aquarel import load_theme

# Adjusting the import path
project_root = Path(os.getcwd()).parents[0]
sys.path.append(os.path.join(project_root, 'core'))
from phenomenologicalModelGeneric import *


# Theme setup for plots
theme = load_theme("boxy_light")
theme.apply()


from phenomenologicalModelGeneric import *

fxa_xb = get_hill(0.1, 9, 0.1, HillType.ACTIVATION)
fxb_xa = get_hill(0.1, 10, 1.5, HillType.INHIBITION)

fxe_xf = get_hill(1, 1, 0.01, HillType.INDEPENDENT)
fxf_xe = get_hill(0.1, 10, 0.01, HillType.INDEPENDENT)

fxi_xj = get_hill(0.1, 9, 0.025, HillType.ACTIVATION)
fxj_xi = get_hill(0.1, 10, 0.01, HillType.INHIBITION)



xa = 0.01
xb = 0.01
xe = 0.01
xf = 0.01
xi = 0.01
xj = 0.01


N0 = np.zeros(6)

N0[0] = xa
N0[1] = xb
N0[2] = xe
N0[3] = xf
N0[4] = xi
N0[5] = xj

epsilon = 0.01


xa_derivDict = {'index':0,
                'growthRate': 0.192,
                'source': [(1, fxb_xa)],
                'sink': [(0, fxa_xb)],
                'interactions':{0:-1,
                                1:0,
                                2:0,
                                3:0,
                                4:0,
                                5:0.0}}  



xb_derivDict = {'index':1,
                'growthRate': 0.978,
                'source': [(0, fxa_xb)],
                'sink': [(1, fxb_xa)],
                'interactions':{0:0,
                                1:-1,
                                2:-0.4,
                                3:0,
                                4:-0.3,
                                5:0.0}}

xe_derivDict = {'index':2,
                'growthRate': 0.921,
                'source': [(3, fxf_xe)],
                'sink': [(2, fxe_xf)],
                'interactions':{0:0,
                                1:-0.9,
                                2:-1,
                                3:0,
                                4:-0.3,
                                5:0.0}}


xf_derivDict = {'index':3,
                'growthRate': 1.19,
                'source': [(2, fxe_xf)],
                'sink': [(3, fxf_xe)],
                'interactions':{0:0,
                                1:-0.5,
                                2:0,
                                3:-1,
                                4:0,
                                5:0.0}}

  
xi_derivDict = {'index':4,
                'growthRate': 0.705,
                'source': [(5, fxj_xi)],
                'sink': [(4, fxi_xj)],
                'interactions':{0:0,
                                1:-0.9,
                                2:-0.4,
                                3:0,
                                4:-1,
                                5:0.0}}

xj_derivDict = {'index':4,
                'growthRate': 0.01,
                'source': [(4, fxi_xj)],
                'sink': [(5, fxj_xi)],
                'interactions':{0:0.1,
                                1:0.1,
                                2:0.1,
                                3:0.1,
                                4:0.1,
                                5:-1}}



dxadt = build_derivative(xa_derivDict)
dxbdt = build_derivative(xb_derivDict)
dxedt = build_derivative(xe_derivDict)
dxfdt = build_derivative(xf_derivDict)
dxidt = build_derivative(xi_derivDict)
dxjdt = build_derivative(xj_derivDict)


def ode(t, N, epsilon):
    
    derivatives = np.zeros(len(N))
    
    derivatives[0] = dxadt(N, epsilon)
    
    derivatives[1] = dxbdt(N, epsilon)
    
    derivatives[2] = dxedt(N, epsilon)
    
    derivatives[3] = dxfdt(N, epsilon)
    
    derivatives[4] = dxidt(N, epsilon)
    
    derivatives[5] = dxjdt(N, epsilon)
    
    return derivatives


def getSS(solObj):
    y = solObj.y
    bh = y[0] + y[1]
    bt = y[2] + y[3]
    ri = y[4] + y[5]
    return np.array([bh[-1],bt[-1],ri[-1]])



def create_heatmap(data, row_labels, x_positions, x_values, xlabel, ylabel, title, fileName = None):
    fig, ax = plt.subplots()
    heatmap = ax.imshow(data, aspect='auto', interpolation = 'bilinear', vmax=1, vmin=0)
    ax.grid(False)

    # Set row labels (y-axis)
    ax.set_yticks(np.arange(data.shape[0]))
    
    yticklabels = ax.set_yticklabels(row_labels)
    yticklabels[-1].set_fontstyle('italic')
    yticklabels[-2].set_fontstyle('italic')
    yticklabels[-3].set_fontstyle('italic')

    # Set specific x-ticks and labels (x-axis)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_values)

    # Set labels and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    # Add grid
    #ax.set_xticks(np.arange(-.5, data.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, data.shape[0], 1), minor=True)
    ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5)
    #ax.tick_params(which="minor", size=0)

    # Add colorbar
    cbar = ax.figure.colorbar(heatmap, ax=ax)
    plt.tight_layout()
    
    if fileName is not None:
        plt.savefig(fileName, transparent=True, dpi=600)

    plt.show()

ss = []
for i in tqdm(np.linspace(0, 0.21, 10000)):
    epsilon = i
    sol = solve_ivp(ode, (0,1000), y0 = N0, t_eval=np.linspace(0,100,10000), args=(epsilon,))
    ss.append(getSS(sol))




ss = np.array(ss)
sum_ss = np.sum(ss,axis=1)
ss = np.array([ss[i]/sum_ss[i] for i in range(len(ss))])

mStates = ['Blautia hydrogenotrophica',
           'Bacteroides thetaiotaomicron',
           'Roseburia intinalis']
tp = [0,200,400,600,800,999]


create_heatmap(ss.T, 
               mStates, 
               tp, 
               np.round(np.linspace(0, 0.21, 6),3), 
               'environmental cue', 
               None, 
               None,
               fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'phenomenological', 'environmtnCue.png'))
