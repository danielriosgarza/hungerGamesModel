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
import matplotlib.colors as mcolors

# Adjusting the import path
project_root = Path(os.getcwd()).parents[0]
sys.path.append(os.path.join(project_root, 'core'))
from phenomenologicalModelGeneric import *


# Theme setup for plots
theme = load_theme("boxy_light")
theme.apply()


from phenomenologicalModelGeneric import *



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



def create_heatmap(data, row_labels, x_positions, x_values, xlabel, ylabel, title, fileName=None):
    fig, ax = plt.subplots()
    
    # Generate a meshgrid for the data points. Add 1 to include the end of the last cell.
    x_edges = np.arange(data.shape[1] + 1) - 0.5
    y_edges = np.arange(data.shape[0] + 1) - 0.5
    
    # Use pcolormesh to plot the heatmap with edges, which naturally includes gaps between cells
    cmap = plt.get_cmap('viridis')  # Choose a colormap
    cmap.set_bad(color='white')  # Set color for NaN or masked values to white
    heatmap = ax.pcolormesh(x_edges, y_edges, np.ma.masked_invalid(data), cmap=cmap, shading='flat')
    
    ax.grid(False)
    
    # Set row labels (y-axis)
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_yticklabels(row_labels)
    
    # Make y-tick labels italic
    for label in ax.get_yticklabels():
        label.set_fontstyle('italic')
    
    # Reverse the y-axis to have the first row at the top
    ax.invert_yaxis()
    
    # Set specific x-ticks and labels (x-axis)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_values)
    
    # Set labels and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    # Add colorbar
    cbar = fig.colorbar(heatmap, ax=ax)
    
    plt.tight_layout()
    
    if fileName is not None:
        plt.savefig(fileName, transparent=True, dpi=600)
    
    plt.show()





##########parameters##############
k_1 = 0.1
h_1 = 9.0
r_1 = 0.025


k_2 = 0.1
h_2 = 10.0
r_2 = 1.5


k_3 = 1.0
h_3 = 1.0
r_3 = 1.495

k_4 = 1.0
h_4 = 1.0
r_4 = 0.001


k_5 = 0.1
h_5 = 9.0
r_5 = 0.021

k_6 = 0.1
h_6 = 10.0
r_6 = 0.214



xa_gr = 0.192
#xa_gr = 0.5
xb_gr = 0.978
#xb_gr = 0.5
xe_gr = 0.921
#xe_gr = 0.5
xf_gr = 1.19
#xf_gr = 0.5
xi_gr = 0.705
#xi_gr = 0.5
xj_gr = 0.01
#xj_gr = 0.5




xaInt = {0:-1,
         1:0,
         2:0,
         3:0,
         4:0,
         5:0.0}

xbInt = {0:0,
         1:-1,
         2:-0.4,
         3:0,
         4:-0.3,
         5:0.0}

xeInt = {0:0,
         1:-0.9,
         2:-1.0,
         3:0,
         4:-0.3,
         5:0.0}


xfInt = {0:0,
         1:-0.4,
         2:0,
         3:-1,
         4:0,
         5:0}


xiInt = {0:0,
         1:-0.9,
         2:-0.4,
         3:0,
         4:-1,
         5:0}

xjInt = {0:0.1,
         1:0.1,
         2:0.1,
         3:0.1,
         4:0.1,
         5:-1}




##########initial conditions##########


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


#############build functions##################

fxa_xb = get_hill(k_1, h_1, r_1, HillType.ACTIVATION)
fxb_xa = get_hill(k_2, h_2, r_2, HillType.INHIBITION)

fxe_xf = get_hill(k_3, h_3, r_3, HillType.INDEPENDENT)
fxf_xe = get_hill(k_4, h_4, r_4, HillType.INDEPENDENT)

fxi_xj = get_hill(k_5, h_5, r_5, HillType.ACTIVATION)
fxj_xi = get_hill(k_6, h_6, r_6, HillType.INHIBITION)


xa_derivDict = {'index':0,
                'growthRate': xa_gr,
                'source': [(1, fxb_xa)],
                'sink': [(0, fxa_xb)],
                'interactions':xaInt}  



xb_derivDict = {'index':1,
                'growthRate': xb_gr,
                'source': [(0, fxa_xb)],
                'sink': [(1, fxb_xa)],
                'interactions':xbInt}

xe_derivDict = {'index':2,
                'growthRate': xe_gr,
                'source': [(3, fxf_xe)],
                'sink': [(2, fxe_xf)],
                'interactions':xeInt}


xf_derivDict = {'index':3,
                'growthRate': xf_gr,
                'source': [(2, fxe_xf)],
                'sink': [(3, fxf_xe)],
                'interactions':xfInt}

  
xi_derivDict = {'index':4,
                'growthRate': xi_gr,
                'source': [(5, fxj_xi)],
                'sink': [(4, fxi_xj)],
                'interactions':xiInt}

xj_derivDict = {'index':4,
                'growthRate': xj_gr,
                'source': [(4, fxi_xj)],
                'sink': [(5, fxj_xi)],
                'interactions':xjInt}



dxadt = build_derivative(xa_derivDict)
dxbdt = build_derivative(xb_derivDict)
dxedt = build_derivative(xe_derivDict)
dxfdt = build_derivative(xf_derivDict)
dxidt = build_derivative(xi_derivDict)
dxjdt = build_derivative(xj_derivDict)


############simulate############
simls = 10000

ss = []
for i in tqdm(np.linspace(0, 0.21, simls)):
    epsilon = i
    sol = solve_ivp(ode, (0,1000), y0 = N0, t_eval=np.linspace(0,100,10000), args=(epsilon,))
    ss.append(getSS(sol))






########plot results################
ss = np.array(ss)
sum_ss = np.sum(ss,axis=1)
ss = np.array([ss[i]/sum_ss[i] for i in range(len(ss))])

mStates = ['Blautia hydrogenotrophica',
           'Bacteroides thetaiotaomicron',
           'Roseburia intinalis']
tp = np.linspace(0, simls,6 )


create_heatmap(ss.T, 
               mStates, 
               tp, 
               np.round(np.linspace(0, 0.21, 6),3), 
               'environmental cue', 
               None, 
               None,
               fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'phenomenological', 'environmtnCue.png'))
