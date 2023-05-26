# -*- coding: utf-8 -*-
"""
Created on Mon May 22 22:21:39 2023

@author: danie
"""
from pathlib import Path
import os
import sys


from pylab import *
import numpy as np
from scipy.integrate import solve_ivp









def getHill(k, h, r, tp = 'inhibition'):
    
    def hill(envP):
        if tp=='inhibition':
            return r*(k**h)/(k**h + envP**h)
        else:
            return r*(envP**h)/(k**h + envP**h)
    
    
    return hill




def dxadt(N):
    
    return N[0] * (mu_xa - N[0] - fxa_xb(epsilon)) + N[1] * fxb_xa(epsilon)


def dxbdt(N):
    
    return N[1] * (mu_xb - N[1] - 0.4*N[2] - 0.3*N[4] - fxb_xa(epsilon)) + N[0]*fxa_xb(epsilon)


def dxedt(N):
    
    return N[2] * (mu_xe - 0.5*N[1] - N[2] - 0.3*N[4] - fxe_xf(epsilon)) + N[3]*fxf_xe(epsilon)


def dxfdt(N):
    
    return N[3] * (mu_xf - N[3] - fxf_xe(epsilon)) + N[2]*fxe_xf(epsilon)


def dxidt(N):
    
    return N[4] * (mu_xi - 0.5 * N[1] - 0.4 * N[2] - N[4] - fxi_xj(epsilon)) + N[5]*fxj_xi(epsilon)


def dxjdt(N):
    
    return N[5] * (mu_xj + 0.1*(N[0] + N[1] + N[2] + N[3] + N[5]) - N[5] - fxi_xj(epsilon)) + N[4]*fxi_xj(epsilon)



def ode(t, N):
    
    derivatives = np.zeros(len(N))
    
    derivatives[0] = dxadt(N)
    
    derivatives[1] = dxbdt(N)
    
    derivatives[2] = dxedt(N)
    
    derivatives[3] = dxfdt(N)
    
    derivatives[4] = dxidt(N)
    
    derivatives[5] = dxjdt(N)
    
    return derivatives
    
    
    

def getSS(solObj):
    y = solObj.y
    bh = y[0] + y[1]
    bt = y[2] + y[3]
    ri = y[4] + y[5]
    return np.array([bh[-1],bt[-1],ri[-1]])
    


#########Parameters#######

epsilon = 0.01
mu_xa = 0.18
mu_xb = 0.79
mu_xe = 0.802
mu_xf = 0.99
mu_xi = 0.8
mu_xj = 0.003
h_xa_xb = 9
h_xb_xa = 10
h_xe_xf = 1
h_xf_xe = 10
h_xi_xj = 9
h_xj_xi = 10
k_xa_xb = 0.1
k_xb_xa = 0.1
k_xe_xf = 0.1
k_xf_xe = 0.1
k_xi_xj = 0.1
k_xj_xi = 0.1
r_xa_xb = 0.1
r_xb_xa = 1.5
r_xe_xf = 0#0.99
r_xf_xe = 0#0.0001
r_xi_xj = 0#0.025
r_xj_xi = 0#0.01

#######Activation functions#########

fxa_xb = getHill(k_xa_xb, h_xa_xb, r_xa_xb, 'activation')
fxb_xa = getHill(k_xb_xa, h_xb_xa, r_xb_xa, 'inhibition')
fxe_xf = getHill(k_xe_xf, h_xe_xf, r_xe_xf, 'activation')
fxf_xe = getHill(k_xf_xe, h_xf_xe, r_xf_xe, 'inhibition')
fxi_xj = getHill(k_xi_xj, h_xi_xj, r_xi_xj, 'activation')
fxj_xi = getHill(k_xj_xi, h_xj_xi, r_xj_xi, 'inhibition')

#########Initial condition##########


xa = 0.01
xb = 0.0
xe = 0.01
xf = 0.0
xi = 0.01
xj = 0.0


N0 = np.zeros(6)

N0[0] = xa
N0[1] = xb
N0[2] = xe
N0[3] = xf
N0[4] = xi
N0[5] = xj

# ss = []
# for i in np.linspace(0, 0.21, 1000):
#     epsilon = i
#     sol = solve_ivp(ode, (0,1000), y0 = N0, t_eval=np.linspace(0,100,10000))
#     ss.append(getSS(sol))




# ss = np.array(ss)
# pcolormesh(ss.T, cmap=cm.coolwarm, vmin = 0.2, vmax = 0.6) 
# colorbar(label='abundance')
# mStates = ['bh', 'bt', 'ri']
# tp = [0,200,400,600,800,999]
# yticks(np.arange(len(mStates))+0.5, mStates)
# xticks(tp, np.round(np.linspace(0, 0.21, 6),3))
# xlabel('environmental cue')
# #savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'twoStatesPhenom.png'), dpi = 600)
# show()

epsilon = 0.1355
sol = solve_ivp(ode, (0,1000), y0 = N0, t_eval=np.linspace(0,100,10000))

y = sol.y
bh = y[0] + y[1]
bt = y[2]+y[3]
ri = y[4]+y[5]
plot(sol.t, bh, color = 'r', label='bh'); plot(sol.t, bt, color = 'orange', label='bt'); plot(sol.t,ri, color='blue', label='ri');legend()
show()