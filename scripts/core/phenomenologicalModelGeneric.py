# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:31:00 2023

@author: danie
"""

import numpy as np


from enum import Enum

class HillType(Enum):
    INHIBITION = 'inhibition'
    ACTIVATION = 'activation'
    INDEPENDENT = 'independent'
    
    

def get_hill(k: float, h: float, r: float, tp: HillType = HillType.INHIBITION):
    """
    Generates a Hill function based on the provided parameters.
    """
    def hill(envP: float) -> float:
        if tp == HillType.INHIBITION:
            return r * (k**h) / (k**h + envP**h)
        elif tp== HillType.ACTIVATION:
            return r * (envP**h) / (k**h + envP**h)
        else:
            return r
            

    return hill


def build_derivative(derivDict: dict):
    '''
    

    Parameters
    ----------
    derivDict : dict
        {
            'index': speciesKindex
            'growthRate': speciesKgrowthRate
            'source' : [(speciesA, hillFuncA), (speciesB, hillFuncB), ...],
            'sink' : [(speciesA, hillFuncA), (speciesB, hillFuncB), ...],
            'interactions' : {'speciesA' : x1,
                              'speciesB' : x2,
                              ...}
            }
    
    'speciesA' and 'speciesB' ... are indexes of the species on a vector of N species
    found in the system

    Returns 
    -------
    the derivative function of a species

    '''
    
    index = derivDict['index']
    growth_rate = derivDict['growthRate']

    def dxdt(N: np.ndarray, epsilon: float) -> float:
        interactions = sum(max(0,N[i]) * derivDict['interactions'][i] for i in range(len(N)))
        sink = sum(t[1](epsilon) * N[index] for t in derivDict['sink'])
        source = sum(max(N[t[0]],0) * t[1](epsilon) for t in derivDict['source'])
        
        return max(0,N[index]) * (growth_rate + interactions) - sink + source

    return dxdt






##########example###########



# fxa_xb = getHill(0.1, 9, 0.1, 'activation')
# fxb_xa = getHill(0.1, 10, 1.5, 'inhibition')

# fxe_xf = getHill(1, 1, 0.99, 'activation')
# fxf_xe = getHill(0.1, 10, 0.0001, 'inhibition')

# fxi_xj = getHill(0.1, 9, 0.025, 'activation')
# fxj_xi = getHill(0.1, 10, 0.01, 'inhibition')



# xa = 0.01
# xb = 0.0
# xe = 0.01
# xf = 0.0
# xi = 0.01
# xj = 0.0


# N0 = np.zeros(6)

# N0[0] = xa
# N0[1] = xb
# N0[2] = xe
# N0[3] = xf
# N0[4] = xi
# N0[5] = xj

# epsilon = 0.01


# xa_derivDict = {'index':0,
#                 'growthRate': 0.18,
#                 'source':[(1, fxb_xa)],
#                 'sink': [(0, fxa_xb)],
#                 'interactions':{0:-1,
#                                 1:0,
#                                 2:0,
#                                 3:0,
#                                 4:0,
#                                 5:0}}  



# xb_derivDict = {'index':1,
#                 'growthRate': 0.79,
#                 'source':[(0, fxa_xb)],
#                 'sink': [(1, fxb_xa)],
#                 'interactions':{0:0,
#                                 1:-1,
#                                 2:-0.4,
#                                 3:0,
#                                 4:-0.3,
#                                 5:0}}

# xe_derivDict = {'index':2,
#                 'growthRate': 0.802,
#                 'source':[(3, fxf_xe)],
#                 'sink': [(2, fxe_xf)],
#                 'interactions':{0:0,
#                                 1:-0.5,
#                                 2:-1,
#                                 3:0,
#                                 4:-0.3,
#                                 5:0}}


# xf_derivDict = {'index':3,
#                 'growthRate': 0.99,
#                 'source':[(2, fxe_xf)],
#                 'sink': [(3, fxf_xe)],
#                 'interactions':{0:0,
#                                 1:0,
#                                 2:0,
#                                 3:-1,
#                                 4:0,
#                                 5:0}}

  
# xi_derivDict = {'index':4,
#                 'growthRate': 0.8,
#                 'source':[(5, fxj_xi)],
#                 'sink': [(4, fxi_xj)],
#                 'interactions':{0:0,
#                                 1:-0.5,
#                                 2:-0.4,
#                                 3:0,
#                                 4:-1,
#                                 5:0}}

# xj_derivDict = {'index':4,
#                 'growthRate': 0.003,
#                 'source':[(4, fxi_xj)],
#                 'sink': [(5, fxj_xi)],
#                 'interactions':{0:0.1,
#                                 1:0.1,
#                                 2:0.1,
#                                 3:0.1,
#                                 4:0.1,
#                                 5:-1}}



# dxadt = buildDerivative(xa_derivDict)
# dxbdt = buildDerivative(xb_derivDict)
# dxedt = buildDerivative(xe_derivDict)
# dxfdt = buildDerivative(xf_derivDict)
# dxidt = buildDerivative(xi_derivDict)
# dxjdt = buildDerivative(xj_derivDict)


# def ode(t, N, epsilon):
    
#     derivatives = np.zeros(len(N))
    
#     derivatives[0] = dxadt(N, epsilon)
    
#     derivatives[1] = dxbdt(N, epsilon)
    
#     derivatives[2] = dxedt(N, epsilon)
    
#     derivatives[3] = dxfdt(N, epsilon)
    
#     derivatives[4] = dxidt(N, epsilon)
    
#     derivatives[5] = dxjdt(N, epsilon)
    
#     return derivatives


# def getSS(solObj):
#     y = solObj.y
#     bh = y[0] + y[1]
#     bt = y[2] + y[3]
#     ri = y[4] + y[5]
#     return np.array([bh[-1],bt[-1],ri[-1]])



# ss = []
# for i in np.linspace(0, 0.21, 1000):
#     epsilon = i
#     sol = solve_ivp(ode, (0,1000), y0 = N0, t_eval=np.linspace(0,100,10000), args=(epsilon,))
#     ss.append(getSS(sol))




# ss = np.array(ss)
# pcolormesh(ss.T, cmap=cm.coolwarm, vmin = 0.2, vmax = 1.2) 
# colorbar(label='abundance')
# mStates = ['bh', 'bt', 'ri']
# tp = [0,200,400,600,800,999]
# yticks(np.arange(len(mStates))+0.5, mStates)
# xticks(tp, np.round(np.linspace(0, 0.21, 6),3))
# xlabel('environmental cue')
# #savefig(os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'twoStatesPhenom.png'), dpi = 600)
# show()
