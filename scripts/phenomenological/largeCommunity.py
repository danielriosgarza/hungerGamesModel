# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 16:29:22 2023

@author: danie
"""

from pathlib import Path
import os
import sys
import sdeint


from pylab import *
import numpy as np
from scipy.integrate import solve_ivp

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'core'))

from phenomenologicalModelGeneric import *


import numpy as np
from scipy.stats import beta

from sklearn.decomposition import PCA
pca = PCA(n_components=2, whiten=1)

import matplotlib.pyplot as plt
from aquarel import load_theme

theme = load_theme("boxy_light")
theme.apply()


def getInteractionV(N, a, b):
    x = beta.rvs(a, b, size=N)
    
    return 2*x - 1





def randomDerivDict(N,a=50, b=51):
    
    
    dicts = []
    
    for i in range(N):
        interactions = getInteractionV(N, a, b)
        interactions[i] = -1
        dicts.append({'index':i,
                      'growthRate': np.random.uniform(),
                      'source':[],
                      'sink': [],
                      'interactions':{k:interactions[k] for k in range(N)}})
        
    return dicts


def ode(N, t, epsilon, derivs):
    derivatives = np.zeros(len(N))
    
    for i in range(len(N)):
        derivatives[i] = max(-N[i], derivs[i](N, epsilon))
    
    return derivatives




N = 50
negInt = 75

def diffusion(N, t):
    return np.diag(N*0.3)


derivDicts = randomDerivDict(N)

fxa_xb = getHill(0.1, 9, 0.1, 'activation')
fxb_xa = getHill(0.1, 10, 1.5, 'inhibition')

d1 = derivDicts[0]
d2 = derivDicts[1]

d1['growthRate'] = 0.18
d1['source'] = [(1, fxb_xa)]
d1['sink'] = [(0, fxa_xb)]
d1['interactions'][1] = 0

ssp = getInteractionV(N, 50, negInt)

ssp[0] = 0 #interaction with its other phenotype
ssp[1] = -1 #self interaction

d2['growthRate'] = 0.79
d2['source'] = [(0, fxa_xb)]
d2['sink'] = [(1, fxb_xa)]

derivDicts[0] = d1.copy()
derivDicts[1] = d2.copy()

for i in range(len(ssp)):
    derivDicts[i]['interactions'][1] = ssp[i]


derivs = [buildDerivative(derivDicts[i]) for i in range(N)]



def getSS(solObj):
    y = solObj.T
    sp1 = y[0] + y[1]
    
    v = [sp1[-1]]
    
    for i in y[2::]:
        v.append(i[-1])
    
    return np.array(v)

def close(func, *args):
    def newfunc(x, t):
        return func(x, t, *args)
    return newfunc

N0 = 0.1*np.ones(N)
#N0[1] = 0.00001


envP = []

ss = []
for i in range(1000):
    
    print(i)
    epsilon = np.random.uniform(low=0, high = 0.25)
    envP.append(epsilon)
    
    # N0[0] = 0.00001
    # N0[1] = 0.1
    
    # if epsilon<0.25/2:
    #     N0[0] = 0.1
    #     N0[1] = 0.00001
    
    args=(epsilon,derivs)
    sol = sdeint.itoint(close(ode, *args), diffusion, N0, np.linspace(0,50,1000))
    ss.append(getSS(sol))

ss = np.array(ss)
pcolormesh(ss.T, cmap=cm.coolwarm) 
colorbar(label='abundance')
plt.show()

pc = pca.fit_transform(ss)

plt.scatter(pc.T[0], pc.T[1], s=10, c=envP, cmap=cm.coolwarm)
colorbar(label='environment cue')
plt.show()


sspA = getInteractionV(10000, 50,51)
sspB = getInteractionV(10000, 50,negInt)


plt.hist(sspA, bins=100, density=True, histtype='step', color='k')
plt.hist(sspA, bins=100, density=True, histtype='bar', color='gray')

plt.hist(sspB, bins=100, density=True, histtype='step', color='red')
plt.hist(sspB, bins=100, density=True, histtype='bar', color='red', alpha=.5)

plt.show()

