# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 15:09:43 2023

@author: danie
"""

import os
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from aquarel import load_theme

theme = load_theme("boxy_light")
theme.apply()
import seaborn as sns

import plotly.io as pio
pio.renderers.default='browser'

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'core'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'db'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))


from mainClasses import *
from parseTable import *
from updateParameters import *
from readModelDB import *
from loadParameters import *


###setup####

def simulateExperiment(group, 
                       experimentLabel, 
                       dbPath, 
                       measuredStates, 
                       combined=False, 
                       intervals=None,
                       starttime = 0,
                       endtime = 120):
    
    
    if intervals is None:
        intervals = [4]*len(measuredStates) 
    ipH_path = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bhbtri_ipH4.tsv')
    
    strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', group)
    
    initialStates = {}
    
    for i,v in enumerate(measuredStates):
        
        initialStates[v] = get_initialState(v, strainSummaryFolder, experimentLabel, combined, intervals[i]) 

    
    db = get_database(dbPath)
    

    wc = createMetabolome(db, 'wc')
    
    
    
    for state in initialStates:
        if ('live' not in state) and (state!='dead') and (state!='pH'):
            wc.metD[state].update(initialStates[state])


    predictpH = getpH(wc.metabolites, ipH_path)
    
    pH =  predictpH(wc.get_concentration())


    wc_f = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    
    
    for state in initialStates:
        if ('live' not in state) and (state!='dead') and (state!='pH'):
            wc_f.metD[state].update(initialStates[state])
    
    

    wc_r = createMetabolome(db, 'wc', pH, pHFunc=predictpH)
    
    
    for state in initialStates:
        if ('live' not in state) and (state!='dead') and (state!='pH'):
            wc_r.metD[state].update(initialStates[state])
    



    species_f = Microbiome({'bh':createBacteria(db, 'bh', 'wc'), 'bt':createBacteria(db, 'bt', 'wc'), 'ri':createBacteria(db, 'ri', 'wc')})
    species_f.subpopD['xa'].count = 0
    species_f.subpopD['xe'].count = 0
    species_f.subpopD['xi'].count = 0
        

    species_r = Microbiome({'bh':createBacteria(db, 'bh', 'wc'), 'bt':createBacteria(db, 'bt', 'wc'), 'ri':createBacteria(db, 'ri', 'wc')})
    species_r.subpopD['xa'].count = 0
    species_r.subpopD['xe'].count = 0
    species_r.subpopD['xi'].count = 0
    
    if 'live' in measuredStates:
        if group == 'bh':
            species_r.subpopD['xa'].count = initialStates['live']
            
        elif group == 'bt':
            species_r.subpopD['xe'].count = initialStates['live']
        
        elif group == 'ri':
            species_r.subpopD['xi'].count = initialStates['live']
            
    if 'live_bh' in measuredStates:
        species_r.subpopD['xa'].count = initialStates['live_bh']
    
    if 'live_bt' in measuredStates:
        species_r.subpopD['xe'].count = initialStates['live_bt']
    
    if 'live_ri' in measuredStates:
        species_r.subpopD['xi'].count = initialStates['live_ri']
    
    


    p1 = Pulse(wc_f, species_f, starttime, endtime, 10000, 0, 0, 0, 0)

    r_species = Reactor(species_r, wc_r, [p1], 60)

    r_species.simulate()
    
    return r_species




def genericSimulation(db):
    ipH_path = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'bhbtri_ipH4.tsv')
        
    wc = createMetabolome(db, 'wc')

    predictpH = getpH(wc.metabolites, ipH_path)

    pH =  predictpH(wc.get_concentration())


    wc_f = createMetabolome(db, 'wc', pH, pHFunc=predictpH)

    wc_r = createMetabolome(db, 'wc', pH, pHFunc=predictpH)


    species_f = Microbiome({'bh':createBacteria(db, 'bh', 'wc'), 'bt':createBacteria(db, 'bt', 'wc'), 'ri':createBacteria(db, 'ri', 'wc')})
    species_f.subpopD['xa'].count = 0
    species_f.subpopD['xe'].count = 0
    species_f.subpopD['xi'].count = 0
        

    species_r = Microbiome({'bh':createBacteria(db, 'bh', 'wc'), 'bt':createBacteria(db, 'bt', 'wc'), 'ri':createBacteria(db, 'ri', 'wc')})


    p1 = Pulse(wc_f, species_f, 0, 120, 10000, 0, 0, 0, 0)

    r_species = Reactor(species_r, wc_r, [p1], 60)

    r_species.simulate()
    
    return r_species





def makeExperimentPlot(group, 
             state, 
             stateType = 'metabolite',
             experiments = ['bhbtri'],
             lables = ['bh3'],
             colors = ['#ff0000'],
             simulObj = [None, None, None],
             alpha=1,
             legend = True,
             ylim = None):
    
    fig, ax = plt.subplots()
    strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', group)
    
    stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
            
    for i,v in enumerate(experiments):
        stateDF = getDFdict(stateTable, state, True)[v]
        sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = colors[i], lw=.666, label=lables[i], alpha=alpha)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})
        
        if simulObj[i] is not None:
            
            if state == 'pH':
                ax.plot(simulObj[i].time_simul, simulObj[i].pH_simul, color = colors[i], label=lables[i] + ' simul', linestyle='--', lw=5)
            
            elif state == 'live':
                
                if group == 'bh':
                    ax.plot(simulObj[i].time_simul, simulObj[i].cellActive_dyn[0], color = colors[i], label=lables[i] + ' simul', linestyle='--', lw=5)
                
                if group == 'bt':
                    ax.plot(simulObj[i].time_simul, simulObj[i].cellActive_dyn[1], color = colors[i], label=lables[i] + ' simul', linestyle='--', lw=5)
                
                if group == 'ri':
                    ax.plot(simulObj[i].time_simul, simulObj[i].cellActive_dyn[2], color = colors[i], label=lables[i] + ' simul', linestyle='--', lw=5)
            
            elif state == 'live_bh':
                
                ax.plot(simulObj[i].time_simul, simulObj[i].cellActive_dyn[0], color = colors[i], label=lables[i] + ' simul', linestyle='--', lw=5)
            
            elif state == 'live_bt':
                
                ax.plot(simulObj[i].time_simul, simulObj[i].cellActive_dyn[1], color = colors[i], label=lables[i] + ' simul', linestyle='--', lw=5)
            
            elif state == 'live_ri':
                
                ax.plot(simulObj[i].time_simul, simulObj[i].cellActive_dyn[2], color = colors[i], label=lables[i] + ' simul', linestyle='--', lw=5)
            
            elif state == 'dead':
                
               
              ax.plot(simulObj[i].time_simul, np.sum(simulObj[i].cellInactive_dyn, axis=0), color = colors[i], label=lables[i] + ' simul', linestyle='--', lw=5)
                
            else:
                
                ax.plot(simulObj[i].time_simul, simulObj[i].met_simul[simulObj[i].metabolome.metabolites.index(state)], color = colors[i], label=lables[i] + ' simul', linestyle='--', lw=5)
                #ax.set_ylim(0, 30)
                
           

    if legend:
        legend_properties = {'size': 18}
        # Set the location of the legend to be outside the plot area
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), prop=legend_properties)
        plt.subplots_adjust(right=0.75)
    else:
        ax.get_legend().remove()

    if stateType =='metabolite':
        ax.set_ylabel(state + ' (mM)', fontsize=26)
        plt.title(state, fontsize=32)
    
    elif stateType =='cells':
        ax.set_ylabel('$10^5$ cells/uL', fontsize=26)
        plt.title(state + ' cells', fontsize=32)
    else:
        ax.set_ylabel('pH', fontsize=26)
        plt.title('pH', fontsize=32)
    ax.set_xlabel('Time (h)', fontsize=26)
    
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])

    
    plt.tight_layout()

    return fig,ax


def makeKineticPlot(x,
                    y,
                    color,
                    legend,
                    xlabel,
                    ylabel,
                    title = None,
                    linestyle = 'o-'):
    '''
    '$10^5$ cells/uL'
    'Time (h)'
    'mM'

    '''
    
    #fig, ax = plt.subplots()
    
    
    plt.plot(x, y, color=color, linestyle = linestyle, lw=3, alpha = 1.0, label=legend)
            
    if legend is not None:
    
        legend_properties = {'size':18}
        plt.legend(fontsize=16, prop=legend_properties, bbox_to_anchor=(1.0, 1.0))
    
    ax = plt.gca()

    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xlabel(xlabel, fontsize=16)
    
    
    if title is not None:
        plt.title(title, fontsize=16)
    
    
    

    
    plt.tight_layout()
    


    


def plot_stacked_bar_charts(data1, labels1, colors1, ylabel1, data2, labels2, colors2, ylabel2, figPath = None):
    """
    Plots two bar charts, one on top of the other, with given data, labels, colors, and y-axis labels.
    """
    # Create figure and subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8)) # Adjusted for two stacked subplots
    
    # Plot first bar chart
    
    for i, (value, color) in enumerate(zip(data1, colors1)):
        ax1.bar(i + 1, value, color=color)
    
    ax1.set_xlim(0, len(labels1) + .5)
    ax1.set_ylabel(ylabel1, fontsize=15)
    ax1.set_xticks([])
    
    # Add legend for the first plot
    patches1 = [mpatches.Patch(color=color, label=label) for label, color in zip(labels1, colors1)]
    ax1.legend(handles=patches1, loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)

    # Plot second bar chart
    
    for i, (value, color) in enumerate(zip(data2, colors2)):
        ax2.bar(i + 1, value, color=color)
    ax2.set_xlim(0, len(labels2) + .5)
    ax2.set_ylabel(ylabel2, fontsize=15)
    ax2.set_xticks([])
    
    # Add legend for the second plot
    patches2 = [mpatches.Patch(color=color, label=label) for label, color in zip(labels2, colors2)]
    ax2.legend(handles=patches2, loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
    
    # Adjust layout
    plt.tight_layout()
    
    if figPath is not None:
        plt.savefig(figpath, dpi = 600)
    # Show plot
    plt.show()
            


# species = 'bh'
# experiments = ['bhbt', 'bhri', 'bhbtri']
# labels = ['bh1', 'bh2', 'bh3']
# colors = ['#00ff26', '#003eff', '#ff0000']


# mod1 = simulateExperiment(species, experiments[0], 'bh_bhbt.txt', 'modelDB_bhbtri_bh1 - Copy.sqlite3', ['live', 'trehalose', 'glucose', 'pyruvate', 'lactate','acetate'])

# mod2 = simulateExperiment(species, experiments[1], 'bh_bhri.txt', 'modelDB_bhbtri_bh2 - Copy.sqlite3', ['live', 'trehalose', 'glucose', 'pyruvate', 'lactate','acetate'])

# mod3 = simulateExperiment(species, experiments[2], 'bh_bhbtri.txt', 'modelDB_bhbtri_bh3 - Copy.sqlite3', ['live', 'trehalose', 'glucose', 'pyruvate', 'lactate','acetate'])


# makeExperimentPlot(species, 'pH', 'pH', ['bhri'], labels, colors)
# makeExperimentPlot(species, 'live', 'cells', ['bhri'], labels, colors)
# makeExperimentPlot(species, 'trehalose', 'metabolite', ['bhri'], labels, colors, simulObj = [mod2])
# makeExperimentPlot(species, 'glucose', 'metabolite', ['bhri'], labels, colors, simulObj = [mod2])
# makeExperimentPlot(species, 'pyruvate', 'metabolite', ['bhri'], labels, colors, simulObj = [mod2])
# makeExperimentPlot(species, 'lactate', 'metabolite', ['bhri'], labels, colors, simulObj = [mod2])
# makeExperimentPlot(species, 'acetate', 'metabolite', ['bhri'], labels, colors, simulObj = [mod2])
