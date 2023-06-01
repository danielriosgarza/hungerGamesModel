# -*- coding: utf-8 -*-
"""
Created on Tue May 30 14:09:47 2023

@author: danie
"""
import os
import sys
from pathlib import Path

from aquarel import load_theme

theme = load_theme("boxy_light")
theme.apply()


import matplotlib.pyplot as plt

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'compare2experiments'))
from general import *


strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'ambr')

stateTable = parseTable(os.path.join(strainSummaryFolder, 'trehalose' + '.tsv'))

vessels = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11']

makeExperimentPlot('ambr', 'trehalose', 'metabolite', ['AMBR7'], ['trehalose'], ['#8900ff']*12)
plt.show()


fig, ax = plt.subplots()
species = 'ambr'
strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', species)

state = 'acetate'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7']
sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = '#003eff', lw=.666, label='acetate', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})

state = 'butyrate'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7']
sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = '#ff00a1', lw=.666, label='butyrate', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})
       
legend_properties = {'size':12}
ax.legend(fontsize=16, prop=legend_properties, bbox_to_anchor=(1.0, 1.0))

ax.set_ylabel(state + ' (mM)', fontsize=16)
plt.title(state, fontsize=16)


plt.tight_layout()

plt.show()


makeExperimentPlot('ambr', 'ri', 'cells', ['AMBR7'], ['AMBR7'], ['#00B8FF']*12)
plt.show()

makeExperimentPlot('ambr', 'bh', 'cells', ['AMBR7'], ['AMBR7'], ['#00B8FF']*12)
plt.show()



fig, ax = plt.subplots()
species = 'ambr'
strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', species)

state = 'ri'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7']
sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = '#003eff', lw=.666, label='Ri', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})

state = 'bh'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7']
sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = '#ff00a1', lw=.666, label='Bh', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})
       
legend_properties = {'size':12}
ax.legend(fontsize=16, prop=legend_properties, bbox_to_anchor=(1.0, 1.0))

ax.set_ylabel('$10^5$ cells/uL', fontsize=16)
plt.title('Bh & Ri', fontsize=16)


plt.tight_layout()

plt.show()

makeExperimentPlot('ambr7p', 'acetate', 'metabolite', ['AMBR7_P'], ['acetate'], ['#003eff'])
plt.show()

makeExperimentPlot('ambr7p', 'butyrate', 'metabolite', ['AMBR7_P'], ['butyrate'], ['#ff00a1'])
plt.show()


fig, ax = plt.subplots()
species = 'ambr7p'
strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', species)

state = 'acetate'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7_P']
sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = '#003eff', lw=.666, label='acetate', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})

state = 'butyrate'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7_P']
sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = '#ff00a1', lw=.666, label='butyrate', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})

plt.vlines(108, 0, 30, color='k', ls='--')
plt.vlines(120, 0, 30, color='k', ls='--')
       
legend_properties = {'size':12}
ax.legend(fontsize=16, prop=legend_properties, bbox_to_anchor=(1.0, 1.0))

ax.set_ylabel(state + ' (mM)', fontsize=16)
plt.title(state, fontsize=16)


plt.tight_layout()

plt.show()

makeExperimentPlot('ambr7p', 'ri', 'cells', ['AMBR7_P'], ['ri'], ['#00B8FF'])
plt.show()

makeExperimentPlot('ambr7p', 'bh', 'cells', ['AMBR7_P'], ['ri'], ['#FF10F0'])
plt.show()

fig, ax = plt.subplots()
species = 'ambr7p'
strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', species)

state = 'ri'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7_P']
#sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = '#003eff', lw=.666, label='Ri', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})

state = 'bh'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7_P']
sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = '#ff00a1', lw=.666, label='Bh', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})
plt.vlines(108, 0, 18, color='k', ls='--')
plt.vlines(120, 0, 18, color='k', ls='--')       
legend_properties = {'size':12}
ax.legend(fontsize=16, prop=legend_properties, bbox_to_anchor=(1.0, 1.0))

ax.set_ylabel('$10^5$ cells/uL', fontsize=16)
plt.title('Bh', fontsize=16)


species = 'ambr'
strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', species)

state = 'ri'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7']
#sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = '', lw=.666, label='Ri_control', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})

state = 'bh'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7']
sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = 'orange', lw=.666, label='Bh_control', alpha=1)#,


plt.tight_layout()

plt.show()




fig, ax = plt.subplots()
species = 'ambr7p'
strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', species)

state = 'acetate'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7_P']
sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = '#003eff', lw=.666, label='acetate', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})

species = 'ambr'
strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', species)

state = 'acetate'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7']
sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = 'purple', lw=.666, label='acetate control', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})

# state = 'butyrate'
# stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
# stateDF = getDFdict(stateTable, state, True)['AMBR7_P']
# sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = '#ff00a1', lw=.666, label='butyrate', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})

plt.vlines(108, 0, 30, color='k', ls='--')
plt.vlines(120, 0, 30, color='k', ls='--')
       
legend_properties = {'size':12}
ax.legend(fontsize=16, prop=legend_properties, bbox_to_anchor=(1.0, 1.0))

ax.set_ylabel(state + ' (mM)', fontsize=16)
plt.title(state, fontsize=16)


plt.tight_layout()

plt.show()


fig, ax = plt.subplots()
species = 'ambr7p'
strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', species)

state = 'butyrate'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7_P']
sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = '#ff00a1', lw=.666, label='butyrate', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})

species = 'ambr'
strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', species)

state = 'butyrate'
stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
stateDF = getDFdict(stateTable, state, True)['AMBR7']
sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = 'purple', lw=.666, label='butyrate control', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})

# state = 'butyrate'
# stateTable = parseTable(os.path.join(strainSummaryFolder, state + '.tsv'))
# stateDF = getDFdict(stateTable, state, True)['AMBR7_P']
# sns.lineplot(x='time', y=state, marker="o", markers=True, data=stateDF, color = '#ff00a1', lw=.666, label='butyrate', alpha=1)#, err_style='bars', err_kws = {'capsize':6, 'fmt':'o'})

plt.vlines(108, 0, 30, color='k', ls='--')
plt.vlines(120, 0, 30, color='k', ls='--')
       
legend_properties = {'size':12}
ax.legend(fontsize=16, prop=legend_properties, bbox_to_anchor=(1.0, 1.0))

ax.set_ylabel(state + ' (mM)', fontsize=16)
plt.title(state, fontsize=16)


plt.tight_layout()

plt.show()