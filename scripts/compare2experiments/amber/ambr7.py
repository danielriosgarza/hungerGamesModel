import os
import sys
from pathlib import Path
figPath = os.path.join(Path(os.getcwd()).parents[2], 'files', 'Figures', 'ambr')


import pandas as pd
from aquarel import load_theme

theme = load_theme("boxy_light")
theme.apply()
import matplotlib.pyplot as plt
import numpy as np

# Assume we have some paired data for each group and condition.
dataCS1 = {
    'Bh_st1': [1.269845537,3.269944015,2.9678026,2.152795559,2.888866196,6.155168623,2.864436871,2.995948778,3.108646167,3.184716984,4.736401188,4.671900025],
    #'Bh_st2': [0.653875033,4.966910072,9.139935421,6.700265642,9.139935421,6.700265642],
    'Bh_st1_shift': [8.178329813,9.211629961,10.52377895,10.04647742,11.12835594,6.946610068,8.876324947,8.957434996,7.80869617,9.577699191,10.78959718,11.01501551],
    'Bt_st1': [7.362715689,7.883947777,9.34436142,9.606885087,9.00274262,18.3827473,9.554987566,14.38988005,9.641181493,7.962850506,9.709456246,8.861140179],
    #'Bt_st2': [18.28725588,18.01829605,9.856078069,12.12778676,9.856078069,12.12778676],
    'Bt_st1_shift': [11.25007155,9.861970123,8.662849726,10.18434076,9.545625748,4.931591588,14.41072622,10.52386681,8.375019927,14.53331247,11.26247945,10.76316734],
    'Ri_st1': [3.060500241,8.526911482,9.843814351,6.975876662,6.7400029,12.58619914,6.337852604,6.020393394,7.992966593,7.630271534,7.08465904,5.908916477],
    #'Ri_st2': [0.495324107,0.390617684,0.719634134,0.054050849,0.719634134,0.054050849]
    'Ri_st1_shift': [9.677690279,5.956884488,5.014326098,3.136599819,3.457476179,4.977637609,3.53377954,6.157676907,6.495434974,6.145248601,4.506173987,4.086962874],
}

dataCS2 = {
    'Bh_before': [9.218419638,9.477715722,9.783836833,6.436661141,23.64454511,6.928172562],
    #'Bh_st1': [8.178329813,9.211629961,10.52377895,10.04647742,11.12835594,6.946610068,8.876324947,8.957434996,7.80869617,9.577699191,10.78959718,11.01501551],
    'Bh_after': [7.076201297,7.967928476,7.218042663,5.801464789,10.27342672,10.02183422],
    #'Bt_before': [18.28725588,18.01829605,9.856078069,12.12778676,9.856078069,12.12778676],
    #'Bt_st1': [11.25007155,9.861970123,8.662849726,10.18434076,9.545625748,4.931591588,14.41072622,10.52386681,8.375019927,14.53331247,11.26247945,10.76316734],
    'Bt_before': [12.50235407,10.54881982,9.062824736,10.01231813,11.5633287,14.09714137],
    'Bt_after': [10.2478823,9.854724403,9.059353725,14.4904099,10.1078602,14.06164337],
    'Ri_before': [9.13224499,8.409204199,6.429643959,6.716487074,4.65650851,1.495002006],
    #'Ri_st1': [9.677690279,5.956884488,5.014326098,3.136599819,3.457476179,4.977637609,3.53377954,6.157676907,6.495434974,6.145248601,4.506173987,4.086962874],
    'Ri_after': [0.715779272,0.425115107,0.261476501,0.452779295,0.336179798,0.063978878]
}

groups = list(dataCS2.keys())
values = list(dataCS2.values())

fig, ax = plt.subplots()

# Define color palette
colors = ['#FF10F0', '#FF10F0', '#ff8300', '#ff8300', '#00B8FF', '#00B8FF']
#colors = ['#FF10F0', '#ff8300', '#00B8FF']

# Create an array for the position of each box. Increase space between different group pairs.
positions = [0, 0.5, 2, 2.5, 4, 4.5]
#positions = [0, 2, 4]

# Plot boxes. Decrease width for thinner boxes.
bp = ax.boxplot(values, positions=positions, notch=False, patch_artist=True, widths=0.4)

# Loop through each box and change color and edge color
for box, color in zip(bp['boxes'], colors):
    box.set(color=color, linewidth=2)
    box.set(facecolor=color) 

# Loop through each median and change color and line width
for median in bp['medians']:
    median.set(color='black', linewidth=3)

ax.set_xticks(positions)
ax.set_xticklabels(groups, fontsize = 16)

plt.xticks(rotation=45)  # Add this line to rotate x-axis labels.
plt.ylabel(r'$10^5$ cells $\mu l^{-1}$', fontsize = 16)  # Add this line to label y-axis.

#plt.title('Composition before and after perturbation')
plt.title('composition shift in chemostat', fontsize = 16)
plt.grid(axis='y')
plt.ylim(0,18.5)
plt.tight_layout()
#plt.savefig(os.path.join(figPath, 'ambr_states_biomass_shift.png'), dpi = 300)

plt.show()




# trehalose
dataCS1 = {
    'trehalose_st1': [0.017496121,0.021628529,0.023543418,0.018147376,0.016457069,0.015207539,0.019566645,0.015115879,0.01690426,0.016607393,0.022788369,0.023299737],
    'trehalose_st2': [0,0,0,0,0,0,0,0,0,0,0,0]
}
groups = list(dataCS1.keys())
values = list(dataCS1.values())
fig, ax = plt.subplots()
colors = ['#8900ff', '#8900ff']
positions = [0, 0.5]
# Plot boxes. Decrease width for thinner boxes.
bp = ax.boxplot(values, positions=positions, notch=False, patch_artist=True, widths=0.4)
# Loop through each box and change color and edge color
for box, color in zip(bp['boxes'], colors):
    box.set(color=color, linewidth=2)
    box.set(facecolor=color) 
    
# Loop through each median and change color and line width
for median in bp['medians']:
    median.set(color='black', linewidth=3)

ax.set_xticks(positions)
ax.set_xticklabels(groups)

plt.xticks(rotation=45, fontsize=16)  # Add this line to rotate x-axis labels.
plt.ylabel(r'mM', fontsize=16)  # Add this line to label y-axis.

plt.title('trehalose', fontsize=16)
plt.grid(axis='y')
plt.tight_layout()
plt.savefig(os.path.join(figPath, 'ambr_states_trehalose.png'), dpi = 300)
plt.show()

# butyrate
dataCS1 = {
    'butyrate_st1': [17.7546871,18.14298889,17.55617059,17.23113546,16.71458793,16.70777307,18.0679005,16.93999824,16.37319457,16.64245222,16.30454622,16.89635301],
    #'butyrate_st2': [12.80619951,13.64017356,14.24910662,11.65572223,10.48205991,11.43286864]
    'butyrate_st1_shift': [12.4748346,13.16629928,12.21266006,11.9914739,12.0336189,14.8915113,12.44713603,12.18398676,11.87508479,12.21627144,11.67489717,11.38928864]
}


dataCS2 = {
    'butyrate_before': [12.80619951,13.64017356,14.24910662,11.65572223,10.48205991,11.43286864],
    'butyrate_after': [10.04084973,9.566914812,8.611403299,8.80122715,8.353366555,7.625320008
]
}

groups = list(dataCS1.keys())
values = list(dataCS1.values())
fig, ax = plt.subplots()
colors = ['#ff00a1', '#ff00a1']
positions = [0, 0.5]
# Plot boxes. Decrease width for thinner boxes.
bp = ax.boxplot(values, positions=positions, notch=False, patch_artist=True, widths=0.4)
# Loop through each box and change color and edge color
for box, color in zip(bp['boxes'], colors):
    box.set(color=color, linewidth=2)
    box.set(facecolor=color) 
    
# Loop through each median and change color and line width
for median in bp['medians']:
    median.set(color='black', linewidth=3)

ax.set_xticks(positions)
ax.set_xticklabels(groups, fontsize=16)

plt.xticks(rotation=45)  # Add this line to rotate x-axis labels.
plt.ylabel(r'mM', fontsize=16)  # Add this line to label y-axis.

plt.title('butyrate shift', fontsize=16)
plt.grid(axis='y')
plt.tight_layout()
plt.savefig(os.path.join(figPath, 'ambr_states_butyrate_shift.png'), dpi = 300)
plt.show()

# acetate
dataCS1 = {
    'acetate_st1' :  [16.62311796,15.03746814,15.2302017,14.53527278,15.7555663,16.5739363,16.32970215,14.76018823,15.8471278,14.35490399,18.38610398,17.67909154],
    #'acetate_st2' : [13.98989495,14.75380638,16.71874716,24.37206401,24.79630629,26.87621888]
    'acetate_st1_shift': [20.90242177,19.87053295,20.42851759,21.51884275,19.57573955,13.87865792,18.3317889,18.60894731,20.12968519,19.07194562,22.00452719,22.47934977]
}
groups = list(dataCS1.keys())
values = list(dataCS1.values())
fig, ax = plt.subplots()
colors = ['#003eff', '#003eff']
positions = [0, 0.5]
# Plot boxes. Decrease width for thinner boxes.
bp = ax.boxplot(values, positions=positions, notch=False, patch_artist=True, widths=0.4)
# Loop through each box and change color and edge color
for box, color in zip(bp['boxes'], colors):
    box.set(color=color, linewidth=2)
    box.set(facecolor=color) 
    
# Loop through each median and change color and line width
for median in bp['medians']:
    median.set(color='black', linewidth=3)

ax.set_xticks(positions)
ax.set_xticklabels(groups, fontsize=16)

plt.xticks(rotation=45)  # Add this line to rotate x-axis labels.
plt.ylabel(r'mM', fontsize=16)  # Add this line to label y-axis.

plt.title('acetate shift', fontsize=16)
plt.grid(axis='y')
plt.tight_layout()
plt.savefig(os.path.join(figPath, 'ambr_states_acetate_shift.png'), dpi = 300)
plt.show()