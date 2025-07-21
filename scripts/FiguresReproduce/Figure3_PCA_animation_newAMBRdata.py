# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 11:09:27 2023

@author: drgarza
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Ellipse
import numpy as np


from aquarel import load_theme

theme = load_theme("boxy_light")
theme.apply()

#from sklearn.manifold import MDS as PCA
from sklearn.covariance import MinCovDet
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
pca = PCA(n_components=2, svd_solver ='randomized', power_iteration_normalizer='QR', iterated_power=10000)





def get_data(fileName):
    
    strainSummaryFolder = os.path.join(Path(os.getcwd()).parents[1], 'files', 'strainSummaries', 'ambr')
    
    
    
    with open(os.path.join(strainSummaryFolder, fileName)) as f:
        labels = f.readline().strip().split('\t')[6::]
        experiments = []
        conditions = []
        description = []
        time = []
        timeCondition = []
        replicate = []
        exclude= [ ]
        data = {i:[] for i in labels if i not in exclude}
        
        for line in f:
            a = line.strip().split('\t')
            experiments.append(a[0])
            conditions.append(a[3])
            description.append(a[2])
            
            for i in range(len(labels)):
                if labels[i] not in exclude:
                    data[labels[i]].append(float(a[i+6]))
            time.append(float(a[4]))
            timeCondition.append(float(a[5]))
            
            ambrExp = a[0].split('_')
            replicate.append(ambrExp[0] + '_' + ambrExp[1] + '_' + a[4])
            
    
    dataM = np.array([np.array(data[i]) for i in labels if i not in exclude]).T
    return labels, experiments, description, conditions, time, timeCondition, replicate, dataM









#############################HeatMap###############
labels, experiments, description, conditions, time, timeCond, replicate, dataM = get_data('ambrSS.txt')


data = dataM.T
fig, ax = plt.subplots()
heatmap = ax.imshow(data, aspect='auto', vmin = 0, vmax =12)
ax.grid(False)
ax.set_yticks(np.arange(data.shape[0]))
yticklabels = ax.set_yticklabels(labels)


# Set specific x-ticks and labels (x-axis)


ax.set_xticks([0.0, 13.0, 22.0, 28.0])
ax.set_xticklabels([0.0, 13.0, 22.0, 28.0])

    # # Set labels and title
    # ax.set_xlabel(xlabel)
    # ax.set_ylabel(ylabel)
    # ax.set_title(title)
    
# Add grid
ax.set_xticks(np.arange(-.5, data.shape[1], 1), minor=True)
ax.set_yticks(np.arange(-.5, data.shape[0], 1), minor=True)
ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5)
ax.tick_params(which="minor", size=0)

# Add colorbar
cbar = ax.figure.colorbar(heatmap, ax=ax)
plt.tight_layout()

# Save the figure
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'composition_heatMap.png')

plt.savefig(fileName, dpi=600)






labels, experiments, description, conditions, time, timeCond, replicate, dataM = get_data('ambrAll.txt')



##################3d species#########################
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
from pathlib import Path

# Assuming 'description', 'dataM', and 'timeCond' are already defined

descCol = []
counted = []
counter = 0
used = []
for i, v in enumerate(description):
    if (v != 'ongoing pH perturbation') and (v != 'ongoing feed perturbation') and (v != 'post pH perturbation'):
        used.append(i)
        if v not in counted:
            descCol.append(counter)
            counted.append(v)
            counter += 1
        else:
            descCol.append(counted.index(v))

descCol = np.array(descCol)
used = np.array(used)
timeCond = np.array(timeCond)

# Add jitter to the data
def add_jitter(arr, scale=0.5):
    return arr + np.random.uniform(-scale, scale, size=arr.shape)

x = add_jitter(dataM.T[-1][used])
y = add_jitter(dataM.T[-2][used])
z = add_jitter(dataM.T[-3][used])

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

cols = ['#00FF00', '#FF007F', '#FFD700']
markers = ['o', '^', 's']  # Circle, triangle_up, square
labels = ['Control', 'Post pH/feed perturbation', 'Post feed perturbation']

for i in range(3):
    idx = descCol == i
    # Larger markers with timeCond mapping
    scatter1 = ax.scatter(x[idx], y[idx], z[idx],
                          s=150,
                          c=timeCond[used][idx],
                          cmap='binary',
                          lw=0,
                          zorder=i,
                          alpha=0.7,
                          marker=markers[i], 
                          vmin=0,
                          vmax=200)
    # Smaller markers with fixed color and marker shape
    scatter2 = ax.scatter(x[idx], y[idx], z[idx],
                          s=40,
                          c=cols[i],
                          marker=markers[i],
                          edgecolors='k',
                          linewidths=0.5,
                          zorder=i+1,
                          alpha=1.0,
                          label=labels[i])

# Labels
ax.set_xlabel("R. intestinalis", fontsize=20)
ax.set_ylabel("B. thetaiotaomicron", fontsize=20)
ax.set_zlabel("B. hydrogenotrophica", fontsize=20)
ax.view_init(elev=45, azim=45)

# Add legend
ax.legend()

# Add colorbar for time information
mappable = plt.cm.ScalarMappable(cmap='binary')
mappable.set_array(np.arange(200))#(timeCond[used])
cbar = plt.colorbar(mappable, ax=ax, pad=0.1)
cbar.set_label('time (h)', fontsize=20)

plt.tight_layout()

fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'multistability', 'ambr_path', 'composition.png')

#plt.savefig(fileName, dpi=600)
plt.show()



# #Background PCA

import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os
from pathlib import Path

# Assuming 'description', 'dataM', 'timeCond', and 'used' are already defined

# Prepare the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(dataM)

pca = PCA(n_components=2)
transformed_data = pca.fit_transform(scaled_data)

x_pc = transformed_data[:, 0][used]
y_pc = transformed_data[:, 1][used]

explained_variance = pca.explained_variance_ratio_

# Define colors, markers, and labels
cols = ['#00FF00', '#FF007F', '#FFD700']  # Bright green, pink, gold
markers = ['o', '^', 's']  # Circle, triangle_up, square
labels = ['Control', 'Post pH/feed perturbation', 'Post feed perturbation']

# Function to add jitter
def add_jitter(arr, scale=0.5):
    return arr + np.random.uniform(-scale, scale, size=arr.shape)

# Add jitter to the PCA components
x_pc_jittered = add_jitter(x_pc)
y_pc_jittered = add_jitter(y_pc)

# Create the plot
fig, ax = plt.subplots(figsize=(8, 6))

for i in range(3):
    idx = descCol == i
    # Larger markers with timeCond mapping
    scatter1 = ax.scatter(x_pc_jittered[idx], y_pc_jittered[idx],
                          s=150,
                          c=timeCond[used][idx],
                          cmap='binary',
                          lw=0,
                          zorder=i,
                          alpha=0.7,
                          marker=markers[i], 
                          vmin=0,
                          vmax=200)
    # Smaller markers with fixed color and marker shape
    scatter2 = ax.scatter(x_pc_jittered[idx], y_pc_jittered[idx],
                          s=40,
                          c=cols[i],
                          marker=markers[i],
                          edgecolors='k',
                          linewidths=0.5,
                          zorder=i+1,
                          alpha=1.0,
                          label=labels[i])
    
    
    

# Labels
ax.set_xlabel(f"PC1 ({explained_variance[0]*100:.2f}% variance explained)", fontsize=20)
ax.set_ylabel(f"PC2 ({explained_variance[1]*100:.2f}% variance explained)", fontsize=20)

# Add legend
ax.legend()

# Add colorbar for time information
mappable = plt.cm.ScalarMappable(cmap='binary')
mappable.set_array(np.arange(200))#(timeCond[used])
cbar = plt.colorbar(mappable, ax=ax, pad=0.1)
cbar.set_label('time (h)', fontsize=20)

plt.tight_layout()

# Save the figure
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'composition_pca.png')

#plt.savefig(fileName, dpi=600)
plt.show()




#################################feed Pertubation##############

labels, experiments, description, conditions, time, timeCond, replicate, dataM = get_data('ambrAll.txt')


controlInd = [i for i in range(len(description)) if description[i]=="control conditions" and timeCond[i]>100]
feedInd = [i for i in range(len(description)) if description[i]=="post feed perturbation" and timeCond[i]>100]

#########cells########

#Bh
bhC = dataM.T[labels.index('Blautia hydrogenotrophica (10^5 cells/uL)')][controlInd]
bhF = dataM.T[labels.index('Blautia hydrogenotrophica (10^5 cells/uL)')][feedInd]


#Bt
btC = dataM.T[labels.index('Bacteroides thetaiotaomicron (10^5 cells/uL)')][controlInd]
btF = dataM.T[labels.index('Bacteroides thetaiotaomicron (10^5 cells/uL)')][feedInd]


#Ri
riC = dataM.T[labels.index('Roseburia intestinalis (10^5 cells/uL)')][controlInd]
riF = dataM.T[labels.index('Roseburia intestinalis (10^5 cells/uL)')][feedInd]

#

# Calculating mean and standard deviation for each group
data_means = [np.mean(bhC), np.mean(bhF), np.mean(btC), np.mean(btF), np.mean(riC), np.mean(riF)]
data_stds = [np.std(bhC), np.std(bhF), np.std(btC), np.std(btF), np.std(riC), np.std(riF)]

# Labels for each condition
labels = [
    'Blautia hydrogenotrophica (Control)', 'Blautia hydrogenotrophica (Feed)',
    'Bacteroides thetaiotaomicron (Control)', 'Bacteroides thetaiotaomicron (Feed)',
    'Roseburia intestinalis (Control)', 'Roseburia intestinalis (Feed)'
]

# Grouping data into pairs for easier side-by-side comparison
# Pairing labels and values for control and feed conditions for each species
group_labels = ['Blautia hydrogenotrophica', 'Bacteroides thetaiotaomicron', 'Roseburia intestinalis']
control_means = [np.mean(bhC), np.mean(btC), np.mean(riC)]
control_stds = [np.std(bhC), np.std(btC), np.std(riC)]
feed_means = [np.mean(bhF), np.mean(btF), np.mean(riF)]
feed_stds = [np.std(bhF), np.std(btF), np.std(riF)]

# Defining the width and positions for each pair of bars
x = np.arange(len(group_labels))  # label locations
width = 0.35  # width of the bars

# Creating the bar plot with side-by-side grouping
plt.figure(figsize=(8, 6))
plt.bar(x - width/2, control_means, width, yerr=control_stds, label='control', capsize=5, color = 'g')
plt.bar(x + width/2, feed_means, width, yerr=feed_stds, label='feed perturbation', capsize=5, color = '#FFD700')

# Adding labels and title
plt.xticks(x, group_labels)
plt.ylabel('Mean Cells (10^5 cells/uL)')
plt.title('Concentrations: Control vs feed perturbation')
plt.legend(title='Condition')
plt.tight_layout()
plt.ylim(0,30)
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'cells.png')

plt.savefig(fileName, dpi=300)


plt.show()


########relevant metabs####
labels, experiments, description, conditions, time, timeCond, replicate, dataM = get_data('ambrAll.txt')
#Bh
#bhC = dataM.T[labels.index('Formate (mM)')][controlInd]
#bhF = dataM.T[labels.index('Formate (mM)')][feedInd]


#Bt
aceC = dataM.T[labels.index('Acetate (mM)')][controlInd]
aceF = dataM.T[labels.index('Acetate (mM)')][feedInd]


#Ri
butC = dataM.T[labels.index('Butyrate (mM)')][controlInd]
butF = dataM.T[labels.index('Butyrate (mM)')][feedInd]

#

# Calculating mean and standard deviation for each group
data_means = [np.mean(aceC), np.mean(aceF), np.mean(butC), np.mean(butF)]
data_stds = [np.std(aceC), np.std(aceF),np.std(butC), np.std(butF)]

# Labels for each condition
labels = [
    'acetate (control)', 'acetate (feed pert.)',
    'butyrate (control)', 'butyrate (feed pert.)'
]

# Grouping data into pairs for easier side-by-side comparison
# Pairing labels and values for control and feed conditions for each species
group_labels = ['Acetate', 'Butyrate']
control_means = [np.mean(aceC), np.mean(butC)]
control_stds = [np.std(aceC), np.std(butC)]
feed_means = [np.mean(aceF), np.mean(butF)]
feed_stds = [np.std(aceF), np.std(butF)]

# Defining the width and positions for each pair of bars
x = np.arange(len(group_labels))  # label locations
width = 0.35  # width of the bars

# Creating the bar plot with side-by-side grouping
plt.figure(figsize=(8, 6))
plt.bar(x - width/2, control_means, width, yerr=control_stds, label='control', capsize=5, color = 'g')
plt.bar(x + width/2, feed_means, width, yerr=feed_stds, label='feed perturbation', capsize=5, color = '#FFD700')

# Adding labels and title
plt.xticks(x, group_labels)
plt.ylabel('mM')
#plt.title('Concentrations: Control vs feed perturbation')
plt.legend(title='Condition')
plt.tight_layout()
plt.ylim(0,30)
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'prodMet.png')

plt.savefig(fileName, dpi=300)

plt.show()

########treh####
labels, experiments, description, conditions, time, timeCond, replicate, dataM = get_data('ambrAll.txt')
#Bh
#bhC = dataM.T[labels.index('Formate (mM)')][controlInd]
#bhF = dataM.T[labels.index('Formate (mM)')][feedInd]


#Bt
treC = dataM.T[labels.index('Trehalose (mM)')][controlInd]
treF = dataM.T[labels.index('Trehalose (mM)')][feedInd]


# Calculating mean and standard deviation for each group
data_means = [np.mean(treC), np.mean(treF)]
data_stds = [np.std(treC), np.std(treF)]

# Labels for each condition
labels = [
    'trehalose (control)', 'trehalose (feed pert.)'
]

# Grouping data into pairs for easier side-by-side comparison
# Pairing labels and values for control and feed conditions for each species
group_labels = ['Trehalose']
control_means = [np.mean(treC)]
control_stds = [np.std(treC)]
feed_means = [np.mean(treF)]
feed_stds = [np.std(treF)]

# Defining the width and positions for each pair of bars
x = np.arange(len(group_labels))  # label locations
width = 0.35  # width of the bars

# Creating the bar plot with side-by-side grouping
plt.figure(figsize=(8, 6))
plt.bar(x - width/2, control_means, width, yerr=control_stds, label='control', capsize=5, color = 'g')
plt.bar(x + width/2, feed_means, width, yerr=feed_stds, label='feed perturbation', capsize=5, color = '#FFD700')

# Adding labels and title
plt.xticks(x, group_labels)
plt.ylabel('mM')
#plt.title('Concentrations: Control vs feed perturbation')
plt.legend(title='Condition')
plt.tight_layout()
plt.ylim(0,0.2)
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures',
                        'multistability', 'ambr_path', 'treh.png')

plt.savefig(fileName, dpi=300)

plt.show()






plt.show()

labels, experiments, description, conditions, time, timeCond, replicate, dataM = get_data('ambrAll.txt')