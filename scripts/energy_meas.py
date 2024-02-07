import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from plotting_prefs import *

root = '../data/power/power-03.11.2023/'

dirs = ['before_window/', 'beamsplit_reflect/']

data_all = pd.DataFrame(columns=['pulse_freq', 'current', 'avg_energy', 'energy_std', 
                        'avg_power', 'Detector Position'])
for dir in dirs:
    files = os.listdir(root + dir)
    #count = 0
    
    for f in files:
        df = pd.read_csv(root + dir + f, sep='\t', skiprows=17)
        # there's typically one entry that contains the mean, std, etc.
        # all other lines have those entries as zero.
        stats = df[~df['Std Deviation'].isna()]

        # in some files, these aren't given, so I'm calculating them myself
        std = np.std(df['Measurement'])
        avg_energy = np.average(df['Measurement'])
        avg_power = np.sum(df['Measurement'])/float(df['Second'].tail(1))

        # [:-4] to ignore the '.csv' file extension
        name_comps = f[:-4].split('_')
        current = int(name_comps[1][:-1])
        freq_khz = int(name_comps[2][:-3])
        
        stats_row = [freq_khz, current, avg_energy, std, avg_power, dir]
        data_all.loc[len(data_all)] = stats_row

"""
Plot the data
"""

data_all = data_all.sort_values('current')

fig, axes = plt.subplots(1,2, figsize=(8,4), sharex=True)

for detect_pos, df in data_all.groupby('Detector Position'):
    for freq_khz, df2 in df.groupby('pulse_freq'):
        ax_ind = dirs.index(detect_pos)
        
        axes[ax_ind].errorbar(df2['current'], df2['avg_energy']*1E6, yerr=df2['energy_std']*1E6,
                   ls='', marker='o', alpha=0.8, capsize=2, label=f'{freq_khz}kHz')
        #axes[ax_ind].plot(df2['current'], df2['avg_energy']*1E6, #yerr=df2['energy_std']*1E6,
        #           ls='-', marker='o', label=f'{freq_khz}kHz')
        #axes[ax_ind].plot(df2['current'], df2['avg_power']*1E6, 'o-', label=f'{freq_khz}kHz',
        #            ms=8) #marker='o', s=16, )

'''fig, ax = plt.subplots(figsize=(6,5))
for freq_khz, df in data_all.groupby('pulse_freq'):
    curr_vals = []
    y_vals = []
    #y_errs = []
    df1 = df[df['Detector Position'] == 'before_window/']
    df2 = df[df['Detector Position'] != 'before_window/']
    for current in range(51, 59):
        print(df2['current'])
        if current in list(df2['current']):
            curr_vals.append(current)
            vals1 = df1[df1['current']==current]#['avg_energy']
            vals2 = df2[df2['current']==current]#['avg_energy']
            #print(vals1['avg_energy'])
            print('****')
            print('vals1', vals1)
            print('*****')
            print('vals2', vals2)
            y_vals.append(vals1.iloc[0]['avg_energy']/vals2.iloc[0]['avg_energy'])
            #y_errs.append(float(np.sqrt(vals1['energy_std']**2 + vals2['energy_std']**2)))
    ax.plot(curr_vals, y_vals, 'o-', ms=8, label=f'{freq_khz}kHz')
    #print('Current values', curr_vals)
    #print(y_vals)
    #print('---')
plt.legend()
plt.xlabel('Current [A]')
plt.ylabel('Ratio')'''
    

titles = ['Measured before chamber window', 'right after beamsplitter']

'''for i, ax in enumerate(axes):
    ax.legend()
    ax.set_xlabel('Current [A]')
    ax.set_ylabel('Energy Per Pulse [$\mu$J]')
    #ax.set_ylabel('Average Power [$\mu$W]')
    ax.title.set_text(titles[i])'''

data_all.to_csv('energy_meas_03Nov.csv', index=False)
plt.tight_layout()
plt.show()