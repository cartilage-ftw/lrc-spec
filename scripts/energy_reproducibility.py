"""
dfsmflksdmkfsdmfmsdkm
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

from plotting_prefs import *

root_dir = '../data/power/09.11.2023/'
subdirs = ['beamsplitter_reflection/', 'before_window/']

#power = pd.read_csv('beamsplitter_reflection/56A_8kHz.csv', sep='\t', 
#                        skiprows=17)


def plot_laser_init():
    path = root_dir + subdirs[0]
    print(path)
    files = os.listdir(path)
    init_dfs = []
    for file in files:
        if '_initialization' in file: #name
            df = pd.read_csv(path + file, sep='\t', skiprows=17)
            init_dfs.append(df)

    fig, ax = plt.subplots(figsize=(6,6))

    start_time = 1E20 # initialized to a very large number
    for df in init_dfs:
        #h, m, s = df['Date']
        time_arr = []
        for i in range(len(df)):
            h, m, s = [float(x) for x in df.iloc[i]['Date'].split(',')[0].split(':')]
            time_arr.append(h*3600 + m*60 + s)
        df['time'] = time_arr
        start_time = min(start_time, np.min(df['time']))
    for df in init_dfs:
        df['start_time'] = df['time'] - start_time
        ax.plot(df['start_time'], df['Measurement']*1E6)
    plt.xlabel('Time [s]')
    plt.ylabel('Energy per pulse [$\mu$J]')
    plt.tight_layout()
    plt.show()


def make_time_arr(df):
    """
    The raw data reports timestamps with date, and hh:mm:ss time
    convert that into a time series starting at 0s

    This method appends a column called 'time_secs' to the passed DataFrame
    which contains that
    """
    h0, m0, s0 = [float(x) for x in df.iloc[0]['Date'].split(',')[0].split(':')]
    start_time = h0*3600 + m0*60 + s0
    time_arr = []
    for i in range(len(df)):
        h, m, s = [float(x) for x in df.iloc[i]['Date'].split(',')[0].split(':')]
        time = (h*3600 + m*60 + s) - start_time
        time_arr.append(time)
    df['time_secs'] = time_arr


def compare_day_var():
    """
    particularly for 8 kHz
    """
    power_nov8 = pd.read_csv('../data/power/08.11.2023_56A_8kHz.csv', sep='\t', 
                        skiprows=17)
    power_nov9 = pd.read_csv('../data/power/09.11.2023/beamsplitter_reflection/09.11.2023_56A_initialization_long_4.csv', sep='\t', 
                        skiprows=17)
    power_nov10 = pd.read_csv('../data/power/10.11.2023_init/10.11.2023_56A_8kHz_stable.csv', sep='\t',
                              skiprows=17)
    fig, axes = plt.subplots(3, 1, sharey=True, figsize=(7,5))
    dates = ['Nov 10', 'Nov 9', 'Nov 8']
    datum_pl = [power_nov10, power_nov9, power_nov8]
    for i, (ax, dat) in enumerate(zip(axes, datum_pl)):
        make_time_arr(dat)
        av = np.average(dat['Measurement'])*1E9
        sigma = np.std(dat['Measurement'])*1E9
        dat['1sig_below'] = dat['Measurement']*1E9 -sigma
        dat['1sig_above'] = dat['Measurement']*1E9+sigma
        ax.plot(dat['time_secs'], dat['Average Value']*1E9, marker='.', markersize=1,
                ls='', label=f'$\mu:{av:.1f}$nJ, $\sigma:{sigma:.1f}$nJ')
        ax.text(x=0, y=620, s=f'{dates[i]}')
        ax.fill_between(dat['time_secs'], y2=dat['1sig_above'],
                    y1=dat['1sig_below'], alpha=0.2)
        ax.legend(loc='upper right')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Energy Per Pulse [nJ]')
    plt.tight_layout()
    plt.show()


'''for subdir in subdirs:
    files = os.listdir(root_dir + subdir)
if '_initialization' not in file_name:'''

if __name__ == '__main__':
    #plot_laser_init()
    compare_day_var()