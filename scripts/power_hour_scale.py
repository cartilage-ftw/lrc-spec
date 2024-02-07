"""
Power measurements carried out for a few hours on Nov 8th, 2023
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from plotting_prefs import *

power = pd.read_csv('../data/power/08.11.2023_56A_8kHz.csv', sep='\t', 
                        skiprows=17)

start_point = power.head(1)['Date']
h0, m0, s0 = [float(x) for x in start_point[0].split(',')[0].split(':')]
# convert into s
start_time = h0*3600 + m0*60 + s0

time_arr = []

for i in range(len(power)):
    h, m, s = [float(x) for x in power.iloc[i]['Date'].split(',')[0].split(':')]
    time_instant = h*3600 + m*60 + s
    time = time_instant - start_time
    time_arr.append(time)

power['time'] = time_arr

power_averaged = power[~power['Average Power'].isna()]
fig, ax = plt.subplots(figsize=(6,6))
ax.plot(power_averaged['time'], power_averaged['Average Value']*1E6, label='mean')

# display uncertainty via shadow/fill
power_averaged['1sig_below'] = (power_averaged['Average Value']-power_averaged['Std Deviation'])*1E6
power_averaged['1sig_above'] = (power_averaged['Average Value']+power_averaged['Std Deviation'])*1E6

ax.fill_between(power_averaged['time'], y2=power_averaged['1sig_above'],
                    y1=power_averaged['1sig_below'], alpha=0.2, label='$1\sigma$')

avg_dev = np.average(power_averaged['Std Deviation'])*1E6
avg_dev_all = np.std(power['Measurement'])*1E6
mean_erg = np.average(power_averaged['Average Value'])*1E6
percent_dev = 100*np.float(avg_dev/mean_erg)
long_term_drift = np.std(power_averaged['Average Value'])*1E6
print(mean_erg, percent_dev)
print('Long long_term_drift:', 100*long_term_drift/mean_erg)
print('avg: ', avg_dev, 'avg_all', avg_dev_all)
ax.legend()
plt.text(x=250, y=0.375, s=f'YAG current: 56A\npulse freq: 8kHz\n' + 
                r'$\sigma$: ' + f'{avg_dev:.2f} $\mu$J ($\pm${percent_dev:.2f}\%)', ha='left')
plt.xlabel('Time [s]')
plt.ylabel('Energy Per Pulse [$\mu$J]')
plt.tight_layout()
plt.show()
