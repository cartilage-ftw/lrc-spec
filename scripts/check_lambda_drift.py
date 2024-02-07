import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from plotting_prefs import *

dat = pd.read_csv('../data/wavelength/30.10.2023, 17.04,  15,798.01845 cm-1.lta', sep='\t',
                    skiprows=80, names=['time', 'wavelength'])
dat2 = pd.read_csv('../data/wavelength/30.10.2023, 17.55,  15,798.01848 cm-1.lta', sep='\t',
                    skiprows=80, names=['time', 'wavelength'])
dat4 = pd.read_csv('../data/wavelength/30.10.2023, 19.34,  15,798.01777 cm-1.lta', sep='\t',
                    skiprows=80, names=['time', 'wavelength'])

dual_monitor = pd.read_csv('../data/wavelength/08.11.2023, 20.23,  long-term recording 1.lta',
                           sep='\t', skiprows=80, names=['time', 'wave1', 'wave2'])

dual_monitor = dual_monitor.replace(-4, np.nan)
dual_monitor = dual_monitor.replace(-3, np.nan)
dat = dat.replace(-4, np.nan)
dat2 = dat2.replace(-4, np.nan)
dat4 = dat4.replace(-4, np.nan)

#convert to seconds
dat['time'] = dat['time']/1000 
dat2['time'] = dat2['time']/1000
dat4['time'] = dat4['time']/1000
dual_monitor['time'] = dual_monitor['time']/1000
#dat.dropna(inplace=True)

c=299_792_458 # m/s
dat['freq'] = (c/dat['wavelength'])#*1E6
dat2['freq'] = c/dat2['wavelength']
dat4['freq'] = c/dat4['wavelength']
dual_monitor['freq1'] = c/dual_monitor['wave1']
dual_monitor['freq2'] = c/dual_monitor['wave2']

#plt.hist(dat2['freq'], bins=100)
fig, axes = plt.subplots(2, 1, sharex=True, figsize=(6,6))

sigma_1 = np.std(dual_monitor['freq1'])*1E3#*1000
sigma_2 = np.std(dual_monitor['freq2'])*1E3#*1000

drate_1 = 1E3*np.abs(min(dual_monitor['freq1']) - max(dual_monitor['freq1']))/60
drate_2 = 1E3*np.abs(min(dual_monitor['freq2'].dropna()) - max(dual_monitor['freq2'].dropna()))/60
print('sigma in wave1:', sigma_1, 'MHz')
print('sigma in wave2', sigma_2, 'MHz' )
print(max(dual_monitor['freq2'].dropna()))

axes[0].plot(dual_monitor['time'], dual_monitor['wave1'], c='crimson',
              label=f'fundamental ($\Delta f={drate_1:.2f}$ MHz/min)')
axes[1].plot(dual_monitor['time'], dual_monitor['wave2'], c='cornflowerblue',
              label=f'second harmonic ($\Delta f={drate_2:.2f}$ MHz/min)')
'''plt.plot(dat['time'], dat['freq'],  label='2 hours')
plt.plot(dat2['time'], dat2['freq'], zorder=-1, label='3 hours')
plt.plot(dat4['time'], dat4['freq'], zorder=-2, label='4.5 hours')'''

for ax in axes:
    ax.legend()
    ax.set_ylabel('Wavelength [nm]')
    ax.ticklabel_format(style='plain', useOffset=False)

axes[1].set_xlabel('Time [s]')
plt.tight_layout()
plt.subplots_adjust(hspace=0)

#plt.tight_layout()
plt.show()