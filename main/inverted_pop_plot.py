import pandas as pd
import io_utils
import numpy as np
import matplotlib.pyplot as plt

nov16_dir = '../data/lrc/16-11-2023 LRC/'
nov16_file = '2023-11-16-17-20-27.csv'#2023-11-16-17-20-27.csv'
nov16_background = '2023-11-16-16-50-32.csv'

nov16_100Hz = io_utils.read_lrc_timeseries(nov16_dir + nov16_file, buncher_cycle_dur=1E-2)
nov16_no_probe = io_utils.read_lrc_timeseries(nov16_dir + nov16_background)

fig, axes = plt.subplots(1,3, figsize=(9,6), sharex=True, sharey=True)
axes[2].hist(nov16_no_probe['cycle_time']*1E6, bins=200, color='royalblue', alpha=0.4, label='Ground State')
axes[2].hist(nov16_no_probe['cycle_time']*1E6, bins=200, fill=False, ec='snow', lw=0.3, alpha=1.0, )

axes[1].hist(nov16_100Hz['cycle_time']*1E6, bins=300, fill=False, ec='k', lw=0.3, alpha=1.0)
axes[1].hist(nov16_100Hz['cycle_time']*1E6, bins=300, color='#ffff00', alpha=0.4, label='Metastable State')

axes[0].hist(nov16_100Hz['cycle_time']*1E6, bins=300, color='#ffff00', alpha=0.4, label='Metastable State')
axes[0].hist(nov16_100Hz['cycle_time']*1E6, bins=300, fill=False, ec='k', lw=0.3, alpha=1.0)
axes[0].hist(nov16_no_probe['cycle_time']*1E6, bins=200, color='royalblue', alpha=0.4, label='Ground State')
axes[0].hist(nov16_no_probe['cycle_time']*1E6, bins=200, fill=False, ec='snow', lw=0.3, alpha=1.0)

axes[1].text(x=310, y=1000, s=r'$F=7/2\to5/2$' + '\n' + r'$\bar{\nu}=28502.46$ cm$^{-1}$')# +
axes[1].text(x=290, y=910, s='Bunching Freq: 100 Hz' + '\n80 shots per bunch')
for ax in axes:
    if ax != axes[0]:
        ax.axvline(x=0.28E3, ymin=0, ymax=1, lw=1.25, c='red', ls='--')
    ax.set_xlabel('Arrival Time [$\mu$s]')
    ax.legend()
plt.xlim(0.14E3, 0.45E3)
plt.tight_layout()
plt.subplots_adjust(wspace=0)
plt.savefig('population_inversion_atd.png', dpi=200)
plt.show()
#print(np.unique(nov16_100Hz['wavenum_req']))
#io_utils.make_spectrum(nov16_100Hz, ms_cut=0.29, filters=[0, 0], save_file=False)
