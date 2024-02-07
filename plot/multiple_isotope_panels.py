import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plotting_prefs import *

ref_freq = 28502.348 * 29.99 # cm-1


def get_freq_offset(wavenum):
    return wavenum*29.99 - ref_freq

def normalize_ms_fraction(df):
    df['Norm MS Fraction'] = (df['MS Fraction'])/np.std(df['MS Fraction']) #  - np.average(df['MS Fraction'])
    df['Freq GHz'] = df['Wavenumber']*29.99
    df['Freq Offset'] = df['Wavenumber']*29.99 - ref_freq

lu176_spec = pd.read_csv("../main/lu176_nov27_hfs.csv")
normalize_ms_fraction(lu176_spec)

lu175_spec = pd.read_csv('../main/lu175_nov27_hfs.csv')
normalize_ms_fraction(lu175_spec)


what = lu175_spec[lu175_spec['Wavenumber'] < 28502.149] # left to the center of the F=6 multiplet in 176
max_lu175 = np.max(what['MS Fraction'])

print('Max fraction', max_lu175)

fig, axes = plt.subplots(2, 1, sharex=True, figsize=(6,6))

axes[0].plot(lu175_spec['Freq Offset'], lu175_spec['MS Fraction'], marker='D', c='dimgray', ms=4,
         mfc='#ff557fff', mec='#ff557fff', label=r'$^{175}\rm{Lu}$')
axes[1].plot(lu176_spec['Freq Offset'], lu176_spec['MS Fraction'], marker='o', c='#767676', ms=4.5,
        mfc='#aaaaff', mec='#aaaaff', label=r'$^{176}\rm{Lu}$')

#ax.axvline(x=get_freq_offset(28502.149), ymin=0, ymax=1, c='darkgray', zorder=-2, ls='--', lw=0.75)
#ax.axhline(y=max_lu175, xmin=0, xmax=1, c='darkgray', zorder=-2, ls='--', lw=0.75)
for ax in axes:
    ax.tick_params(which='major', axis='both', length=6, labelsize=13)
#ax.text(x=-14, y=26.5, s=f'${get_freq_offset(28502.149):.1f}$ GHz', c='#aaaaff')
#ax.text(x=-27, y=2.5, s=f'${max_lu175:.2f}$', c='#ff557fff')
    ax.set_xlabel('Frequency Detuning [GHz]', fontsize=13)
    ax.axvline(x=0, ymin=0, ymax=1, lw=0.75, zorder=-1, ls='--',c='dimgray')
    ax.set_ylabel('Scaled MS Fraction [\%]', fontsize=13)
    ax.legend(fontsize=13, loc='upper left')
plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.savefig('isotope_shift.png', dpi=300)
plt.show()