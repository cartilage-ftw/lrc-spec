"""
Done on the day the external trigger was set up.
"""

import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt

from plotting_prefs import *

dat = pd.read_csv('../data/power/31-10-2023.tsv', sep='\t',header=2)

dat['pow_corrected'] = dat['Power'] - 1.0
freq_color_dict = {1: 'deeppink', }
ls_dict = {10: '--', 1: '-'}
col_dict = {'I': 'cornflowerblue', 'E':'deeppink'}
trig_dict= {'I': 'internal', 'E': 'external'}
for trigger, df1 in dat.groupby('Trigger'):
        for freq, df2 in df1.groupby('Frequency'):
            plt.errorbar(df2['Current'], df2['pow_corrected'], yerr=df2['Power_err'],
                         ls=ls_dict[freq], c=col_dict[trigger], marker='o', barsabove=True,
                         capsize=3, label=f'{freq} kHZ, {trig_dict[trigger]}')
plt.legend()
plt.xlabel('Current [A]')
plt.ylabel('Power [mW]')
plt.tight_layout()
plt.show()