import numpy as np
import pandas as pd

import io_utils
import ion_cloud_extent
import matplotlib.pyplot as plt

from line_fitting import fit_gaussian
from scipy.fft import fft, fftfreq


dir_path = '../data/lrc/16-11-2023 LRC/'
bg1_file = '2023-11-16-16-53-31.csv'
bg2_file = '2023-11-16-16-50-32.csv'

dir_path_nov22 = '../data/lrc/Ion Cloud 22-11-2023/'

bg_nov22 = io_utils.read_lrc_timeseries(dir_path_nov22 + ion_cloud_extent.clean_bg_nov22 + '.csv')
#bg2 = io_utils.read_lrc_timeseries(dir_path + bg1_file)
bg1 = io_utils.read_lrc_timeseries(dir_path + bg2_file)

y_vals, bins = np.histogram(bg_nov22['cycle_time']*1E3, bins='stone')
bin_centers = (bins[:-1] + bins[1:])/2

bin_width = bins[1] - bins[0]
fit = fit_gaussian(bg_nov22)
print(fit.fit_report())
plt.show()


plt.show()