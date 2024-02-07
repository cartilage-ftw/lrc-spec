"""
To see what line profiles describe Kim's Lu measurements best
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from lmfit.models import LinearModel, GaussianModel, VoigtModel

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
plt.rcParams['figure.dpi'] = 150


lu_data = pd.read_csv('../data/lu/dat', sep='\t', names=['wavenumber', 'x_err', 'y', 'y_err'])

# Model
offset = LinearModel(prefix='base_')
# hyperfine multiplets
hf1 = VoigtModel(prefix='F1_')
hf2 = VoigtModel(prefix='F2_')
hf3 = VoigtModel(prefix='F3_')

pars = hf1.make_params(center=28502.379, amplitude=0.8, sigma=0.085)#gamma=0.2, 
pars.update(hf2.make_params(center=28503.01, amplitude=0.876,sigma=0.106)) # gamma=0.25, 
pars.update(hf3.make_params(center=28503.72, amplitude=0.94,sigma=0.136))#  gamma=0.32, 
pars.update(offset.make_params(slope=0, intercept=0.04))

pars['base_intercept'].set(max=0.1)
# composite model
composite = hf1 + hf2 + hf3 + offset
# fit results
fit = composite.fit(lu_data['y'], pars, x=lu_data['wavenumber'])

print(dir(fit.result))
print('Chi Square:', fit.result.chisqr)

fig, axes = plt.subplots(2,1, figsize=(5,6), height_ratios=[4,1], sharex=True)

for param, val in fit.values.items():
    print(param, val)

axes[0].plot(lu_data['wavenumber'], lu_data['y'], marker='o', ls='', c='deeppink')
#axes[0].plot(lu_data['wavenumber'], fit.best_fit, c='dimgray')
finer_grid = np.linspace(min(lu_data['wavenumber']), max(lu_data['wavenumber']), 1000)
axes[0].plot(finer_grid, fit.eval(x=finer_grid), c='dimgray')


# plot residuals
axes[1].plot(lu_data['wavenumber'], fit.result.residual, c='pink')
axes[1].set_xlabel("Wavenumber [cm$^{-1}$]")
axes[1].ticklabel_format(style='plain', useOffset=False)
plt.tight_layout()
plt.show()

'''first_peak = lu_data[lu_data['wavenumber'] < 28502.55]
g1 = hf1 + offset
refit = g1.fit(first_peak['y'], pars, x=first_peak['wavenumber'])
axes[0].plot(finer_grid, refit.eval(x=finer_grid), ls='--')'''