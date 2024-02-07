import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from lmfit.models import GaussianModel


plt.rcParams['figure.dpi'] = 150
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']

ra_decay = pd.read_csv('../data/2023-10-06-16-34-09.csv',
                       # I forgot which column was what, but only the first is relevant
                    names=['volts', 'col2', 'col3', 'col4'])

# volt (=energy, related by a scale) of each detection
count_pos = ra_decay['volts']

calib_factor = 1

# currently assumes 1V <-> 1 MeV (which may not be true, but that can be calibrated.)
ra_decay['energy'] = ra_decay['volts']* calib_factor

print("Number of counts detected:", len(count_pos))

ra_counts = ra_decay[ra_decay['volts'] > 5.0]
ra_counts = ra_counts[ra_counts['volts'] < 5.38]
print(np.sort(ra_counts['volts']))
print('Counts coming from Ra decays:', ra_counts)


bin_width = 0.0025 # in volts

num_bins = len(np.arange(min(count_pos), max(count_pos)+bin_width, bin_width))
vals, bin_edges = np.histogram(count_pos, bins=num_bins)

print(bin_edges[-2:])
print(max(count_pos)+bin_width)
print(max(count_pos))
#plt.hist(count_pos, bins=bin_edges, ec='k', fc='w') #fc='w', ec='k',

ra_peaks = [5.044, 5.136, 5.198, 5.2964, 5.3258]
ra_known = [5.43, 5.54, 5.61, 5.72, 5.75]

calib_fit = np.polyfit(ra_peaks, ra_known, deg=1)
print('x and constant coefficients', calib_fit)
print(ra_known)
print(np.array(ra_peaks)*calib_fit[0] + calib_fit[1])
print(np.array(ra_peaks) + 0.4)

plt.plot(ra_peaks, ra_known, 'o-')
#calibrate(ra_decay, bin_width)

fig, ax = plt.subplots(figsize=(6,6))
#ax.minorticks_on()
ax.tick_params(axis='both', which='major', direction='out', length=6)
#ax.tick_params(axis='both', which='minor', direction='out', length=4)
plt.stairs(vals, bin_edges*calib_fit[0]+calib_fit[1], ec='k', lw=0.5) #
plt.xlabel('Energy [MeV]')
plt.ylabel(f'Counts [{1000*bin_width*calib_fit[0]:.2f}' +  ' keV$^{-1}$]')
'''plt.text(5.6, 110, s='$^{223}$Ra', ha='center', va='bottom')

ymins = np.array([2.7, 12.9, 29.4, 94/1.41, 14.5])/100 + 0.05
for i, x in enumerate(ra_known):
    plt.axvline(x, ymin=ymins[i], ymax=1.05/1.41, lw=0.25, c='r')'''
plt.tight_layout()
plt.show()




'''def calibrate(data, bin_width):
    vals, bin_edges = np.histogram(data['volts'],
            bins=len(np.arange(min(data['volts']), max(data['volts']+bin_width), bin_width)))
    bin_means = (bin_edges[:-1] + bin_edges[1:])/2
    # now (bin_means, vals) are an (x,y) tuple list that we can fit Gaussians to.
    
    guess_means = [5.04, 5.133, 5.197, 5.293, 5.328]
    guess_amp = [4, 28, 62, 169, 36.3]

    gauss_multi = [GaussianModel(prefix=f'g{i}_') for i in range(5)]

    pars = gauss_multi[0].make_params()
    for i in range(1,5):
        pars.update(gauss_multi[i].make_params())
    for i in range(5):
        pars[f'g{i}_center'].set(guess_means[i])
        pars[f'g{i}_amplitude'].set(guess_amp[i])
        pars[f'g{i}_sigma'].set(0.02/2.355) # 0.02 was the rough FWHM I inspected
    composite_model = gauss_multi[0] + gauss_multi[1] + gauss_multi[2] + gauss_multi[3] + gauss_multi[4]

    left_ind = 9999
    right_ind = 0
    for i, val in enumerate(bin_means):
        if val > 5.0:
            left_ind = min(left_ind, i)
        if val > 5.372:
            right_ind = i
            break
    vals_to_fit = vals[left_ind:right_ind]
    bin_means_fit = bin_means[left_ind:right_ind]

    #init = composite_model.eval(vals_to_fit, x=bin_means_fit)
    fit = composite_model.fit(vals_to_fit, pars, x=bin_means_fit)
    plt.scatter(bin_means_fit, vals_to_fit, fc='deeppink', ec='k', lw=0.75, s=6)
    plt.plot(bin_means_fit, fit.best_fit, c='k')

    print(fit.params)
    plt.show()'''