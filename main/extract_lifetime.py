import io_utils
import saturation_curve

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lmfit import Model
from scipy.constants import c, h

pulse_duration = 8.3E-9 # ns
scaling_factor = 7.92 
spot_size = 2.0E-6 # 2 mm^2


"""
file containing summary of linewidth measurements
"""
lw_summary = pd.read_csv("../data/derived/lu175_dec12_data.csv", skiprows=1, sep='\t')

def add_intensity_col(df, laser_freq=8000):
    """
    Using the info about bunching frequency, laser operation frequency, etc. compute the "intensity"
    of laser radiation, taking into account spot size, etc.
    """
    od_filters_trans = []
    for i in range(len(df)):
        filter1, filter2 = str(df.iloc[i]['od_filters']).split(', ')
        # TODO: the get_transm method returns a percentage not fraction
        # This is not a safe way to do things.
        od_filters_trans.append(saturation_curve.get_transm(float(filter1), float(filter2))/100)
    df['Transm'] = od_filters_trans
    df['Number of Shots'] = laser_freq/df['bunching_freq']
    df['Energy Per Pulse'] = 1E-9*df['Transm']*scaling_factor*df['per_pulse_energy']
    df['Intensity'] = (df['Energy Per Pulse']*df['Number of Shots']/(df['bunching_freq']*spot_size))



def saturated_linewidth(intensity, tau, g1, g2, wavelength):
    i_sat = (g1/g2)*np.pi*h*c/((wavelength**3)*tau)
    gamma_sat = (1/tau)*(1 + intensity/i_sat)**(1/2)
    # convert to GHz
    return gamma_sat*h


def extract_tau(df):
    saturation_model = Model(saturated_linewidth)

    gl = 2*(7/2) + 1
    gu_1 = 2*(5/2) + 1
    gu_2 = 2*(7/2) + 1
    gu_3 = 2*(9/2) + 1

    lambda1 = (1/28502.36)*1E7
    lambda2 = (1/28503.02)*1E7
    lambda3 = (1/28503.7)*1E7

    line1_pars = saturation_model.make_params(tau=80E-8, g1=gl, g2=gu_1, wavelength=lambda1)
    line2_pars = saturation_model.make_params(tau=80E-8, g1=gl, g2=gu_2, wavelength=lambda2)
    line3_pars = saturation_model.make_params(tau=80E-8, g1=gl, g2=gu_3, wavelength=lambda3)
    for pars in [line1_pars, line2_pars, line3_pars]:
        pars['g1'].set(vary=False)
        pars['g2'].set(vary=False)
        pars['wavelength'].set(vary=False)

    df1 = df.dropna(subset='gamma1')
    df2 = df.dropna(subset='gamma2')
    df3 = df.dropna(subset='gamma3')
    line1_fit = saturation_model.fit(df1['gamma1'], line1_pars, intensity=df1['Intensity'])
    line2_fit = saturation_model.fit(df2['gamma2'], line2_pars, intensity=df2['Intensity'])
    line3_fit = saturation_model.fit(df3['gamma3'], line3_pars, intensity=df3['Intensity'])

    for fit in [line1_fit, line2_fit, line3_fit]:
        print(fit.fit_report())

# compute intensities
add_intensity_col(df=lw_summary)
print(lw_summary)
extract_tau(df=lw_summary)

lw_summary = lw_summary.sort_values('Intensity')
fig, axes = plt.subplots(3,1,figsize=(6,7))
axes[0].errorbar(lw_summary['Intensity'], lw_summary['gamma1'], yerr=lw_summary['gamma1_err'],
                 marker='o', capsize=2, c='deeppink', label='$F=5/2$')

axes[1].errorbar(lw_summary['Intensity'], lw_summary['gamma2'], yerr=lw_summary['gamma2_err'],
                 marker='o', capsize=2, c='gold', label='$F=7/2$')

axes[2].errorbar(lw_summary['Intensity'], lw_summary['gamma3'], yerr=lw_summary['gamma3_err'],
                 marker='o', capsize=2, c='cornflowerblue', label='$F=9/2$')
for ax in axes:
    ax.set_xlabel('Intensity')
    ax.set_ylabel('$\gamma_s$ [GHz]')
    ax.legend(loc='upper left')
plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.show()


"""
Previous version of the script
"""
"""
measured = pd.read_csv('../data/derived/lu175_summary.csv', sep='\t')
measured.dropna(subset=['OD Filters'], inplace=True)

print('THESE ESFLKASMLK *******')
print(measured['OD Filters'])
print('******')

for i in range(len(measured)):
    filter1, filter2 = str(measured.iloc[i]['OD Filters']).split(', ')
    od_filters_trans.append(saturation_curve.get_transm(float(filter1), float(filter2)))
measured['Transm'] = od_filters_trans

measured['Energy Per Pulse'] = 1E-9*measured['Transm']*scaling_factor*measured['Laser Pulse Energy']
measured['Intensity'] = (measured['Energy Per Pulse']*measured['Number of Shots']/(spot_size*pulse_duration))

print(measured['Intensity'].head())

fig, ax = plt.subplots(figsize=(6,6))

ax.scatter(measured['Intensity'], measured['Gamma2'], label='$F=7/2$', c='deeppink', marker='o', s=12)
ax.scatter(measured['Intensity'], measured['Gamma3'], label='$F=9/2$', c='cornflowerblue', marker='D', s=12)
"""