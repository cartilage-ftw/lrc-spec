import matplotlib.pyplot as plt
import numpy as np


def estimate_mspop(data, ms_cut=0.275):
    total_ions = len(data)
    ms_ions = len(data[data['cycle_time']*1E3 < ms_cut])
    return ms_ions/total_ions


def plot_spectrum(spec_data, species='Lu$^+$', date_obs='', alt_text=''):
    plt.plot(spec_data['Wavenumber'], spec_data['MetastablePop'], c='gray')
    plt.errorbar(spec_data['Wavenumber'], spec_data['MetastablePop'], marker='o', ms=4, capsize=2, ls='',
                c='cornflowerblue', xerr=spec_data['Wavenumber_err'], label=species)
    
    x_min = np.min(spec_data['Wavenumber'])
    y_max = np.max(spec_data['MetastablePop'])

    pulse_details = 'Energy: 0.303$\mu$J/pulse'
    plt.text(x_min, 0.9*y_max, s=f'Date Obs: {date_obs}' + 
                    '\n'+f'Buncher Freq: {alt_text}' +'\n' + f'{pulse_details}')
    plt.xlabel('Wavenumber [cm$^{-1}$]')
    plt.legend()
    plt.tight_layout()
    plt.show()