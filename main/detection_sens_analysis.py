import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import io_utils


def show_decreasing_ct_effect():
    fig, axes = plt.subplots(4, 3, figsize=(12,9), sharex='col')

    spectra_files = ['2023-12-22-17-21-27_spectrum', # 2000 cts per step
                     '2023-12-22-18-49-46_spectrum', # 500
                     '2023-12-22-18-56-11_spectrum', # 100
                     '2023-12-22-19-00-58_spectrum'] # 10
    
    off_resonant_atds = ['2023-12-22-17-21-27_atd_28501.62cm-1',
                         '2023-12-22-18-49-46_atd_28501.62cm-1',
                         '2023-12-22-18-56-11_atd_28501.59cm-1',
                         '2023-12-22-19-00-58_atd_28501.62cm-1']
    
    resonant_atds = ['2023-12-22-17-21-27_atd_28503.01cm-1',
                     '2023-12-22-18-49-46_atd_28503.02cm-1',
                     '2023-12-22-18-56-11_atd_28503.01cm-1',
                     '2023-12-22-19-00-58_atd_28503.05cm-1']
    
    counts = [2000, 500, 100, 10]
    # plot the spectra
    for i, spec_file in enumerate(spectra_files):
        spec_data = pd.read_csv('./' + spec_file + '.csv')
        axes[i,0].plot(spec_data['Wavenumber'], spec_data['MS Fraction'], c='gray')
        axes[i,0].errorbar(spec_data['Wavenumber'], spec_data['MS Fraction'], marker='o', ms=6,
                         xerr=spec_data['Wavenumber_err'], ls='', c='#aaaaff', capsize=2)
        axes[i,0].minorticks_on()
        axes[i,0].set_ylabel("MS Fraction \%", fontsize=13)
        axes[i,0].set_xlabel("Wavenumber [cm$^{-1}$]", fontsize=13)
        axes[i,0].ticklabel_format(useOffset=False) # suppress scientific notation with offset
        axes[i,0].tick_params(axis='both', which='major', length=6, labelsize=13)
        axes[i,0].tick_params(axis='both', which='minor', length=3)
        axes[i,0].set_xticks([28502.0, 28503.0, 28504.0])

        max_y = np.max(spec_data['MS Fraction'])
        axes[i,0].text(x=28501.55, y=0.9*max_y, s=f'{counts[i]} counts/' + r'$\bar{\nu}$', fontsize=13)
    
    # now draw the off-resonant ATDs
    for i, off_res_file in enumerate(off_resonant_atds):
        wavenum = off_res_file.split('_')[-1][:-4]
        print('Wavenumber', wavenum)
        atd_data = pd.read_csv('./' + off_res_file + '.csv')
        weights, edges = np.histogram(atd_data['cycle_time']*1E3, bins='fd')
        axes[i, 1].hist(atd_data['cycle_time']*1E3, color='#ffff00', bins='fd', alpha=0.3)
        axes[i, 1].hist(atd_data['cycle_time']*1E3, fill=False, lw=0.3, bins='fd')
        #axes[i, 1].step(edges[1:], bins, c='red')
        axes[i, 1].set_xlabel('Arrival Time [ms]', fontsize=13)
        axes[i, 1].set_xlim(.150, .305)
        axes[i, 1].set_ylabel("Counts", fontsize=13)
        axes[i, 1].text(x=0.155, y=0.9*np.max(weights), s=wavenum + " cm$^{-1}$", fontsize=13)
        axes[i, 0].axvline(x=float(wavenum), ymin=0, ymax=1, lw=2, c='#ffff00', alpha=0.3)
        axes[i, 1].axvline(x=0.194, ymin=0, ymax=1, c='red', ls='--', zorder=-1, label='MS Cut')
        axes[i,1].legend(loc='upper right')
        axes[i,1].tick_params(axis='both', which='major', length=6, labelsize=13)

    for i, on_res_file in enumerate(resonant_atds):
        wavenum = on_res_file.split('_')[-1][:-4]
        print('Wavenumber', wavenum)
        atd_data = pd.read_csv('./' + on_res_file + '.csv')
        if i==3:
            # i had to adjust bins manually for this one
            axes[i, 2].hist(atd_data['cycle_time']*1E3, bins=100, color='pink', alpha=0.3)
            #axes[i, 1].hist(atd_data['cycle_time']*1E3, fill=False, lw=0.3, bins=100) # idk why 100 bins works
            axes[i, 2].text(x=0.155, y=0.9*7, s=wavenum + " cm$^{-1}$", fontsize=13)
            axes[i, 2].hist(atd_data['cycle_time']*1E3, fill=False, ec='red', lw=0.3, bins=100)
        else:
            weights, edges = np.histogram(atd_data['cycle_time']*1E3, bins='fd')
            axes[i, 2].hist(atd_data['cycle_time']*1E3, bins='fd', color='pink', alpha=0.3)
            #axes[i, 1].hist(atd_data['cycle_time']*1E3, fill=False, lw=0.3, bins='fd',
            #                 zorder=-1, label='Resonant')
            axes[i, 2].hist(atd_data['cycle_time']*1E3, fill=False, ec='red', lw=0.3, bins='fd')
            axes[i, 2].text(x=0.155, y=0.9*np.max(weights), s=wavenum + " cm$^{-1}$", fontsize=13)
        axes[i, 2].set_xlabel('Arrival Time [ms]', fontsize=13)
        axes[i, 2].set_xlim(.150, .305)
        axes[i, 2].set_ylabel("Counts", fontsize=13)
        axes[i, 0].axvline(x=float(wavenum), ymin=0, ymax=1, lw=2, c='pink', alpha=0.3)
        axes[i, 2].axvline(x=0.194, ymin=0, ymax=1, c='red', ls='--', zorder=-1, label='MS Cut')

        axes[i,2].tick_params(axis='both', which='major', length=6, labelsize=13)
        axes[i,2].tick_params(axis='both', which='minor', length=3)
    plt.tight_layout()
    #plt.subplots_adjust(hspace=0)
    plt.savefig('effect_decreasing_count_per_step.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    show_decreasing_ct_effect()