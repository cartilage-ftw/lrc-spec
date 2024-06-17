import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import level_pop

from scipy.stats import bootstrap
from line_fitting import fit_two_gaussians, fit_gaussian

pd.set_option("display.precision", 8) # show up to 8 decimal places
data_nov = '../data/lrc/2023-11-09-18-04-07.csv'

plt.rcParams['figure.dpi'] = 150
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']


class ATDistribution:
    """
    Defines an object that lets you calculate the MS fraction
    If it's unobvious why I defined a separate class for one method,
    it was to facilitate uncertainty analysis with the "bootstrap" method
    the callable "ms_fraction" was needed
    """
    def __init__(self, gs_cutoff):
        self.gs_cutoff = gs_cutoff

    def ms_fraction(self, data):
        ms_counts = np.sum(np.where(data < self.gs_cutoff, 1, 0))#len(data[data['cycle_time'] < self.gs_cutoff])
        ms_frac = ms_counts/len(data)
        return ms_frac


def calculate_bootstrap_error(data, ms_cut):
    """
    TODO: document this properly
    """
    atd = data['cycle_time']*1E3
    # the scipy method requires you to put the data in a sequence
    atd_sample = (atd, ) # thus, I turned it into a tuple
    atd_distrib = ATDistribution(gs_cutoff=ms_cut)
    res = bootstrap(atd_sample, atd_distrib.ms_fraction)
    print("Error on MS Fraction", res.standard_error)
    return res


def read_lrc_timeseries(file_name, buncher_delay=0, buncher_cycle_dur=1E-3, discard_garbage=True):
    """
    Keyword arguments:
        file_name --- name of the data file
        buncher_delay --- in seconds
        buncher_cycle_dur --- in secs

    Data acquisition clock frequency is 40MHz, which corresponds to 25 ns ticks.
    """

    print("Reading:", file_name)
    dat = pd.read_csv(file_name, names=['Mass Voltage', 'Ticks', 'wavenum_obs', 'wavenum_req'])
    
    print("Number of counts recorded:", len(dat))
    garbage = dat[dat['wavenum_obs'] < 0.0]
    print(len(garbage), 'wavenumber recordings with no meaning')
    print('Unique garbage wavenumber values', set(garbage['wavenum_obs']))
    # When does it report -3333333.333 and when does it -2500000.0?
    if discard_garbage == True:
        useful = dat[dat['wavenum_obs'] >= 0.]
    else:
        useful = dat
    useful['time_s'] = useful['Ticks']/40E6 # tick rate is 40 MHz, divide by that to get time in secs
    useful['cycle_time'] = (useful['time_s'] % buncher_cycle_dur) + buncher_delay

    #print("Duration of operation:", useful.iloc[-1]['time_s'])
    print(useful.head())
    return useful

"""
heights, bin_edges = np.histogram(data['cycle_time'], bins=1000)
    bin_centers = (bin_edges[1:] + bin_edges[:-1])/2
    #plt.plot(bin_centers, heights)#,# marker='D', mfc='w', label='$^{175}$Lu$^+$ [Nov 9]')
    plt.step(data['Ticks']/40E3, data['cycle_time'])
    data.to_csv('useful.csv', index=False)
    #plt.hist(data['cycle_time']*1E3, bins=1000)
    plt.xlabel("Arrival Time [ms]")
    plt.tight_layout()
    plt.show()
"""

lu175_hfs_wavenums = [
    (28502.0, 28502.667),
    (28502.7, 28503.3),
    (28503.35, 28504.2)
]


def get_peak_heights(wavenums, ms_percents, wavenum_regions=lu175_hfs_wavenums):
    heights = []
    for i, (wav_min, wav_max) in enumerate(wavenum_regions):
        #print('Min wav', wav_min)
        wav_arr = np.array(wavenums)
        indices = np.where((wav_arr > wav_min) & (wav_arr < wav_max))[0]
        #print('condition is', cond)
        max_count = np.max(np.array(ms_percents)[indices])
        #print(f'Max height of peak {i+1}:', max_count)
        heights.append(max_count)
    return heights


def make_spectrum(data, ms_cut = 0.29, filters=[0, 0], transm_percent=100., file_name='spectrum.csv', save_file=True):
    # truncate steps with garbage wavenum reading (filled with -33333.3 or -25000)
    data_all = data[data['wavenum_obs'] > 0.]
    #plt.hist(data_all['cycle_time'], bins='stone')
    fig, ax = plt.subplots(figsize=(6,6))
    ms_percents = []
    wav_avgs = []
    wav_errs = []
    ms_errors_percent = []
    for wavenum_step, data_step in data.groupby('wavenum_req'):
        wave_obs_av = np.average(data_step['wavenum_obs'])
        wave_obs_std = np.std(data_step['wavenum_obs'])
        arrival_times = data_step['cycle_time']*1E3 # converted to millisec
        atd = ATDistribution(gs_cutoff=ms_cut)
        atd_sample = (arrival_times, )
        print("Running bootstrap for step", wavenum_step)
        #bootstrap_res = bootstrap(atd_sample, atd.ms_fraction)
        ms_percents.append(100*atd.ms_fraction(arrival_times))
        #ms_errors_percent.append(100*bootstrap_res.standard_error)
        wav_avgs.append(wave_obs_av)
        wav_errs.append(wave_obs_std)
    
    laser_text = 'Laser: 8kHz, 56A'
    #ax.errorbar(wav_avgs, xerr=wav_errs, y=ms_percents, yerr=ms_errors_percent, capsize=2, marker='o', markersize=6,
    #            c='deeppink', ls='', label=laser_text)
    #ax.plot(wav_avgs, ms_percents, c='dimgray')
    plt.xlabel('Wavenumber [cm-1]')
    plt.ylabel('Metastable Pop \%')
    plt.ticklabel_format(style='plain', useOffset=False)
    plt.tight_layout()
    
    #ax.text(x=np.min(data_all['wavenum_obs']), y=0.95*np.max(ms_percents),
    #                s=f'OD={filters[0]:.1f}+{filters[1]:.1f}'+'\n'+f'Transmission: {transm_percent}\%')
    #plt.show()
    spec_data = pd.DataFrame(data={'Wavenumber':wav_avgs, 'Wavenumber_err':wav_errs,
                     'MS Fraction':ms_percents})#, 'MS Error':ms_errors_percent})
    if save_file == True:
        spec_data.to_csv(file_name, index=False)
    return spec_data#get_peak_heights(wav_avgs, ms_percents)


if __name__ == "__main__":
    file_names = ['2023-11-16-17-41-06.csv', '2023-11-16-18-03-22.csv', '2023-11-16-18-17-39.csv',
                    '2023-11-16-18-45-27.csv', '2023-11-16-18-55-57.csv']
    dir_path = '../data/lrc/16-11-2023 LRC/'
    buncher_freqs = ['100 Hz', '200 Hz', '500 Hz', '1kHz (Run 1)', '1kHz (Run 2)']
    #data = read_lrc_timeseries(data_nov)
    for file, bunch_freq in zip(file_names, buncher_freqs):
        data = read_lrc_timeseries(dir_path + file)
        count = 0
        wav_avgs = []
        wav_stds = []
        ms_pops = []
        wav_offsets = []
        for wavenum, wavenum_data in data.groupby('wavenum_req'):
            count += 1
            wave_obs_av = np.average(wavenum_data['wavenum_obs'])
            wave_obs_std = np.std(wavenum_data['wavenum_obs'])
            wav_offsets.append(wave_obs_av - wavenum*2)
            wav_stds.append(wave_obs_std)
            #plt.scatter(wavenum_data['cycle_time']*1E3, wavenum_data['wavenum_obs'], label=f'Requested: {wavenum*2} cm-1')
            plt.ticklabel_format(style='plain', useOffset=False)
            bins = np.histogram_bin_edges(wavenum_data['cycle_time']*1E3, bins='fd')
            bin_centers = (bins[:-1] + bins[1:])/2
            #fit = fit_two_gaussians(wavenum_data)
            #if fit.params['ground_amplitude'].stderr > 100:
            #    fit = fit_gaussian(wavenum_data)
            #plt.plot(bin_centers, fit.eval(x=bin_centers))
            ms_percent = 100*level_pop.estimate_mspop(wavenum_data)
            ms_pops.append(ms_percent)
            wav_avgs.append(wave_obs_av)
            #print(f'{ms_percent}% pop in metastable at', wave_obs_av)
            '''plt.hist(wavenum_data['cycle_time']*1E3, bins=bins, label=f'{wave_obs_av:.3f} $\pm$ {wave_obs_std:.3f} cm$-1$ ')
            plt.legend()
            plt.xlabel("Arrival Time [ms]")
            plt.tight_layout()
            plt.show()'''
        spectrum_data = pd.DataFrame(data={'Wavenumber': wav_avgs, 'Wavenumber_err':wav_stds,
                        'MetastablePop': ms_pops})
        level_pop.plot_spectrum(spectrum_data, species='$^{175}$Lu$^+$',
                                 alt_text=bunch_freq, date_obs='Nov 16, 2023')
        #print('Typical wavenum offset', np.average(wav_offsets), 'cm-1; scatter:', np.std(wav_offsets))
        #print('Typical wavenum jitter', np.average(wav_stds), 'cm-1')
        print(count, "measured wavenumber steps")
    