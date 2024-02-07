import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


import io_utils, level_pop, line_fitting

"""
A dictionary of which OD value translates to what percentage transmission (for 351 nm, UV)
"""
density_transmit_dict = {0.0 :100,
                         0.1: 82.78, 
                         0.3: 54.77,
                         0.4: 43.12,
                         0.5: 35.27,
                         0.6: 29.91,
                         1.0: 13.03,
                         2.0: 1.77,
                         3.0: 0.24,
                         4.0: 0.03}


def get_transm(filter1, filter2):
    filter1 = float(f'{filter1:.1f}')
    filter2 = float(f'{filter2:.1f}')
    return (density_transmit_dict[filter1]*density_transmit_dict[filter2])/100


file_filter_dict_1kHz = {
    '2023-11-20-17-42-46': [0, 0], # no filters
    '2023-11-20-17-45-02': [1.0, 0],
    '2023-11-20-17-47-10': [0.6, 0],
    '2023-11-20-17-49-23': [0.3, 0],
    '2023-11-20-17-51-32': [0.1, 0],
    '2023-11-20-18-14-10': [0, 0.4], # up to here the scan had 1000 counts per wavenumber step
    '2023-11-20-18-16-11': [0, 0.4], # here onwards, 2000 counts per wavenumber
    '2023-11-20-18-18-48': [1.0, 0.4],
    '2023-11-20-18-21-34': [0.6, 0.4],
    '2023-11-20-18-23-47': [0.3, 0.4],
    # there should've been one measurement for 0.1 + 0.4 
    # 2023-11-20-18-26-04. Maybe it's in the SHE folder and I didn't copy properly
}

file_filter_dict_8kHz = {
    '2023-11-20-17-13-39': [2.0, 0],
    '2023-11-20-17-17-02': [1.0, 0],
    '2023-11-20-17-28-15': [0.6, 0],
    '2023-11-20-17-30-37': [0.3, 0],
    '2023-11-20-17-33-03': [0.1, 0],
    '2023-11-20-17-06-40': [0, 0]
}

file_filter_dict_8kHz_21 = {
    '2023-11-21-17-07-14': [0, 0.4],
    '2023-11-21-17-12-42': [2.0, 0.4],
    '2023-11-21-17-14-36': [2.0, 0.4],
    '2023-11-21-17-19-41': [2.0, 0.5],
    '2023-11-21-17-22-36': [0, 0.5],
    '2023-11-21-17-31-10': [0, 3.0],
    '2023-11-21-17-34-13': [0.1, 3.0],
    '2023-11-21-17-36-52': [0.3, 3.0],
    '2023-11-21-17-39-50': [0.6, 3.0],
    '2023-11-21-17-42-53': [1.0, 3.0],
    '2023-11-21-17-45-32': [1.0, 3.0],
    '2023-11-21-17-50-44': [2.0, 3.0]
}

file_filter_dict_8kHz_nov24 = {
    '': 1
}


def make_saturation_curve(files_and_filters, root_dir, per_pulse_energy, pulse_energy_std=0):
    transmissions = []
    heights = []
    for file_name, filters in files_and_filters.items():
        dat = io_utils.read_lrc_timeseries(root_dir + file_name + '.csv')
        
        # number of laser shots a wavelength step has seen on average
        # it can be done more precisely, separately for each step

        total_steps = len(dat.groupby('wavenum_req')) 
        total_bunches = 1E3 * dat.iloc[-1]['time_s']
        bunches_per_step = total_bunches/total_steps
        num_shots_per_step = bunches_per_step*8 #each bunch saw exactly 8 shots
        
        transmission = get_transm(*filters) # there are 2 filters
        #print(f'Transmission percent for {filters}: {transmission}')
        peak_heights = io_utils.make_spectrum(dat, filters=filters, transm_percent= transmission,
                                              ms_cut=0.302)
        transmissions.append(transmission)
        heights.append(peak_heights)
    fig, axes = plt.subplots(3, 1,sharex=True, figsize=(5,7))
    for i in range(3):
        #print(heights)
        trans = np.array(heights).transpose()
        #print("Transpose of arr is", trans)
        axes[i].scatter(transmissions, trans[i], marker='d', c='cornflowerblue')
        #axes[i].plot(transmissions, trans[i], ls='-', c='gray', zorder=-1)
        axes[i].set_ylabel(f'Peak {i+1}')
        dat = pd.DataFrame(data={
            'Energy Per Pulse [nJ]':(np.array(transmissions)/100)*per_pulse_energy*1E9,
            'Energy_err [nJ]': 1E9*pulse_energy_std/np.sqrt(num_shots_per_step),
            'Peak Amplitude': trans[i],
            'Amp_err':np.zeros(len(transmissions))
            })
        #dat.to_csv(f'saturation_data_8kHz_peak{i+1}_nov21.csv', index=False)
    axes[-1].set_xlabel(f'Transmission \% [of ${per_pulse_energy*1E6}\mu$J]')
    
    plt.subplots_adjust(hspace=0)
    plt.show()


energy_20Nov = pd.read_csv('../data/power/20.11.2023_56A_8kHz_no_od.csv', sep='\t', skiprows=17)
avg_energy_20nov = np.average(energy_20Nov['Measurement'])
std_energy_20nov = np.std(energy_20Nov['Measurement'])
energy_21Nov = pd.read_csv('../data/power/21.11.2023_56A_8kHz_2.csv', sep='\t', skiprows=17)
avg_energy_21nov = np.average(energy_21Nov['Measurement'])
std_energy_21nov = np.std(energy_21Nov['Measurement'])

energy_scaling_8kHz = 2.638/0.333

saturation_dir = '../data/saturation/'
saturation_dat = ['saturation_data_8kHz_peak1_nov21', 'saturation_data_8kHz_peak2_nov21',
                  'saturation_data_8kHz_peak3_nov21']


if __name__ == '__main__':
    '''fig, axes = plt.subplots(3, 1, figsize=(8,6))
    for i in range(len(saturation_dat)):
        df = pd.read_csv(saturation_dir + saturation_dat[i] + '.csv')
        axes[i].errorbar(df['Energy Per Pulse [nJ]'], df['Energy_err [nJ]'],
                          xerr=df['Energy_err [nJ]'], ls='', marker='d', ms=6)
    plt.show()'''
    '''make_saturation_curve(file_filter_dict_1kHz, root_dir='../data/lrc/20-11-2023/',
                            per_pulse_energy=72E-9)'''
    '''lu176 = io_utils.read_lrc_timeseries('../data/lrc/27-11-2023/' + '2023-11-27-20-04-53.csv') # '2023-11-27-18-59-41.csv')
    io_utils.make_spectrum(lu176, ms_cut=0.302, filters=[0, 4.0], transm_percent=get_transm(*[0, 4.0]),
                           file_name='lu175_nov27_hfs.csv')'''
    '''make_saturation_curve(file_filter_dict_8kHz, root_dir='../data/lrc/20-11-2023/',
                            per_pulse_energy=avg_energy_20nov*energy_scaling_8kHz,
                            pulse_energy_std=std_energy_20nov*energy_scaling_8kHz)'''

    make_saturation_curve(file_filter_dict_8kHz_21, root_dir='../data/lrc/21-11-2023/',
                            per_pulse_energy=avg_energy_21nov*energy_scaling_8kHz,
                            pulse_energy_std=std_energy_21nov*energy_scaling_8kHz)