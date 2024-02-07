import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import io_utils
import level_pop
import line_fitting

file_name_pos_dict_nov16 = {'2023-11-16-16-35-35' : 4.0,
              '2023-11-16-16-36-27': 4.0,
              '2023-11-16-16-39-48': 3.0,
              '2023-11-16-16-42-05': 3.0,
              '2023-11-16-16-44-36': 2.0,
              '2023-11-16-16-46-21': 1.0,
              '2023-11-16-16-47-36': 0.0,
              '2023-11-16-16-56-36': 4.0, 
              '2023-11-16-16-59-16': 5.0, 
              '2023-11-16-17-02-23': 6.0,
              '2023-11-16-17-05-35': 7.0
              }

file_name_pos_1kHz_nodelay = {
        '2023-11-22-15-23-31': 4.0,
        '2023-11-22-15-41-34': 4.0,
        '2023-11-22-15-44-38': 3.0,
        '2023-11-22-15-46-27': 2.0,
        '2023-11-22-15-48-49': 1.0,
        '2023-11-22-15-50-50': 0.0,
        '2023-11-22-15-53-13': 4.0,
        '2023-11-22-15-54-50': 5.0,
        '2023-11-22-16-03-57': 6.0,
}

fileName_pos_1kHz_delayed = {
    '2023-11-22-16-05-37': 6.0,
    '2023-11-22-16-08-02': 7.0,
    '2023-11-22-16-09-52': 5.0,
    '2023-11-22-16-11-40': 4.0,
    '2023-11-22-16-15-22': 4.5,
    '2023-11-22-16-20-26': 3.5,
    '2023-11-22-16-22-07': 3.0,
    '2023-11-22-16-24-22': 2.0,
    '2023-11-22-16-31-56': 1.0,
    # i forgot to copy 2023-11-22-16-29-01.csv, for 0.0 mm reading
}

file_name_pos_2kHz_extended = {
    '2023-11-22-16-39-11': 1.0,
    '2023-11-22-16-41-22': 2.0,
    '2023-11-22-16-43-11': 3.0,
    '2023-11-22-16-46-34': 3.5,
    '2023-11-22-16-51-09': 4.0,
    '2023-11-22-17-06-13': 4.5,
    '2023-11-22-17-08-24': 5.0,
    '2023-11-22-17-10-13': 6.0
}

file_name_pos_2kHz_compact = {
    '2023-11-22-18-38-51': 4.0,
    '2023-11-22-18-41-23': 4.5,
    '2023-11-22-18-43-52': 5.0,
    '2023-11-22-18-45-32': 5.5,
    '2023-11-22-18-47-13': 6.0,
    '2023-11-22-18-49-17': 6.5,
    '2023-11-22-18-51-51': 4.0,
    '2023-11-22-18-52-54': 4.0,
    '2023-11-22-18-54-55': 3.5,
    '2023-11-22-18-56-24': 3.0,
    '2023-11-22-18-57-45': 2.5,
    '2023-11-22-18-59-19': 2.0,
    '2023-11-22-19-00-55': 1.5,
    '2023-11-22-19-02-37': 1.0,
    '2023-11-22-19-04-00': 0.5,
    '2023-11-22-19-05-47': 4.0
}


filename_pos_10kHz_nov27 = {
    '2023-11-27-18-11-47': 4.0,
    '2023-11-27-18-15-51': 3.5,
    '2023-11-27-18-16-53': 3.5, # 10,000 counts in each file from this one onwards
    '2023-11-27-18-18-12': 3.0,
    '2023-11-27-18-19-15': 2.5,
    '2023-11-27-18-20-25': 2.0,
    '2023-11-27-18-21-37': 4.0,
    '2023-11-27-18-22-28': 4.5,
    '2023-11-27-18-23-25': 5.0,
    '2023-11-27-18-24-26': 5.5,
    '2023-11-27-18-25-15': 6.0,
    '2023-11-27-18-26-37': 4.0,
    '2023-11-27-18-27-55': 3.75
}
clean_bg_nov22 = '2023-11-22-16-33-14'
# with the probe laser blocked; testing if the ablation itself produced any background
# metastable ions
clean_bg = ['2023-11-16-16-50-32', '2023-11-16-16-53-31']

def plot_cloud_extent(file_name_pos_dict, data_dir, ms_cut=0.275, export_file_name=''):
    fig, ax = plt.subplots(figsize=(7,5))
    wavenums = []
    wavenum_stds = []
    #max = 0
    ms_counts = []
    positions = []
    for file_name, pos in file_name_pos_dict.items():
        datum = io_utils.read_lrc_timeseries(data_dir + file_name + '.csv', discard_garbage=False)
        #datum = datum[datum['wavenum_obs'] > 0.0]
        ms_count = len(datum[datum['cycle_time']*1E3 < ms_cut])
        #max = np.max([ms_count, max])
        wavnum_avg = np.average(datum[datum['wavenum_obs'] > 0.]['wavenum_obs'])
        wavenums.append(wavnum_avg)
        wavnum_std = np.std(datum[datum['wavenum_obs'] > 0.]['wavenum_obs'])
        wavenum_stds.append(wavnum_std)
        total_count = len(datum)
        ms_fraction = 100*(ms_count/total_count)
        ms_counts.append(ms_fraction)
        positions.append(pos)
        #ax.errorbar(pos, ms_fraction, xerr=0.08, capsize=2, marker='o', ms=4,
        #            label=f'${wavnum_avg:.3f} \pm {wavnum_std:.3f}$')

    wavenum_text = r'$\bar{\nu}$: ' +f'${np.average(wavenums):.3f}\pm{np.average(wavenum_stds):.3f}$ ' + 'cm$^{-1}$'
    ax.errorbar(positions, ms_counts, xerr=0.08, capsize=3, marker='o', #mec='dimgray', mew=0.4,
                c='deeppink', ms=4, ls='', label=wavenum_text)

    gauss_fit = line_fitting.fit_two_gaussians(ms_counts, positions)
    x = np.linspace(0, 7.5, 100)
    fit_y = gauss_fit.eval(x=x)

    ax.plot(x, fit_y, c='gray', zorder=-1, label=f"FWHM: {gauss_fit.params['ground_fwhm'].value:.2f} mm")
    ax.minorticks_on()
    ax.tick_params(which='major', axis='both', length=8, direction='in', labelsize=13)
    ax.tick_params(which='minor', axis='both', length=4, direction='in', labelsize=13)
    #ax.text(x=0.0, y=0.95*np.max(ms_counts), s=)
    plt.legend(loc='upper left', fontsize=13)
    plt.ylabel('Excited Population', fontsize=13)
    plt.xlabel('Slider Position [mm]', fontsize=13)
    plt.tight_layout()
    df = pd.DataFrame(data={'Excited Population': ms_counts, 'Slider Position [mm]':positions})
    df.to_csv(export_file_name, index=False)
    #plt.savefig('../figures/ion_cloud_extent.png', dpi=300)
    plt.show()

if __name__ == "__main__": 
    '''plot_cloud_extent(file_name_pos_dict_nov16, '../data/lrc/16-11-2023 LRC/', ms_cut=0.275,
                            export_file_name='cloud_8kHz_extended.csv')

    plot_cloud_extent(file_name_pos_1kHz_nodelay, '../data/lrc/Ion Cloud 22-11-2023/', ms_cut=0.3125,
                            export_file_name='cloud_1kHz_nodelay.csv')

    plot_cloud_extent(fileName_pos_1kHz_delayed, '../data/lrc/Ion Cloud 22-11-2023/', ms_cut=0.3125,
                            export_file_name='cloud_1kHz_delayed.csv')'''

'''    plot_cloud_extent(file_name_pos_2kHz_extended, '../data/lrc/Ion Cloud 22-11-2023/', ms_cut=0.3125,
                            export_file_name='cloud_2kHz_extende.csv')

    plot_cloud_extent(file_name_pos_2kHz_compact, '../data/lrc/Ion Cloud 22-11-2023/', ms_cut=0.3125,
                            export_file_name='cloud_2kHz_extende.csv')'''

'''plot_cloud_extent(filename_pos_10kHz_nov27, '../data/lrc/27-11-2023/', ms_cut=0.29,
                        export_file_name='cloud_nov27_compact.csv')'''
