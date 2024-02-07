import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import saturation_curve, io_utils

from scipy import signal

# laser operated at 8 kHz, with OD filters 2.0 + 3.0 (one after the other)
lu176_nov24_bunch100 = ['2023-11-24-18-53-00', '2023-11-24-18-36-50']

background_nov24 = ['2023-11-24-16-16-03',
                    '2023-11-24-16-19-46',
                    '2023-11-24-16-24-10',
                    '2023-11-24-16-30-53']

lu175_nov24_bunch1k_dict = {
    (0.3, 3.0) : ['2023-11-24-16-36-48',
                    '2023-11-24-16-42-30',
                    '2023-11-24-16-53-52'],
    (0.6, 3.0): ['2023-11-24-17-04-18',
                 '2023-11-24-17-07-43',
                 '2023-11-24-17-10-52'],
    (1.0, 3.0): ['2023-11-24-17-40-23',
                '2023-11-24-17-44-26',
                '2023-11-24-17-48-38']
        }

lu175_nov24_bunch100_dict = {
    (1.0, 3.0): ['2023-11-24-17-53-04',
                 '2023-11-24-17-59-27',
                 '2023-11-24-18-05-48'],
    (2.0, 3.0): ['2023-11-24-18-12-27',
                 '2023-11-24-18-18-12',
                 '2023-11-24-18-23-45']
        }

filters_used = [2.0, 3.0]

output_dir = '../data/spectra/'


'''def compare_wav_coverage(wav_coverages):
    pairwise_diffs = []
    for i in range(1, len(wav_coverages)):
        pairwise_diffs.append(wav_coverages[i-1] - wav_coverages[i])
    '''

def get_wavenum_lag(spec1, spec2):
    """
    Returns the number of indices spec 1 is lagging by, relative to spec 2
    Negative values mean it's leading instead of lagging
    """
    corr = signal.correlate(spec1['Wavenumber'], spec2['Wavenumber'])
    lags = signal.correlation_lags(spec1['Wavenumber'].size, spec2['Wavenumber'].size,
                     mode='full')
    return lags[np.argmax(corr)]


def stack_spectra(file_names, root_dir, density_filters=[0,0]):
    spectra_dfs = []
    wav_coverages = []
    wav_steps = []
    for f in file_names:
        if '.csv' not in f:
            f = f + '.csv'
        arrtime_dist = io_utils.read_lrc_timeseries(root_dir + f, discard_garbage=True)
        # the above method call also discards readings with negative wavenumber values
        spec = io_utils.make_spectrum(arrtime_dist, ms_cut=0.29, filters=density_filters,
                     transm_percent=saturation_curve.get_transm(*density_filters),
                       file_name=output_dir + f[:-4] + '_spectrum.csv')
        spectra_dfs.append(spec)
        wav_coverages.append(np.max(spec['Wavenumber']) - np.min(spec['Wavenumber']))
        wav_steps.append(len(spec['Wavenumber']))
    
    # which spectrum out of these  has the widest wavenumber coverage? 
    widest_spec_index = np.argmax(wav_coverages) # index of that

    # WARNING: don't concatenate array[-2:] + array[:-2] blindly
    # This way the sliced off initial wavenumbers may appear at the end, and the spectrum
    # will be non-sensical. You have to sacrifice and truncate

    # take this as reference, and add shifted/corrected spectra on top of this.
    final_spectrum = spectra_dfs[np.argmin(wav_steps)].to_numpy()
    # assume all bad wavenumbers are those on the red
    for i in range(len(spectra_dfs)):
        if i != np.argmin(wav_steps): # except this one
            spec_i = spectra_dfs[i].to_numpy()
            delta = len(spec_i[:,0]) - np.min(wav_steps)
            final_spectrum[:,0] += spec_i[delta:, 0] # 
            final_spectrum[:,1] = np.sqrt(final_spectrum[:,1]**2 + spec_i[delta:, 1]**2)
            final_spectrum[:,2] += spec_i[delta:, 2]
    '''for i in range(len(spectra_dfs)):
        print(final_spectrum[:5,])
        if i != widest_spec_index:
            lag = get_wavenum_lag(spec1=spectra_dfs[i], spec2=spectra_dfs[widest_spec_index])
            spec_i_arr = spectra_dfs[i].to_numpy()
            # we have to truncate points that aren't averaged.
            last_ind = len(spec_i_arr[lag:,0])
            print('Lag is', lag)
            #final_spectrum[last_ind:,] = np.nan
            # column 0 is wavenum, 1 is wavnum std, 2 is MS frac
            final_spectrum[:last_ind,0] += spec_i_arr[lag:,0] # sum wavenumbers (to be averaged)
            final_spectrum[:last_ind,1] = np.sqrt(final_spectrum[:last_ind,1]**2 + spec_i_arr[lag:,1]**2)
            final_spectrum[:last_ind,2] += spec_i_arr[lag:,2] '''
    final_df = pd.DataFrame(data={'Wavenumber': final_spectrum[:,0]/len(spectra_dfs),
                                  'Wavenumber_err': final_spectrum[:,1],
                                  'MS Fraction': final_spectrum[:,2]/len(spectra_dfs)})
    fig, ax = plt.subplots(figsize=(6,6))
    ax.plot(final_df['Wavenumber'], final_df['MS Fraction'], c='#707070', marker='o', mfc='#aaaaffff', mec='#aaaaffff')
    ax.text(x=np.min(final_df['Wavenumber']), y=0.95*np.max(final_df['MS Fraction']),
                    s=f'OD={density_filters[0]:.1f}+{density_filters[1]:.1f}'+
                    '\n'+f'Transmission: {saturation_curve.get_transm(*density_filters)}\%')
    plt.xlabel('Wavenumber [cm-1]')
    plt.ylabel('Metastable Pop \%')
    plt.ticklabel_format(style='plain', useOffset=False)

    final_df.to_csv('../data/spectra/' + f[:-4] + '_stacked_spec.csv', index=False)
    plt.tight_layout()
    plt.show()
    print('Warning: some wavenumbers may not be overlapping. While averaging, the non-overlapping numbers')
    print("will be stripped off")



if __name__ == "__main__":
    root = '../data/lrc/24-11-2023/'
    stack_spectra(lu176_nov24_bunch100, root_dir='../data/lrc/24-11-2023/', density_filters=filters_used)
    for dict in [lu175_nov24_bunch1k_dict, lu175_nov24_bunch100_dict]:
        for density_filters, files in dict.items():
            stack_spectra(files, root_dir=root, density_filters=density_filters)
        #break
    #stack_spectra(background_nov24, root)
    '''rng = np.random.default_rng()
    x = rng.standard_normal(1000)
    y = np.concatenate([rng.standard_normal(100), x])
    correlation = signal.correlate(x, y, mode="full")
    lags = signal.correlation_lags(x.size, y.size, mode="full")
    lag = lags[np.argmax(correlation)]
    print(lag)
    print(lags)
    print(correlation)
    print(y)
    fig, ax = plt.subplots()
    plt.plot(np.linspace(0, len(x), len(x)), x, label='x')
    plt.plot(np.linspace(0, len(y), len(y)), y+2)
    plt.show()'''


    """
    fig, ax = plt.subplots(figsize=(6,6))
            '''ax.plot(spectra_dfs[i].iloc[-lag:]['Wavenumber'], 30+spectra_dfs[i].iloc[-lag:]['MS Fraction'], marker='o', ls='--',
                        label='Upon correlation')'''
            ax.plot(spectra_dfs[widest_spec_index]['Wavenumber'], 15+spectra_dfs[widest_spec_index]['MS Fraction'],
                     'ko-', label='template')
            ax.plot(spectra_dfs[i]['Wavenumber'], spectra_dfs[i]['MS Fraction'], 'mo-',
                        label='without correlation')
            ax.legend()
            plt.show()
    """