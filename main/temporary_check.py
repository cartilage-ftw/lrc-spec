import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import io_utils
import saturation_curve

file_name = '2023-12-12-17-56-09.csv'
'''file_names = [
        '2023-12-21-16-12-35', # 0.258
        '2023-12-21-16-25-28',
        '2023-12-21-16-30-16',
        '2023-12-21-16-34-45',
        '2023-12-21-16-40-38', # 0.234
        '2023-12-21-16-46-10', 
        '2023-12-21-16-52-52', #0.215
        '2023-12-21-17-00-15',
    ]'''
#root_dir = '../data/lrc/21-12-2023/'
root_dir = '../data/lrc/22-12-2023/'
file_names = [
    # ms cut at 0.194 ms would be fine
    '2023-12-22-17-21-27', # 2000 cts/step hfs scan
    '2023-12-22-18-49-46', # 500 cts/step hfs scan
    '2023-12-22-18-56-11', # 100 cts/step hfs scan; done with larger delM, to make each step slow
                            # this is important, to provide the motor enough time to properly adjust wavelength
                            # between steps
    '2023-12-22-19-00-58', # 10 cts/step, hfs scan; even larger del M
    '2023-12-22-19-09-16', # 10 cts/step with probe laser blocked (a background scan)
]
for file in file_names:
    dat = io_utils.read_lrc_timeseries(root_dir + file+'.csv', buncher_delay=0, 
                                    buncher_cycle_dur=0.005, discard_garbage=True)
    i = 0
    for wavenum_step, wavenum_data in dat.groupby('wavenum_req'):
        mean_obs_wavenum = np.average(wavenum_data['wavenum_obs'])
        if mean_obs_wavenum > 28503.0 or i == 1:
            fig, ax = plt.subplots(figsize=(6,6))
            
            plt.hist(wavenum_data['cycle_time'], label=wavenum_step, bins='fd')
            #plt.savefig('lol.png', dpi=300)
            plt.legend()
            plt.tight_layout()
            plt.show()
            wavenum_data.to_csv(f'{file}_atd_{mean_obs_wavenum:.2f}cm-1.csv')
            break
        i += 1


'''bunching_freq_dict = {
    '2023-12-12-17-19-49': 50,
    '2023-12-12-16-44-54': 100, # with OD = [0, 4.0]
    '2023-12-12-16-52-32': 100,
    '2023-12-12-17-00-34': 200,
    '2023-12-12-17-05-35': 400,
    #'2023-12-12-17-08-53': 500, # i actually don't know what this is
    '2023-12-12-17-13-48': 500,
    '2023-12-12-17-16-47': 1000,
}

intensity_variation_dict_50hz = {
    '2023-12-12-17-19-49': [0.0, 4.0],
    '2023-12-12-17-32-37': [0.3, 4.0],
    '2023-12-12-17-44-18': [0.6, 4.0],
    '2023-12-12-17-56-09': [1.0, 4.0],
    '2023-12-12-18-08-19': [2.0, 4.0]
}

root_dir = '../data/lrc/12-12-2023/'

if __name__ == '__main__':
    fig, axes = plt.subplots(len(bunching_freq_dict.keys()), 1, figsize=(6,8))
    for i, (file_name, filters) in enumerate(intensity_variation_dict_50hz.items()):
        dat = io_utils.read_lrc_timeseries(root_dir + file_name + '.csv', buncher_delay=0,
                                     buncher_cycle_dur=1/50, discard_garbage=True)
        spec = io_utils.make_spectrum(dat, ms_cut=0.295, filters=filters,
                     transm_percent=saturation_curve.get_transm(*filters),
                     file_name=file_name+'_spectrum.csv')'''