import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from lmfit.models import GaussianModel

root_dir = '../data/bandwidth/15-01-2024/'

plt.rcParams['figure.dpi'] = 150
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']


horiz_meas_dict = {
    0.0: '2024-01-15-15-20-29',
    0.25: '2024-01-15-15-23-26',
    0.5: '2024-01-15-15-25-12',
    0.75: '2024-01-15-15-27-03',
    1.0: '2024-01-15-15-29-13',
    1.25: '2024-01-15-15-35-17',
    1.5: '2024-01-15-15-37-30',
    1.75: '2024-01-15-15-38-09',
    2.0: '2024-01-15-15-38-46',
    2.25: '2024-01-15-15-39-50',
    2.45: '2024-01-15-15-40-27',
    2.50: '2024-01-15-15-41-16',
    2.75: '2024-01-15-15-41-55',
    3.0: '2024-01-15-15-42-37',
    3.25: '2024-01-15-15-43-40',
    3.5: '2024-01-15-15-44-27',
    3.75: '2024-01-15-15-45-12',
    4.0: '2024-01-15-15-46-08',
    4.25: '2024-01-15-15-47-50',
    4.5: '2024-01-15-15-48-51',
    4.75: '2024-01-15-15-49-54',
    5.0: '2024-01-15-15-50-30',
    5.25: '2024-01-15-15-51-15',
    5.5: '2024-01-15-15-52-20',
    5.75: '2024-01-15-15-53-03',
    5.95: '2024-01-15-15-53-46',
    6.0: '2024-01-15-15-54-49',
    6.25: '2024-01-15-15-55-44',
    6.50: '2024-01-15-15-56-32',
    6.75: '2024-01-15-15-57-15',
    7.0: '2024-01-15-15-58-27',
    7.25: '2024-01-15-15-59-13',
    7.5: '2024-01-15-16-00-22',
    7.75: '2024-01-15-16-01-18',
    8.0: '2024-01-15-16-02-05',
    8.25: '2024-01-15-16-03-25',
    8.5: '2024-01-15-16-04-07',
    8.75: '2024-01-15-16-05-11',
    9.0: '2024-01-15-16-05-52',
    9.25: '2024-01-15-16-06-30',
    9.5: '2024-01-15-16-07-19'
}

vert_pos_dict = {
    4.0: '2024-01-15-17-33-09',
    4.12: '2024-01-15-17-34-11',
    4.25: '2024-01-15-17-35-16',
    4.38: '2024-01-15-17-36-11',
    4.50: '2024-01-15-17-37-13',
    4.62: '2024-01-15-17-38-14',
    4.75: '2024-01-15-17-39-04',
    4.88: '2024-01-15-17-39-50',
    5.00: '2024-01-15-17-40-31',
    5.12: '2024-01-15-17-41-09',
    5.25: '2024-01-15-17-42-24'
}

new_vert_strip = {
    21.0: '2024-02-02-20-10-28',
    # 21.25 was empty
    21.5: '2024-02-02-20-13-05',
    21.63: '2024-02-02-20-17-51',
    21.75: '2024-02-02-20-14-42',
    21.88: '2024-02-02-20-16-43',
    22.0: '2024-02-02-20-15-40',
    #22.25 also empty '2024-02-02-20-18-59
    22.38: '2024-02-02-20-21-25',
    22.5: '2024-02-02-20-20-13',
    22.62: '2024-02-02-20-22-36',
    22.75: '2024-02-02-20-23-36',
    22.88: '2024-02-02-20-24-57',
    23.00: '2024-02-02-20-25-56',
    23.13: '2024-02-02-20-27-05',
    23.25: '2024-02-02-20-28-01',
    23.38: '2024-02-02-20-29-08',
    23.50: '2024-02-02-20-29-59',
    23.63: '2024-02-02-20-30-53',
    23.75: '2024-02-02-20-31-55',
    23.88: '2024-02-02-20-33-04',
    24.00: '2024-02-02-20-34-03',
    24.12: '2024-02-02-20-35-05',
    24.25: '2024-02-02-20-36-01',
    24.37: '2024-02-02-20-37-10',
    #24.50: ''
}

new_vert_positions = []
new_vert_counts = []
for pos, file_name in new_vert_strip.items():
    new_vert_positions.append(pos)
    df = pd.read_csv('../data/bandwidth/02-02-2024/' + file_name + '.csv', sep=',')
    new_vert_counts.append(len(df))

new_vert_positions += [24.50, 24.62, 24.75, 24.88, 25.00, 25.12, 25.25, 25.37,
                      25.50, 25.62, 25.75, 25.88, 26.00, 26.12, 26.25, 26.37,
                      26.50, 26.62, 26.75, 26.88, 27.00, 27.12, 27.25, 27.50,
                      27.75, 28.00, 28.50, 29.00, 29.50, 29.25]

new_vert_files = ['2024-02-02-20-38-35.csv',
                    '2024-02-02-20-39-36.csv',
                    '2024-02-02-20-40-27.csv',
                    '2024-02-02-20-41-32.csv',
                    '2024-02-02-20-42-36.csv',
                    '2024-02-02-20-43-39.csv',
                    '2024-02-02-20-44-36.csv',
                    '2024-02-02-20-45-36.csv',
                    '2024-02-02-20-46-35.csv',
                    '2024-02-02-20-47-56.csv',
                    '2024-02-02-20-48-56.csv',
                    '2024-02-02-20-49-58.csv',
                    '2024-02-02-20-51-23.csv',
                    '2024-02-02-20-52-35.csv',
                    '2024-02-02-20-53-32.csv',
                    '2024-02-02-20-54-46.csv',
                    '2024-02-02-20-55-52.csv',
                    '2024-02-02-20-57-16.csv',
                    '2024-02-02-20-58-28.csv',
                    '2024-02-02-20-59-50.csv',
                    '2024-02-02-21-00-47.csv',
                    '2024-02-02-21-01-41.csv',
                    '2024-02-02-21-02-42.csv',
                    '2024-02-02-21-03-53.csv',
                    '2024-02-02-21-04-52.csv',
                    '2024-02-02-21-05-54.csv',
                    '2024-02-02-21-06-54.csv',
                    '2024-02-02-21-07-55.csv',
                    '2024-02-02-21-08-55.csv',
                    '2024-02-02-21-09-57.csv']

for file_name in new_vert_files:
    df = pd.read_csv('../data/bandwidth/02-02-2024/' + file_name, sep=',')
    new_vert_counts.append(len(df))

date_2nd_feb = ['Feb 2'] * len(new_vert_counts)

vert_pos_df = pd.DataFrame(data={'Position': new_vert_positions,
                                 'Counts': new_vert_counts, 'Date': date_2nd_feb})

feb4_pos = [28.25, 28.50, 29.00, 29.50, 30.00, 30.25, 29.25, 29.75, 30.50, 31.00, 30.62, 30.50, 30.75,
            31.25, 31.50, 31.75, 32.00, 32.25, 32.50, 32.12, 32.75, 33.00, 33.25,
            33.50, 33.75, 34.00, 34.25, 34.50]

feb4_files = ['2024-02-04-20-04-01.csv'
	,'2024-02-04-20-06-12.csv'
	,'2024-02-04-20-07-46.csv'
	,'2024-02-04-20-09-09.csv'
	,'2024-02-04-20-10-16.csv'
	,'2024-02-04-20-11-27.csv'
	,'2024-02-04-20-12-53.csv'
	,'2024-02-04-20-14-00.csv'
	,'2024-02-04-20-15-12.csv'
	,'2024-02-04-20-16-34.csv'
	,'2024-02-04-20-18-08.csv'
	,'2024-02-04-20-19-24.csv'
	,'2024-02-04-20-20-36.csv'
	,'2024-02-04-20-21-51.csv'
	,'2024-02-04-20-23-41.csv'
	,'2024-02-04-20-24-48.csv'
	,'2024-02-04-20-25-49.csv'
	,'2024-02-04-20-26-56.csv'
	,'2024-02-04-20-27-57.csv'
	,'2024-02-04-20-28-59.csv'
	,'2024-02-04-20-30-00.csv'
	,'2024-02-04-20-31-36.csv'
	,'2024-02-04-20-32-59.csv'
	,'2024-02-04-20-34-00.csv'
	,'2024-02-04-20-35-50.csv'
	,'2024-02-04-20-36-54.csv'
	,'2024-02-04-20-37-56.csv'
	,'2024-02-04-20-38-47.csv']


feb8_positions = [30.50, 30.70, 30.90, 31.10, 31.30, 31.50,
                  31.63, 31.75, 31.88, 32.00, 32.25, 32.50,
                   
                    30.50, 30.75, 31.00, 
                    
                  30.25, 30.00, 29.75, 29.50, 29.25, 29.00,
                  28.75, 28.50, 28.25, 28.00, 27.75, 27.50, 
                  27.25, 27.00, 26.75, 26.50, 26.25, 26.00,
                  25.75, 25.50, 25.25, 25.00, 24.75, 24.50,
                  24.25, 24.00, 23.75, 23.50, 23.25, 23.00, 22.75, 
                  22.50, 22.25, 22.00, 21.75, 21.50, 21.25]

print('Feb 8 measurements', len(feb8_positions), 'positions')
feb8_files = ['2024-02-08-17-13-47.csv'
	,'2024-02-08-17-16-00.csv'
	,'2024-02-08-17-17-01.csv'
	,'2024-02-08-17-18-29.csv'
	,'2024-02-08-17-19-30.csv'
	,'2024-02-08-17-20-34.csv'
	,'2024-02-08-17-24-57.csv'
	,'2024-02-08-17-26-01.csv'
	,'2024-02-08-17-26-59.csv'
	,'2024-02-08-17-30-03.csv'
	,'2024-02-08-17-31-37.csv'
	,'2024-02-08-17-34-24.csv'
	,'2024-02-08-17-39-34.csv'
	,'2024-02-08-17-41-10.csv'
	,'2024-02-08-17-42-24.csv'
	,'2024-02-08-17-44-32.csv'
	,'2024-02-08-17-46-13.csv'
	,'2024-02-08-17-49-06.csv'
	,'2024-02-08-17-50-23.csv'
	,'2024-02-08-17-51-23.csv'
	,'2024-02-08-17-52-53.csv'
	,'2024-02-08-17-53-59.csv'
	,'2024-02-08-17-55-25.csv'
	,'2024-02-08-17-58-50.csv'
	,'2024-02-08-18-00-26.csv'
	,'2024-02-08-18-30-11.csv'
	,'2024-02-08-18-32-54.csv'
	,'2024-02-08-18-33-59.csv'
	,'2024-02-08-18-35-12.csv'
	,'2024-02-08-18-36-21.csv'
	,'2024-02-08-18-37-47.csv'
	,'2024-02-08-18-42-45.csv'
	,'2024-02-08-18-43-44.csv'
	,'2024-02-08-18-45-06.csv'
	,'2024-02-08-18-46-04.csv'
	,'2024-02-08-18-47-13.csv'
	,'2024-02-08-18-48-11.csv'
	,'2024-02-08-18-49-42.csv'
	,'2024-02-08-18-51-02.csv'
	,'2024-02-08-18-52-15.csv'
	,'2024-02-08-18-53-28.csv'
	,'2024-02-08-18-54-37.csv'
	,'2024-02-08-18-55-47.csv'
	,'2024-02-08-18-57-13.csv'
	,'2024-02-08-18-59-15.csv'
	,'2024-02-08-19-00-19.csv'
	,'2024-02-08-19-01-21.csv'
	,'2024-02-08-19-02-33.csv'
	,'2024-02-08-19-03-31.csv'
	,'2024-02-08-19-05-00.csv'
	,'2024-02-08-19-06-14.csv'
	,'2024-02-08-19-07-49.csv']

print(len(feb8_files), 'files')
feb8_counts = []

for file_name in feb8_files:
    try:
        df = pd.read_csv('../data/bandwidth/08-02-2024/' + file_name, sep=',')
        feb8_counts.append(len(df))
    except:
        print('Error reading', file_name, ': Potentially empty, or not found')
        feb8_counts.append(0)


feb4_counts = []
for file_name in feb4_files:
    try:
        df = pd.read_csv('../data/bandwidth/04-02-2024/' + file_name, sep=',')
        feb4_counts.append(len(df))
    except:
        feb4_counts.append(0)
feb4_date = ['Feb 4'] * len(feb4_files)

vert_pos_feb4 = pd.DataFrame(data={'Position': feb4_pos, 'Counts': feb4_counts, 
                                'Date': feb4_date})

vert_pos_df = pd.concat([vert_pos_df, vert_pos_feb4])

'''for i in range(len(vert_pos_df)):
    print(vert_pos_df.iloc[i])'''
def fit_and_plot(df, region, ax):
    #print(df)
    data = df[(df['Position'] > region[0]) & (df['Position'] < region[1])]
    line_model = GaussianModel()
    line_pars = line_model.make_params(sigma=4., center=np.average(data['Position']),
                     amplitude=np.max(data['Counts']))
    line_pars['center'].set(min=np.min(df["Position"]), max=np.max(df['Position']))
    
    #print(len(data), 'items in there')
    line_fit = line_model.fit(data['Counts'], line_pars, x=data['Position'])
    
    line_region = np.linspace(region[0], region[1], 100)
    line_fitted_y = line_fit.eval(x=line_region)
    print(f"Line Center:{line_fit.params['center']} ; FWHM: {line_fit.params['fwhm']}")
    ax.plot(line_region, line_fitted_y, c='gray')
    ax.axvline(x=line_fit.params['center'].value, ymin=0, ymax=1, lw=0.5, c='gray', ls='--')


fig, axes = plt.subplots(2, 1, figsize=(8,9), sharex=True)

ax = axes[0]
for date, meas_df in vert_pos_df.groupby('Date'):
    ax.scatter(meas_df['Position'], meas_df['Counts'], label=date)

axes[1].scatter(feb8_positions, feb8_counts, c='green', label='Feb 8')

fit_and_plot(vert_pos_df, [21.4, 22.1], ax)
fit_and_plot(vert_pos_df, [22.3, 23.14], ax)
fit_and_plot(vert_pos_df, [23.2, 24.3], ax)
fit_and_plot(vert_pos_df, [24.5, 25.8], ax)
fit_and_plot(vert_pos_df, [27.4, 29.01], ax)
fit_and_plot(vert_pos_df, [29.36, 31.1], ax)
fit_and_plot(vert_pos_df, [30.9, 32.16], ax)

'''def compute_bandwidth(width_r2, D2, screen_dist=49, v_fsr=0.6645, ref_ind=1.4757, wavelength=351E-7):
    bandwidth= width_r2*D2/(2*(ref_ind**2)*screen_dist**2*wavelength)
    print('The bandwidth of the laser is', bandwidth)


positions = []
counts = []

for pos, file_name in horiz_meas_dict.items():
    df = pd.read_csv(root_dir + file_name + '.csv')
    cts = len(df)
    positions.append(pos)
    counts.append(cts)

fig, axes = plt.subplots(1, 2, figsize=(8,6))

vert_positions = []
counts_vert = []

for pos, file_name in vert_pos_dict.items():
    df = pd.read_csv(root_dir + file_name + '.csv')
    cts = len(df)
    vert_positions.append(pos)
    counts_vert.append(cts)

axes[0].plot(positions, counts, marker='o')
axes[1].plot(vert_positions, counts_vert, marker='o')'''
# width of first ring is FWHM ~= 0.415 mm

feb8_df = pd.DataFrame(data={'Position': feb8_positions, 'Counts':feb8_counts})

print("--------- Fitting to Feb 8 data -----------")
fit_and_plot(feb8_df, [21.1, 21.84], axes[1])
fit_and_plot(feb8_df, [21.84, 22.89], axes[1])
fit_and_plot(feb8_df, [22.82, 24.10], axes[1])
fit_and_plot(feb8_df, [24.12, 25.87], axes[1])

fit_and_plot(feb8_df, [28.46, 29.87], axes[1])
fit_and_plot(feb8_df, [30.2, 30.91], axes[1])
fit_and_plot(feb8_df, [31.2, 32.0], axes[1])

for _ax in axes:
    _ax.axhline(y=0, xmin=0, xmax=1, ls='--', c='dimgray', zorder=-1, lw=0.5)
    _ax.legend()

axes[1].set_xlabel('Vertical Height [mm]', fontsize=13)
axes[0].set_ylabel('Counts', fontsize=13)
axes[0].set_title('Fabry-Perot Pattern')
plt.tight_layout()
plt.savefig('fabry_perot_new.png', dpi=300)
plt.show()