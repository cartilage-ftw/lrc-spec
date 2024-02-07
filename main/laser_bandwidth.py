import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from lmfit.models import LorentzianModel, GaussianModel, VoigtModel
import lmfit

plt.rcParams['figure.dpi'] = 150
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']

root_dir = '../data/bandwidth/16-01-2024/'

"""
Data files
"""

# for this, horizontal was kept at 0 mm (this needs to be remembered)
vertical_scan_dict = {
    1.75: '2024-01-16-17-22-00',
    2.0: '2024-01-16-17-20-22',
    2.25: '2024-01-16-17-19-26',
    2.5: '2024-01-16-17-18-14',
    2.75: '2024-01-16-17-15-39',
    3.0: '2024-01-16-17-14-41',
    3.25: '2024-01-16-17-13-31',
    3.5: '2024-01-16-17-11-24',

    2.0: '2024-01-16-17-24-41',
    1.75: '2024-01-16-17-26-25',

    1.5: 0.0, # the method needs to check the datatype. If it's a float, then use the value, instead of
    1.25: 0.0, # counting from the list
    1.0: 0.,
    0.75: 0.,

    0.0: '2024-01-16-17-31-05',
    0.25: '2024-01-16-17-32-30',
    0.5: '2024-01-16-17-33-33',
    0.33: '2024-01-16-17-34-26',
    0.12: '2024-01-16-17-35-32',
    -0.12: '2024-01-16-17-36-40',
    -0.25: '2024-01-16-17-37-38',
    -0.37: '2024-01-16-17-39-18'
}


# horiz scan with vertical at 4 mm
horiz_scan_4mm_vert = {
    0.: '2024-01-16-17-41-43',
    0.25: '2024-01-16-17-43-31',
    0.5: '2024-01-16-17-44-40',
    0.75: '2024-01-16-17-46-25',
    0.75: '2024-01-16-17-47-36',
    1.0: '2024-01-16-17-48-38',
    1.25: '2024-01-16-17-50-08',
    1.5: '2024-01-16-17-51-17',
    1.75: '2024-01-16-17-52-16',
    2.0: '2024-01-16-17-53-13',
    2.25: '2024-01-16-17-54-28',
    2.5: '2024-01-16-17-55-20',
    2.75: '2024-01-16-17-56-20',
    3.0: '2024-01-16-17-57-37',
    #3.25: '2024-01-16-17-58-52',
    #3.25: '2024-01-16-18-00-57',
    3.5: '2024-01-16-18-02-04',
    3.75: '2024-01-16-18-02-59',
    4.0: '2024-01-16-18-04-09',
    4.25: '2024-01-16-18-05-16',
    4.5: 0.0,
    4.75: 0.,
    5.0: 0.,
    5.25: 0.,
    5.5: '2024-01-16-18-07-00',
    5.75: '2024-01-16-18-08-21',
    5.88: '2024-01-16-18-09-39',
    6.0: '2024-01-16-18-10-41',
    6.12: '2024-01-16-18-11-53',
    6.25: '2024-01-16-18-12-56',
    6.37: '2024-01-16-18-13-54',
    #6.50: '2024-01-16-18-15-03', last 2 probably didn't move the slider
    #6.63: '2024-01-16-18-16-10',
}

# vertical at 3 mm
horiz_scan_3mm_vert = {
    -0.75: 0.,
    -0.25: '2024-01-16-18-22-19',
    0.0: '2024-01-16-18-23-22',
    0.25: '2024-01-16-18-24-39',
    0.50: '2024-01-16-18-25-50',
    0.75: '2024-01-16-18-26-54',
    1.00: '2024-01-16-18-28-01',
    1.25: '2024-01-16-18-29-23',
    1.50 : 0.,
    1.75: 0.,
    2.0: 0.,
    2.25: 0.,
    2.5: 0.,
    2.75: 0.,
    3.0: 0.,
    3.25: 0.,
    3.5: 0.,
    3.75: 0.,
    4.0: 0.,
    4.25: 0.,
    4.5: 0.,
    4.75: '2024-01-16-18-32-04',
    5.0: '2024-01-16-18-33-13',
    5.25: '2024-01-16-18-34-16',
    5.5: '2024-01-16-18-35-16',
    5.75: '2024-01-16-18-36-24',
    6.0: '2024-01-16-18-37-27',    
}


def get_count_from_file(file_name, flatten=False):
    if type(file_name) != str:
        return file_name # probably a float or an int
    
    df = []
    try:
        df = pd.read_csv(root_dir + file_name + '.csv',
                 names=['Voltage', 'Time Tick', 'Irrelevant1', 'Irrelevant2'])
    

    except:
        print("Error reading from ", file_name)

    if flatten == True:
        pass
    else: 
        return len(df)
    '''tick0 = df.iloc[0]['Time Tick']
    df['Abs Time'] = (df['Time Tick'].to_numpy() - tick0)/40_000_000
        #print(df.columns)
    counts, bin_edges = np.histogram(df['Abs Time'], bins='fd')
    count = get_normalized_counts(counts, bin_edges)
    plt.hist(df['Abs Time'], bins='fd')
    plt.show()
        #return count'''
    

"""
Function to fit the laser intensity fluctiation
"""
def line_inten_model(t, a0, a1, omega, phi):
    return a0 + a1*omega*(np.cos(omega*t + phi)**2)


def model_residual(params, t, recorded_counts):
    par_vals = params.valuesdict()
    a0 = params['a0']
    a1 = params['a1']
    omega = params['omega']
    phi = params['phi']

    return (line_inten_model(t, a0, a1, omega, phi) - recorded_counts)


def get_normalized_counts(counts, bin_edges):
    bin_means = (bin_edges[:-1] + bin_edges[1:])/2

    if len(bin_means) < 4:
        print('Too few constraints; skipping normalization')
        return counts
    '''print('Trying to fit')
    print(bin_means)
    print(counts)
    print(bins)'''
    guess_a1 = np.average(counts)
    guess_a0 = abs(np.max(counts) - guess_a1)
    inten_model = Model(line_inten_model, independent_vars=['t'])
    pars = inten_model.make_params(a0=guess_a0, a1=guess_a1, omega=30/6.28, phi=0)
    pars['phi'].set(min=-np.pi, max=np.pi)
    pars['a0'].set(min=0.5*guess_a0)
    pars['a1'].set(min=0.5*guess_a1)
    print(pars)
    #fit_params = create_params(a0 = guess_a0, a1=10, omega=40_000, phi=0)
    fit = inten_model.fit(counts, pars, t=bin_means)
    #fit = minimize(model_residual, fit_params, args=(bin_means,), kws={'recorded_counts': counts})
    print(fit.fit_report())


"""
Plot the 2D image
"""

vert_pos = []
vert_counts = []

horiz_3mm_pos = []
horiz_3mm_counts = []

horiz_4mm_pos = []
horiz_4mm_counts = []

for pos, file_name in vertical_scan_dict.items():
    counts = 0
    if type(file_name) != str:
        counts = file_name
    else:
        counts = get_count_from_file(file_name)
    vert_pos.append(pos)
    vert_counts.append(counts)

for pos, file_name in horiz_scan_3mm_vert.items():
    counts = 0
    if type(file_name) != str:
        counts = file_name
    else:
        counts = get_count_from_file(file_name)
    horiz_3mm_pos.append(pos)
    horiz_3mm_counts.append(counts)

for pos, file_name in horiz_scan_4mm_vert.items():
    counts = 0
    if type(file_name) != str:
        counts = file_name
    else:
        counts = get_count_from_file(file_name)
    horiz_4mm_pos.append(pos)
    horiz_4mm_counts.append(counts)


fig, axes = plt.subplots(2,2, figsize=(8,8))

image = axes[0,0].scatter(np.zeros(len(vert_pos)), vert_pos, c=vert_counts,
                           marker='s', cmap='Blues', label='vertical scan')
image_3mm = axes[0,0].scatter(horiz_3mm_pos, 3*np.ones(len(horiz_3mm_pos)), c=horiz_3mm_counts,
                           marker='s', cmap='Blues', label='horiz scan 3mm vert')
image_4mm = axes[0,0].scatter(horiz_4mm_pos, 4*np.ones(len(horiz_4mm_pos)), c=horiz_4mm_counts,
                           marker='s', cmap='Blues', label='horiz scan 3 mm vert')

axes[0,0].set_ylabel('Vertical Pos [mm]')
axes[0,0].set_xlabel('Horizontal Pos [mm]')

#cbar = plt.colorbar(mappable=image, label='Counts')

axes[0,1].scatter(vert_counts, vert_pos, marker='o', s=16, c='pink',)
axes[0,1].text(x=-0.25, y=0.9*np.max(vert_counts),  s='Horizontal: $0$ mm [fixed]')
axes[0,1].set_ylabel('Vertical Pos [mm]')

axes[1,0].scatter(horiz_3mm_pos, horiz_3mm_counts, marker='o', c='#aaaaff', s=16)
axes[1,0].text(x=-0.25, y=0.9*np.max(horiz_3mm_counts), s='Vertical: $3$ mm [fixed]')
axes[1,0].set_xlabel('Horizontal Pos [mm]')

axes[1,1].scatter(horiz_4mm_pos, horiz_4mm_counts, marker='o', c='#aaaaff', s=16)
axes[1,1].text(x=-0.25, y=0.9*np.max(horiz_4mm_counts), s='Vertical: $4$ mm [fixed]')
axes[1,1].set_xlabel('Horizontal Pos [mm]')



"""
Fit some widths
"""
vertical_df = pd.DataFrame(data={'Position': vert_pos, 'Counts': vert_counts})

vert_line1_data = vertical_df[(vertical_df['Position'] < 3.55) & (vertical_df['Position'] > 1.45)]
line2_vert = GaussianModel()
line2_pars = line2_vert.make_params(sigma=2., amplitude=50_000)
line2_fit = line2_vert.fit(vert_line1_data['Counts'], line2_pars, x=vert_line1_data['Position'])
line2_region = np.linspace(3.55, 1.45, 100)
line2_fitted_y = line2_fit.eval(x=line2_region)

vert_line1_data = vertical_df[(vertical_df['Position'] < 0.5) & (vertical_df['Position'] > -0.4)]
line1_vert = GaussianModel()
line1_pars = line1_vert.make_params(sigma=2., amplitude=50_000)
line1_fit = line1_vert.fit(vert_line1_data['Counts'], line1_pars, x=vert_line1_data['Position'])
line1_region = np.linspace(0.5, -0.4, 100)
line1_fitted_y = line1_fit.eval(x=line1_region)
axes[0, 1].plot(line1_fitted_y, line1_region, c='red')
axes[0, 1].plot(line2_fitted_y, line2_region, c='red')
axes[0, 1].axhline(y=line1_fit.params['center'].value, xmin=0, xmax=1, ls='--', c='red', lw=0.5)
axes[0, 1].axhline(y=line2_fit.params['center'].value, xmin=0, xmax=1, ls='--', c='red', lw=0.5)

print('***Vertical strip ****')
print(f"Line 1 Center:{line1_fit.params['center']} ; FWHM: {line1_fit.params['fwhm']}")
print(f"Line 2 Center:{line2_fit.params['center']} ; FWHM: {line2_fit.params['fwhm']}")
print('*** End of Verticaal strip /****')


def fit_and_plot(df, region, ax):
    #print(df)
    data = df[(df['Position'] > region[0]) & (df['Position'] < region[1])]
    line_model = GaussianModel()
    line_pars = line_model.make_params(sigma=4., center=np.average(data['Position']),
                     amplitude=np.max(df['Counts']))
    
    #print(len(data), 'items in there')
    line_fit = line_model.fit(data['Counts'], line_pars, x=data['Position'])
    
    line_region = np.linspace(region[0], region[1], 100)
    line_fitted_y = line_fit.eval(x=line_region)
    print(f"Line Center:{line_fit.params['center']} ; FWHM: {line_fit.params['fwhm']}")
    ax.plot(line_region, line_fitted_y)
    ax.axvline(x=line_fit.params['center'].value, ymin=0, ymax=1, lw=0.5, ls='--')

"""
Horizontal strip
"""
horiz_3mm = pd.DataFrame(data={'Position': horiz_3mm_pos, 'Counts': horiz_3mm_counts})
fit_and_plot(horiz_3mm, [-0.76, 1.26], axes[1,0])
fit_and_plot(horiz_3mm, [4.5, 6.01], axes[1,0])

horiz_4mm = pd.DataFrame(data={'Position': horiz_4mm_pos, 'Counts': horiz_4mm_counts})
fit_and_plot(horiz_4mm, [2.24, 4.26], axes[1,1])
fit_and_plot(horiz_4mm, [5.49, 6.26], axes[1,1])


""" Fitting a circle into D2
"""

def circle_residual(params, x, y):
    h = params['h'].value
    k = params['k'].value
    a = params['a'].value
    b = params['b'].value
    #cos_phi = np.cos(params['phi'].value)
    res_h = ((x-h)/a)**2
    res_k = ((y-k)/b)**2
    res_r = res_h + res_k - 1**2
    return np.array([res_h, res_k, res_r])


circle2_params = lmfit.Parameters()
circle2_params.add('h', value=0.25, max=-0.01)
circle2_params.add('k', value=0.5,  min=3.0)
circle2_params.add('a', value=5)
circle2_params.add('b', value=5)

x = np.array([0.0, 5.45, 6.19])
y = np.array([-0.06, 3.0, 4.0])
fitted_circ = lmfit.minimize(circle_residual, circle2_params, args=(x,y))


print(lmfit.fit_report(fitted_circ))

h = fitted_circ.params['h'].value
k = fitted_circ.params['k'].value
a = fitted_circ.params['a'].value
b = fitted_circ.params['b'].value

circ = matplotlib.patches.Ellipse((h,k), width=a, height=b, fc='none', ec='r')
axes[0,0].add_patch(circ)

#axes[0,0].set_xlim(-7, 7)
#axes[0,0].set_ylim(-7, 7)
plt.subplots_adjust(top=0.95, left=0.1, right=0.95)
plt.tight_layout()
plt.savefig('laser_fabry_perot.png', dpi=300)
plt.show()