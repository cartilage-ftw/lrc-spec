import matplotlib.pyplot as plt
import numpy as np

from lmfit.models import Model, ExponentialModel, GaussianModel, VoigtModel




'''exp_mod = ExponentialModel(prefix='exp_')
pars = exp_mod.guess(y, x=x)'''

def chopped_decaying_exp(x, decay_rate, amplitude, x_cut):
    """
    An exponential function that starts at x > x_cut (and is zero before that)
    """
    y = amplitude*np.exp(-decay_rate*(x-x_cut)) # a normal exponential
    # cut-off region before x_cut (set to 0)
    y[np.where(x <= x_cut)] = 0
    return y



def lin_offset(x, offset):
    return offset


def fit_separate_voigts(spec_data, pos1, pos2, pos3, gauss_sigma, fit_region_width=0.35,
                         vary_gauss_component=True):
    """
    
    Keyword arguments
        fit_region_width -- mask out stuff outside pos +/- fit_region_width while fitting
        an individual line
    """
    v1 = VoigtModel(prefix='v1_')
    v2 = VoigtModel(prefix='v2_')
    v3 = VoigtModel(prefix='v3_')
    #total = v1 + v2 + v3 
    #pars = None
    
    fit_results = []
    for index, (v, pos) in enumerate(zip([v1, v2, v3], [pos1, pos2, pos3])):
        spec_region = spec_data[(spec_data['Wavenumber'] < (pos + fit_region_width))
                                & (spec_data['Wavenumber'] > (pos - fit_region_width))]
        v_pars = v.guess(spec_data['MS Fraction'], x=spec_data['Wavenumber'], center=pos)#
        v_pars[f'v{index+1}_' + 'gamma'].set(value=1, min=0., max=1, vary=True)
        v_pars[f'v{index+1}_' + 'amplitude'].set(min=0., value=20, max=100.)
        v_pars[f'v{index+1}_' + 'center'].set(min=np.min(spec_region['Wavenumber']),
                                max=np.max(spec_region['Wavenumber']))
        # fix the Gaussian sigma
        v_pars['v{0}_sigma'.format(index+1)].set(value=gauss_sigma, vary=vary_gauss_component)

        normed_ms = spec_region['MS Fraction']/100
        epsilon = 0.05 # a smoothing parameter in case a 0.0 value is encountered
        weights = 1/((normed_ms + epsilon) * (1-normed_ms))
        v_fit = v.fit(spec_region['MS Fraction'], v_pars, x=spec_region["Wavenumber"],
                            weights=weights)
                      #method='ampgo')
        fit_results.append(v_fit)
    return fit_results


def fit_triple_voigt(spec_data, pos1, pos2, pos3, gauss_sigma, vary_gauss_component=True):
    v1 = VoigtModel(prefix='v1_')
    v2 = VoigtModel(prefix='v2_')
    v3 = VoigtModel(prefix='v3_')
    #y_offset = Model(lin_offset)
    model = v1 + v2 + v3 #+ y_offset

    pars = v1.guess(spec_data['MS Fraction'], x=spec_data['Wavenumber'], center=pos1)#.set(value=pos1)#guess(,)
    pars.update(v2.guess(spec_data['MS Fraction'],  x=spec_data['Wavenumber'],center=pos2))
    pars.update(v3.guess(spec_data['MS Fraction'], x=spec_data['Wavenumber'],center=pos3))
    #pars.update(y_offset.make_params(offset=0))
    #pars['offset'].set(min=0.0)
    for v in ['v1_', 'v2_', 'v3_']:
        pars[v + 'gamma'].set(value=1, min=0.0, max=0.1, vary=True)
        pars[v + 'amplitude'].set(min=0., value=20, max=100.)
        pars[v + 'center'].set(min=np.min(spec_data['Wavenumber']),
                                max=np.max(spec_data['Wavenumber']))
        # passing vary=False makes the Lorentzian gamma independent of the Gaussian sigma
        pars[v + 'sigma'].set(value=gauss_sigma, vary=vary_gauss_component)
    
    # need to make an initial guess
    #print(pars)
    # the guess() method doesn't work for composite models. I need to guess manually
    #init = model.guess(spec_data['MS Fraction'], x=spec_data['Wavenumber'])
    # now fit
    
    normed_ms = spec_data['MS Fraction']/100
    epsilon = 0.05
    weights = (epsilon + normed_ms)*(1-normed_ms)#np.ones(len(normed_ms))#1/np.sqrt() ##
    return model.fit(spec_data['MS Fraction'], pars, x=spec_data['Wavenumber'], weights=weights)#, method='emcee')
    #print(fit.fit_report())


def fit_two_gaussians(count, pos):#df
    #counts, time_bins = np.histogram(df['cycle_time']*1E3, bins='fd')
    #bin_means = (time_bins[:-1] + time_bins[1:])/2

    gauss_g = GaussianModel(prefix='ground_')
    gauss_m = GaussianModel(prefix='metastable_')
    #pars = gauss_g.make_params(amplitude=30, center=0.30, sigma=0.14)
    #pars.update(gauss_m.make_params(amplitude=10, center=0.23, sigma=0.14))
    #pars['ground_center'].set(min=0.25, max=0.35)
    
    pars = gauss_g.make_params(amplitude=30, center=4, sigma=0.5)
    #pars.update(gauss_m.make_params(amplitude=6, center=7, sigma=0.2))
    
    offset = Model(lin_offset)
    pars.update(offset.make_params(offset=2))
    #pars['metastable_center'].set(min=5, max=10)
    #pars['metastable_amplitude'].set(min=0, max=100)
    #pars['metastable_sigma'].set(min=0.1, max=10)
    '''for x in ['ground_', 'metastable_']:
        pars[x + 'amplitude'].set(min=0.)
        pars[x + 'sigma'].set(max=1.)'''
    total = gauss_g + offset # + gauss_m 
    fit = total.fit(count, pars, x=pos)#x=bin_means)
    print(fit.fit_report())
    return fit


def fit_gaussian(df):
    df['cycle_time_ms'] = df['cycle_time']*1E3

    counts, time_bins = np.histogram(df['cycle_time_ms'], bins='fd', density=True)
    bin_means = (time_bins[:-1] + time_bins[1:])/2
    gauss1 = GaussianModel(prefix='g1_')
    # initialize parameters
    pars = gauss1.make_params(amplitude=30, center=0.30, sigma=0.14)

    '''exp_part = Model(chopped_decaying_exp)
    pars.update(exp_part.make_params(x_cut = 0.30, amplitude=50, decay_rate=10))
    pars['x_cut'].set(min=0.34)'''
    #pars['amplitude'].set(min=2.0)
    #total = gauss1 + exp_part

    fit = gauss1.fit(counts, pars, x=bin_means)

    plt.hist(df['cycle_time_ms'], density=True, bins=time_bins)
    plt.plot(bin_means, fit.best_fit)
    plt.legend()
    plt.show()

    return fit


'''print('\n**->now with updates**\n')
for p in pars:
    print(p)'''

'''
pars['g1_sigma'].set(value=15, min=3)
pars['g1_center'].set(value=105, min=75, max=125)
pars['g1_amplitude'].set(value=2000, min=10)
'''

'''gauss2 = GaussianModel(prefix='g2_')
pars.update(gauss2.make_params())

pars['g2_center'].set(value=155)#, min=125, max=175)
pars['g2_sigma'].set(value=15)#, min=3)
pars['g2_amplitude'].set(value=2000)#, min=10)

mod = gauss1 + gauss2 + exp_mod
print('type of mod', type(mod))
init = mod.eval(pars, x=x)
out = mod.fit(y, pars, x=x)

#print(dir(out.params))
print(out.params)
print('----')
print(type(out.params))
print('****')
print(out.fit_report())
#print(out.fit_report(min_correl=0.5))

fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8))
axes[0].plot(x, y)
axes[0].plot(x, init, '--', label='initial fit')
axes[0].plot(x, out.best_fit, '-', label='best fit')
axes[0].legend()

comps = out.eval_components(x=x)
axes[1].plot(x, y)
axes[1].plot(x, comps['g1_'], '--', label='Gaussian component 1')
axes[1].plot(x, comps['g2_'], '--', label='Gaussian component 2')
axes[1].plot(x, comps['exp_'], '--', label='Exponential component')
axes[1].legend()

plt.show()'''