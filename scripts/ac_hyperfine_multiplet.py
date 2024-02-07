import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lmfit.models import GaussianModel, LorentzianModel, VoigtModel

plt.rcParams['figure.dpi'] = 150
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']

# energy of the level in absence of splitting
ref_energy = 0.

# I have averaged the FSCC and CI + MBPT estimates
# energy values are in MHz
A = (-2097.62 - 2808.948)/2 # for singlet ^1P_1 (6d7p)
#B =  #

# for triplet ^3P_1 (7s7p)
# A = ()
# B = 

i_spin = float(eval(input("Enter nuclear spin I:\n")))
j_spin = float(eval(input("Enter electronic spin J:\n")))

f_spins = np.arange(i_spin - j_spin, i_spin+ j_spin + 1, 1)


lu175_data = pd.read_csv('../data/lu/dat', sep='\t',
                          names=['wavenumber', 'x_err', 'y', 'y_err'])
lu175_data['E'] = lu175_data['wavenumber'] * 29.979
lu175_data['del_E'] = lu175_data['E'] - 28503.149*29.979
def compute_levels(I, J, F_list, A, B) -> np.array:
    """
    A, B are the hyperfine constants.
    """
    dE_list = []
    for F in F_list:
        K = F*(F+1) - I*(I+1) - J*(J+1)
        dE = (1/2)*A*K #+ B*( ((3/2)*K*(K+1) - 2*I*(I+1)*J(J+1))/
                         #       (2*I*(2*I-1)*2*J*(2*J-1)) )
        dE_list.append(dE)
        print(f'F:{F}, E = {dE}')
    return np.array(dE_list)


def draw_multiplets(F_list:list, energy_levels:np.array, sigma:float):
    #while drawing a synthetic spectrum, draw the prediction in the energy range [min_E, max_E]
    min_E = -9999.
    max_E = 9999.

    intensity = None # reinitialize later on
    components = []
    for i, (F, dE) in enumerate(zip(F_list, energy_levels)):
        # energy
        min_E = min(energy_levels - 4*sigma)
        max_E = max(energy_levels + 4*sigma)
        line = VoigtModel(prefix=f'F{i+1}_')
        pars = line.make_params()
        pars[f'F{i+1}_center'].set(dE)
        # TODO: change this to intensity derived from Wigner n-j symbols
        pars[f'F{i+1}_amplitude'].set(1)
        pars[f'F{i+1}_sigma'].set(sigma)

        E_range = np.linspace(min_E, max_E, 1000)
        func = line.eval(pars, x=E_range)
        components.append(func)
        if intensity is None:
            intensity = func
        else:
            intensity = intensity + func # this works if they are both numpy arrays

    norm = max(intensity)
    fig, ax = plt.subplots(figsize=(6,6))
    ax.plot(E_range, intensity/norm, label='total', c='dimgray')
    for F, comp in zip(F_list, components):
        frac_val = str(F.as_integer_ratio()[0]) + '/' + str(F.as_integer_ratio()[1])
        ax.plot(E_range, np.array(comp)/norm, ls='--', label=f'F={frac_val}')

    ax.scatter(lu175_data['del_E'], lu175_data['y'], label='$^{175}$Lu$^+$')
    ax.set_ylim(bottom=0.)
    ax.set_xlabel('$\Delta3/ E$ [GHz]')
    ax.set_ylabel('Intensity')
    plt.legend()
    plt.tight_layout()
    plt.show()


# soon, turn that into a dictionary of levels and intensities
# for now, equal intensity is assumed
levels = compute_levels(i_spin, j_spin, f_spins, A, B=0.)

# sigma of line is given in GHz
draw_multiplets(f_spins, levels/1000, sigma=2.59) # divide by 1000 to convert to GHz