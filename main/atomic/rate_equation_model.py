import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd

#import saturation_curve

from scipy.integrate import quad, solve_ivp
from scipy.constants import hbar, c, pi

plt.rcParams['figure.dpi'] = 150
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Serif']



spot_area = 2E-6 # 2 mm^2

pulse_energy = 300E-9*7.68 # nJ per pulse


wav_12 = 351E-9 # nm
omega_12 = 2*np.pi*c/wav_12

buncher_freq = 1000
laser_freq = 8_000 # 8 kHz

laser_pulse_fwhm = 8E-9 # 8 ns

pulse_delay = 0#400E-9 # 400 ns
laser_bandwidth = 2.5E9 # GHz
laser_omega_sigma = 2*np.pi*laser_bandwidth/2.355



def gaussian_pulse(t, fwhm=laser_pulse_fwhm):
    sigma = fwhm/2.355
    return np.exp(-(t-pulse_delay)**2/(2*sigma**2))



def integrate_lineshape(omega, A21, freq_mode_hop):
    """
    Convolves the Gaussian and Lorentzian parts
    """
    gamma_linewidth = (1/2*np.pi)*(A21 + freq_mode_hop) # 1.4 GHz mode hopping
    return 1/(np.sqrt(2*pi)*(laser_omega_sigma/2.355)) \
                    * np.exp(-(omega - omega_12)**2/(2*laser_omega_sigma**2)) \
                    * (gamma_linewidth/(2*pi)) / ((omega-omega_12)**2 + gamma_linewidth**2/4)

def S(wav, A21):
    omega = 2*pi*c/wav
    integ = quad(integrate_lineshape,
                     omega_12 - 2*laser_bandwidth,
                     omega_12 + 2*laser_bandwidth,
                     args=(A21, 1.4E9))
    s = ((wav**2*pulse_energy/spot_area)/(2*laser_pulse_fwhm*hbar*omega_12)) * integ[0]
    #print('value of s', s)
    return s


def f(t):
    """
    Defines the laser pulses
    """
    num_pulses = laser_freq*(1/buncher_freq)
    # time between two successive pulses
    pulse_sep = 1/laser_freq
    # the pulse takes some extra time between the buncher trigger, etc. to reach
    # if enough pulses have already passed and the bunch has probably been injected
    # into the drift tube
    return np.where( (t+pulse_delay) // pulse_sep < num_pulses, # if enough pulses have not passed
                        gaussian_pulse(t % pulse_sep), # give a pulse
                        0) # otherwise truncate

def RHS(t, rho_vec, A21, A23, A31, transmission, sat_parameter):
    """
    t: time
    rho_vec: (rho_1, rho_2, rho_3)
    """
    d_rho1_dt = A21*rho_vec[1] + A31*rho_vec[2] \
                         - (1/2)*(rho_vec[0] - rho_vec[1])*A21*f(t)*transmission*sat_parameter
    d_rho2_dt = -(A21 + A23)*rho_vec[1] \
                        + (1/2)*(rho_vec[0] - rho_vec[1])*A21*f(t)*transmission*sat_parameter
    d_rho3_dt = A23*rho_vec[1] - (A31 + 0)*rho_vec[2] # could include a damping constant here
    return d_rho1_dt, d_rho2_dt, d_rho3_dt

"""
Saturation Data
Measurements from Nov 20 and 21
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
    return (density_transmit_dict[filter1]*density_transmit_dict[filter2])/1E4

# combines measurements from Nov 20 and 21
sat_data = pd.read_csv('./lu175_nov_saturation.csv', sep='\t', skiprows=1)


transmissions = []
for i in range(len(sat_data)):
    filter1 = sat_data.iloc[i]['Filter 1']
    filter2 = sat_data.iloc[i]['Filter 2']
    transmissions.append(get_transm(filter1, filter2))
sat_data['Transmission'] = transmissions

print(sat_data)
#print(sat_data['Transmission'])


def solve_populations(t_eval, rho_0, A21, transmission):
    print('Starting to solve!')
    ti = time.time()
    # compute frequency dependent saturation parameter (this doesn't need to be re-evaluated constantly)
    sat_parameter = S(omega_12, A21) 
    print('Took', time.time() - ti, 'secs to compute S(omega, A21)')
    rho_sol = solve_ivp(RHS, [0, 1E-3], y0=rho_0, t_eval=t_eval, method='Radau',# jac='array_like',
                         args=(A21, A23, A31, transmission, sat_parameter))
    tf = time.time()
    print(f'Finished! Took {tf-ti} seconds!')

    print(rho_sol.y)
    print('Final level populations:', [rho_sol.y[:,-1]])
    return rho_sol

if __name__ == "__main__":

    A21 = 1.25E7
    A23 = 1.25E7
    A31 = 2.5E-6

    # time interval during which to evaluate
    t_eval = np.arange(0, 1E-3, 1E-9)
    # initial populations
    rho_0 = np.array([1.0, 0.0, 0.0])

    test_A21s = [1.25E7]
    test_transmissions = [0.000718, 0.007632, 0.0177, 0.1303, 0.3527, 1.0] # 0.001987, 
    final_pops_dict = {}
    for A21_i in test_A21s:
        pops_list = []
        for transm in test_transmissions:
            rho_sol = solve_populations(t_eval, rho_0, A21_i, transm)
            pops_list.append(rho_sol.y[:,-1])
            fig2, axes2 = plt.subplots(3, 1)
            for i in range(3):
                axes2[i].plot(rho_sol.t, rho_sol.y[i], label=rf'$|{i+1}\rangle$ population')
            #plt.show()
            #break
        final_pops_dict[A21_i] = pops_list
        break


    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(6,6))
    for i in range(len(axes)):
        axes[i].plot(sat_data['Transmission']*300E-3*(7.92),
                     sat_data[f'Peak {i+1}'], marker='x', c='deeppink', lw=1, mew=0.4, mec='k',
                     label='Nov $^{175}Lu^+$')
        #axes[i].plot(rho_sol.t, rho_sol.y[i], label=rf'$|{i+1}\rangle$ population')
        axes[i].set_ylabel(rf'$\rho_{i+1}$')#(f'Peak {i+1}')##

        for A21_i, populations in final_pops_dict.items():
            print(populations)
            axes[i].plot(np.array(test_transmissions)*300E-3*7.92, np.array(populations)[:, 2],
                         c='cornflowerblue', label='$A_{12}=$' + f"${A21_i}$", marker='o')
        axes[i].legend()
    #axes[0].plot(t_eval, 0.1*f(t_eval)+0.9, c='pink', label='Laser Pulse', zorder=-1)
    axes[2].set_xlabel('Pulse Energy [$\mu$J]')

    plt.tight_layout()
    plt.savefig('saturation_fit.png', dpi=300)
    plt.show()