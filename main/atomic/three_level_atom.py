from astropy.constants import c, hbar, k_B, m_p
from scipy.integrate import odeint, quad, solve_ivp
from scipy.special import voigt_profile
from functools import lru_cache

import numpy as np
import astropy.units as u
import scipy.stats
import matplotlib.pyplot as plt

"""
The important thing is to use astropy.units as it's easy to get a factor of \hbar or c wrong
and forgetting to convert from freq to angular freq, etc. Also, the differential equations
are "stiff" for a pulsed laser (i.e. when the intensity changes as a step function
with each pulse). You can use an implicit solver (e.g. 'LSODA' in scipy)
The number of points it needs to sample at can be reduced dramatically 

Note that the relation betweeen A_21 and B_21 coefficients depends on the definition
of laser intensity/energy density chosen.
"""

ENERGY_PER_PULSE = 2.4 * u.nJ
# TODO: OD filters

LASER_PULSE_DUR = 8 * u.ns
LASER_SPOT_SIZE = 1 * u.mm**2

LASER_REP_FREQ = 8_000 * u.Hz
RESONANCE_FREQ = (np.array([28503.01]) / u.cm).to('Hz', equivalencies=u.spectral())
OMEGA_12 = RESONANCE_FREQ * u.cycle
LASER_FWHM_GAUSS = 4.6 * u.GHz

# TODO: this is just for testing.
OMEGA_L = OMEGA_12 
T = 300 * u.K
M_175 = 175 * m_p # Lu-175

DOPPLER_WIDTH_SIGMA = ((OMEGA_12 / c) * np.sqrt(k_B * T / M_175))
NUCLEAR_SPIN = 0. #7/2
# to convert freq to angular freq, multiply by u.cycle
# for the other way, use .to('Hz', equivalencies=[(u.cy/u.s, u.Hz)])
# see https://github.com/astropy/astropy/issues/6032


A21 = 1.25E7 #/ u.s
A23 = 1.88E7 #/ u.s
A31 = 1E-2 #/ u.s # NOTE: this is random value 


p_buncher = 2.5E-2 * u.mbar
p_DT = 3 * u.mbar

#@lru_cache
def pulse_func(t, tau_pulse=LASER_PULSE_DUR, f_rep=LASER_REP_FREQ):
    """
    A function that mimics a laser pulse of a certain duration
    While satisfying the property that area under the curve for each rectangular pulse is 1
    """
    # separation between each pulse in ns
    pulse_sep = (1/f_rep).to('s').value
    
    height = 1/tau_pulse.to('s').value # such that height*width = 1
    #print("Evaluating pulse function")
    return np.where( (t % pulse_sep) < tau_pulse.to('s').value, height, 0)

#@lru_cache
def laser_energy_density(omega, omega_L):
    """
    Centered at laser wavelength (and varies at each point of the laser scan).
    Independent of the atomic line frequency.

    Defined such that it includes the Gaussian lineshape centered at omega_L, with
    FWHM coming from the measured laser bandwidth
    the integral of this energy density over a pulse duration should contain
    """
    # convert to 1 sigma from FWHM; omega_L is in angular units, so u.cycle
    sigma = (LASER_FWHM_GAUSS / 2.355) * u.cycle
    integrated_energy_density = ENERGY_PER_PULSE*LASER_REP_FREQ.value/(c*LASER_SPOT_SIZE)
    #print("Evaluating Laser energy density")
    return integrated_energy_density.to('J s / cm3').value * scipy.stats.norm(loc=omega_L, scale=sigma).pdf(omega)

#@lru_cache
def atomic_lineshape(omega, omega_res, collision_rates, sigma_doppler=DOPPLER_WIDTH_SIGMA):
    """
    Contains the atomic part only. 
    Includes doppler broadening and saturated atomic lineprofile (due to intensity broadening)
     centered at the atom's rest frequency. Does not depend on laser frequency
    Doppler broadening at ~300K is ~340 MHz (1 sigma)
    """
    # NOTE: the collision broadening should include both elastic and inelastic collision rates
    # Here, what I've included is only the inelastic part.
    #gamma = 2*collision_rates
    #d = # diameter of atoms
    #mean_collision_time = (1/(n_He * d**2))*np.sqrt(M_175/(16*np.pi * k_B *T))
    #gamma = 2/mean_collision_time  
    
    #omega *= u.cycle / u.s
    #print("Dimeinsions of omega_res", omega_res)
    #phi = lambda omega: voigt_profile(omega-omega_res, sigma_doppler, gamma)
    #gamma = 1
    #print("evaluating atomic lineshape")
    phi = scipy.stats.norm(loc=omega_res, scale=sigma_doppler)
    return phi.pdf(omega)#.value

#@lru_cache
def S_factor(t, omega_12, omega_L, collision_rates):
    #print("Evaluating S-factor")
    #print("Dimensions of laser energy density", laser_energy_density(1, 1))
    combined_lineshape = lambda omega: (atomic_lineshape(omega, omega_12, collision_rates)*laser_energy_density(omega, omega_L))
    units = u.J * u.s / u.cm**3
    #print("Integrating lineshape")
    int_lineshape = quad(combined_lineshape, 0, np.inf) 
    #print("Successfully computed the lineshape and S_factor!")
    omega_12 = omega_12.to('Hz', equivalencies=[(u.cy/u.s, u.Hz)])
    return (((np.pi**2 * c**3)/(hbar*omega_12**3)) * pulse_func(t) * int_lineshape * units).value

def rate_equations(t_mus, rho_vec, A21, A31, A23, alpha_31, g2, g1, omega_12, omega_L):
    """
    The rate equations are as derived in eqn. 3.7 to 3.9 in Aayush's master thesis.
    """
    rho_1, rho_2, rho_3 = rho_vec
    t = t_mus/1e15
    if t > 1e-8:
        pass
    #print("Evaluating at t=",t)
    #print("Beginning evaluation of the derivatives")
    # absorption - stimulated emission (multiplied by energy density, integrated lineshape function etc.)
    abs_stim_term  = (rho_2 - (g2/g1)*rho_1)*A21*pulse_func(t) #*np.sum(S_factor(t, omega_12, omega_L, alpha_31)))
    # If there are multiple hyperfine multiplets, the location of the resonance omega_12 is an array,
    # and each of those has an associated S-factor. The sum over all these S-factors is taken to take into account
    # pumping neighbouring hyperfine levels at the same time.
    #print("Unit of S-factor:", S_factor(t, omega_12, omega_L, alpha_31))
    #print("No crashes until evaluation of absorption/stimulated emission")
    # the derivatives
    
    drho1_dt = A21*rho_2 + A31*rho_3 \
                + abs_stim_term \
                #+ alpha_31*rho_3 # collisional de-excitation from metastable
    drho2_dt = -(A21 + A23)*rho_2 - abs_stim_term # de-excitation from 2->3 has been neglected
    #print("Successfully evaluated first two derivatives")
    drho3_dt =  A23*rho_2 - (A31 )*rho_3#+ alpha_31
    #print("The value of A31*rho_3", A31*rho_3)
    if drho3_dt < -1:
        pass
        #print("drho_3 too large:", drho3_dt)
        #print("at t", t%(1/LASER_REP_FREQ.to('Hz').value))
    #print("No issues including collision rates\n Finished calculating rate eq derivatives")
    return [drho1_dt/1e15, drho2_dt/1e15, drho3_dt/1e15]

NUM_PULSES = 100 + 0.001/50
tau = LASER_PULSE_DUR.to('s').value
f = LASER_REP_FREQ.to('Hz').value
events = lambda t, rho, *args: np.sin((t-tau/2)*np.pi*(f))
events.terminal = False
events.direction = 1

t_pulses = np.arange(NUM_PULSES) / LASER_REP_FREQ.value

F_u = 1 + NUCLEAR_SPIN
F_l = 0 + NUCLEAR_SPIN

g_u = 2*F_u + 1
g_l = 2*F_l + 1

# assumed as per Kim et al (2024)
rate_coefficient = 1E-11 * u.cm**3 / u.s # 
# P is the pressure (in DT or buncher, wherever you want to use this)
inelastic_collision_rate = lambda P: (rate_coefficient * P/(k_B * T)).to("1/s")

rho_initial = [1.0, 0.0, 0.0]
rho_array = []
t_array = []
rho0 = rho_initial
print("Initializing solution of rate equations")
for t_start, t_end in zip(t_pulses[:-1], t_pulses[1:]):
    rho = solve_ivp(rate_equations,
                    t_span=(LASER_PULSE_DUR.to('s').value/10_000*1e15,
                        (1/LASER_REP_FREQ.to('Hz').value - LASER_PULSE_DUR.to('s').value/10_000)*1e15), # t_span
                    y0=rho0,
                    method='BDF',
                    args=(A21, A31, A23,
                    inelastic_collision_rate(p_buncher).to('1/s').value, # NOTE: this will evaluate the collisions in the bunhcer only
                    g_u,
                    g_l,
                    OMEGA_12,
                    OMEGA_L),
                    atol = 1e-10,
                    rtol=1e-10,
                    first_step = 1e-0,
                    #max_step=1e9,
                     #events=events,
                     #full_output=True,
                     #t_eval=events
            )
    rho0 = rho.y[:,-1]
    rho_array.append(rho.y)
    t_array.extend(rho.t/1e15+t_start)
print("Succesfully integrated rate equations for", NUM_PULSES, "laser pulses")
rho1, rho2, rho3 = np.hstack(rho_array)
t = np.array(t_array)
fig, axes = plt.subplots(3, 1, figsize=(6,8), sharex=True)

#rate_equation_model
axes[0].plot(t*1E6, rho1, c='mediumpurple', ls='--', label=r'$\rho_1$')
axes[1].plot(t*1E6, rho2, c='mediumpurple', ls='--', label=r'$\rho_2$')
axes[2].plot(t*1E6, rho3, c='mediumpurple', ls='--', label=r'$\rho_3$')
#axes[0].plot(t*1E3, rho1, c='mediumpurple', ls='--', label=r'$\rho_1$')
#axes[1].plot(t*1E3, rho2, c='mediumpurple', ls='--', label=r'$\rho_2$')
#axes[2].plot(t*1E3, rho3, c='mediumpurple', ls='--', label=r'$\rho_3$')
axes[2].plot(t*1E6, pulse_func(t)/1.25E8, 'tab:green', alpha=0.4, label='laser')

axes[0].set_ylabel(r'$\rho_1$')
axes[1].set_ylabel(r'$\rho_2$')
axes[2].set_ylabel(r'$\rho_3$')
for ax in axes:
    ax.set_xlabel('Time [ms]')
    ax.grid()
    ax.legend()
    #ax.set_xscale('log')

plt.tight_layout()
plt.scatter(t*1E6, np.zeros_like(t), alpha=0.1)
#plt.scatter(rho.t_events[0]*1E3, np.zeros_like(rho.t_events[0]), marker='x')
#plt.xlim(-4e8, 4e-8)
plt.show()
plt.savefig('three_level_solution.png', dpi=300)