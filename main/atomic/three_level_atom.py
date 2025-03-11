from astropy.constants import c, hbar, k_B, m_p
from scipy.integrate import odeint, quad, solve_ivp
from scipy.special import voigt_profile
from functools import lru_cache
from sympy.physics.wigner import wigner_6j

import numpy as np
import astropy.units as u
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.formatter.useoffset'] = False

"""
The important thing is to use astropy.units as it's easy to get a factor of \hbar or c wrong
and forgetting to convert from freq to angular freq, etc. Also, the differential equations
are "stiff" for a pulsed laser (i.e. when the intensity changes as a step function
with each pulse). You can use an implicit solver (e.g. 'LSODA' in scipy)
The number of points it needs to sample at can be reduced dramatically 

Note that the relation betweeen A_21 and B_21 coefficients depends on the definition
of laser intensity/energy density chosen.
"""

ENERGY_PER_PULSE = 2.4E-9 * u.J / 1e6
# TODO: OD filters

BUNCHER_FREQ = 1000 * u.Hz
LASER_PULSE_DUR = 8 * u.ns
LASER_SPOT_SIZE = 1.5 * u.mm**2

LASER_REP_FREQ = 8000 * u.Hz
RESONANCE_FREQ = (np.array([28502.3, 28503.06, 28503.7]) / u.cm).to('Hz', equivalencies=u.spectral())
OMEGA_12 = RESONANCE_FREQ * u.cycle
LASER_FWHM_GAUSS = 4.6E9 * u.Hz # NOTE: u.GHz wasn't doing the right thing

T = 300 * u.K
M_175 = 175 * m_p # Lu-175

DOPPLER_WIDTH_SIGMA = ((OMEGA_12 / c) * np.sqrt(k_B * T / M_175))
NUCLEAR_SPIN = 7/2 
# to convert freq to angular freq, multiply by u.cycle
# to get Hz back WITHOUT the 2pi factor, use .to('Hz', equivalencies=[(u.cy/u.s, u.Hz)])
# else, equivalencies=u.dimensionless_angles()
# see https://github.com/astropy/astropy/issues/6032


A21 = 1.25E7 #/ u.s
A23 = 1.88E7 #/ u.s
A31 = 1E-2 #/ u.s # NOTE: this is random value 


p_buncher = 2.5E-2 * u.mbar
p_DT = 3 * u.mbar

def A_FF(A_J, J_u, J_l, F_u, F_l, I):
    """
    The Einstein A_21 coefficient for a hyperfine transition from F' -> F
    in terms of A_21 coefficient for the fine structure transition.
    This takes care of the Clebsch-Gordon coefficients/Wigner nj symbols 
    """
    return (2*F_l + 1)*(2*F_u + 1) * wigner_6j(J_l, I, F_l, F_u, J_l, J_u)**2 * A_J

#@lru_cache
def pulse_func(t: float, tau_pulse=LASER_PULSE_DUR, f_rep=LASER_REP_FREQ):
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
    if isinstance(omega, u.Quantity):
        omega = omega.to('Hz', equivalencies=[(u.cy/u.s, u.Hz)]).value
    if isinstance(omega_L, u.Quantity):
        omega_L = omega_L.to('Hz', equivalencies=[(u.cy/u.s, u.Hz)]).value
    # convert to 1 sigma from FWHM
    sigma = (LASER_FWHM_GAUSS / 2.355).to('Hz').value
    # the ions never stay longer than 1 buncher cycle duration. The amount of energy deposited
    # is limited by that. TODO: Decide what to do here
    num_pulses = (1/LASER_REP_FREQ).value# /np.min(1., BUNCHER_FREQ)
    integrated_energy_density = ENERGY_PER_PULSE*num_pulses/(c*LASER_SPOT_SIZE)
    return integrated_energy_density * np.exp(-(omega - omega_L)**2 / (2*sigma**2))

print("Note: Collisional broadening has been ignored.")
#@lru_cache
def atomic_lineshape(omega, omega_21, collision_rates, sigma_doppler):
    """
    Contains the atomic part only. 
    Includes doppler broadening and saturated atomic lineprofile (due to intensity broadening)
     centered at the atom's rest frequency. Does not depend on laser frequency
    Doppler broadening at ~300K is ~340 MHz (1 sigma)
    """
    if isinstance(omega, u.Quantity):
        omega = omega.to('Hz', equivalencies=[(u.cy/u.s, u.Hz)]).value
    if isinstance(omega_21, u.Quantity):
        omega_21 = omega_21.to('Hz', equivalencies=[(u.cy/u.s, u.Hz)]).value
    if isinstance(sigma_doppler, u.Quantity):
        sigma_doppler = sigma_doppler.to("Hz").value
    #sigma_doppler = sigma_doppler.to('Hz', equivalencies=[(u.cy/u.s, u.Hz)]).value
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
    phi = scipy.stats.norm(loc=omega_21, scale=sigma_doppler)
    return phi.pdf(omega)#.value

@lru_cache
def S_factor(omega_12, omega_L, doppler_sigma, collision_rates):
    """
    The arguments of this method had to be dimensionless to be cache-able
    """
    combined_lineshape = lambda omega: (atomic_lineshape(omega, omega_12, collision_rates, doppler_sigma)\
                                        *laser_energy_density(omega, omega_L))#.value
    units = combined_lineshape(omega_12).unit
    # for integrating with scipy
    dimless_lineshape = lambda omega: combined_lineshape(omega).value
    # NOTE: the integral shouldn't be done from 0, np.inf because the numerical step size would not
    # resolve the lineshape at all.
    int_lineshape = quad(dimless_lineshape,
                        np.min([omega_12, omega_L]) - 5*LASER_FWHM_GAUSS.value,
                        np.max([omega_12, omega_L]) + 5*LASER_FWHM_GAUSS.value
                    )[0]
    #print("Value of integrated lineshape", int_lineshape)
    omega_12 *= 2*np.pi #.to('Hz', equivalencies=u.dimensionless_angles())
    # NOTE: the 2*pi factor in the angular frequency must be correctly included, otherwise
    # the value will be wrong by 8pi^3. I first tried doing this dynamically
    # But then the cache would not work because the arguments were not dimensionless
    # so a unit in 'Hz' is passed to the arguments, without the 2pi factor
    return (((np.pi**2 * c**3)/(hbar*omega_12**3)) * int_lineshape * units).si.value


'''omega = np.linspace(OMEGA_12.value - 5*LASER_FWHM_GAUSS.value, OMEGA_12.value + 5*LASER_FWHM_GAUSS.value, 100) * u.cycle * u.Hz
#print("Value of S-factor", S_factor(1, OMEGA_12, OMEGA_L, A31))
S = [S_factor(OMEGA_12, o, A31) for o in omega]
fig, ax = plt.subplots(figsize=(6,6))
half_max = np.max(S)/2
ax.plot((omega-OMEGA_12)/1E9, S)
ax.axhline(half_max, xmin=0, xmax=1, c='gray', ls='--')
ax.set_xlabel("Frequency detuning [GHz]")
ax.set_ylabel("S-factor [dimensionless]")
print("Width of it", 2.355*np.std(S)/1E9)
plt.show()'''


def rate_equations(t_mus, rho_vec, A21, A31, A23, alpha_31, g2, g1, omega_12, omega_L):
    """
    The rate equations are as derived in eqn. 3.7 to 3.9 in Aayush's master thesis.
    """
    rho_1, rho_2, rho_3 = rho_vec
    t = t_mus#/1e12#/1e6#/1e15
    #print("Evaluating at t=",t)
    #print("Beginning evaluation of the derivatives")
    omega_12 = omega_12.to('Hz', equivalencies=[(u.cy/u.s, u.Hz)]).value
    omega_L = omega_L.to('Hz', equivalencies=[(u.cy/u.s, u.Hz)]).value
    S_factors = 0.
    # If there are multiple hyperfine multiplets, the location of the resonance omega_12 is an array,
    # and each of those has an associated S-factor. The sum over all these S-factors is taken to take into account
    # pumping neighbouring hyperfine levels at the same time.
    for omega_F, d_omega in zip(omega_12, DOPPLER_WIDTH_SIGMA):
        S_factors += S_factor(omega_F, omega_L, alpha_31, d_omega)
    # absorption - stimulated emission (multiplied by energy density, integrated lineshape function etc.)
    abs_stim_term = (rho_2 - (g2/g1)*rho_1)*A21*pulse_func(t) + S_factors

    if np.abs(abs_stim_term) > 1e10: # if the value is too large, it's already saturated
        r = np.rint(np.log10(np.abs(abs_stim_term))) - 10
        # for numerical safety, cap it here. It's not going to change the solution as A21+A23 ~10^7 << 10^10
        abs_stim_term /= 10**r
        #exit()
    #print("Unit of S-factor:", S_factor(t, omega_12, omega_L, alpha_31))
    #print("No crashes until evaluation of absorption/stimulated emission")
    # the derivatives
    
    #print("Value of absorption term", abs_stim_term)
    drho1_dt = A21*rho_2 + A31*rho_3 \
                + abs_stim_term \
                #+ alpha_31*rho_3 # collisional de-excitation from metastable
    #print("Value of derivative 1", drho1_dt)
    drho2_dt = -(A21 + A23)*rho_2 - abs_stim_term # de-excitation from 2->3 has been neglected
    #print("Successfully evaluated first two derivatives")
    drho3_dt =  A23*rho_2 - (A31 )*rho_3#+ alpha_31
    #print("The value of A31*rho_3", A31*rho_3)
    if drho3_dt < -1:
        pass
        #print("drho_3 too large:", drho3_dt)
        #print("at t", t%(1/LASER_REP_FREQ.to('Hz').value))
    #print("No issues including collision rates\n Finished calculating rate eq derivatives")
    return np.array([drho1_dt, drho2_dt, drho3_dt])#/1e12

NUM_PULSES = round((LASER_REP_FREQ / BUNCHER_FREQ).value) #+ 0.001/50
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
print("Initializing solution of rate equations")


def calculate_populations(omega_L):
    # initial level populations.
    rho0 = [1.0, 0.0, 0.0]
    # and variables to store the full solution
    rho_array = np.array([rho0])
    t_array = [0.]
    for n in range(NUM_PULSES):
        t_dur = 1*LASER_PULSE_DUR.to('s').value
        h = t_dur/10_000 # a very small quantity
        time_between_pulses = 1/LASER_REP_FREQ.to('Hz').value
        # slice the time into parts with and without the pulse
        t_with_pulse = np.array([
                            h, t_dur + h,
                        ]) + n*time_between_pulses
        t_no_pulse =  np.array([
                            t_dur + h, time_between_pulses + h,
            ]) + n*time_between_pulses
        #max_step_sizes = [1e2, 1e6]

        # now integrate the rate equations for each of these slices
        for i, (start, end) in enumerate([t_with_pulse, t_no_pulse]):
            #events = lambda t, rho, *args: np.sin((t-(end-start-2*h)/2)*np.pi*(f))
            #events = lambda t, rho, *args: np.sin(np.pi*(t-start)/(end-start)) - 1e-8
            #events.terminal = False
            #events.direction = 0
            rho = solve_ivp(rate_equations,
                        t_span=(start, end), # t_span
                        y0=rho0,
                        method='LSODA',
                        args=(A21, A31, A23,
                        inelastic_collision_rate(p_buncher).to('1/s').value, # NOTE: this will evaluate the collisions in the bunhcer only
                        g_u,
                        g_l,
                        OMEGA_12,
                        omega_L),
                        atol = 1e-7,
                        #rtol=1e-7,
                        #init_step=1e0,
                        #max_step=1e12/1e6,
                        #events=events,
                        #max_step=max_step_sizes[i],
                )
            rho0 = rho.y[:,-1]
            #print(rho.y)
            rho_array = np.append(rho_array, rho.y.T, axis=0)
            '''fig, ax = plt.subplots()
            ax.plot(rho.t/1e3, rho.y.T, label=['rho1', 'rho2', 'rho3'])
            ax.plot(rho.t/1e3, pulse_func(rho.t/1e12)/1.25E8, 'k', ls='--', label='laser')
            if rho.t_events is not None:
                ax.scatter(rho.t_events[0], np.zeros_like(rho.t_events[0]), marker='x')
            ax.scatter(rho.t/1e3, np.zeros_like(rho.t), alpha=0.1)
            print(rho_array.shape)
            print("Values of t in the solution", rho.t)
            #print("")
            ax.legend()'''
            t_array = np.append(t_array, rho.t, axis=0)
            #plt.show()
            #break

    print("Succesfully integrated rate equations for", NUM_PULSES, "laser pulses at wavenumber",
                         (omega_L/u.cycle).to('1/cm', equivalencies=u.spectral()))
    #rho1, rho2, rho3 = np.hstack(rho_array)
    return np.array(t_array), rho_array
    #print("Values in t_array", t/1E-6)



def plot_optical_pumping(rho_array, t):
    fig, axes = plt.subplots(3, 1, figsize=(6,8), sharex=True)
    #rate_equation_model
    #print("The values of rho_array", rho_array)
    axes[0].plot(t*1E6, rho_array[:,0], c='mediumpurple', ls='--', label=r'$\rho_1$')
    axes[1].plot(t*1E6, rho_array[:,1], c='mediumpurple', ls='--', label=r'$\rho_2$')
    axes[2].plot(t*1E6, rho_array[:,2], c='mediumpurple', ls='--', label=r'$\rho_3$')
    #axes[0].plot(t*1E3, rho1, c='mediumpurple', ls='--', label=r'$\rho_1$')
    #axes[1].plot(t*1E3, rho2, c='mediumpurple', ls='--', label=r'$\rho_2$')
    #axes[2].plot(t*1E3, rho3, c='mediumpurple', ls='--', label=r'$\rho_3$')

    axes[0].set_ylabel(r'$\rho_1$')
    axes[1].set_ylabel(r'$\rho_2$')
    axes[2].set_ylabel(r'$\rho_3$')

    t_all = np.linspace(0, NUM_PULSES*(1/LASER_REP_FREQ.to('Hz').value), 1000000)
    axes[2].plot(t_all*1E6, pulse_func(t_all)/1.25E8, 'tab:pink',label='laser')
    axes[2].plot(t*1E6, pulse_func(t)/1.25E8, 'tab:green', label='laser')
    for ax in axes:
        ax.set_xlabel('Time [$\mu$s]')
        ax.grid()
        ax.legend()
        #ax.set_xscale('log')

    plt.tight_layout()
    plt.scatter(t*1E6, np.zeros_like(t), alpha=0.1)
    #plt.scatter(rho.t_events[0]/1e12*1e6, np.zeros_like(rho.t_events[0]), marker='x')
    #plt.xlim(-4e8, 4e-8)
    plt.show()
    plt.savefig('three_level_solution.png', dpi=300)


def omega_L_cm(om):
    """
    Take laser frequency in 1/cm and convert to Hz
    """
    return om.to('Hz', equivalencies=u.spectral()) * u.cycle


if __name__ == "__main__":
    # TODO: this is just for testing.
    #OMEGA_L = 28503.01 / u.cm
    #t, rho_array = calculate_populations(omega_L_cm(OMEGA_L))
    laser_scan_wavenum = np.arange(28500, 28510, 0.5) / u.cm
    # metastable population as a function of laser frequency
    rho3_scan = []
    for omega_laser in omega_L_cm(laser_scan_wavenum):
        t, rho_array = calculate_populations(omega_laser)
        #print(rho_array.shape)
        rho3_scan.append(rho_array[-1, -1])
        #plot_optical_pumping(rho_array, t)
        #break
    #rho_array, t = calculate_populations(omega_L)
    #plot_optical_pumping(rho_array, t)
    fig, ax = plt.subplots(figsize=(6,6))
    ax.plot(laser_scan_wavenum, rho3_scan)
    plt.show()