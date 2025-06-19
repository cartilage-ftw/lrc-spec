from astropy.constants import c, hbar, k_B, m_p
from scipy.integrate import quad, solve_ivp
from scipy.special import voigt_profile
from functools import lru_cache
import scipy.special
from sympy.physics.wigner import wigner_6j
from pathos.multiprocessing import ProcessingPool

import multiprocessing as mp
import time
import numpy as np
import astropy.units as u
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.formatter.useoffset'] = False
plt.rcParams['figure.dpi'] = 150
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
plt.rcParams['font.size'] = 13

# to later use all CPU cores using `pathos`
num_cores = mp.cpu_count()

"""
The important thing is to use astropy.units as it's easy to get a factor of \hbar or c wrong
and forgetting to convert from freq to angular freq, etc. Also, the differential equations
are "stiff" for a pulsed laser (i.e. when the intensity changes as a step function
with each pulse). You can use an implicit solver (e.g. 'LSODA' in scipy)
The number of points it needs to sample at can be reduced dramatically 

Note that the relation betweeen A_21 and B_21 coefficients depends on the definition
of laser intensity/energy density chosen.
"""

# NOTE: Something is off; I need to suppress the intensity by *quite a lot* to not get flat top profiles
ENERGY_PER_PULSE = 2.4E-9 * u.J / 1e10
# TODO: OD filters

BUNCHER_FREQ = 1000 * u.Hz
LASER_PULSE_DUR = 8E-9 * u.s
LASER_SPOT_SIZE = 1.5 * u.mm**2

LASER_REP_FREQ = 8000 * u.Hz
RESONANCE_FREQ = (np.array([28502.3, 28503.01, 28503.7]) / u.cm).to('Hz', equivalencies=u.spectral())
OMEGA_12 = RESONANCE_FREQ * u.cycle
LASER_FWHM_GAUSS = 4.6E9 * u.Hz # NOTE: u.GHz wasn't doing the right thing

# TODO: this is just for testing.
#OMEGA_L = OMEGA_12 #- (20E9 * u.Hz * u.cycle)
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
p_DT = 5 * u.mbar

@lru_cache(maxsize=None)
def A_FF(A_J, J_u, J_l, F_u, F_l, I):
    """
    The Einstein A_21 coefficient for a hyperfine transition from F' -> F
    in terms of A_21 coefficient for the fine structure transition.
    This takes care of the Clebsch-Gordon coefficients/Wigner nj symbols 
    """
    #TODO: double check the prefactors.
    return (2*F_l + 1)*(2*F_u + 1) * wigner_6j(J_l, I, F_l, F_u, 1, J_u)**2 * A_J

@lru_cache(maxsize=None)
def pulse_func(t: float, tau_pulse=LASER_PULSE_DUR.to("s").value, f_rep=LASER_REP_FREQ.to("Hz").value,
                buncher_cycle_duration=1/BUNCHER_FREQ.to("Hz").value):
    """
    A function that mimics a laser pulse of a certain duration
    While satisfying the property that area under the curve for each rectangular pulse is 1
    """
    # separation between each pulse in ns
    pulse_sep = (1/f_rep)#.to('s').value
    # such that height*width = 1
    height = 1/tau_pulse#.to('s').value 
    # flag whether this time 't' is within the buncher cycle duration, or after
    # (there's no radiation after leaving buncher)
    t_within_buncher = np.where(t < buncher_cycle_duration, 1, 0)
    return t_within_buncher * np.where(((t % pulse_sep) < tau_pulse), height, 0)

@lru_cache(maxsize=None)
def laser_energy_density(omega, omega_L, per_pulse_energy):
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
    laser_sigma = (LASER_FWHM_GAUSS / 2.355).to('Hz').value
    integrated_energy_density = per_pulse_energy*LASER_REP_FREQ.value/(c*LASER_SPOT_SIZE)
    return integrated_energy_density * np.exp(-(omega - omega_L)**2 / (2*laser_sigma**2))

print("Note: Collisional broadening has been ignored.")
@lru_cache(maxsize=None)
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
        sigma_doppler = sigma_doppler.to('Hz', equivalencies=[(u.cy/u.s, u.Hz)]).value
    # NOTE: the collision broadening should include both elastic and inelastic collision rates
    # Here, what I've included is only the inelastic part.
    # TODO: either estimate the elastic collision rate from interaction potential, or "fit" it
    gamma = 1e9 # 2*collision_rates
    #print("Inelastic collision rates", collision_rates)
    #d = # diameter of atoms
    #mean_collision_time = (1/(n_He * d**2))*np.sqrt(M_175/(16*np.pi * k_B *T))
    #gamma = 2/mean_collision_time  
    
    #omega *= u.cycle / u.s
    #print("Dimeinsions of omega_res", omega_res)
    #phi = lambda omega: voigt_profile(omega-omega_res, sigma_doppler, gamma)
    #gamma = 1
    #print("evaluating atomic lineshape")
    #sigma_doppler = 100_000_000
    #phi = scipy.stats.norm(loc=omega_21, scale=sigma_doppler)
    #return np.exp(-(omega-omega_21)/(2*(sigma_doppler**2)))
    return scipy.special.voigt_profile(omega-omega_21, sigma_doppler, gamma)#.value

@lru_cache(maxsize=None)
def S_factor(omega_12, omega_L, energy_per_pulse, collision_rates, doppler_sigma):
    """
    The arguments of this method had to be dimensionless to be cache-able
    """
    combined_lineshape = lambda omega: (atomic_lineshape(omega, omega_12, collision_rates, doppler_sigma)\
                                        *laser_energy_density(omega, omega_L, energy_per_pulse))#.value
    units = combined_lineshape(omega_12).unit
    # for integrating with scipy
    dimless_lineshape = lambda omega: combined_lineshape(omega).value
    # NOTE: the integral shouldn't be done from 0, np.inf because the numerical step size would not
    # resolve the lineshape at all.
    #print(omega_12, omega_L, LASER_FWHM_GAUSS)
    #print(dimless_lineshape(864e14))
    #lower = 
    #upper = np.max([omega_12, omega_L])
    '''omega = np.linspace(np.min([omega_12, omega_L]) - 5*LASER_FWHM_GAUSS.value,
                        np.max([omega_12, omega_L]) + 5*LASER_FWHM_GAUSS.value,
                        100000000)
    lineshape = dimless_lineshape(omega)'''
    #print(lineshape)
    #int_lineshape = np.trapz(omega, lineshape)
    int_lineshape = quad(dimless_lineshape,
                        np.min([omega_12, omega_L]) - 5*LASER_FWHM_GAUSS.value,
                        np.max([omega_12, omega_L]) + 5*LASER_FWHM_GAUSS.value,
                        #method='tanh-sinh',
                        #maxdegree=100,
                        #n=10000000,
                        #tol=1e-12,
                        #rtol=1e-12
                        #n=1
                        #limit=10000000,
                        #rtol=1e-10,
                        #tol=1e-10
                        #epsrel=1e-12,
                        #epsabs=1e-12
                    )[0]
    #print("Value of integrated lineshape", int_lineshape)
    omega_12 *= 2*np.pi #.to('Hz', equivalencies=u.dimensionless_angles())
    # NOTE: the 2*pi factor in the angular frequency must be correctly included, otherwise
    # the value will be wrong by 8pi^3. I first tried doing this dynamically
    # But then the cache would not work because the arguments were not dimensionless
    # so a unit in 'Hz' is passed to the arguments, without the 2pi factor
    return int_lineshape * (((np.pi**2 * c**3)/(hbar*omega_12**3)) * units).si.value


'''omega = np.linspace(OMEGA_12.value - 5*LASER_FWHM_GAUSS.value, OMEGA_12.value + 5*LASER_FWHM_GAUSS.value, 100) * u.cycle * u.Hz
#print("Value of S-factor", S_factor(1, OMEGA_12, OMEGA_L, A31))
S = []
for o in omega:
    S.append(S_factor(OMEGA_12[1].value, o[0].value, ENERGY_PER_PULSE, A31, DOPPLER_WIDTH_SIGMA[1].value))
fig, ax = plt.subplots(figsize=(6,6))
half_max = np.max(S)/2
ax.plot((omega-OMEGA_12)/1E9, S)
ax.axhline(half_max, xmin=0, xmax=1, c='gray', ls='--')
ax.set_xlabel("Frequency detuning [GHz]")
ax.set_ylabel("S-factor [dimensionless]")
print("Width of it", 2.355*np.std(S)/1E9)
plt.show()'''


def rate_equations(t, rho_vec, A21, A31, A23, energy_per_pulse, alpha_31,
                    g2, g1, F2, F1, I, J_u, J_l,
                    omega_12, omega_L):
    """
    The rate equations are as derived in eqn. 3.7 to 3.9 in Aayush's master thesis.
    """
    rho_1, rho_2, rho_3 = rho_vec
    omega_12 = omega_12.to('Hz', equivalencies=[(u.cy/u.s, u.Hz)]).value
    omega_L = omega_L.to('Hz', equivalencies=[(u.cy/u.s, u.Hz)]).value
    doppler_widths = DOPPLER_WIDTH_SIGMA.to('Hz', equivalencies=[(u.cy/u.s, u.Hz)]).value

    # absorption - stimulated emission (multiplied by energy density, integrated lineshape function etc.)
    abs_stim_term = 0.
    # If there are multiple hyperfine multiplets, the location of the resonance omega_12 is an array,
    # and each of those has an associated S-factor. The sum over all these S-factors is taken to take into account
    # pumping neighbouring hyperfine levels at the same time.
    #A_sum = 0.
    for i in range(len(omega_12)):
        abs_stim_term += (rho_2 - (g2[i]/g1[i])*rho_1) \
                            *S_factor(omega_12[i], omega_L, energy_per_pulse, alpha_31, doppler_widths[i]) \
                            * A21 \
                            *pulse_func(t) 
        #print(f"Is A(F'->F) = A21 for {F2[i]} -> {F1[i]}",
        #                A21, A_FF(A21, J_u, J_l, F2[i], F1[i], I))
        #A_sum += A_FF(A21, J_u, J_l, F2[i], F1[i], I)
    #print("Sum of A_FF and A_21", A_sum, A21)
    #exit
    if np.abs(abs_stim_term) > 1e10: # if the value is too large, it's already saturated
        r = np.rint(np.log10(np.abs(abs_stim_term))) - 10
        abs_stim_term /= 10**r
        # for numerical safety, cap it here. It's not going to change the solution as A21+A23 ~10^7 << 10^10

    #print("No crashes until evaluation of absorption/stimulated emission")

    # the derivatives in the rate equation
    drho1_dt = A21*rho_2 + A31*rho_3 \
                + abs_stim_term \
                + alpha_31*rho_3 # collisional de-excitation from metastable
    drho2_dt = -(A21 + A23)*rho_2 - abs_stim_term # de-excitation from 2->3 has been neglected
    #print("Successfully evaluated first two derivatives")
    drho3_dt =  A23*rho_2 - (A31 + alpha_31)*rho_3
    return np.array([drho1_dt, drho2_dt, drho3_dt])#/1e12

NUM_PULSES = round((LASER_REP_FREQ/BUNCHER_FREQ).value) #+ 0.001/50
tau = LASER_PULSE_DUR.to('s').value
f = LASER_REP_FREQ.to('Hz').value
events = lambda t, rho, *args: np.sin((t-tau/2)*np.pi*(f))
events.terminal = False
events.direction = 1

t_pulses = np.arange(NUM_PULSES) / LASER_REP_FREQ.value

J_u = 1
J_l = 0
F_u = np.array([-1, 0, 1]) + NUCLEAR_SPIN
F_l = np.zeros_like(F_u) + NUCLEAR_SPIN

g_u = 2*F_u + 1
g_l = 2*F_l + 1

# assumed as per Kim et al (2024)
rate_coefficient = 2E-14 * u.cm**3 / u.s # 
# P is the pressure (in DT or buncher, wherever you want to use this)
inelastic_collision_rate = lambda P: (rate_coefficient * P/(k_B * T)).to("1/s")


def calc_populations(omega_L, energy_per_pulse, arrival_time_gs):
    # initial level populations.
    rho0 = [1.0, 0.0, 0.0]
    # and variables to store the full solution
    rho_array = np.array([rho0])
    t_array = [0.]
    print("Initializing solution of rate equations")
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
                        energy_per_pulse,
                        inelastic_collision_rate(p_buncher).to('1/s').value, # NOTE: this will evaluate the collisions in the bunhcer only
                        g_u,
                        g_l, F_u, F_l, NUCLEAR_SPIN, J_u, J_l, 
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
        
    # now, the "optical pumping" is done. Let's evolve the system in the drift tube
    # where pressure is higher, until they arrive at the detector
    t_buncher = (1/BUNCHER_FREQ).to('s').value
    if isinstance(arrival_time_gs, u.Quantity):
        arrival_time_gs = arrival_time_gs.to('s').value
    rho = solve_ivp(rate_equations,
                        t_span=(t_buncher, t_buncher + arrival_time_gs), # t_span
                        y0=rho0,
                        method='LSODA',
                        args=(A21, A31, A23,
                        energy_per_pulse,
                        # use the collision rate appropriate for the DT
                        inelastic_collision_rate(p_DT).to('1/s').value,
                        g_u,
                        g_l, F_u, F_l, NUCLEAR_SPIN, J_u, J_l, 
                        OMEGA_12,
                        omega_L),
                        atol = 1e-7,
                )
    rho_array = np.append(rho_array, rho.y.T, axis=0)
    t_array = np.append(t_array, rho.t, axis=0)

    # and NOW we're done!
    print("Successfully integrated rate equations for", NUM_PULSES, "laser pulses at wavenum",
          (omega_L/u.cycle).to("1/cm", equivalencies=u.spectral()), "for pulse energy", energy_per_pulse)
    
    t = np.array(t_array)
    return t, rho_array
    #print("Values in t_array", t/1E-6)


def plot_populations(t, rho_array):
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
    max_ms = np.max(rho_array[:,2])
    #axes[2].plot(t_all*1E6, pulse_func(t_all)/1.25E8*max_ms, 'tab:pink',label='laser')
    axes[2].plot(t*1E6, pulse_func(t)/1.25E8*max_ms, 'tab:green', label='laser', alpha=0.5)

    # shade the part of the time spend in the drift tube
    t_buncher = (1/BUNCHER_FREQ).to('us').value # the plot is in microseconds
    t_DT = np.linspace(t_buncher, t[-1]*1E6, 100)
    axes[1].text(x=t_DT[50], y=np.max(rho_array[:,1])/2, s=r"\textsc{Drift Tube}", c='dimgray',
                     va='center', ha='center')
    
    # some characteristics
    axes[0].set_title(
        # the He-Lu^+ inelastic collision rate coefficient (at 300K)
        r"$\alpha_{31}=$" + f'{rate_coefficient.to("cm3 / s").value:.1e}' +" cm$^{3}/$s"
        + ",\quad" # 
        + "$p(\mathrm{DT})=" + f"{p_DT.to('mbar').value:.1f}" + "$ mbar")
    for i, ax in enumerate(axes):
        shade_ymin = np.zeros_like(t_DT)
        shade_ymax = np.max(rho_array[:,i])*np.ones_like(t_DT)

        ax.fill_between(t_DT, shade_ymin, shade_ymax, alpha=0.5, fc='silver')
        ax.axvline(x=t_buncher, ymin=0, ymax=1, c='gray', ls='--', lw=1)
        ax.set_xlabel('Time [$\mu$s]')
        ax.grid()
        ax.legend()
        #ax.set_xscale('log')

    plt.tight_layout()
    plt.scatter(t*1E6, np.zeros_like(t), alpha=0.1)
    #plt.scatter(rho.t_events[0]/1e12*1e6, np.zeros_like(rho.t_events[0]), marker='x')
    #plt.xlim(-4e8, 4e-8)
    plt.savefig('three_level_solution.png', dpi=300)
    plt.show()


def plot_spectra(energies, wavenums, ms_all):
    fig, ax = plt.subplots()
    plateau = np.max(ms_all)

    #alphas= np.linspace(0.6, 0.8, len(energies))
    colors = mpl.colormaps['cividis'](np.linspace(0, 1.0, len(energies)))
    #print("The plateau occured at", plateau)
    #print("Length of energies", len(energies))
    for i in range(len(energies)):
        #print("Was this hit?")
        ax.plot(wavenums, ms_all[i,:], c=colors[i], lw=1., alpha=0.8,
                label=f"{energies[i].to('nJ').value*1E9:.1e}" + r" $\times 10^{-9}$ nJ/pulse")
    ax.axhline(y=plateau, xmin=0, xmax=1, ls='--', c='dimgray', lw=1.)

    ax.text(x=28501.52, y=plateau*1.03,
            s='Plateau due to collisional de-excitation',
            c='k',
            ha='left',
            # va='top', ha='right'
        )
    ax.fill_between(wavenums.to('1/cm').value, plateau*1.02, plateau*0.98, fc='silver', alpha=0.3, zorder=-1)
    # some characteristics
    ax.set_title(
        # the He-Lu^+ inelastic collision rate coefficient (at 300K)
        r"$\alpha_{31}=$" + f'{rate_coefficient.to("cm3 / s").value:.1e}' +" cm$^{3}/$s"
        + ",\quad" # 
        + "$p(\mathrm{DT})=" + f"{p_DT.to('mbar').value:.1f}" + "$ mbar")
    # Shrink current axis's height by 10% on the bottom
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0 + box.height * 0.1,
    #                box.width, box.height * 0.9])
    #ax.legend(bbox_to_anchor=(0.05, -0.1), loc="upper left", ncol=3)
    ax.set_ylim(top=plateau*1.1)
    ax.set_xlabel("Wavenumber [cm$^{-1}$]")
    ax.set_ylabel("Metastable Fraction")
    plt.tight_layout()
    plt.savefig("theoretical_spectra_5mbar.png", dpi=300)
    plt.show()

ARRIVAL_TIME_GS = 450 * u.us # microseconds
wavenum_range = np.linspace(28501.5, 28504.5, 80) / u.cm
omega_laser = lambda wavenum: wavenum.to('Hz', equivalencies=u.spectral()) * u.cycle
#omega_L = OMEGA_12 #- (15E9 * u.cycle * u.Hz)

energies = ENERGY_PER_PULSE*np.logspace(-4, -1.5, 6)
ms_final = []
rho_all = []
#t_all = []

grid = [(E, omega_L) for E in energies for omega_L in omega_laser(wavenum_range)]

t_i = time.time()

def calc_metastable_pop(energy_per_pulse, omega_L):
    t, rho_array = calc_populations(omega_L, energy_per_pulse, ARRIVAL_TIME_GS)
    return rho_array[-1, 2]

with ProcessingPool(num_cores) as pool:
    calc_results = pool.map(lambda grid_pars: calc_metastable_pop(*grid_pars), grid)

metastable_pop = np.array(calc_results).reshape((len(energies), len(wavenum_range)))

'''for per_pulse_energy in energies:
    ms_this_wavenum = []
    for omega_L in omega_laser(wavenum_range):
        t, rho_array = calc_populations(omega_L, per_pulse_energy, ARRIVAL_TIME_GS)
        ms_this_wavenum.append(rho_array[-1, 2])
    ms_final.append(ms_this_wavenum)
    rho_all.append(rho_array)'''
    #t_all.append(t)
    #plot_populations(t, rho_array)

print("Everything took", time.time() - t_i, "seconds")
#t_all = np.array(t_all)
ms_final = metastable_pop#np.array(ms_final)

plot_spectra(energies, wavenum_range, ms_final)