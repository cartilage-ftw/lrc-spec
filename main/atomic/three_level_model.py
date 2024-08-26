"""
    Three level rate equation model
    author: Aayush
"""

from scipy.integrate import odeint, quad, solve_ivp
from scipy.constants import c, hbar, pi, k # boltzmann constant
import numpy as np
import matplotlib.pyplot as plt


plt.rcParams['figure.dpi'] = 150
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
plt.rcParams['font.size'] = 13
plt.rcParams['axes.formatter.useoffset'] = False


E = 2.4E-9 # joules

LASER_PULSE_DUR = 8E-9 # in secs
LASER_REP_FREQ = 8_000 # Hz
LASER_SPOT_SIZE = 1E-6 # 1 mm^2

LASER_FREQ = 28503.01 * c/100 # in Hz
LASER_ANG_FREQ = LASER_FREQ * 2*pi

LASER_GAUSS_FWHM = 4.6E9


def pulse_func(t, tau_pulse=LASER_PULSE_DUR, f_rep=LASER_REP_FREQ):
    """
    A function that mimics a laser pulse of a certain duration
    While satisfying the property that area under the curve for each rectangular pulse is 1
    """
    # separation between each pulse in ns
    pulse_sep = 1/f_rep
    
    height = 1/tau_pulse # such that height*width = 1
    return np.where( (t % pulse_sep) < tau_pulse, height, 0)


def integrate_laser_profile(omega_L:float, f_FWHM:float = 4.6E9) -> float:
    """Mimics the laser line profile
    Args:
        omega_L (float): The central (angular) freq of the laser while operating
        f_FWHM (float): the Gaussian FWHM [GHz] of the laser profile, in frequency units (not angular)
    Returns:
        float: Value of the profile evaluated at that point
    """
    # convert into angular units
    omega_FWHM = 2*pi * f_FWHM
    # FWHM to sigma
    sigma = omega_FWHM/2.355

    area_under_gaussian = sigma * np.sqrt(2*pi)
    return area_under_gaussian


INTEGRATED_LASER_PROFILE = integrate_laser_profile(omega_L=LASER_ANG_FREQ, f_FWHM=LASER_GAUSS_FWHM)

def u(#omega:float, omega_L:float,
      t:float):
    return (E/(LASER_SPOT_SIZE*LASER_PULSE_DUR))*INTEGRATED_LASER_PROFILE*pulse_func(t)


def calc_collision_rate(alpha, pressure, T):
    # assume ideal gas pV = NkT => N/V = p/kT
    p = pressure*100 # convert mbar to pascals
    num_density = pressure/(k*T) # per cubic metres
    # NOTE: Remember to convert from m^-3 to cm^-3
    alpha_cm3 = alpha*1E-6
    return alpha_cm3*num_density
    

def atomic_lineshape(omega_12):
    return 1


step = 0
def rate_equation_model(rho_vec, t, A21, A31, A23, alpha_31, g2, g1, omega_12):
    global step
    step += 1
    total_steps = int(duration/step_size)
    if step % 10000 == 0:
        print(f'Calculating step {step}, {100*step/total_steps:.3f}% done')
    # for the sake of readability, I expand these
    rho_1, rho_2, rho_3 = rho_vec
    
    wavelength = c/(omega_12/(2*pi))
    
    abs_stim_term = (rho_2 - (g2/g1)*rho_1) * (A21*wavelength**2/(2*hbar*omega_12)) \
                                                *u(t)*atomic_lineshape(omega_12)
    drho1_dt = A21*rho_2 + A31*rho_3 \
                        + abs_stim_term \
                        + alpha_31*rho_3
    drho2_dt= -(A23 + A21)*rho_2 - abs_stim_term# ignore quenching from 2->3
    drho3_dt = A23*rho_2 -A31*rho_3 - alpha_31*rho_3
    #print("Called!")
    return [drho1_dt, drho2_dt, drho3_dt]


duration = 1.00001E-3
step_size=1E-9

t = np.linspace(0, duration, int(duration/step_size))

laser_pulses = pulse_func(t)
J_u = 1
J_l = 0
alpha_31 = 1E-11
collision_rate_buncher = calc_collision_rate(alpha_31, pressure=5E-2, T=300)# 

rho1 = np.empty_like(t)
rho2 = np.empty_like(t)
rho3 = np.empty_like(t)

rho0 = [1.0, 0.0, 0.0]
rho1[0] = rho0[0]
rho2[0] = rho0[1]
rho3[0] = rho0[2]

for i in range(1,len(t)):
    t_span = [t[i-1], t[i]]
    rho = odeint(rate_equation_model, rho0, t_span, args=(
                          1.21E7, # = A21
                          3.6E-3, # = A31 NOTE: random value; assumes 1 hour lifetime of metastable
                          1.88E7, # = A23
                          collision_rate_buncher, 
                          2*J_u + 1, # = g2
                          2*J_l + 1, # = g1
                          LASER_ANG_FREQ # = omega_12 NOTE: atomic resonance and laser freq need not be the same
            ))
    rho1[i] = rho[1][0]
    rho2[i] = rho[1][1]
    rho3[i] = rho[1][2]
    rho0 = rho[1]


'''result = solve_ivp(rate_equation_model,
                          t_span=[0, duration],
                          t_eval=t,
                          y0=[1.0, 0.0, 0.0],
                          method='BDF', 
                          args=(
                          1.21E7, # = A21
                          3.6E-3, # = A31 NOTE: random value; assumes 1 hour lifetime of metastable
                          1.88E7, # = A23
                          collision_rate_buncher, 
                          2*J_u + 1, # = g2
                          2*J_l + 1, # = g1
                          LASER_ANG_FREQ # = omega_12 NOTE: atomic resonance and laser freq need not be the same
                          ))'''


fig, axes = plt.subplots(3, 1, figsize=(6,8))


axes[0].plot(t*1E3, rho1, ls='--', label=r'$\rho_1$')
axes[1].plot(t*1E3, rho2, label=r'$\rho_2$')
axes[2].plot(t*1E3, rho3, label=r'$\rho_3$')
axes[2].plot(t*1E3, laser_pulses/1.25E8, c='orange', alpha=0.4, label='laser')

axes[0].set_ylabel(r'$\rho_1$')
axes[1].set_ylabel(r'$\rho_2$')
axes[2].set_ylabel(r'$\rho_3$')
for ax in axes:
    ax.set_xlabel('Time [ms]')
    ax.legend()

plt.tight_layout()
plt.savefig('three_level_solution.png', dpi=300)