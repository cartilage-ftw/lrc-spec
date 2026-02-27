import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
import time
import matplotlib as mpl
from scipy.special import voigt_profile
from astropy.constants import c, h, k_B

mpl.rcParams['axes.formatter.useoffset'] = False
plt.rcParams['figure.dpi'] = 150
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
plt.rcParams['font.size'] = 13

# Constants in SI units
C_SI = c.si.value
H_SI = h.si.value
KB_SI = k_B.si.value

# Experimental parameters for laser, buncher etc.
LASER_PULSE_ENERGY = 2.4e-9  # J (2.4 nJ)
LASER_SPOT_SIZE = 1.5e-6     # m^2 (1.5 mm^2)
LASER_PULSE_DUR = 8e-9       # s
LASER_BANDWIDTH_FWHM = 4.5e9 # Hz (4.5 GHz)
LASER_SIGMA = LASER_BANDWIDTH_FWHM / 2.355

# Timing
LASER_REP_FREQ = 8000.0  # Hz
BUNCHER_FREQ = 100.0 # Hz
NUM_PULSES = int(LASER_REP_FREQ / BUNCHER_FREQ)
T_CYCLE = 1.0 / LASER_REP_FREQ
T_OFF = T_CYCLE - LASER_PULSE_DUR

# Atomic Parameters (Lu+)
M_175 = 175 * 1.66e-27 # kg
TEMP = 300.0 # K
A21 = 1.25e7
A23 = 1.88e7
A31 = 1.0e-2 # Negligible compared to collisions
ALPHA_31 = 8e-14 * 1e-6 # cm3/s -> m3/s
P_BUNCHER = 2.5e-2 * 100 # mbar -> Pa

# Resonance Lines (cm-1)
RESONANCES_CM = np.array([28502.3, 28503.01, 28503.7])
# Hyperfine Degeneracies (2F + 1)
# Note: Ensure these match the specific transitions you are modeling
g_u = np.array([6, 8, 10]) # Example: F_u = 5/2, 7/2, 9/2 -> 2F+1 = 6, 8, 10
g_l = np.array([8, 8, 8])  # Example: F_l = 7/2 -> 8

def get_rate_matrix(omega_laser_hz, per_pulse_energy, pressure_pa):
    """
    Constructs the Rate Matrix M.
    This version uses the analytic Voigt convolution for speed and stability.
    """
    # 1. ATOMIC FREQUENCIES (SI Units)
    # FIX 1: Multiply wavenumber by 100 to get m^-1
    nu0_list = RESONANCES_CM * 100 * C_SI 
    
    # 2. INTENSITY PREFACTOR
    # Mean Intensity J_bar prefactor: E / (Area * Time * 4pi)
    # This has units of [W/m^2] integrated intensity equivalent
    J_prefactor = per_pulse_energy / (4 * np.pi * LASER_SPOT_SIZE * LASER_PULSE_DUR)
    
    # 3. CALCULATE RATES
    W_12_total = 0.0
    W_21_total = 0.0
    
    # Doppler Sigma (Atomic)
    # sigma_D = (nu0 / c) * sqrt(kT / m)
    sigma_doppler = (nu0_list[0] / C_SI) * np.sqrt(KB_SI * TEMP / M_175)
    
    # Combined Sigma for Convolution (Laser + Doppler)
    sigma_total = np.sqrt(LASER_SIGMA**2 + sigma_doppler**2)
    gamma = (A21 + A23) / 2 # Approx natural linewidth (HWHM)
    gamma += 1e9 # broadening due to everything else.
    for i, nu0 in enumerate(nu0_list):
        # Einstein Coefficient Ratio (Specific Intensity definition)
        # B = A / (2 h nu^3 / c^2)
        ratio_A_B = (2 * H_SI * nu0**3) / (C_SI**2)
        B21 = A21 / ratio_A_B
        
        # B. Mean Intensity Overlap
        # Voigt profile at detuning (nu_laser - nu_atom)
        # Note: voigt_profile is normalized to integral=1 over frequency
        detuning = omega_laser_hz - nu0
        overlap = voigt_profile(detuning, sigma_total, gamma)
        
        # J_bar_line [W/m^2 * s] (Spectral Intensity effective)
        # Units check: [W/m^2] * [1/Hz] = [J/m^2]
        J_bar_spectral = J_prefactor * overlap
        
        # Rate = B * J_bar
        # Units check: [m^2/J] * [J/m^2] = [s^-1]
        rate_stim = B21 * J_bar_spectral
        
        W_21_total += rate_stim
        W_12_total += rate_stim * (g_u[i] / g_l[i])

    # COLLISIONAL RATE
    # rate = alpha * density; density = P / kT
    coll_rate = ALPHA_31 * (pressure_pa / (KB_SI * TEMP))
    
    # BUILD MATRICES
    # M_off (No Laser)
    M_off = np.zeros((3,3))
    M_off[0, 1] = A21
    M_off[0, 2] = A31 + coll_rate
    M_off[1, 1] = -(A21 + A23)
    M_off[2, 1] = A23
    M_off[2, 2] = -(A31 + coll_rate)

    # M_on (Laser)
    M_on = M_off.copy()
    M_on[0, 0] -= W_12_total
    M_on[0, 1] += W_21_total
    M_on[1, 0] += W_12_total
    M_on[1, 1] -= W_21_total
    
    return M_on, M_off

def solve_system(wavenum_grid, per_pulse_energy):
    # Vectorized loop for speed
    ms_populations = []
    
    # Pre-compute time steps
    # We assume constant pressure/rates throughout the pulse train for simplicity here
    # If pressure changes (DT vs Buncher), we construct matrices accordingly
    
    # It is faster to loop over frequency than to vectorize the linear algebra
    # because matrix exponentiation is not natively vectorized in scipy.
    
    tau = LASER_PULSE_DUR
    
    for wn in wavenum_grid:
        # Convert wavenumber to Frequency (Correctly!)
        freq_hz = wn * 100 * C_SI 
        
        M_on, M_off = get_rate_matrix(freq_hz, per_pulse_energy, P_BUNCHER)
        
        # Propagators
        U_on = scipy.linalg.expm(M_on * tau)
        U_off = scipy.linalg.expm(M_off * T_OFF)
        U_cycle = U_off @ U_on
        
        # Total evolution for N pulses (Binary Exponentiation)
        U_total = np.linalg.matrix_power(U_cycle, NUM_PULSES)
        
        # Initial State: |1>
        rho_0 = np.array([1.0, 0.0, 0.0])
        rho_final = U_total @ rho_0
        
        ms_populations.append(rho_final[2])
        
    return np.array(ms_populations)

if __name__ == "__main__":
    # Scan Range
    wavenums = np.linspace(28501.5, 28504.5, 300)
    
    t0 = time.time()
    ms_pop = solve_system(wavenums, LASER_PULSE_ENERGY)
    t1 = time.time()
    ms_pop_enhanced = solve_system(wavenums, 10*LASER_PULSE_ENERGY)
    
    print(f"Computed {len(wavenums)} points in {t1-t0:.3f} seconds.")
    
    # Plot
    plt.plot(wavenums, ms_pop, lw=1.5, label=f'E={LASER_PULSE_ENERGY*1e9} nJ')
    plt.plot(wavenums, ms_pop_enhanced, lw=1.5, label=f'E={10*LASER_PULSE_ENERGY*1e9} nJ',
             c='gray')
    #plt.title(f"Metastable Pop. (E={LASER_PULSE_ENERGY*1e9} nJ)")
    plt.xlabel("Wavenumber [cm$^{-1}$]")
    plt.ylabel("Population $\\rho_3$")
    #plt.grid(True, alpha=0.3)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()