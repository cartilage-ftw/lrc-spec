import numpy as np
import pandas as pd
import astropy.units as u
from astropy.constants import m_p, h, c, hbar
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from scipy.special import spherical_jn, spherical_yn, lpmv
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.rcParams.update({
    'figure.dpi':150,
    'text.usetex': True,
    'font.family':'serif',
    'font.serif': ['Computer Modern Serif'],
    'font.size': 13,
})

hbar = hbar.si.value
h = h.si.value
angstrom_to_m = 1e-10          # Å → m
cm1_to_J      = (h * c * 100).si.value  # cm⁻¹ → J


def load_potential(path):
    if path.endswith('.dat'):
        df = pd.read_csv(path, sep=',',  skiprows=2,
                         names=['empty', 'r_A', 'V0_cm1', 'Vp_cm1', 'Vm_cm1'])
        r = df['r_A'].values * angstrom_to_m
        V = df['V0_cm1'].values * cm1_to_J
    else:
        data = np.loadtxt(path, delimiter=',', skiprows=1)
        r = data[:, 0] * angstrom_to_m
        V = data[:, 1] * cm1_to_J
    sort_idx = np.argsort(r)
    print("Reading", path)
    print(r)
    print(V)
    #exit()
    return interp1d(r[sort_idx], V[sort_idx], kind='cubic', fill_value='extrapolate')


def numerov_step(u_im1, u_i, Q_im1, Q_i, Q_ip1, dr):
    return (2*u_i*(1 - dr**2 * Q_i/12) - u_im1*(1 + dr**2 * Q_im1/12)) / (1 + dr**2 * Q_ip1/12)

def integrate_outward(r, E, mu, V_func, l):
    dr = r[1] - r[0]
    N = len(r)
    u = np.zeros(N)
    u[0] = r[0]**(l+1);  u[1] = r[1]**(l+1)
    Q = (2*mu/hbar**2)*(E - V_func(r)) - l*(l+1)/r**2
    for i in range(1, N-1):
        u[i+1] = numerov_step(u[i-1], u[i], Q[i-1], Q[i], Q[i+1], dr)
    return u

def integrate_inward(r, E, mu, V_func, l):
    dr = r[1] - r[0]
    N = len(r)
    u = np.zeros(N)
    kappa = np.sqrt(2*mu*abs(E))/hbar
    u[-1] = 0.0
    u[-2] = np.exp(-kappa * r[-2])
    Q = (2*mu/hbar**2)*(E - V_func(r)) - l*(l+1)/r**2
    for i in range(N-2, 0, -1):
        u[i-1] = numerov_step(u[i+1], u[i], Q[i+1], Q[i], Q[i-1], dr)
    return u


def shoot(E, r, mu, V_func, l):
    u_out = integrate_outward(r, E, mu, V_func, l)
    u_in  = integrate_inward(r, E, mu, V_func, l)
    m = len(r)//2
    dr = r[1] - r[0]
    dout = (u_out[m+1] - u_out[m-1])/(2*dr) / u_out[m]
    din  = (u_in[m+1]  - u_in[m-1]) /(2*dr) / u_in[m]
    return dout - din

def find_bound_states(r, mu, V_func, l, E_min, E_max, nroots=1):
    roots = []
    E_br = np.linspace(E_min, E_max, 200)
    s0 = shoot(E_br[0], r, mu, V_func, l)
    for i in range(1, len(E_br)):
        s1 = shoot(E_br[i], r, mu, V_func, l)
        if s0*s1 < 0 and len(roots)<nroots:
            sol = root_scalar(shoot, args=(r, mu, V_func, l),
                              bracket=[E_br[i-1], E_br[i]], maxiter=50)
            roots.append(sol.root)
        s0 = s1
    return roots


def numerov_scatter(u0, u1, r, k2, mu, V_func, l):
    dr = r[1] - r[0]
    N = len(r); u = np.zeros(N)
    u[0], u[1] = u0, u1
    Q = (2*mu/hbar**2)*(k2*hbar**2/(2*mu) - V_func(r)) - l*(l+1)/r**2
    for i in range(1, N-1):
        u[i+1] = numerov_step(u[i-1], u[i], Q[i-1], Q[i], Q[i+1], dr)
    return u

def compute_phase_shift(r, u, k, l):
    dr = r[1] - r[0]; idx = -2
    rm, um = r[idx], u[idx]; up, um1 = u[idx+1], u[idx-1]
    uprime = (up - um1)/(2*dr)
    jl  = spherical_jn(l, k*rm); nl  = spherical_yn(l, k*rm)
    jl_p= spherical_jn(l, k*rm, derivative=True); nl_p= spherical_yn(l, k*rm, derivative=True)
    num = um*jl_p - uprime*jl; den = uprime*nl - um*nl_p
    return np.arctan2(num, den)

def differential_cross_section(thetas, k, deltas):
    cos_t = np.cos(thetas); f = np.zeros_like(thetas, dtype=complex)
    for l, delta in enumerate(deltas):
        P_l = lpmv(0, l, cos_t)
        f += (2*l+1)*np.exp(1j*delta)*np.sin(delta)*P_l
    return np.abs(f)**2/k**2

def total_cross_section(k, deltas):
    l = np.arange(len(deltas))
    return 4*np.pi/k**2 * np.sum((2*l+1)*np.sin(deltas)**2)

# =============================================================================
# Main execution
# =============================================================================
# potentials
V_gs = load_potential('../structure/Archiv2/Lu+inGS.csv')
V_ms = load_potential('../structure/Archiv2/SO3D1_.dat')
# reduced mass
m_Lu = 175 * m_p.si.value; m_He = 4 * m_p.si.value
mu = (m_Lu * m_He)/(m_Lu + m_He)
# grids
r = np.linspace(0.5e-10, 50e-10, 20000)
E_min, E_max = -5000*cm1_to_J, 0.0

# find bound states for l=0
roots_gs = find_bound_states(r, mu, V_gs, l=0, E_min=E_min, E_max=E_max, nroots=1)
roots_ms = find_bound_states(r, mu, V_ms, l=0, E_min=E_min, E_max=E_max, nroots=1)
print("Bound energies (cm-1) GS:", np.array(roots_gs)/cm1_to_J)
print("Bound energies (cm-1) MS:", np.array(roots_ms)/cm1_to_J)

# scattering parameters
E_grid = np.array([100, 500, 1000])  # in cm-1
E_SI   = E_grid * cm1_to_J
k_grid = np.sqrt(2*mu*E_SI)/hbar
l_max  = 10

# plot setup
colors = mpl.colormaps['coolwarm'](np.linspace(0,1,len(E_grid)))
fig, axes = plt.subplots(1,2, figsize=(16,6))

for j, V_func in enumerate([V_gs, V_ms]):
    sigma_tots = []
    for i, E_cm in enumerate(E_grid):
        k = k_grid[i]; E_SI_val = E_SI[i]
        # compute phase shifts
        deltas = [compute_phase_shift(r,
                  numerov_scatter(r[0]**(l+1), r[1]**(l+1), r, E_SI_val, mu, V_func, l),
                  k, l)
                  for l in range(l_max+1)]
        # differential
        thetas = np.linspace(0, np.pi, 180)
        dsdΩ = differential_cross_section(thetas, k, deltas)
        # total
        σ_tot = total_cross_section(k, deltas)
        sigma_tots.append(σ_tot)
        axes[0].plot(np.rad2deg(thetas), dsdΩ,
                     label=f"{['GS','MS'][j]} E={E_cm} cm" + "$^{-1}$",
                     c=colors[i], ls=['-','--'][j])
    # total vs energy
    axes[1].plot(E_grid, sigma_tots, label=['GS','MS'][j],
                 ls=['-','--'][j])

for ax in axes:
    ax.set_yscale('log')
axes[0].set_xlabel(r'$\theta$ (degrees)')
axes[0].set_ylabel(r'$d\sigma/d\Omega$')
axes[1].set_xlabel(r'$E$ (cm$^{-1}$)')
axes[1].set_ylabel(r'$\sigma_{\rm tot}$')
for ax in axes:
    ax.legend(ncols=2, loc='upper right')
plt.tight_layout()
plt.savefig("cross_section_correct_units_QM.png")
plt.show()
