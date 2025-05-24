import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
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


def load_potential(csv_file):
    if '.dat' in csv_file:
        data = pd.read_csv(csv_file,sep=',', skiprows=5,
                                names=['empty', 'r', 'V_Omega0', 'V_Omega+1', 'V_Omega-1'])
        #r_grid = data['r'].to_numpy()
        r_grid = data['r'].to_numpy()
        V_grid = data['V_Omega0'].to_numpy()
        #print(r_grid)
        #exit()
    else:
        data = np.loadtxt(csv_file, delimiter=',', skiprows=1)
        r_grid = data[:, 0]
        V_grid = data[:, 1]
        #print(r_grid)
        #exit()
    return interp1d(r_grid, V_grid, kind='cubic', fill_value='extrapolate')


def numerov(u0, u1, r, k, mu, V_func, l):
    dr = r[1] - r[0]
    u = np.zeros_like(r)
    u[0], u[1] = u0, u1
    Q = k**2 - 2*mu*V_func(r) - l*(l+1)/r**2

    for i in range(1, len(r)-1):
        k_im1, k_i, k_ip1 = Q[i-1], Q[i], Q[i+1]
        u[i+1] = (2*u[i]*(1 - dr**2 * k_i/12) - u[i-1]*(1 + dr**2 * k_im1/12)) \
                 / (1 + dr**2 * k_ip1/12)
    return u


def compute_phase_shift(r, u, k, l):
    dr = r[1] - r[0]
    # match near the end
    idx = -2
    rm, um = r[idx], u[idx]
    up, um1 = u[idx+1], u[idx-1]
    uprime = (up - um1) / (2*dr)

    jl  = spherical_jn(l, k*rm)
    nl  = spherical_yn(l, k*rm)
    jl_p = spherical_jn(l, k*rm, derivative=True)
    nl_p = spherical_yn(l, k*rm, derivative=True)
    num = um * jl_p - uprime * jl
    den = uprime * nl - um * nl_p
    return np.arctan2(num, den)


def differential_cross_section(thetas, k, deltas):
    cos_t = np.cos(thetas)
    f = np.zeros_like(thetas, dtype=complex)
    for l, delta in enumerate(deltas):
        P_l = lpmv(0, l, cos_t)
        f += (2*l+1) * np.exp(1j*delta) * np.sin(delta) * P_l
    return np.abs(f)**2 / k**2

def total_cross_section(k, deltas):
    l = np.arange(len(deltas))
    return 4*np.pi/k**2 * np.sum((2*l+1) * np.sin(deltas)**2)


mu = 0.5           # reduced mass 
E_grid  = np.array([0.1, 1.0, 100, 500])           # energy
k_grid  = np.sqrt(2*mu*E_grid)
l_max = 10         # truncate partial wave sum
r_min, r_max = 1e-5, 50.0
N = 50000
r = np.linspace(r_min, r_max, N)

# Load potential from 'potential.csv'
V_func = load_potential('../structure/Archiv2/Lu+inGS.csv')
V_func_ms = load_potential("../structure/Archiv2/SO3D1_.dat")
def compute_phase_shifts(V_func, k):
    # Compute phase shifts
    deltas = []
    for l in range(l_max+1):
        # initial behavior u ~ r^(l+1)
        u0 = r_min**(l+1)
        u1 = r[1]**(l+1)
        u = numerov(u0, u1, r, k, mu, V_func, l)
        deltas.append(compute_phase_shift(r, u, k, l))
    return deltas


colors = mpl.colormaps['coolwarm'](np.linspace(0,1,len(E_grid)))

fig, axes = plt.subplots(1,2, figsize=(16,6))
for j, V_ in enumerate([V_func, V_func_ms ]):
    sigma_tots = []
    for i, k in enumerate(k_grid):
        deltas = compute_phase_shifts(V_, k)
        #print(deltas)
        #exit()
        # Differential cross section
        thetas = np.linspace(0, np.pi, 180)
        #print(thetas)
        #exit()
        dsdOmega = differential_cross_section(thetas, k, deltas)
        #print(dsdOmega)
        #exit()
        # Total cross section
        sigma_tot = total_cross_section(k, deltas)
        sigma_tots.append(sigma_tot)
        E = k**2/(2*mu)
        axes[0].plot(np.rad2deg(thetas), np.abs(dsdOmega), label=['GS', 'MS'][j] + f' Energy: {E:0.1f}', c=colors[i],
                ls=['-', '--'][j])
    axes[1].plot(E_grid, sigma_tots, label=['GS', 'MS'][j], ls=['-', '--'][j])
        #break
# Display results
print(f"Total elastic cross section σ_el = {sigma_tot:.6f} (units)")

#plt.plot(thetas * 180/np.pi, dsdOmega)
for ax in axes:
    ax.legend(ncols=2, loc='upper right')
    ax.set_yscale('log')
axes[0].set_xlabel(r'Scattering angle $\theta$ (degrees)')
#plt.yscale("log")
axes[0].set_ylabel(r'Differential Cross Section')
axes[1].set_xlabel("Energy")
axes[1].set_ylabel("$\sigma_{tot}$")
axes[0].set_title(r'$d\sigma/d\Omega$ - Differential Elastic Cross Section')
plt.tight_layout()
plt.savefig("cross_section_QM.png", dpi=150)
plt.show()
