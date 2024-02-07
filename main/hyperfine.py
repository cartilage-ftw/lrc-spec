import numpy as np
from lmfit import Model


def calculate_hfs_energy(F, A_hfs, B_hfs, E0, I, J):
    K = F*(F+1) - I*(I+1) - J*(J+1)
    delE = A_hfs*K/2 + B_hfs*((3/2)*K*(K+1) - 2*I*(I+1)*J*(J+1))/(2*I*(2*I-1)*2*J*(2*J-1))
        #B_hfs*( ((3/2)*K*(K+1) - 2*I*(I+1)*J*(J+1)) / (2*I*(2*I-1)*2*J*(2*J-1)) )
    #print('Delta E', delE)
    return delE + E0


def calculate_hfs_diff(F, A_hfs, B_hfs, I, J):
    F1 = F
    F2 = F-1
    delE1 = calculate_hfs_energy(A_hfs, B_hfs, I, J, F1)
    delE2 = calculate_hfs_energy(A_hfs, B_hfs, I, J, F2)
    delta = delE1 - delE2
    print('Caculated energy diff:', delta/1000, 'GHz')
    return delta




#deltaEs = [0.66 * 30000, 0.68* 30000] # MHz
Fs = np.array([5/2, 7/2, 9/2])
Fs_176 = np.array([6,7,8])

Es_nov24 = np.array([28502.36, 28503.02, 28503.7])#*30_000
E_nov27_176 = np.array([28502.149, 28503.050, 28503.908])
#Es_nov27 = np.array([28502.348, 28503.008, 28503.697])
hfs_model = Model(calculate_hfs_energy)

E0_175 = np.sum(((2*Fs+1)*Es_nov24)/(np.sum(2*Fs+1)))
E0_176 = np.sum(((2*Fs_176+1)*E_nov27_176)/np.sum(2*Fs_176+1))
pars = hfs_model.make_params(A_hfs=3000, B_hfs=-1000, E0=E0_175, I=7/2, J=1)
pars['I'].set(vary=False)
pars['J'].set(vary=False)
#pars['E0'].set(vary=False)

pars_176 = hfs_model.make_params(A_hfs=3000, B_hfs=-1000, E0=E0_176, I=7, J=1)
pars_176['I'].set(vary=False)
pars_176['J'].set(vary=False)
#pars_176['E0'].set(vary=False)
print('Initialized with parameters', pars_176)



fit_176 = hfs_model.fit(E_nov27_176, pars_176, F=Fs_176)
fit = hfs_model.fit(Es_nov24, pars, F=Fs)#, method='emcee')
print(fit.fit_report())
print(fit_176.fit_report())

print('Calculated E0 for 175', E0_175)
print('Calcualted E0 for 176', E0_176)