import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import io_utils



plt.rcParams['axes.formatter.useoffset'] = False

def gauss(x, A, mu, sigma):
	y= A*np.exp(-(x-mu)**2/(2*(sigma**2)))
	return y


def gaussian_spectrum(E_list: iter, F_list: iter, E_range: np.array, freq_fwhm=4_600):
	"""

	E_list: list of the energy gaps between the hyperfine multiplets
	E_range: to evaluate the spectrum in
	freq_fwhm: Laser Bandwidth in terms of frequency (FWHM).

	A Gaussian line profile is assumed.
	For deciding line intensities, proportionality to (2F'+1) since the lower level is taken to be unsplitted
	"""

	intensities = 2*np.array(F_list)+1
	norm_int = 100*intensities/np.max(intensities)
	laser_sigma_cm = (freq_fwhm/2.355)/29_999 # bandwidth in cm-1
	#laser_sigma_cm)
	#print(norm_int)
	signal_y = np.zeros(len(E_range))
	#print(signal_y)
	for line_pos, intensity in zip(E_list, norm_int):
		print("LINE POS:", line_pos)
		signal_y += gauss(E_range, intensity, line_pos, laser_sigma_cm)
		#print(signal_y)
	return signal_y, intensities

def calc_hfs_energy(A_hfs, B_hfs, I, J, F_list):
	"""
	Calculates energy splitting of upper level
	"""
	# total spins F of different multiplets
	K = F_list*(F_list+1) - I*(I+1) - J*(J+1)
	delE = A_hfs*K/2 + B_hfs*((3/2)*K*(K+1) - 2*I*(I+1)*J*(J+1))/(2*I*(2*I-1)*2*J*(2*J-1))
	return delE

J = 1
I = 3/2
J_list = []
if J > 0:
	for j in range(-J,J+1):
		J_list.append(j)
	J_list = np.sort(np.array(J_list))
F_list = I + J_list

E0_1P1 = 22180.533
E0_3P1 = 29250.366

# HFS splitted energies of the 
ac227_1Po1_E = E0_1P1 + calc_hfs_energy(A_hfs=-1689, B_hfs=352, I=3/2, J=1, F_list=F_list)/29_999
ac227_3Po1_E = E0_3P1 + calc_hfs_energy(A_hfs=8712, B_hfs=-874, I=3/2, J=1, F_list=F_list)/29_999
lu175_3Po1_E = 28503.138 + calc_hfs_energy(A_hfs=4952, B_hfs=-1962, I=7/2, J=1,
											 F_list=np.array([5/2, 7/2, 9/2]))/29_999

print(ac227_1Po1_E)
E_grid_1P1 = np.linspace(E0_1P1-0.5, E0_1P1+0.5, 1000)
E_grid_3P1 = np.linspace(E0_3P1-2, E0_3P1+2, 1000)
lu175_grid = np.linspace(28501, 28505, 300)

ac227_1Po1_spec, intensities_1Po1 = gaussian_spectrum(ac227_1Po1_E, F_list, E_grid_1P1, freq_fwhm=2_000)
ac227_3Po1_spec, intensities_3Po1 = gaussian_spectrum(ac227_3Po1_E, F_list, E_grid_3P1, freq_fwhm=4_500)

lu175_3Po1_spec, intensities_lu175 = gaussian_spectrum(lu175_3Po1_E, np.array([5/2, 7/2, 9/2]), lu175_grid)
fig, axes = plt.subplots(3, 1, figsize=(5,9))

axes[0].plot(E_grid_1P1, ac227_1Po1_spec, c='dimgray', lw=1, label=r'$^{227}$Ac$^+$ $^1\textrm{P}^o_1$')
axes[0].text(x=22180.73, y=85, s=r'$\Delta f=2.0$ GHz')
axes[1].plot(E_grid_3P1, ac227_3Po1_spec, c='dimgray', lw=1, label=r'$^{227}$Ac$^+$ $^3\textrm{P}^o_1$')
axes[1].text(x=29251.2, y=85, s=r'$\Delta f=4.5$ GHz')
axes[2].plot(lu175_grid, lu175_3Po1_spec, c='dimgray', lw=1, label=r'$^{175}$Lu$^+$ $^3\textrm{P}^o_1$')


def lorentzian(x_range, mean, height, y_max, fwhm=10):
    # 10 MHz FWHM assumed
    gamma = fwhm/29_972 # MHz into cm-1
    print("Rescaling with height", height)
    y_values = (1/np.pi)*( (gamma/2)/ ( (x_range-mean)**2 + (gamma/2)**2))
    # normalize the distribution
    norm_values = y_values/np.max(y_values)
    # rescale to have the highest peak being 100%
    rescaled_values = (100/y_max)*norm_values
    return height * rescaled_values

lu175_meas = pd.read_csv('./spectra_and_atds/2023-11-27-20-04-53_spectrum.csv')
print(lu175_meas.head())
# also display measured spectrum for Lu-175
axes[2].errorbar(lu175_meas['Wavenumber'], lu175_meas['MS Fraction']*2.5, yerr=lu175_meas['MS Error']*2.5, capsize=2,
				  #c='#aaaaff',
      			c='mediumpurple', marker='o', ms=4, ls='', label='Measured')

for ax, lines, intensities in zip(axes, [ac227_1Po1_E, ac227_3Po1_E, lu175_3Po1_E],
                                        [intensities_1Po1, intensities_3Po1, intensities_lu175]):
	for line_pos, intensity in zip(lines, intensities):
		# +/- 30 MHz around the line
		x_range = np.linspace(line_pos-0.15, line_pos+0.15, 100)
		ax.plot(x_range, lorentzian(x_range, mean=line_pos, height=intensity, y_max=np.max(intensities)),
          				ls='--', lw=1, zorder=1, c='deeppink')
		#ax.axvline(x=line, ymin=0, ymax=1, c='deeppink', ls='--', lw=1, zorder=-1)
	ax.legend(loc='upper left')
	

plt.tight_layout()
plt.savefig('actinium_mock_spectra.png')
plt.show()