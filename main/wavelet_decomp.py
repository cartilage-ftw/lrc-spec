import numpy as np
import matplotlib.pyplot as plt
import scipy

import pywt
#from wavelet import WaveletTransform, getExponent


'''t = np.linspace(-1, 1, 200, endpoint=False)
sig  = np.cos(2 * np.pi * 7 * t) + np.real(np.exp(-7*(t-0.4)**2)*np.exp(1j*2*np.pi*2*(t-0.4)))
widths = np.arange(1, 31)
cwtmatr, freqs = pywt.cwt(sig, widths, 'mexh')
plt.imshow(cwtmatr, extent=[-1, 1, 1, 31], cmap='PRGn', aspect='auto',
			vmax=abs(cwtmatr).max(), vmin=-abs(cwtmatr).max())  # doctest: +SKIP
plt.show()'''


signal1 = scipy.stats.norm(loc=5.5, scale=4.5)
amp1 = 200.
signal2 = scipy.stats.norm(loc=6.8, scale=0.5)
amp2 = 25

x = np.linspace(0, 11, 1000)
#x = np.arange(8)
noise = np.random.uniform(size=len(x))

component1 = amp1*signal1.pdf(x)
component2 = amp2*signal2.pdf(x)
y = (component1 + component2) + noise


print('End haha ******')

#inv_app_coeffs = pywt.upcoef('a', coeffs[0], wavelet='sym2', level = len(coeffs)-1, take=len(y))
#inv_det_coeffs = [pywt.upcoef('d', coeffs[i], wavelet='sym2', take=len(y), level=len(coeffs)-i)
#                  			for i in range(1, len(coeffs))]

mra_coeffs = pywt.mra(y, wavelet='coif2', transform='dwt')

#inv_all = pywt.waverec(coeffs, 'db5')

fig, axes = plt.subplots(len(mra_coeffs)+1, 2, figsize=(8,12), sharex=True)

axes[0, 0].step(x, y, label='$f(x)$')
axes[0, 1].step(x, component2, label='component 2', c='cornflowerblue')
axes[0, 1].step(x, component1, label='component 1', c='pink')

#axes[1, 0].step(x, mra_coeffs[0])

coeff_sum = np.zeros(shape=(len(mra_coeffs[0]),))

#axes[1, 1].step(x, coeff_sum)

def get_label(i):
	if i > 1:
		return '$c_{D_' + f'{i+1}' + '}$'
	else: 
		return '$c_A$'
	

for i, det_coeffs in enumerate(mra_coeffs):
	coeff_sum += mra_coeffs[i]
	axes[i+1, 0].step(x, mra_coeffs[i], label=get_label(i))
	axes[i+1, 1].step(x, coeff_sum, label=get_label(i))

for i, j in np.ndindex(axes.shape):
	axes[i,j].legend(loc='upper left')


plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.savefig('wavelet_decomp.png', dpi=300)
plt.show()