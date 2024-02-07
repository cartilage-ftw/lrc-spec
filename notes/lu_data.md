## New Fit to November Measurements

### Saturated Linewidths

For the $F=5/2$ multiplet

| $\gamma_s$| Per Pulse Energy (reflection) | Conversion Factor | Density Filter |  Transmission $\%$ | Final Assumed Energy Density |
| --- | --- | --- | --- | --- | --- |
| | | | | | |


For the $F=7/2$ multiplet

For the $F=9/2$ multiplet


Summary of derived quantities

**Note**: this currently only includes uncertainty due to the fit to the line, not systematic uncertainties in wavelength measurement ($\sim 0.002$ cm-1) and calibration offset (~4th decimal place) and the near $0.001$ cm-1 dispersion in $\bar{\nu}$ for each scan step, while recording the hfs spectrum.

* for the $^1S_0 \to ^3P^o_1$ transition in Lu
    * center of gravity for the two isotopes
        * $^{175}\rm{Lu}$: $ 28503.138 \pm 0.003$ cm-1
        *   $^{176}\rm{Lu}$: $28503.113 \pm 0.003$ cm-1
    *    $\Delta \nu^{175,176} = \nu^{176} - \nu^{175} = 0.025 \pm 0.003$ cm-1

#### TODO:
* Check for hyperfine anomalies
* Take $\Delta \nu^{A,A'} = K_{MS}\left(\frac{1}{M_A} - \frac{1}{M_{A'}}\right) + F\delta \langle r^2\rangle^{AA'}$ and estimate the difference in nuclear charge radii

| Species | $A_{hfs}$ [$10^{-3}\ \rm{cm}^{-1}$] | $B_{hfs}$ [$10^{-3}\ \rm{cm}^{-1}$] | Study |
| --- | --- | --- | --- |
| $^{175}\rm{Lu}$ | $165.16 $ | $-65.56$ | this work |
| $^{175}\rm{Lu}$ | $165.5 \pm 0.3$ | $-62.4\pm 1.0$ | [Den Hartog _et al._ (2020), _ApJS_](https://ui.adsabs.harvard.edu/abs/2020ApJS..248...10D/abstract) |
| |
| $^{176}\rm{Lu}$ | $116.6$ | $-86.8$ | this work | 
| $^{176}\rm{Lu}$ |


Quadrupole moments

| Species | $Q_s$ [b] | Study | method |
| --- | --- | --- | --- |
| $^{175}\rm{Lu}$ |  |  this work |
| $^{175}\rm{Lu}$ | $+3.49 \pm 0.02$ | [Dey _et al._ (1979)]() |
| $^{175}\rm{Lu}$ | $+3.415 \pm 0.034$ | [Haiduke _et al._ (2007), _Chem. Phys. Lett_](https://www.sciencedirect.com/science/article/pii/S000926140700961X) | molecular method |
| $^{175}\rm{Lu}$ | | |
| |
| $^{176}\rm{Lu}$ | | this work |
|  $^{176}\rm{Lu}$ | $+4.818 \pm 0.0048 $ | [Haiduke _et al._ (2007), _Chem. Phys. Lett_](https://www.sciencedirect.com/science/article/pii/S000926140700961X)|  |
| $^{176}\rm{Lu}$ | $+4.92 \pm 0.03$ | | ratio of $B_{hfs}$ factors|

### Nov 27, 2023


$^{176}$ Lu ![Voigt fits to the Nov 27 data for Lu](../figures/lu176_fit_nov27.png)
|F | Wavenumber [cm-1] | $\Gamma$ [GHz] | Height |
|---|---|---|---|
|$6$|$28502.149 \pm 0.002$ | $2.12 \pm 0.09$  | - |
|$7$|$28503.050 \pm 0.003$| $2.36 \pm 0.10$ | - |
|$8$| $28503.908 \pm 0.003$ | $2.19 \pm 0.09$ | - |

--- 

$^{175}$ Lu
![](../figures/fitted_hfs_multiplet.png)


|F | Wavenumber [cm-1] | $\Gamma$ [GHz] | Height |
|---|---|---|---|
|$5/2$|28502.348 | 1.65 | 5.76 |
|$7/2$|28503.008 | 1.58 | 7.37 |
|$9/2$| 28503.697 | 2.41 | 10.95 |


### Nov 24, 2023 

|F | Wavenumber [cm-1] | $\Gamma$ [GHz] | Height |
|---|---|---|---|
|$5/2$|28502.36 | 2.84 | 6.62 |
|$7/2$|28503.02 | 2.58 | 7.89 |
|$9/2$| 28503.70 | 3.38 | 13.17 |
## Old Fits to data by Kim

Assuming three Voigt profiles, and a linear offset (~0.08)

| | Wavenumber [cm-1] | Height | FWHM [GHz]| FWHM [cm-1]  |Gaussian $\sigma$ [cm-1]| Lorentzian $\Gamma$ [cm-1] |
| --- | --- | --- | --- | --- | --- | --- |
|1. |28502.38 | 0.65 | 6.1|0.203 | 0.056 | 0.056 |
|2. |28503.03 | 0.78 | 6.9 |0.231 | 0.064 | 0.064 | 
|3. |28503.72| 0.88 | 7.9 |0.264 | 0.073 | 0.073

Note that the Gaussian $\sigma$ and Lorentzian $\Gamma$ were constrained to be equal by the fitting routine. Also, if you translate the FWHM to a purely Gaussian $\sigma =$ FWHM $/ 2.355$ (worth comparing with laser profile of 1.8 GHz at $470$ $nm$ to notice additional broadening)
| $\sigma$ [GHz] |
| --- |
| 2.59 |
| 2.94 |
| 3.36 |

For a transition with lifetime $\tau=10^{-8}$ s, and laser with $\lambda=500$ nm, the saturation intensity should be about $I_{sat} =2×10^2$ $W/m^2$ and natural linewidth some 16 MHz.

If a laser is operating at $1$ $mW$, for a spot size of $100\mu m$, that corresponds to an intensity of $\sim 10^5 W/m^2$, which is $~1000$ times $I_{sat}$.

 if saturated linewidth goes as $\gamma_s \sim \gamma\sqrt{I/I_{sat}}$, then the 16 MHz line becomes 500 MHz ($\Delta\nu_{\rm{FWHM}}$).