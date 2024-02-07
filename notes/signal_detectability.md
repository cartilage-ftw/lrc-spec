All scans were homogeneously processed by making an MS cut at $0.194$ ms (decided by visual inspection). The operation conditions were kept constant for all of them
* p(DT) = 0.25 mbar, p(buncher) = , 
* Laser operation at 8 kHz, 56A YAG current. See *** file for pulse energy density measurements.
* 200 Hz bunching frequency

Data files that were used:
* 2000 counts/step: `2023-12-22-17-21-27.csv`
* 500 cts/step: `2023-12-22-18-49-46`
* 100 cts/step: `2023-12-22-18-56-11`
* 10 cts/step: `2023-12-22-19-00-58`
	* In the lab notebook, you will find multiple scans at this setting. I used this specific one because for it, I increased $\Delta M$ sufficiently high to make the count rate per second very low. This is important, as the dye laser motor needs adequate time for switching wavelength between successive steps, and steps of 10 counts might otherwise be very fast.
	* A "background" signal scan with probe laser blocked was taken: `2023-12-22-19-09-16`

The figure below shows spectra extracted from all four scans, and the arrival time distribution for the counts that came at certain wavelength steps. 
The resonant one corresponds to the peak of the $F=7/2\to 7/2$ hyperfine multiplet, while off-resonant corresponds to the first or second step in the scan (see spectra in the left panels for a visual idea).

![](../figures/effect_decreasing_count_per_step.png)

This is extremely helpful in discovering signal regions when beam time is limited.

In absence of prior knowledge of mobilities in the different electronic states, one must rely on signal processing methods to discover faint signals.