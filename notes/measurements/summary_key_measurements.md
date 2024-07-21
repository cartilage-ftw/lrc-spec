# Summary of Key Measurements [Nov 2023-Sep 2024]

the digital notebook of Aayush's measurements (this is a selection; see physical notebook for complete list).

**TODO**: For some measurements, I haven't listed QMS settings $m$ and $\delta m$ that was chosen, etc. Refer to the stored file configurations and write those down as well. They are important to know to make sure there isn't contamination from other species.

### Nov 16, 2023: Effect of number of laser pulses per bunch

p(stop cell) = $64.01$ mbar
p(drift tube) = $2.90$ mbar
p(buncher) = $4.66 \times 10^{-2}$ mbar
Laser operated at 8 kHz (and 56A YAG current)
No optical density filters were used. Unfortunately, I also don't see a measurement of the per pulse energy (but the reflection from the beamsplitter may have been 330 nJ, from memory)

**Note:** These were not complete hfs scans, but measurements taken at $\sim 28503.3$ cm$^{-1}$  

| Buncher Freq | # of shots | File Name           |
| ------------ | ---------- | ------------------- |
| 1 kHz        | 8          | 2023-11-16-17-14-00 |
| 500 Hz       | 16         | 2023-11-16-17-16-23 |
| 250 Hz       | 32         | 2023-11-16-17-17-54 |
| 200 Hz       | 40         | 2023-11-16-17-19-17 |
| 100 Hz       | 80         | 2023-11-16-17-20-27 |
Near complete population inverse was achieved in the 100 Hz case.

**Complete HFS scans** were also done; listed below
Due to no OD filters being used, the lines were extremely power broadened (still not flat top)

| Buncher Freq | # of shots | File Name           |
| ------------ | ---------- | ------------------- |
| 100 Hz       | 80         | 2023-11-16-17-41-06 |
| 200 Hz       | 40         | 2023-11-16-18-03-22 |
| 500 Hz       | 16         | 2023-11-16-18-17-39 |
| 1 kHz        | 8          | 2023-11-16-18-45-27 |
|              |            | 2023-11-16-18-55-57 |
**Note**: There are 2 separate scans for the 1 kHz case. Each with 2000 counts per wavenum scan step.

### Nov 20, 2023: Saturation curve using OD filters
p(DT) = 3.27 mbar
p(buncher) = $5.09 \times 10^-2$ mbar
p(stop. cell) = $63.82$ mbar

Laser operation: 8 kHz (56A current for pump laser)
Buncher: 1 kHz

**Energy measurement**: 20.11.2023_56A_8kHz_no_od.csv

In the following table, I scratched the first measurement, and the second should be preferred. I adjusted the position of the last mirror before the glass window (by sliding the linear stage, to position 4.0 mm)

| OD Filter | File Name               |
| --------- | ----------------------- |
| 0         | ~~2023-11-20-16-59-07~~ |
| 0         | 2023-11-20-17-06-40     |
| 2.0       | 2023-11-20-17-13-39     |
| 1.0       | 2023-11-20-17-17-02     |
| 0.6       | 2023-11-20-17-28-15     |
| 0.3       | 2023-11-20-17-30-37     |
| 0.1       | 2023-11-20-17-33-03     |

Laser frequency changed from 8 kHz to 1 kHz.
**Energy measurement at 1kHz**: 20.11.2023_56A_1kHz.csv

| OD Filter | File Name           |
| --------- | ------------------- |
| 0         | 2023-11-20-17-42-46 |
| 1.0       | 2023-11-20-17-45-02 |
| 0.6       | 2023-11-20-17-47-10 |
| 0.3       | ... - 17-49-23      |
| 0.1       | ... 17-51-32        |
| 0.4       | ... 18-14-10        |

Meausrements done with better statistics (with 2000 cts/step)
 Energy re-measured to check stability: **20.11.2023_56A_1kHz_2.csv**

| OD Filter | File Name           |
| --------- | ------------------- |
| 0.4       | 2023-11-20-18-16-11 |
| 1.0 + 0.4 | 2023-11-20-18-18-48 |
| 0.6 + 0.4 | 2023-11-20-18-21-34 |
| 0.3 + 0.4 | ... -18-23-47       |
| 0.1 + 0.4 | ... -18-26-04       |
|           |                     |

### Nov 22, 2023: 
* We learnt that the time delay between the DAQ trigger sent to the YAG laser, and the arrival of the pulse at the window (as measured via a PIN diode) is about $500$ $\mathrm{ns}$ 
* Increasing the trap depth from 1V to 2V makes the spatial confinement of the ions narrower.

### Nov 24, 2023: Several (3 or so) scans at same configuration

drift  tube gas inlet rate $5.6$ mbar L/s
$p(DT) = 3.13$ mbar
p(buncher) = $5.2\times 10^{-2}$ mbar
p(stop. cell) = 63.54 mbar

* I started the measurements with $\delta m$ initially at $3.9$ (assuming I needed a high mass resolution). To compensate for the decreased count rate, I decided to increase ablation laser intensity too much. This resulted in substantial metastable population coming even without the involvement of the resonant probe laser. 
  **Therefore, the background level was quite high**!
  **Note**: I believe we had the Yb sample in the sample holder as well. And it may well be that what I thought was contamination from Lu II metastable, was just the ground state of Yb II (which we know very well that it is the case; e.g. [Kim et al. 2024](https://arxiv.org/abs/2404.05383))
This persisted even when I lowered the intensity again, and reduced $\delta m \to 3.5$ 

$f_{buncher} = 1$ kHz
* Background scan (probe laser blocked)
	*`2023-11-24-16-16-03`
	`2023-11-24-16-19-46`
	`2023-11-24-16-24-10`
	`2023-11-24-16-30-53`
with Laser Operation at 8 kHz, 56A, $f_{buncher} = 1$ kHz, buncher delay=0

**Note**: Yes, it's true that an individual scan was repeated many times. Each file contains a complete hyperfine scan of the

| Isotope        | $f_{buncher}$ [Hz] | Optical Density Filters | File Name                                                               |
| -------------- | ------------------ | ----------------------- | ----------------------------------------------------------------------- |
| $^{175}$Lu     | $1,000$            | 0.3 + 3.0               | `2023-11-24-16-36-48`<br> `2023-11-24-16-42-30`<br> `...-16-53-52`      |
| :              | $1,000$            | 0.6 + 3.0               | `2023-11-24-17-04-18`<br>`2023-11-24-17-07-43`<br>`2023-11-24-17-10-52` |
| :              | $1,000$            | 1.0 + 3.0               | `2023-11-24-17-40-23`<br>`2023-11-24-17-44-26`<br>`2023-11-24-17-48-38` |
| :              | $100$              | 1.0 + 3.0               | `2023-11-24-17-54-04`<br>`2023-11-24-17-59-27`<br>`2023-11-24-18-05-48` |
| $^{175}$Lu     | $100$              | 2.0 + 3.0               | `2023-11-24-18-12-27`<br>`2023-11-24-18-18-12`<br>`2023-11-24-18-23-45` |
| **$^{176}$Lu** | $100$              | 2.0 + 3.0               | `2023-11-24-18-36-50`<br>`2023-11-24-18-53-00`                          |





