import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import constants
import saturation_curve

"""
Some "utility" methods, as I decided to call them, for performing certain calculations
"""


def lifetime_to_isat(lifetime, wavelength):
    """
    Estimates I_sat for a level with lifetime \tau
    """
    return np.pi*constants.h*constants.c/(3*lifetime*(wavelength**3))

def inten_to_power(inten, spot_area):
    return inten * spot_area


def pulseenergy_for_power(power, num_pulses, pulse_dur):
    return pulse_dur*power/num_pulses


def sat_linewidth(A_ki, energy, sat_energy):
    # the saturation "pulse energy" is used for calculation here
    return A_ki*np.sqrt(1 + energy/sat_energy)


# per pulse energy at 8 kHz, 56 A
#pp_energy = 2.4E-6 # 2.4 mu J

transm_filters = [[2.0, 4.0], [1.0, 4.0], [0.6, 4.0], [0.3, 4.0],
                  [0, 4.0], [0.6, 3.0], [0.3, 3.0], [0, 3.0], [2.0, 0],
                  [1.0, 0.5], [1.0, 0], [0.6, 0], [0.3, 0], [0, 0]]


def make_theor_saturation(A_ki, pulse_energy, sat_energy, filters):
    """
    Note that once again, instead of intensity/power this directly deals with
     per pulse energy
    """
    energies = []
    linewidths = []
    for filter_comb in filters:
        percent = saturation_curve.get_transm(*filter_comb)
        energies.append(erg)
        linewidths.append(sat_linewidth(A_ki, erg, sat_energy))
    fig, ax = plt.subplots(figsize=(6,6))
    ax.plot(100*np.array(energies)/pulse_energy, np.array(linewidths)/1E9, marker='o', ms=7, c='gray', mfc='cornflowerblue', mec='cornflowerblue',
                ls='-', label=r'$A_{ki}=1.25\times 10^{7}$ s$^{-1}$')
    ax.set_xlabel('Transmission [\%] ')
    ax.set_ylabel('$\gamma_s (E)$ [GHz]')
    #ax.set_title(r'Predicted linewidth $\gamma_s$ for the Lu $^1S_0 \to ^3P_1$ transition probed with a laser operated at 8 kHz, 56A ')
    plt.tight_layout()
    plt.legend()
    plt.show()


if __name__ == '__main__':
    A_ki = 1.25E7
    wavelength = 1E-2 * (1/28502)
    print('Saturation intensity for state with lifetime', 1/A_ki)
    int = lifetime_to_isat(1/A_ki, wavelength)
    print('is ', int, 'Watt per m^2 per s')
    spot_area = 2.0E-6 # 2 mm^2
    power = inten_to_power(int, spot_area)
    print("For a spot area of 2 squared mm, this corresponds to a power of", power*1E6, 'microWatts')
    print('And for 80 pulses of 8 ns each')
    sat_pulse_energy = pulseenergy_for_power(power, 80, 8E-9)
    print(sat_pulse_energy*1E9, 'nJ is needed')
    print('When providing 1nJ, we can expect a linewidth of')
    print(sat_linewidth(A_ki, 0.72E-9, sat_pulse_energy)/1E9, 'GHz')

    make_theor_saturation(A_ki, pulse_energy=2.4E-6, sat_energy=sat_pulse_energy,
                           filters=transm_filters)