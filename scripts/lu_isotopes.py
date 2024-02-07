import pandas as pd
import matplotlib.pyplot as plt

lu175 = pd.read_csv('../data/lu/175Lu.csv', sep=',', names=['x', 'y'])
lu176 = pd.read_csv('../data/lu/176Lu.csv', sep=',', names=['x', 'y'])

lu175['wavenumber'] = 2*(lu175['x'] + 14250)
lu176['wavenumber'] = 2*(lu176['x'] + 14250)

lu_data = pd.read_csv('../data/lu/dat', sep='\t', names=['wavenumber', 'x_err', 'y', 'y_err'])
fig, ax = plt.subplots(figsize=(6,6))

ax.plot(lu175['wavenumber'], lu175['y'], label='$^{175}$Lu$^+$')
ax.plot(lu176['wavenumber'], lu176['y'], label='$^{176}$Lu$^+$')
ax.plot(lu_data['wavenumber'], lu_data['y']*100, label='old $^{175}$Lu$^+$')
plt.legend()
plt.show()