"""
Makes conservation plot for the transport-only test
"""

import matplotlib.pyplot as plt
from tomplot import set_tomplot_style
import pandas as pd
import numpy as np
from matplotlib import ticker

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
base_results_dir = '/data/users/abrown/data/pdc_idealised_paper/transport'
plot_dir = '.'
plot_name = f'{plot_dir}/conservation.png'
test_cases = ['stretchy_conv_adv', 'stretchy_conv']
labels = ['Advective', 'Conservative']
linestyles = ['--', '-']
xlabel = r'$t \ / $ s'
ylabel = r'$(M(t) - M(0)) / M(0)$'
column_indices = {'timestep': 6, 'mass': 7}
dt = 4.0

# Create figure here
set_tomplot_style()
fig, ax = plt.subplots(1, 1, figsize=(6,5))

for test_case, label_stem, linestyle in zip(test_cases, labels, linestyles):

    file_name = f'{base_results_dir}/{test_case}/raw_data/output.log'

    # ------------------------------------------------------------------------ #
    # Read in raw data as a "CSV"
    # ------------------------------------------------------------------------ #
    # key point is skipinitialspace=True, which means that not every space is a delimiter (only bunches of spaces)
    raw_data = pd.read_csv(file_name, header=None, sep=' ', skipinitialspace=True, usecols=column_indices.values())

    # Rename columns to human-readable values
    # First make inverse dictionary of columns
    column_titles = {}
    for key, value in column_indices.items():
        column_titles[value] = key
    raw_data = raw_data.rename(columns=column_titles)

    # ------------------------------------------------------------------------ #
    # Clean up data
    # ------------------------------------------------------------------------ #

    raw_data['mass'] = raw_data['mass'].astype(float)
    raw_mass = raw_data['mass'].values

    # Extract values, keep first but otherwise remove all odd values
    mass = np.zeros(int(len(raw_mass) / 2)+1)
    mass[0] = raw_mass[0]
    mass[1:] = raw_mass[1::2]
    mass = (mass - mass[0]) / mass[0]  # Normalise
    time = np.array([dt * i for i in range(len(mass))])

    # ------------------------------------------------------------------------ #
    # Plot data
    # ------------------------------------------------------------------------ #

    data_range = np.max(mass) - np.min(mass)
    label = f'{label_stem}, range = {data_range:1.1e}'

    ax.plot(time, mass, linestyle=linestyle, color='black', label=label)

ax.legend()

ax.set_xlabel(xlabel, labelpad=-10)
ax.set_ylabel(ylabel, labelpad=-30)
ax.set_xticks([0, 2000])
ax.set_yticks([-2e-4, 0, 2e-4])
ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1g}"))

print(f'Plotting to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()
