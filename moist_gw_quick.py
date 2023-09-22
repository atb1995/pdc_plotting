"""
A tomplot example, for making a quick single plot from LFRic data.
"""

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, tomplot_contours, tomplot_cmap,
                     plot_contoured_field, add_colorbar_ax,
                     tomplot_field_title, extract_lfric_vertical_slice)

def add_theta_e(dict_in, coords_Z):
    """
    Takes name of file and returns cube for theta_e variable.
    """
    Lv = 2.501e6
    cp = 1005.0
    exner = dict_in['exner']
    theta = dict_in['theta']
    mr_v = dict_in['m_v']
    exner_wth = np.zeros(np.shape(theta))
    # Get height data
    nz = np.shape(exner)[1]
    nz_full = np.shape(theta)[1]

    print(np.shape(exner))
    zmin = 0.0
    zmax = 10000.

    z1d_full = np.linspace(zmin, zmax, nz_full)
    z1d_half = 0.5*(z1d_full[1:] + z1d_full[0:nz_full-1])
    print(np.shape(z1d_half))

    for j in range(1,nz_full-1):   # Loop through columns
        # Map from W3 to Wtheta
        weight_denom = z1d_half[j] - z1d_half[j-1]
        weight_upper = z1d_full[j] - z1d_half[j-1]
        weight_lower = z1d_half[j] - z1d_full[j]
        # Initial data has different structure to time evolving data
        exner_wth[:, j] = weight_upper / weight_denom *exner[:, j] + weight_lower / weight_denom *exner[:, j]

    # Bottom
    weight_denom = z1d_half[1] - z1d_half[0]
    weight_upper = z1d_full[0] - z1d_half[0]
    weight_lower = z1d_half[1] - z1d_full[0]
    exner_wth[:, 0] = weight_upper / weight_denom *exner[:, 1] + weight_lower / weight_denom *exner[:, 0]
    # Top
    weight_denom = z1d_half[nz_full-2] - z1d_half[nz_full-3]
    weight_upper = z1d_full[nz_full-1] - z1d_half[nz_full-3]
    weight_lower = z1d_half[nz_full-2] - z1d_full[nz_full-1]
    exner_wth[:, nz_full-1] = weight_upper / weight_denom *exner[:, nz_full-2] + weight_lower / weight_denom *exner[:, nz_full-3]

    # Theta_e
    Temp = theta * exner_wth
    exp_arg = Lv * mr_v / (cp * Temp)
    theta_e = theta * np.exp(exp_arg)
    dict_in['theta_e'] = theta_e


# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
# When copying this example these should not be relative to this file
results_dir = '/hpc/scratch/d02/abrown/cylc-run/r40541_pdc_main_branch-gungho_model-meto-xc40-pdc_idealised/share/output/intel_64-bit_fast-debug/gungho_model/moist_gravity_wave/BiP200x4-1500x500_dt-3p6/results'
plot_dir = f'{abspath(dirname(__file__))}'
results_file_name = f'{results_dir}/lfric_diag.nc'
height_file_name = f'{results_dir}/lfric_diag.nc'
background_file_name = f'{results_dir}/lfric_initial.nc'
plot_name = f'{plot_dir}/moist_gw.png'
# ---------------------------------------------------------------------------- #
# Things that should be altered based on the plot
# ---------------------------------------------------------------------------- #
field_names = ['exner','theta','m_v']
colour_scheme = 'Blues'
time_idx = -1
field_label = r'$\theta $ K'
contour_method = 'contour'
# ---------------------------------------------------------------------------- #
# Things that are likely the same for all plots
# ---------------------------------------------------------------------------- #
set_tomplot_style()
data_file = Dataset(results_file_name, 'r')
background_data_file = Dataset(background_file_name, 'r')
# height_file = Dataset(height_file_name, 'r')
# ---------------------------------------------------------------------------- #
# Data extraction
# ---------------------------------------------------------------------------- #
background_dict= {}
data_dict = {}
for field_name in field_names:
    field_data, coords_X, _, coords_Z = \
        extract_lfric_vertical_slice(data_file, field_name, time_idx,
                                    slice_along='y', slice_at=249.99998)
    background_field_data, coords_X, _, coords_Z = \
        extract_lfric_vertical_slice(background_data_file, field_name, time_idx,
                                    slice_along='y', slice_at=249.99998)
    background_dict[field_name] = background_field_data
    data_dict[field_name] = field_data
time = data_file['time'][time_idx]

#breakpoint()

# ---------------------------------------------------------------------------- #
# Data manipulation
# ---------------------------------------------------------------------------- #
print(data_dict)
add_theta_e(data_dict, coords_Z)
add_theta_e(background_dict, coords_Z)


data_dict['theta_e_pert'] = data_dict['theta_e'] - background_dict['theta_e'][0]

field_data = data_dict['theta_e_pert']
field_name='theta_e_pert'
# print(np.linalg.norm(field_data - 0.5)/np.linalg.norm(field_data))
#breakpoint()
# ---------------------------------------------------------------------------- #
# Plot data
# ---------------------------------------------------------------------------- #
fig, ax = plt.subplots(1, 1)
contours = tomplot_contours(field_data)
cmap, lines = tomplot_cmap(contours, colour_scheme)
cf, _ = plot_contoured_field(ax, coords_X, coords_Z, field_data, contour_method,
                             contours, cmap=cmap, line_contours=lines)
add_colorbar_ax(ax, cf, field_label)
tomplot_field_title(ax, f't = {time:.1f}', minmax=True, field_data=field_data)
# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()

