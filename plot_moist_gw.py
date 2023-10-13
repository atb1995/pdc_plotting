"""
Error value and final field plot for the moist gravity wave test case.
"""

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from os.path import abspath, dirname
from scipy.interpolate import griddata
from tomplot import (set_tomplot_style, tomplot_contours, tomplot_cmap,
                     plot_contoured_field, add_colorbar_fig, plot_convergence,
                     tomplot_field_title, extract_lfric_vertical_slice,
                     tomplot_legend_ax)

# ============================================================================ #
# Routine to compute theta_e from prognostic variables
# ============================================================================ #
def calculate_theta_e(data_dict):
    """Takes a dictionary of data (variables are keys) and returns theta_e."""

    Lv = 2.501e6
    cp = 1005.0

    exner = data_dict['exner']
    theta = data_dict['theta']
    mr_v = data_dict['m_v']

    # Need to convert Exner to Wtheta
    exner_wth = np.zeros_like(theta)

    # Assume uniform extrusion
    exner_wth[:,0] = 2.0*exner[:,0] - exner[:,1]
    exner_wth[:,-1] = 2.0*exner[:,-1] - exner[:,-2]
    for j in range(1, np.shape(theta)[1]-1):
        exner_wth[:,j] = 0.5*(exner[:,j]+exner[:,j-1])

    # Finally, calculate theta_e
    T = theta * exner_wth
    exp_arg = Lv * mr_v / (cp * T)
    theta_e = theta * np.exp(exp_arg)

    return theta_e

# ============================================================================ #
# Routines to regrid theta_e to calculate error
# ============================================================================ #

def regrid(old_field_data, old_coords_X, old_coords_Y, new_coords_X, new_coords_Y):
    new_coords = (new_coords_X, new_coords_Y)
    old_coords = (old_coords_X.flatten(), old_coords_Y.flatten())
    old_field_data = old_field_data.flatten()

    nearest_field_data = griddata(old_coords, old_field_data, new_coords, method='nearest')
    new_field_data = griddata(old_coords, old_field_data, new_coords, method='linear')
    new_field_data[np.isnan(new_field_data)] = nearest_field_data[np.isnan(new_field_data)]
    return new_field_data

def calculate_error(field, true_field):
    """Computes the L2 error in a field, relative to some 'true' value."""
    diff_data = field - true_field
    return np.linalg.norm(diff_data) / np.linalg.norm(true_field)

# ============================================================================ #
# Plotting script
# ============================================================================ #

# Specify all of the different meshes / data files that we need
resolution_pairs = [  # (nx_dyn, nx_phys), grouped by nx_dyn
[(600, 1200), (600, 600), (600, 300), (600, 200), (600, 150), (600, 120)],
[(400, 1200), (400, 800), (400, 400), (400, 200), (400, 100), (400, 80), # (400, 50)  # This point makes plot look messy!
],
[(300, 1200), (300, 900), (300, 600), (300, 300), (300, 150), (300, 100), (300, 75), (300, 60)],
[(200, 800), (200, 600), (200, 400), (200, 200), (200, 100), (200, 50)],
[(150, 150)],
[(120, 120)],
[(100, 600), (100, 400), (100, 300), (100, 200), (100, 100), (100, 50)],
[(75, 75)],
[(60, 60)]
]

dx_dict = {1200: 250, 900: 333, 800: 375, 600: 500, 400: 750, 300: 1000,
           200: 1500, 150: 2000, 120: 2500, 100: 3000, 75: 4000, 50: 5000}

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
base_results_dir = '/data/users/abrown/data/pdc_idealised_paper/moist_gw'
plot_dir = '.'
plot_name = f'{plot_dir}/moist_gw.png'

# Create figure here
set_tomplot_style()
fig, axarray = plt.subplots(1, 2, figsize=(12,5))

# Settings needed for both plots
field_names = ['exner','theta','m_v']  # Get coords from m_v so this is last
time_idx = -1

# ============================================================================ #
# Convergence plot
# ============================================================================ #
# ---------------------------------------------------------------------------- #
# Read "true" data
# ---------------------------------------------------------------------------- #
results_dir = f'{base_results_dir}/D1200_P1200'
results_file_name = f'{results_dir}/lfric_diag.nc'
height_file_name = f'{results_dir}/lfric_initial.nc'

data_file = Dataset(results_file_name, 'r')
height_file = Dataset(height_file_name, 'r')

data_dict = {}

for field_name in field_names:
    field_data, true_coords_X, _, true_coords_Z = \
        extract_lfric_vertical_slice(data_file, field_name, time_idx,
                                     slice_along='y', slice_at=250.0,
                                     height_dataset=height_file)
    data_dict[field_name] = field_data

true_theta_e_fine = calculate_theta_e(data_dict)

# ---------------------------------------------------------------------------- #
# Calculate errors
# ---------------------------------------------------------------------------- #

all_errors = []
all_dx_phys = []

for res_pairs in resolution_pairs:
    errors = []
    dx_phys = []

    for i, res_pair in enumerate(res_pairs):
        print(f'Computing error for {res_pair}')

        results_dir = f'{base_results_dir}/D{res_pair[0]}_P{res_pair[1]}'
        results_file_name = f'{results_dir}/lfric_diag.nc'
        height_file_name = f'{results_dir}/lfric_initial.nc'

        data_file = Dataset(results_file_name, 'r')
        height_file = Dataset(height_file_name, 'r')

        data_dict = {}

        for field_name in field_names:
            field_data, coords_X, _, coords_Z = \
                extract_lfric_vertical_slice(data_file, field_name, time_idx,
                                            slice_along='y', slice_at=250.0,
                                            height_dataset=height_file)
            data_dict[field_name] = field_data

        theta_e = calculate_theta_e(data_dict)

        if i == 0:
            true_theta_e = regrid(true_theta_e_fine, true_coords_X, true_coords_Z,
                                  coords_X, coords_Z)
        nx_phys = res_pair[1]
        error = calculate_error(theta_e, true_theta_e)
        dx_phys.append(300000. / nx_phys)
        errors.append(error)

    all_errors.append(errors)
    all_dx_phys.append(dx_phys)

# ---------------------------------------------------------------------------- #
# Arrange errors into useful way for convergence plot
# ---------------------------------------------------------------------------- #
cdfp_errors = []
cdfp_dx_phys = []
non_pdc_errors = []
non_pdc_dx_phys = []
fdcp_errors = []
fdcp_dx_phys = []

for errors, dx_phys_values, res_pairs in zip(all_errors, all_dx_phys, resolution_pairs):
    for error, dx_phys, res_pair in zip(errors, dx_phys_values, res_pairs):
        if res_pair[0] > res_pair[1]:
            fdcp_errors.append(error)
            fdcp_dx_phys.append(dx_phys)
        elif res_pair[0] == res_pair[1]:
            non_pdc_errors.append(error)
            non_pdc_dx_phys.append(dx_phys)
        else:
            cdfp_errors.append(error)
            cdfp_dx_phys.append(dx_phys)

# ---------------------------------------------------------------------------- #
# Convergence plot
# ---------------------------------------------------------------------------- #
cdfp_label = r'$\Delta x_{dyn} > \Delta x_{phys}$'
non_pdc_label = r'$\Delta x_{dyn} = \Delta x_{phys}$'
fdcp_label = r'$\Delta x_{dyn} < \Delta x_{phys}$'

# Plot lines to join up values
for errors, dx_phys_values in zip(all_errors, all_dx_phys):
    plot_convergence(axarray[0], dx_phys_values, errors, color='black',
                     marker='', linestyle_between_points='--',
                     log_by='data', best_fit=False)

# Plot points based on CDFP/FDCP/non-PDC
markersize=8
plot_convergence(axarray[0], cdfp_dx_phys, cdfp_errors, label=cdfp_label,
                 marker='o', color='red', log_by='data', best_fit=False,
                 markersize=markersize)
plot_convergence(axarray[0], non_pdc_dx_phys, non_pdc_errors, label=non_pdc_label,
                 marker='^', color='black', log_by='data', best_fit=False,
                 markersize=markersize)
plot_convergence(axarray[0], fdcp_dx_phys, fdcp_errors, label=fdcp_label,
                 marker='s', color='blue', log_by='data', best_fit=False,
                 markersize=markersize)

axarray[0].legend(ncols=3, loc='lower left', bbox_to_anchor=(-0.1, 1.05, 0.6, 0.2),
                  columnspacing=0.25, handletextpad=0.2)

axarray[0].set_xlabel(r'$\log(\Delta x_{phys} \ /$ m)', labelpad=-10)
axarray[0].set_ylabel(r'$\log(||\theta_e - \theta_e^{true}|| / || \theta_e^{true}||)$', labelpad=-10)
xlims = [6, 8]
ylims = [-17, -14]
axarray[0].set_xticks(xlims)
axarray[0].set_xticklabels(xlims)
axarray[0].set_yticks(ylims)
axarray[0].set_yticklabels(ylims)
axarray[0].grid()

# ============================================================================ #
# Field contour plot
# ============================================================================ #
results_dir = f'{base_results_dir}/D1200_P1200'
results_file_name = f'{results_dir}/lfric_diag.nc'
height_file_name = f'{results_dir}/lfric_initial.nc'
background_file_name = f'{results_dir}/lfric_initial.nc'
# ---------------------------------------------------------------------------- #
# Things that should be altered based on the plot
# ---------------------------------------------------------------------------- #
colour_scheme = 'RdBu_r'
field_label = r"$\theta'_e \ /$ K"
contour_method = 'contour'
# ---------------------------------------------------------------------------- #
# Things that are likely the same for all plots
# ---------------------------------------------------------------------------- #
data_file = Dataset(results_file_name, 'r')
background_data_file = Dataset(background_file_name, 'r')
height_file = Dataset(height_file_name, 'r')
# ---------------------------------------------------------------------------- #
# Data extraction
# ---------------------------------------------------------------------------- #
background_dict = {}
data_dict = {}

for field_name in field_names:
    field_data, coords_X, _, coords_Z = \
        extract_lfric_vertical_slice(data_file, field_name, time_idx,
                                     slice_along='y', slice_at=250.0,
                                     height_dataset=height_file)
    background_field_data, _, _, _ = \
        extract_lfric_vertical_slice(background_data_file, field_name, time_idx,
                                     slice_along='y', slice_at=250.0,
                                     height_dataset=height_file)
    background_dict[field_name] = background_field_data
    data_dict[field_name] = field_data

# ---------------------------------------------------------------------------- #
# Data manipulation
# ---------------------------------------------------------------------------- #
theta_e_final = calculate_theta_e(data_dict)
theta_e_init = calculate_theta_e(background_dict)

# Subtract first column of initial valuse as the background
field_data = theta_e_final - theta_e_init[0]

coords_X /= 1000.
coords_Z /= 1000.

# ---------------------------------------------------------------------------- #
# Plot data
# ---------------------------------------------------------------------------- #
contours = np.linspace(-0.003, 0.003, 10)
cmap, lines = tomplot_cmap(contours, colour_scheme)
cf, _ = plot_contoured_field(axarray[1], coords_X, coords_Z, field_data, contour_method,
                             contours, cmap=cmap, line_contours=lines)
add_colorbar_fig(fig, cf, field_label, cbar_format='.3f')

# ---------------------------------------------------------------------------- #
# Labels
# ---------------------------------------------------------------------------- #
axarray[1].set_xlabel(r'$x \ /$ km', labelpad=-10)
axarray[1].set_ylabel(r'$z \ /$ km', labelpad=-10)
xlims = [-150, 150]
ylims = [0, 10]
axarray[1].set_xlim(xlims)
axarray[1].set_xticks(xlims)
axarray[1].set_xticklabels(xlims)
axarray[1].set_ylim(ylims)
axarray[1].set_yticks(ylims)
axarray[1].set_yticklabels(ylims)

fig.subplots_adjust(wspace=0.15)

# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()

