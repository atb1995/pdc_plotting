"""
Error value and final field plot for the moist gravity wave test case.
"""

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from os.path import abspath, dirname
from scipy.interpolate import griddata
from tomplot import (set_tomplot_style, tomplot_contours, tomplot_cmap,
                     plot_contoured_field, add_colorbar_fig,
                     tomplot_field_title, extract_lfric_field,
                     extract_lfric_coords, plot_field_quivers,
                     regrid_horizontal_slice)

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
results_dir = '/data/users/abrown/data/pdc_idealised_paper/transport/stretchy_conv'
plot_dir = '.'
plot_name = f'{plot_dir}/coarse_transport.png'

# Create figure here
set_tomplot_style()
fig, axarray = plt.subplots(1, 3, sharey='row', figsize=(12,5))

# Settings needed for both plots
field_name = 'w3_aerosol'
wind_names = ['u_in_w2h', 'v_in_w2h']
time_idxs = [None, 'halfway', -1]
titles = [r'$t=0$', r'$t=\tau/2$', r'$t=\tau$']
contours = np.linspace(10, 80, 13)
colour_scheme = 'OrRd'
field_label = r"$\bar{a}_Y \ /$ g kg$^{-1}$"
contour_method = 'tricontour'
file_names = [f'{results_dir}/lfric_initial.nc', f'{results_dir}/lfric_diag.nc', f'{results_dir}/lfric_diag.nc']

for i, (ax, file_name, title, time_idx) in \
        enumerate(zip(axarray, file_names, titles, time_idxs)):

    data_file = Dataset(file_name, 'r')
    if time_idx == 'halfway':
        time_idx = int((len(data_file['time_instant'][:])+1)/2 - 1)

    # ------------------------------------------------------------------------ #
    # Scalar data extraction
    # ------------------------------------------------------------------------ #
    field_data = extract_lfric_field(data_file, field_name, time_idx=time_idx, level=0)
    coords_X, coords_Y = extract_lfric_coords(data_file, field_name)

    # Change units
    field_data *= 1000.0

    # ------------------------------------------------------------------------ #
    # Plot scalar data
    # ------------------------------------------------------------------------ #
    cmap, lines = tomplot_cmap(contours, colour_scheme)
    cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data, contour_method,
                                 contours, cmap=cmap, line_contours=lines)

    # ------------------------------------------------------------------------ #
    # Labels
    # ------------------------------------------------------------------------ #

    if i == 0:
        ylims = [-90, 90]
        ax.set_ylabel(r'$\phi \ / $ deg', labelpad=-10)
        ax.set_ylim(ylims)
        ax.set_yticks(ylims)
        ax.set_yticklabels(ylims)

    xlims = [-180, 180]
    ax.set_xlabel(r'$\lambda \ /$ deg', labelpad=-10)
    ax.set_xlim(xlims)
    ax.set_xticks(xlims)
    ax.set_xticklabels(xlims)
    ax.set_title(title)

    # ------------------------------------------------------------------------ #
    # Vector data extraction
    # ------------------------------------------------------------------------ #
    field_data_X = extract_lfric_field(data_file, wind_names[0], time_idx=time_idx, level=0)
    field_data_Y = extract_lfric_field(data_file, wind_names[1], time_idx=time_idx, level=0)
    old_coords_X, old_coords_Y = extract_lfric_coords(data_file, wind_names[0])

    # Need to regrid to lon-lat grid to get quivers that look good on sphere
    new_X_1d = np.linspace(-180, 180, 60)
    new_Y_1d = np.linspace(-90, 90, 60)
    new_coords_X, new_coords_Y = np.meshgrid(new_X_1d, new_Y_1d, indexing='ij')
    field_data_X = regrid_horizontal_slice(new_coords_X, new_coords_Y,
                                           old_coords_X, old_coords_Y,
                                           field_data_X, periodic_fix='sphere')
    field_data_Y = regrid_horizontal_slice(new_coords_X, new_coords_Y,
                                           old_coords_X, old_coords_Y,
                                           field_data_Y, periodic_fix='sphere')

    # ------------------------------------------------------------------------ #
    # Plot vector data
    # ------------------------------------------------------------------------ #
    _ = plot_field_quivers(ax, new_coords_X, new_coords_Y, field_data_X, field_data_Y,
                           magnitude_filter=5.0, spatial_filter_step=4, spatial_filter_offset=1)


add_colorbar_fig(fig, cf, field_label, cbar_format='.0f')

fig.subplots_adjust(wspace=0.15)

# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()

