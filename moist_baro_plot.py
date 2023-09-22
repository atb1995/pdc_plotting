"""
A tomplot example, for making a quick single plot from LFRic data.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, tomplot_contours, tomplot_cmap,
                     plot_contoured_field, add_colorbar_fig,
                     tomplot_field_title, extract_lfric_field, extract_lfric_coords)

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
# When copying this example these should not be relative to this file
results_dir = '/hpc/scratch/d02/abrown/cylc-run/r40541_pdc_main_branch-gungho_model-meto-xc40-pdc_idealised/share/output/intel_64-bit_fast-debug/gungho_model'
plot_dir = f'.'
results_file_names = [f'{results_dir}/moist_baroclinic_orog/C96_MG_dt-900p0/results/lfric_diag.nc', f'{results_dir}/moist_baroclinic_orog_pdc/C48_cdfp_dt-900p0/results/lfric_diag.nc' \
                      ,f'{results_dir}/moist_baroclinic_orog/C48_MG_dt-900p0/results/lfric_diag.nc',f'{results_dir}/moist_baroclinic_orog_pdc/C48_fdcp_dt-900p0/results/lfric_diag.nc', \
                      '/hpc/scratch/d02/abrown/cylc-run/r40541_pdc_main_branch-gungho_model-meto-xc40-moist_baro_C24/share/output/intel_64-bit_fast-debug/gungho_model/moist_baroclinic_orog/C24_MG_dt-900p0/results/lfric_diag.nc']
plot_name = f'{plot_dir}/moist_baroclinic_orog.png'
res_title = ['C96', "D48 P96", "C48", "D48 P24", "C24"]

#------------------------------------------------------------------------------#
# Create axes                                                                  #
#------------------------------------------------------------------------------#

fig, ax = plt.subplots(5, 2, figsize=(8.27, 11.69), sharex='col', sharey='row')

for i in range(5):
    # ---------------------------------------------------------------------------- #
    # Things that should be altered based on the plot
    # ---------------------------------------------------------------------------- #
    field_names = ['theta', 'exner', 'm_cl']
    colour_schemes = ['OrRd', 'OrRd','Blues' ] 
    time_idx = -1
    field_labels = [r'$\theta \ /$ K',r'$m_{cl} \ / $ kg kg$^{-1}$',]
    levels= [2,0,9]
    contour_method = 'tricontour'
    # ---------------------------------------------------------------------------- #
    # Things that are likely the same for all plots
    # ---------------------------------------------------------------------------- #
    set_tomplot_style()
    data_file = Dataset(results_file_names[i], 'r')

    # ---------------------------------------------------------------------------- #
    # Data extraction
    # ---------------------------------------------------------------------------- #
    field_data = \
        extract_lfric_field(data_file, field_names[0], time_idx, levels[0])
    coords_X, coords_Y= \
        extract_lfric_coords(data_file, field_names[0])
    time = data_file['time'][time_idx]
    # ---------------------------------------------------------------------------- #
    # Data manipulation
    # ---------------------------------------------------------------------------- #
    df = pd.DataFrame({'coords_X': coords_X,
                       'coords_Y': coords_Y,
                       'field_data': field_data})
    coord_filter_x = [-150, 180] 
    coord_filter_y = [0, 90]
    if coord_filter_x is not None:
        df = df[(df['coords_X'] > coord_filter_x[0]) & (df['coords_X'] < coord_filter_x[1])]
    if coord_filter_y is not None:
        df = df[(df['coords_Y'] > coord_filter_y[0]) & (df['coords_Y'] < coord_filter_y[1])]

    coords_X = df['coords_X'].values
    coords_Y = df['coords_Y'].values
    field_data = df['field_data'].values
    # ---------------------------------------------------------------------------- #
    # Plot data
    # ---------------------------------------------------------------------------- #

    contours = tomplot_contours(field_data, min_num_bins=6)
    cmap, lines = tomplot_cmap(contours, colour_schemes[0], cmap_rescale_type= 'bottom')
    cf1, _ = plot_contoured_field(ax[i,0], coords_X, coords_Y, field_data, contour_method,
                                contours, cmap=cmap, plot_contour_lines=False)
    
    ax[i,0].set_xticks([-140,0,140]) 
    ax[i,0].set_xticklabels([-140,0,140], fontsize=12)
    ax[i,1].set_xticks([-140,0,140]) 
    ax[i,1].set_xticklabels([-140,0,140], fontsize=12)
    
    ax[i,1].set_yticks([10,80]) 
    ax[i,1].set_yticklabels([10,80], fontsize=12)
    ax[i,0].set_yticks([10,80]) 
    ax[i,0].set_yticklabels([10,80], fontsize=12)
    # ---------------------------------------------------------------------------- #
    # Data extraction
    # ---------------------------------------------------------------------------- #
    field_data = \
        extract_lfric_field(data_file, field_names[1], time_idx, levels[1])
    coords_X, coords_Y= \
        extract_lfric_coords(data_file, field_names[1])
    time = data_file['time'][time_idx]
    # ---------------------------------------------------------------------------- #
    # Data manipulation
    # ---------------------------------------------------------------------------- #
    df = pd.DataFrame({'coords_X': coords_X,
                    'coords_Y': coords_Y,
                    'field_data': field_data})
    coord_filter_x = [-150, 180] 
    coord_filter_y = [0, 90]
    if coord_filter_x is not None:
        df = df[(df['coords_X'] > coord_filter_x[0]) & (df['coords_X'] < coord_filter_x[1])]
    if coord_filter_y is not None:
        df = df[(df['coords_Y'] > coord_filter_y[0]) & (df['coords_Y'] < coord_filter_y[1])]

    coords_X = df['coords_X'].values
    coords_Y = df['coords_Y'].values
    field_data = df['field_data'].values
    # ---------------------------------------------------------------------------- #
    # Plot data
    # ---------------------------------------------------------------------------- #

    contours = tomplot_contours(field_data, min_num_bins=10)
    cmap, lines = tomplot_cmap(contours, colour_schemes[1])
    cf, _ = plot_contoured_field(ax[i,0], coords_X, coords_Y, field_data, contour_method,
                                contours, cmap=cmap, plot_filled_contours=False)

    ax[i,0].set_ylabel('Latitude', fontsize=15)
    tomplot_field_title(ax[i,0], res_title[i])
    # ---------------------------------------------------------------------------- #
    # Data extraction
    # ---------------------------------------------------------------------------- #
    field_data = \
        extract_lfric_field(data_file, field_names[2], time_idx, levels[2])
    coords_X, coords_Y= \
        extract_lfric_coords(data_file, field_names[2])
    time = data_file['time'][time_idx]
    # ---------------------------------------------------------------------------- #
    # Data manipulation
    # ---------------------------------------------------------------------------- #
    df = pd.DataFrame({'coords_X': coords_X,
                    'coords_Y': coords_Y,
                    'field_data': field_data})
    coord_filter_x = [-150, 180] 
    coord_filter_y = [0, 90]
    if coord_filter_x is not None:
        df = df[(df['coords_X'] > coord_filter_x[0]) & (df['coords_X'] < coord_filter_x[1])]
    if coord_filter_y is not None:
        df = df[(df['coords_Y'] > coord_filter_y[0]) & (df['coords_Y'] < coord_filter_y[1])]

    coords_X = df['coords_X'].values
    coords_Y = df['coords_Y'].values
    field_data = df['field_data'].values
    # ---------------------------------------------------------------------------- #
    # Plot data
    # ---------------------------------------------------------------------------- #

    contours = tomplot_contours(field_data, min_num_bins=4)
    cmap, lines = tomplot_cmap(contours, colour_schemes[2])
    cf, _ = plot_contoured_field(ax[i,1], coords_X, coords_Y, field_data, contour_method,
                                contours, cmap=cmap,
                                line_contours=lines)

    tomplot_field_title(ax[i,1], res_title[i])


ax[-1,0].set_xlabel('Longitude', fontsize=15)
ax[-1,1].set_xlabel('Longitude', fontsize=15)

# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
set_tomplot_style(fontsize=16, family='serif', usetex=True)
plt.subplots_adjust(left=0.1,
                    bottom=0.05,
                    right=0.9,
                    top=0.9,
                    wspace=0.2,
                    hspace=0.4)
add_colorbar_fig(fig, cf, field_labels[1], location='bottom', ax_idxs=[9],cbar_ticks=[0,0.005], cbar_format=".4g")
add_colorbar_fig(fig, cf1, field_labels[0], location='bottom', ax_idxs=[8],cbar_ticks=[240,320], cbar_format=".3g")
plt.subplots_adjust(left=0.1,
                    bottom=0.05,
                    right=0.9,
                    top=0.9,
                    wspace=0.2,
                    hspace=0.4)
fig.savefig(plot_name, bbox_inches='tight')
plt.close()


