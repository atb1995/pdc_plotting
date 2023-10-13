"""
A tomplot example, for making a quick single plot from LFRic data.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, tomplot_contours, tomplot_cmap,
                     plot_contoured_field, add_colorbar_fig,
                     tomplot_field_title, extract_lfric_field, extract_lfric_coords)
def get_surface_pressure(exner_l0, exner_l1, height_wth_l0, height_w3_l0, height_w3_l1):
    
    # Convert to hPa
    rd = 287.05
    p0 = 100000.0
    kappa = rd/1005.0
    pressure_l0 = 0.01*exner_l0**(1.0/kappa) * p0
    pressure_l1 = 0.01*exner_l1**(1.0/kappa) * p0

    # Vertically interpolate pressure to surface
    pressure = pressure_l0 + (height_wth_l0 - height_w3_l0)*(pressure_l0 - pressure_l1)/(height_w3_l0 - height_w3_l1)

    return pressure

def get_surface_temperature(theta_l0, exner_l0, exner_l1, height_wth_l0, height_w3_l0, height_w3_l1):
    
    # Vertically interpolate exner to surface
    exner = exner_l0 + (height_wth_l0 - height_w3_l0)*(exner_l0 - exner_l1)/(height_w3_l0 - height_w3_l1)

    # Calculate temperature
    temperature = exner*theta_l0

    return temperature

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
results_dir = '/data/users/abrown/data/pdc_idealised_paper/moist_baroclinic'
plot_dir = '.'
res_title = ['D96_P96', 'D48_P96', 'D48_P48', 'D48_P24', 'D24_P24']
results_file_names = [f'{results_dir}/{res_name}/lfric_diag.nc' for res_name in res_title]
height_file_names = [f'{results_dir}/{res_name}/lfric_initial.nc' for res_name in res_title]
plot_name = f'{plot_dir}/moist_baroclinic_orog.png'

#------------------------------------------------------------------------------#
# Create axes                                                                  #
#------------------------------------------------------------------------------#
set_tomplot_style()
fig, ax = plt.subplots(5, 2, figsize=(12, 12), sharex='col', sharey='row')

for i in range(5):
    
    # ---------------------------------------------------------------------------- #
    # Things that should be altered based on the plot
    # ---------------------------------------------------------------------------- #
    field_names = ['theta', 'exner', 'm_cl']
    colour_schemes = ['OrRd', 'OrRd','Blues' ] 
    time_idx = -1
    field_labels = [r'T / K',r'$m_{cl} \ / $ kg kg$^{-1}$',]
    contour_method = 'tricontour'
    
    # ---------------------------------------------------------------------------- #
    # Things that are likely the same for all plots
    # ---------------------------------------------------------------------------- #
    data_file = Dataset(results_file_names[i], 'r')
    height_file = Dataset(height_file_names[i], 'r')

    # ---------------------------------------------------------------------------- #
    # Data extraction
    # ---------------------------------------------------------------------------- #
    theta_l0= \
        extract_lfric_field(data_file, field_names[0], time_idx, 0)
    exner_l0= \
        extract_lfric_field(data_file, field_names[1], time_idx, 0)
    exner_l1 = \
        extract_lfric_field(data_file, field_names[1], time_idx, 1)
    height_w3_l0 = \
        extract_lfric_field(height_file, 'height_w3', time_idx, 0)
    height_w3_l1 = \
        extract_lfric_field(height_file, 'height_w3', time_idx, 1)
    height_wth_l0 = \
        extract_lfric_field(height_file, 'height_wth', time_idx, 0)
    field_data = get_surface_temperature(theta_l0, exner_l0, exner_l1, height_wth_l0, height_w3_l0, height_w3_l1)
    
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
    contours = np.linspace(220, 330, 12)
    cmap, lines = tomplot_cmap(contours, colour_schemes[0], cmap_rescale_type= 'bottom')
    cf1, _ = plot_contoured_field(ax[i,0], coords_X, coords_Y, field_data, contour_method,
                                contours, cmap=cmap, plot_contour_lines=False)
    
    xticks = [-120, 0, 120]
    yticks = [10, 80]
    ax[i,0].set_xticks(xticks) 
    ax[i,0].set_xticklabels(xticks)
    ax[i,1].set_xticks(xticks) 
    ax[i,1].set_xticklabels(xticks)
    
    ax[i,1].set_yticks(yticks) 
    ax[i,1].set_yticklabels(yticks)
    ax[i,0].set_yticks(yticks) 
    ax[i,0].set_yticklabels(yticks)
    
    # ---------------------------------------------------------------------------- #
    # Data extraction
    # ---------------------------------------------------------------------------- #
    field_data = get_surface_pressure(exner_l0, exner_l1, height_wth_l0, height_w3_l0, height_w3_l1)
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
    contours = np.linspace(900, 1020, 9)
    cmap, lines = tomplot_cmap(contours, colour_schemes[1])
    cf, _ = plot_contoured_field(ax[i,0], coords_X, coords_Y, field_data, contour_method,
                                contours, cmap=cmap, plot_filled_contours=False)

    ax[i,0].set_ylabel(r'$\phi \ / $ deg', labelpad=-10)
    tomplot_field_title(ax[i,0], res_title[i].replace('_', ' '))

    # ---------------------------------------------------------------------------- #
    # Data extraction
    # ---------------------------------------------------------------------------- #
    field_data = \
        extract_lfric_field(data_file, field_names[2], time_idx, 9)
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
    contours = np.linspace(0.0, 0.005, 6)
    cmap, lines = tomplot_cmap(contours, colour_schemes[2])
    cf, _ = plot_contoured_field(ax[i,1], coords_X, coords_Y, field_data, contour_method,
                                contours, cmap=cmap,
                                line_contours=lines)

    tomplot_field_title(ax[i,1], res_title[i].replace('_', ' '))


ax[-1,0].set_xlabel(r'$\lambda \ / $ deg')
ax[-1,1].set_xlabel(r'$\lambda \ / $ deg')

# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
cbar_padding = 0.15

# Pad width before adding cbars
plt.subplots_adjust(wspace=0.1)

add_colorbar_fig(fig, cf, field_labels[1], location='bottom', ax_idxs=[9],
                 cbar_ticks=[0,0.005], cbar_format=".4g", cbar_padding=cbar_padding)
add_colorbar_fig(fig, cf1, field_labels[0], location='bottom', ax_idxs=[8],
                 cbar_ticks=[220,330], cbar_format=".3g", cbar_padding=cbar_padding)
plt.subplots_adjust(bottom=0.05,
                    hspace=0.2)
fig.savefig(plot_name, bbox_inches='tight')
plt.close()


