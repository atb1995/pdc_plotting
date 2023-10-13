"""
A tomplot example, for making a quick single plot from LFRic data.
"""

import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, tomplot_contours, tomplot_cmap,
                     plot_contoured_field, add_colorbar_ax, add_colorbar_fig,
                    extract_lfric_field, tomplot_field_title, extract_lfric_vertical_slice,
                    extract_lfric_coords)
def get_zonal_average(data_file, time_idx, field_name, full):
    if (len(data_file[field_name].dimensions) == 2):
        num_levels = np.shape(data_file[field_name])[0]
    elif (len(data_file[field_name].dimensions) == 3):
        num_levels = np.shape(data_file[field_name])[1]
    else:
        raise RuntimeError('extract_lfric_vertical_slice: cannot work with '
                            + 'data of this shape')
    levels = range(num_levels)
    nx = 720
    ny = 360 
    coords_X, coords_Y= \
            extract_lfric_coords(data_file, field_name)
    field_data=np.zeros((ny,nx,num_levels))
    xmin = np.amin(coords_X)
    xmax = np.amax(coords_X)
    ymin = np.amin(coords_Y)
    ymax = np.amax(coords_Y)
    # Generate a regular grid to interpolate the data.
    xi = np.linspace(xmin, xmax, nx)
    yi = np.linspace(ymin, ymax, ny)

    xf, yf = np.meshgrid(xi, yi)

    # dcmip Stretched grid
    mu = 15.
    n_full=31
    zf = np.zeros([n_full])
    lid=30
    for k in range(n_full):
        zf[k] = lid * (np.sqrt(mu*(float(k)/float(n_full-1))**2 + 1.) - 1.)/(np.sqrt(mu+1.) - 1)
    zh = 0.5*(zf[1:] + zf[0:n_full-1])
    # ------------------------------------------------------------------------ #
    # Loop through levels and build up data
    # ------------------------------------------------------------------------ #
    for lev_idx, level in enumerate(levels):
        # Extract field data for this level
        field_data_level = extract_lfric_field(data_file, field_name,
                                            time_idx=time_idx, level=level)
        fi = griddata((coords_X, coords_Y), field_data_level, (xf, yf), method='linear')
        fi_n = griddata((coords_X, coords_Y), field_data_level, (xf, yf), method='nearest')
        fi[np.isnan(fi)] = fi_n[np.isnan(fi)]
        field_data[:,:, level] = field_data[:,:,level] + fi
    mean_field_data = np.mean(field_data,axis=1)
    if(full):
        z = zf
    else:
        z = zh
    lat, height = np.meshgrid(yi, z)

    return mean_field_data.T,lat, height
# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
# When copying this example these should not be relative to this file
results_dir = '/data/users/abrown/data/pdc_idealised_paper/held_suarez'
plot_dir = '.'
res_title = ['D96_P96', 'D48_P96', 'D48_P48', 'D48_P24', 'D24_P24']
results_file_names = [f'{results_dir}/{res_name}/lfric_averaged.nc' for res_name in res_title]
plot_name = f'{plot_dir}/held_suarez.png'
# ---------------------------------------------------------------------------- #
# Things that should be altered based on the plot
# ---------------------------------------------------------------------------- #
field_names = ['theta', 'u_in_w2h']
colour_schemes = ['OrRd','RdBu_r']
time_idx = -1
field_labels = [r'$\theta \ / $ K', r'$u \ / $ m s$^{-1}$']
contour_method = 'contour'

set_tomplot_style()
fig, ax = plt.subplots(5, 2, figsize=(12, 12), sharex='col', sharey='row')

for i in range(5):
    # ---------------------------------------------------------------------------- #
    # Things that are likely the same for all plots
    # ---------------------------------------------------------------------------- #
    data_file = Dataset(results_file_names[i], 'r')
    # ---------------------------------------------------------------------------- #
    # Data extraction
    # ---------------------------------------------------------------------------- #

    print(f'Compute zonal average for {field_names[0]}, {res_title[i]}')
    field_data, coords_Y, coords_Z = get_zonal_average(data_file, time_idx, field_names[0],full=True)
    time = data_file['time'][time_idx]
    # ---------------------------------------------------------------------------- #
    # Plot data
    # ---------------------------------------------------------------------------- #
    contours = tomplot_contours(field_data)
    cmap, lines = tomplot_cmap(contours, colour_schemes[0])
    cf1, _ = plot_contoured_field(ax[i, 0], coords_Y, coords_Z, field_data, contour_method,
                                contours, cmap=cmap, line_contours=lines)
    ax[i,0].set_ylabel(r'$z \ /$ km', labelpad=-10)
    tomplot_field_title(ax[i,0], res_title[i].replace('_', ' '))

    yticks = [0, 30]
    ax[i,0].set_yticks(yticks) 
    ax[i,0].set_yticklabels(yticks)
    ax[i,1].set_yticks(yticks) 
    ax[i,1].set_yticklabels(yticks)

    # ---------------------------------------------------------------------------- #
    # Data extraction
    # ---------------------------------------------------------------------------- #

    print(f'Compute zonal average for {field_names[1]}, {res_title[i]}')
    field_data, coords_Y, coords_Z = get_zonal_average(data_file, time_idx, field_names[1],full=False)
    time = data_file['time'][time_idx]
    # ---------------------------------------------------------------------------- #
    # Plot data
    # ---------------------------------------------------------------------------- #
    contours = tomplot_contours(field_data, divergent_flag=True)
    cmap, lines = tomplot_cmap(contours, colour_schemes[1], remove_contour=0)
    cf, _ = plot_contoured_field(ax[i,1], coords_Y, coords_Z, field_data, contour_method,
                                contours, cmap=cmap, line_contours=lines)
    tomplot_field_title(ax[i,1], res_title[i].replace('_', ' '))
# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
ax[-1,0].set_xlabel(r'$\phi \ / $ deg')
ax[-1,1].set_xlabel(r'$\phi \ / $ deg')
print(f'Saving figure to {plot_name}')
# Pad width before adding cbars
plt.subplots_adjust(wspace=0.1)
cbar_padding = 0.15
add_colorbar_fig(fig, cf, field_labels[1], location='bottom', ax_idxs=[9],
                 cbar_ticks=[-40,40], cbar_format=".2g", cbar_padding=cbar_padding)
add_colorbar_fig(fig, cf1, field_labels[0], location='bottom', ax_idxs=[8],
                 cbar_ticks=[200,800], cbar_format=".3g", cbar_padding=cbar_padding)
plt.subplots_adjust(bottom=0.05,
                    hspace=0.2)
fig.savefig(plot_name, bbox_inches='tight')
plt.close()


