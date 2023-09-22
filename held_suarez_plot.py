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
        z=zf
    else:
        z=zh
    lat, height = np.meshgrid(yi, z)
    #breakpoint()
    return mean_field_data.T,lat, height
# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
# When copying this example these should not be relative to this file
plot_dir = f'{abspath(dirname(__file__))}'
results_file_names = ['/data/users/abrown/data/held_suarez/lfric_averaged_C96.nc', \
                      '/data/users/abrown/data/held_suarez_pdc/lfric_averaged_C48_cdfp.nc', \
                      '/data/users/abrown/data/held_suarez/lfric_averaged_C48.nc', \
                      '/data/users/abrown/data/held_suarez_pdc/lfric_averaged_C48_fdcp.nc', \
                      '/data/users/abrown/data/held_suarez/lfric_averaged_C24.nc']
plot_name = f'{plot_dir}/held_suarez.png'
# ---------------------------------------------------------------------------- #
# Things that should be altered based on the plot
# ---------------------------------------------------------------------------- #
field_names = ['theta', 'u_in_w2h']
colour_schemes = ['OrRd','RdBu_r']
time_idx = -1
field_labels = [r'$\theta \ / $ K', r'$u \ / $ m s$^{-1}$']
contour_method = 'contour'
res_title = ["C96", "D48 P96", "C48", "D48 P24", "C24"]

fig, ax = plt.subplots(5, 2, figsize=(8.27, 11.69), sharex='col', sharey='row')

for i in range(5):
    # ---------------------------------------------------------------------------- #
    # Things that are likely the same for all plots
    # ---------------------------------------------------------------------------- #
    set_tomplot_style()
    data_file = Dataset(results_file_names[i], 'r')
    # ---------------------------------------------------------------------------- #
    # Data extraction
    # ---------------------------------------------------------------------------- #

    field_data, coords_Y, coords_Z = get_zonal_average(data_file, time_idx, field_names[0],full=True)
    time = data_file['time'][time_idx]
    # ---------------------------------------------------------------------------- #
    # Plot data
    # ---------------------------------------------------------------------------- #
    contours = tomplot_contours(field_data)
    cmap, lines = tomplot_cmap(contours, colour_schemes[0])
    cf1, _ = plot_contoured_field(ax[i, 0], coords_Y, coords_Z, field_data, contour_method,
                                contours, cmap=cmap, line_contours=lines)
    ax[i,0].set_ylabel(r'$z \ /$ km', fontsize=16)
    tomplot_field_title(ax[i,0], res_title[i])
    
    ax[i,0].set_yticks([0,30]) 
    ax[i,0].set_yticklabels([0,30], fontsize=12)
    ax[i,1].set_yticks([0,30]) 
    ax[i,1].set_yticklabels([0,30], fontsize=12)


    # ---------------------------------------------------------------------------- #
    # Data extraction
    # ---------------------------------------------------------------------------- #

    field_data, coords_Y, coords_Z = get_zonal_average(data_file, time_idx, field_names[1],full=False)
    time = data_file['time'][time_idx]
    # ---------------------------------------------------------------------------- #
    # Plot data
    # ---------------------------------------------------------------------------- #
    contours = tomplot_contours(field_data, divergent_flag=True)
    cmap, lines = tomplot_cmap(contours, colour_schemes[1], remove_contour=0)
    cf, _ = plot_contoured_field(ax[i,1], coords_Y, coords_Z, field_data, contour_method,
                                contours, cmap=cmap, line_contours=lines)
    tomplot_field_title(ax[i,1], res_title[i])
# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
ax[-1,0].set_xlabel('Latitude', fontsize=15)
ax[-1,1].set_xlabel('Latitude', fontsize=15)
print(f'Saving figure to {plot_name}')
set_tomplot_style(fontsize=16, family='serif', usetex=True)
plt.subplots_adjust(left=0.1,
                    bottom=0.05,
                    right=0.9,
                    top=0.9,
                    wspace=0.2,
                    hspace=0.4)
add_colorbar_fig(fig, cf, field_labels[1], location='bottom', ax_idxs=[9],cbar_ticks=[-40,40], cbar_format=".2g")
add_colorbar_fig(fig, cf1, field_labels[0], location='bottom', ax_idxs=[8],cbar_ticks=[200,800], cbar_format=".3g")
plt.subplots_adjust(left=0.1,
                    bottom=0.05,
                    right=0.9,
                    top=0.9,
                    wspace=0.2,
                    hspace=0.4)
fig.savefig(plot_name, bbox_inches='tight')
plt.close()


