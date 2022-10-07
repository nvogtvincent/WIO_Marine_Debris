#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualise marine debris results
@author: noam
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import cmasher as cmr
import xarray as xr
from matplotlib import ticker
from osgeo import gdal, osr
from skimage.measure import block_reduce
from matplotlib.gridspec import GridSpec
from sys import argv
from scipy.ndimage import gaussian_filter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import rasterio

# PARAMETERS
param = {'grid_res': 1.0,                          # Grid resolution in degrees
         'lon_range': [20, 130],                   # Longitude range for output
         'lat_range': [-40, 30],                   # Latitude range for output

         # Analysis parameters
         'us_d': argv[1],    # Sinking timescale (days)
         'ub_d': argv[2],      # Beaching timescale (days)

         # Time range
         'y0'  : 1993,
         'y1'  : 2012,

         # Physics
         'mode': argv[3],

         # Source/sink time
         'time': 'sink',

         # Sink sites
         'sites': np.array([13,14,15,16,17,18]),

         # Names
         'sink': 'Seychelles Plateau',
         'class': argv[4],

         # Plot ship tracks
         'tracks': True,

         # Export
         'export': True}

# DIRECTORIES
dirs = {'script': os.path.dirname(os.path.realpath(__file__)),
        'data': os.path.dirname(os.path.realpath(__file__)) + '/../POSTPROC/',
        'fisheries': os.path.dirname(os.path.realpath(__file__)) + '/../FISHERIES/DATA/PROC/',
        'shipping': os.path.dirname(os.path.realpath(__file__)) + '/../SHIPPING/',
        'fig': os.path.dirname(os.path.realpath(__file__)) + '/../FIGURES/'}

# FILE HANDLES
array_str = np.array2string(param['sites'], separator='-').translate({ord(i): None for i in '[]'})
fh = {'flux': dirs['data'] + '/marine_clim_flux_' + param['mode'] + '_' + array_str + '_s' + str(param['us_d']) + '_b' + str(param['ub_d']) + '.nc',
      'debris': dirs['fisheries'] + 'IOTC_monclim.nc',
      'shipping': dirs['shipping'] + 'global_shipping.tif'}

n_years = param['y1']-param['y0']+1
##############################################################################
# LOAD DATA                                                                  #
##############################################################################

fmatrix = xr.open_dataarray(fh['flux'])

# Sum over all sink months
fmatrix = fmatrix.sum(dim='sink_month')

# Find proportion of released debris that arrived at site by dividing by the
# total number released from each cell, i.e. 36*4*12 per year
fmatrix = fmatrix/(1728*n_years)

# Now degrade resolution by factor 12 (-> 1deg)
fmatrix12 = fmatrix.coarsen(latitude=12, longitude=12).mean()

# Open debris input functions (from fisheries)
debris = xr.open_dataset(fh['debris'])

# Open ship tracks
with rasterio.open(fh['shipping']) as src:
    sf = 10
    img = block_reduce(src.read(1), block_size=(sf, sf), func=np.sum)
    height = img.shape[0]
    width = img.shape[1]
    cols, rows = np.meshgrid(np.arange(0, width*10, sf), np.arange(0, height*10, sf))
    xs, ys = rasterio.transform.xy(src.transform, rows, cols)
    lons= np.array(xs)
    lats = np.array(ys)

    # Apply gaussian filter to improve contourf quality
    img = gaussian_filter(img, sigma=1)

##############################################################################
# PLOT                                                                       #
##############################################################################

f = plt.figure(constrained_layout=True, figsize=(27, 11))
gs = GridSpec(2, 4, figure=f, width_ratios=[2, 0.03, 0.98, 0.03], height_ratios=[1, 1])
ax = []
ax.append(f.add_subplot(gs[:, 0], projection = ccrs.PlateCarree())) # Main figure (flux probability)
ax.append(f.add_subplot(gs[:, 1])) # Colorbar for main figure
ax.append(f.add_subplot(gs[0, 2], projection = ccrs.PlateCarree())) # Fisheries 1
ax.append(f.add_subplot(gs[1, 2], projection = ccrs.PlateCarree())) # Fisheries 2
ax.append(f.add_subplot(gs[:, 3])) # Colorbar for fisheries

gl = []
hist = []

land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='k',
                                        facecolor='w',
                                        zorder=1)

land_10k = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='w',
                                        facecolor='k',
                                        zorder=1)

n_lon = param['lon_range'][1] - param['lon_range'][0]
n_lat = param['lat_range'][1] - param['lat_range'][0]

lon_bnd1 = np.linspace(param['lon_range'][0], param['lon_range'][1], num=int(n_lon)+1)
lat_bnd1 = np.linspace(param['lat_range'][0], param['lat_range'][1], num=int(n_lat)+1)

for i in [0, 2, 3]:
    ax[i].set_aspect(1)
    ax[i].set_facecolor('k')

# Total flux
hist.append(ax[0].pcolormesh(lon_bnd1, lat_bnd1,
                             fmatrix12.sum(dim='source_month'), cmap=cmr.sunburst_r,
                             norm=colors.LogNorm(vmin=1e-5, vmax=1e-2),
                             transform=ccrs.PlateCarree()))

ax[0].add_feature(land_10k)
gl.append(ax[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.5, color='white', linestyle='-', zorder=11))
gl[0].xlocator = mticker.FixedLocator(np.arange(-210, 210, 10))
gl[0].ylocator = mticker.FixedLocator(np.arange(-90, 120, 10))
gl[0].xlabels_top = False
gl[0].ylabels_right = False
gl[0].ylabel_style = {'size': 20}
gl[0].xlabel_style = {'size': 20}
ax[0].text(21, -39, 'Likelihood of ' + param['class'] + ' debris beaching at ' + param['sink'], fontsize=28, color='k', fontweight='bold')
ax[0].set_xlim([20, 120])
ax[0].set_ylim([-40, 30])

# Add an overlay with ship tracks
if param['tracks']:
    thresh = 1e7
    img = np.ma.masked_where(img < thresh, img)
    sig_plot = ax[0].contourf(lons, lats, img,
                              levels=np.array([thresh, np.max(img)*2]), colors='none', hatches=['\\\\\\'])


cb0 = plt.colorbar(hist[0], cax=ax[1], fraction=0.05)
cb0.set_label('Mass fraction', size=24)
ax[1].tick_params(axis='y', labelsize=22)

fishing_cmap = cmr.torch_r

# Input from purse seiners
ps_input = debris.ps_effort.transpose('lat', 'lon', 'month').rename({'lat': 'latitude',
                                                                     'lon': 'longitude',
                                                                     'month': 'source_month'})
ps_input.coords['longitude'] = fmatrix12.longitude.values # Have to reset this due to float issue
ps_input.coords['latitude'] = fmatrix12.latitude.values # Have to reset this due to float issue

ps_flux = (ps_input*fmatrix12).sum(dim='source_month')
hist.append(ax[2].pcolormesh(lon_bnd1, lat_bnd1,
                             ps_flux/ps_flux.sum(),
                             cmap=fishing_cmap, norm=colors.LogNorm(vmin=1e-4, vmax=1e-1),
                             transform=ccrs.PlateCarree()))
ax[2].add_feature(land_10k)
ax[2].text(21, -39, param['class'] + ' purse-seine debris', fontsize=22, color='k', fontweight='bold', zorder=12)
ax[2].set_facecolor('w')
ax[2].set_extent([20, 120, -40, 30], crs=ccrs.PlateCarree())

# Input from drifting and set longlines
ll_input = debris.ll_effort.transpose('lat', 'lon', 'month').rename({'lat': 'latitude',
                                                                     'lon': 'longitude',
                                                                     'month': 'source_month'})
ll_input.coords['longitude'] = fmatrix12.longitude.values # Have to reset this due to float issue
ll_input.coords['latitude'] = fmatrix12.latitude.values # Have to reset this due to float issue

ll_flux = (ll_input*fmatrix12).sum(dim='source_month')
ll_flux = ll_flux.coarsen(longitude=5, latitude=5).sum() # Re-degrade

hist.append(ax[3].pcolormesh(lon_bnd1[::5], lat_bnd1[::5],
                             ll_flux/ll_flux.sum(),
                             cmap=fishing_cmap, norm=colors.LogNorm(vmin=1e-4, vmax=1e-1),
                             transform=ccrs.PlateCarree()))
ax[3].add_feature(land_10k)
ax[3].text(21, -39, param['class'] + ' longline debris', fontsize=22, color='k', fontweight='bold', zorder=12)
ax[3].set_facecolor('w')
ax[3].set_extent([20, 120, -40, 30], crs=ccrs.PlateCarree())

cb1 = plt.colorbar(hist[1], cax=ax[4])
cb1.set_label('Normalised risk from fishery', size=24)
ax[4].tick_params(axis='y', labelsize=22)

for ii, i in enumerate([2, 3]):
    gl.append(ax[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                              linewidth=0.5, color='k', linestyle='-', zorder=11))
    gl[ii+1].xlocator = mticker.FixedLocator(np.arange(0, 240, 20))
    gl[ii+1].ylocator = mticker.FixedLocator(np.arange(-80, 80, 20))
    gl[ii+1].ylabel_style = {'size': 20}
    gl[ii+1].xlabel_style = {'size': 20}
    gl[ii+1].xlabels_top = False
    gl[ii+1].ylabels_right = False
    gl[ii+1].ylabels_left = True
    if ii == 1:
        gl[ii+1].xlabels_bottom = True
    else:
        gl[ii+1].xlabels_bottom = False

if param['export']:
    plt.savefig(dirs['fig'] + 'marine_sources_' + param['mode'] + '_' + np.array2string(param['sites'], separator='-') + '_s' + str(param['us_d']) + '_b' + str(param['ub_d']) + '.pdf', bbox_inches='tight', dpi=300)








