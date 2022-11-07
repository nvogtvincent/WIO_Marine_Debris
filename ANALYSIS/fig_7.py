#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Carry out EOF analysis on MD data
@author: noam
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cmasher as cmr
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xskillscore as xs
from matplotlib.gridspec import GridSpec
from sys import argv
from scipy import signal
from scipy.fft import rfft, rfftfreq

# PARAMETERS
param = {'degrade': 12,                            # Degradation factor
         'lon_range': [20, 130],                   # Longitude range for output
         'lat_range': [-40, 30],                   # Latitude range for output

         # Analysis parameters
         'us_d': int(argv[1]),    # Sinking timescale (days)
         'ub_d': int(argv[2]),     # Beaching timescale (days)

         # Physics
         'physics': argv[3],

         # Source/sink time
         'time': argv[4],

         # Sink sđites
         'sites': np.array([13,14,15,16,17,18]),

         # Set significance threshold (for log transform)
         'sig_thresh': 1e-9}

try:
    param['name'] = argv[5] + ' '
except:
    param['name'] = ''

# DIRECTORIES
dirs = {'script': os.path.dirname(os.path.realpath(__file__)),
        'data': os.path.dirname(os.path.realpath(__file__)) + '/../MATRICES/',
        'ref': os.path.dirname(os.path.realpath(__file__)) + '/../REFERENCE/',
        'fig': os.path.dirname(os.path.realpath(__file__)) + '/../FIGURES/'}

# FILE HANDLES
array_str = np.array2string(param['sites'], separator='-').translate({ord(i): None for i in '[]'})
fh = {'flux': dirs['data'] + '/marine_' + param['time'] + '_flux_' + param['physics'] + '_' + array_str + '_s' + str(param['us_d']) + '_b' + str(param['ub_d']) + '.nc'}
time_var = 'sink_time' if param['time'] == 'sink' else 'source_time'

##############################################################################
# LOAD DATA                                                                  #
##############################################################################

# Only use 1995-2012 for sink (2-year spinup, NOT valid for us > 1 year)
if param['time'] == 'sink':
    fmatrix = xr.open_dataarray(fh['flux'])[:, :, 24:240]
else:
    fmatrix = xr.open_dataarray(fh['flux'])[:]

# Degrade resolution
fmatrix = fmatrix.coarsen(dim={'longitude': param['degrade'], 'latitude': param['degrade'], time_var: 1}).sum()
fmatrix = fmatrix.transpose(time_var, 'longitude', 'latitude')

lons, lats = np.meshgrid(fmatrix.coords['longitude'], fmatrix.coords['latitude'])

# Log transform data:
# We release 36*4*(12^2) particles per degree grid cell (except for edge cases)
# so let's set everything less than 0.1% of this to 'negligible'.
sig_thresh = 36*4*param['degrade']*param['degrade']*param['sig_thresh']
fmatrix = fmatrix.where(fmatrix > sig_thresh)
fmatrix = fmatrix.fillna(sig_thresh)
fmatrix = np.log10(fmatrix)

# Remove mean
fmatrix = fmatrix - fmatrix.mean(dim=time_var)
fmatrix.values = signal.detrend(fmatrix.values, axis=0, type='linear')

##############################################################################
# CARRY OUT FOURIER TRANSFORM                                                #
##############################################################################
omega = 2*np.pi

yf = rfft(fmatrix.values, axis=0)
xf = rfftfreq(len(fmatrix.coords[time_var]), 1/12)

seas_freq_idx = np.where(xf == 1)
phase = -np.angle(yf[seas_freq_idx, :, :])[0, 0, :, :]

ts = np.arange(0, len(fmatrix.coords[time_var])/12, 1/12) + 1/24
ts_arr = np.zeros_like(fmatrix.values)
phase_arr = np.zeros_like(fmatrix.values)

ts_arr[:] = ts[:, np.newaxis, np.newaxis]
phase_arr[:] = phase[np.newaxis, :, :]

seas_clim = np.cos((omega*ts_arr) - phase_arr)
seas_clim = xr.DataArray(seas_clim, coords={time_var: fmatrix.coords[time_var],
                                            'longitude': fmatrix.coords['longitude'],
                                            'latitude': fmatrix.coords['latitude']})

##############################################################################
# CALCULATE CORRELATION                                                      #
##############################################################################

corr = xs.pearson_r(fmatrix, seas_clim, dim=time_var)
corr_p_eff = xs.pearson_r_eff_p_value(fmatrix, seas_clim, dim=time_var)

##############################################################################
# PLOT                                                                       #
##############################################################################

f = plt.figure(figsize=(20, 11.9), constrained_layout=True)
gs = GridSpec(1, 2, figure=f, width_ratios=[0.99, 0.01])
ax = []
ax.append(f.add_subplot(gs[0, 0], projection = ccrs.PlateCarree())) # Main figure (eof)
ax.append(f.add_subplot(gs[0, 1])) # Colorbar for main figure

land_10m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='white',
                                        facecolor='black',
                                        zorder=1)
ax[0].set_aspect(1)

# Plot correlations
corr_plot = ax[0].contourf(fmatrix.coords['longitude'], fmatrix.coords['latitude'], corr.T,
                           cmap=cmr.fall_r, transform=ccrs.PlateCarree(), extend='neither',
                           levels=np.linspace(0, 1, 11), rasterized=True)
sig_plot1 = ax[0].contourf(fmatrix.coords['longitude'], fmatrix.coords['latitude'], corr_p_eff.where((corr_p_eff > 0.01)*(corr_p_eff <= 0.05)).T,
                            levels=np.array([0.01, 0.05]), hatches=['/'], colors='none')
sig_plot2 = ax[0].contourf(fmatrix.coords['longitude'], fmatrix.coords['latitude'], corr_p_eff.where(corr_p_eff <= 0.01).T,
                            levels=np.array([0, 0.01]), hatches=['.'], colors='none')

ax[0].add_feature(land_10m)
gl = ax[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.5, color='white', linestyle='--', zorder=11)
gl.xlocator = mticker.FixedLocator(np.arange(-210, 210, 10))
gl.ylocator = mticker.FixedLocator(np.arange(-90, 120, 10))
gl.xlabels_top = False
gl.ylabels_right = False
gl.ylabel_style = {'size': 24}
gl.xlabel_style = {'size': 24}

title = 'Seasonal correlation'
ax[0].text(22, -38, title, fontsize=32, color='k', fontweight='bold')
ax[0].set_xlim([20.5, 129.5])
ax[0].set_ylim([-39.5, 29.5])

if param['time'] == 'sink':
    ax[0].set_title(param['name'] + ' (Beaching time)', fontsize=32, color='k', fontweight='bold')
else:
    ax[0].set_title(param['name'] + ' (Source time)', fontsize=32, color='k', fontweight='bold')

cb0 = plt.colorbar(corr_plot, cax=ax[1], fraction=0.1)
ax[1].tick_params(axis='y', labelsize=24)
cb0.set_label('Correlation', size=28)
cb0.set_ticks(np.linspace(0, 1, 6))

# Save
plt.savefig(dirs['fig'] + 'Seasonal_correlation_' + param['time'] + '_' + param['physics'] + '_' + array_str + '_s' + str(param['us_d']) + '_b' + str(param['ub_d']) + '.pdf', bbox_inches='tight', dpi=300)

##############################################################################
# PLOT PHASE                                                                 #
##############################################################################

# # Convert phase to peak month (0 -> Jan, 1 -> Feb, -1 -> Dec...)
phase *= 6/np.pi
phase[phase < 0] = phase[phase < 0] + 12

# Plot

f = plt.figure(figsize=(21, 11.9), constrained_layout=True)
gs = GridSpec(1, 2, figure=f, width_ratios=[0.99, 0.01])
ax = []
ax.append(f.add_subplot(gs[0, 0], projection = ccrs.PlateCarree())) # Main figure (eof)
ax.append(f.add_subplot(gs[0, 1])) # Colorbar for main figure

land_10m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='white',
                                        facecolor='black',
                                        zorder=1)
ax[0].set_aspect(1)

# Plot phase offset
corr_plot = ax[0].pcolormesh(fmatrix.coords['longitude'], fmatrix.coords['latitude'], phase.T,
                             cmap=cmr.infinity_s, transform=ccrs.PlateCarree(), rasterized=True)
sig_plot1 = ax[0].contourf(fmatrix.coords['longitude'], fmatrix.coords['latitude'], corr_p_eff.where((corr_p_eff > 0.01)*(corr_p_eff <= 0.05)).T,
                            levels=np.array([0.01, 0.05]), hatches=['/'], colors='none')
sig_plot2 = ax[0].contourf(fmatrix.coords['longitude'], fmatrix.coords['latitude'], corr_p_eff.where(corr_p_eff <= 0.01).T,
                            levels=np.array([0, 0.01]), hatches=['.'], colors='none')

ax[0].add_feature(land_10m)
gl = ax[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.5, color='white', linestyle='--', zorder=11)
gl.xlocator = mticker.FixedLocator(np.arange(-210, 210, 10))
gl.ylocator = mticker.FixedLocator(np.arange(-90, 120, 10))
gl.xlabels_top = False
gl.ylabels_right = False
gl.ylabel_style = {'size': 24}
gl.xlabel_style = {'size': 24}

title = 'Seasonal cycle peak'
ax[0].text(22, -38, title, fontsize=32, color='k', fontweight='bold')
ax[0].set_xlim([20.5, 129.5])
ax[0].set_ylim([-39.5, 29.5])

if param['time'] == 'sink':
    ax[0].set_title(param['name'] + ' beaching time', fontsize=36, color='k', fontweight='bold')
else:
    ax[0].set_title(param['name'] + ' source time', fontsize=36, color='k', fontweight='bold')

cb0 = plt.colorbar(corr_plot, cax=ax[1], fraction=0.1)
ax[1].tick_params(axis='y', labelsize=24)
cb0.set_label('Peak month', size=28)
cb0.set_ticks(np.linspace(0.5, 11.5, 12))
cb0.set_ticklabels(['January', 'February', 'March', 'April', 'May', 'June',
                    'July', 'August', 'September', 'October', 'November', 'December'])
cb0.ax.invert_yaxis()
# Save
plt.savefig(dirs['fig'] + 'Seasonal_phase_' + param['time'] + '_' + param['physics'] + '_' + array_str + '_s' + str(param['us_d']) + '_b' + str(param['ub_d']) + '.pdf', bbox_inches='tight', dpi=300)
