#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 15:34:50 2022
Estimate fisheries-related debris beaching rate at Seychelles taking into account input
distribution
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe
import cmasher as cmr
import xarray as xr
from matplotlib.ticker import FormatStrFormatter
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from tqdm import tqdm
import geopandas as gpd
import pandas as pd
from sys import argv


# PARAMETERS
param = {# Analysis parameters
         'us_d': ['30', '90', '360', '1800'],
         'ub_d': 30,

         # Physics
         'mode': ['CS0', 'CS1', 'CS3', 'CS5'],

         # Sink sites
         'sites': np.array([13,14,15,16,17,18]),
         'label': 'Seychelles_Plateau',
        }

# DIRECTORIES
dirs = {'script': os.path.dirname(os.path.realpath(__file__)),
        'data': os.path.dirname(os.path.realpath(__file__)) + '/../MATRICES/',
        'fisheries': os.path.dirname(os.path.realpath(__file__)) + '/../FISHERIES/DATA/PROC/',
        'shipping': os.path.dirname(os.path.realpath(__file__)) + '/../SHIPPING/',
        'dfad': os.path.dirname(os.path.realpath(__file__)) + '/../GDP_DATA/dFADs/',
        'fig': os.path.dirname(os.path.realpath(__file__)) + '/../FIGURES/'}

# FILE HANDLES
array_str = np.array2string(param['sites'], separator='-').translate({ord(i): None for i in '[]'})
fh = {'fisheries': dirs['fisheries'] + 'IOTC_monclim.nc',
      'dfad': dirs['dfad'] +'dfad_deploy.gpkg',}

###############################################################################
# LOAD DATA ###################################################################
###############################################################################
# Get grid
with xr.open_dataset(fh['fisheries']) as file:
    lon_bnd = file.coords['lon_bnd'].values
    lat_bnd = file.coords['lat_bnd'].values

# Grid dFADs
dfad_data = gpd.read_file(fh['dfad'])
dfad_data['pt_date'] = pd.to_datetime(dfad_data['pt_date'], format='%Y-%m-%dT%H:%M:%S')
dfad_data['source_month'] = dfad_data['pt_date'].dt.month
dfad_grid = np.zeros((len(lat_bnd)-1, len(lon_bnd)-1, 12, 12), dtype=np.float64)

# Grid dFAD deployments
for month in range(12):
    # i.e. broadcast source month to all sink months
    source_month_hist = np.histogram2d(dfad_data.loc[dfad_data['source_month'] == month+1]['geometry'].x.values,
                                       dfad_data.loc[dfad_data['source_month'] == month+1]['geometry'].y.values,
                                       bins=[lon_bnd, lat_bnd],
                                       density=True)[0].T
    dfad_grid[:, :, month, :] = np.repeat(source_month_hist[:, :, np.newaxis], 12, axis=2)

# Get grid
with xr.open_dataset(fh['fisheries']) as file:
    lon_bnd = file.coords['lon_bnd'].values
    lat_bnd = file.coords['lat_bnd'].values
    fisheries_data = {'Longline': np.repeat(np.transpose(file.ll_effort.values, axes=(1, 2, 0))[:, :, :, np.newaxis], 12, axis=3),
                      'Purse seine': np.repeat(np.transpose(file.ps_effort.values, axes=(1, 2, 0))[:, :, :, np.newaxis], 12, axis=3),
                      'dFAD': dfad_grid}

time_series = {'Longline': {'Class A': [], 'Class B': [], 'Class C': []},
               'Purse seine': {'Class A': [], 'Class B': [], 'Class C': []},
               'dFAD': []}

# PLOT
f = plt.figure(constrained_layout=True, figsize=(20, 20))
gs = GridSpec(3, 1, figure=f)
ax = []
ax.append(f.add_subplot(gs[0, 0]))
ax.append(f.add_subplot(gs[1, 0]))
ax.append(f.add_subplot(gs[2, 0]))
linewidth_list = [2, 2, 2, 2]
marker_list = ['o', 'v', '^', 's']
color_list = ['darkred', 'darkorange', 'teal', 'olivedrab']

# Now loop through trajectory files
for fishery_i, fishery in enumerate(['Longline', 'Purse seine', 'dFAD']):

    max_flux = 0
    ax[fishery_i].set_xlim([1-0.5, 12+0.5])
    ax[fishery_i].set_xlim([1-0.5, 12+0.5])

    if fishery == 'dFAD':
        fh_flux = dirs['data'] + '/marine_clim_flux_C0_' + array_str + '_s1800_b20.nc'

        with xr.open_dataarray(fh_flux) as file:
            raw_flux = file.coarsen({'latitude': 12, 'longitude': 12, 'source_month': 1, 'sink_month': 1}).sum()
            source_normalised_flux = fisheries_data[fishery]*raw_flux

        # Now average spatially and across source months
        fishery_time_series = source_normalised_flux.sum(dim=['latitude', 'longitude', 'source_month'])
        time_series['dFAD'] = fishery_time_series/np.sum(fishery_time_series)
        max_flux = np.max([np.max(time_series['dFAD']), max_flux])

        ax[fishery_i].plot(np.arange(1, 13), time_series['dFAD'], ls='-',
                           color='k', label='dFAD', marker='P',
                           linewidth=2, markersize=15)
    else:
        for class_i, debris_class in enumerate(['Class A', 'Class B', 'Class C', 'Class D']):
            fh_flux = dirs['data'] + '/marine_clim_flux_' + param['mode'][class_i] + '_' + array_str + '_s' + str(param['us_d'][class_i]) + '_b30.nc'

            with xr.open_dataarray(fh_flux) as file:
                raw_flux = file.coarsen({'latitude': 12, 'longitude': 12, 'source_month': 1, 'sink_month': 1}).sum()
                source_normalised_flux = fisheries_data[fishery]*raw_flux

            # Now average spatially and across source months
            fishery_time_series = source_normalised_flux.sum(dim=['latitude', 'longitude', 'source_month'])
            time_series[fishery][debris_class] = fishery_time_series/np.sum(fishery_time_series)
            max_flux = np.max([np.max(time_series[fishery][debris_class]), max_flux])

            ax[fishery_i].plot(np.arange(1, 13), time_series[fishery][debris_class], ls='-',
                               color=color_list[class_i],
                               label=debris_class, marker=marker_list[class_i],
                               linewidth=linewidth_list[class_i], markersize=15)

    # Plot
    max_flux *= 1.1

    ax[fishery_i].fill_between([0, 2.5], [max_flux, max_flux],
                               color='none', hatch='/', edgecolor='k', linewidth=0)
    ax[fishery_i].fill_between([11.5, 12.5], [max_flux, max_flux],
                               color='none', hatch='/', edgecolor='k', linewidth=0)
    ax[fishery_i].fill_between([5.5, 8.5], [max_flux, max_flux],
                               color='none', hatch='\\', edgecolor='k', linewidth=0)
    ax[fishery_i].spines['top'].set_visible(False)
    ax[fishery_i].spines['right'].set_visible(False)
    ax[fishery_i].set_yticks(np.arange(0,1,0.1))
    ax[fishery_i].set_ylim([0, max_flux])
    ax[fishery_i].tick_params(axis='y', labelsize=22)
    ax[fishery_i].text(0.6, max_flux*0.95, fishery, ha='left', va='top',
                       fontsize=36, path_effects=[pe.withStroke(linewidth=6, foreground='white')])

    if fishery_i == 0:
        # ax[fishery_i].set_ylabel('Normalised beaching rate', fontsize=26)
        legend = ax[fishery_i].legend(loc="lower right", fontsize=26, fancybox=False, framealpha=1)
        legend.get_frame().set_linewidth(0)

    if fishery_i == 2:
        ax[fishery_i].set_xticks(np.arange(1, 13, 1))
        ax[fishery_i].set_xlabel('Beaching month', fontsize=32)
        ax[fishery_i].tick_params(axis='x', labelsize=26)
        ax[fishery_i].set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul',
                                       'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
        # ax[fishery_i].set_ylabel('Normalised beaching rate', fontsize=26)
        legend = ax[fishery_i].legend(loc="lower right", fontsize=26, fancybox=False, framealpha=1)
        legend.get_frame().set_linewidth(0)
    else:
        ax[fishery_i].set_xticks([])

plt.savefig(dirs['fig'] + param['label'] + '_FisheriesRate.pdf' , dpi=300, bbox_inches="tight")
