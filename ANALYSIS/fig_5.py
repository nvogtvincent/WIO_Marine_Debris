#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analyse time-series of marine debris accumulation at sites of interest
@author: noam
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import cmasher as cmr
import xarray as xr
from matplotlib.ticker import FormatStrFormatter

# PARAMETERS
param = {# Analysis parameters
         'us_d': 30.0,    # Sinking timescale (days)
         'ub_d': 30.0,      # Beaching timescale (days)

         # Sink
         'sink': ['Aldabra', 'Praslin'],
         'sink_name': ['Aldabra', 'Praslin'],
         'linestyle': ['-', '--'],

         # Physics
         'mode': '0000',

         # CMAP
         'cmap': cmr.guppy_r,
         'write_cmap': True, # Whether to write cmap data (good w/ 100/0010)
         'n_source': 10,

         # Export
         'name': 'Class A',
         'export': True}

# DIRECTORIES
dirs = {'script': os.path.dirname(os.path.realpath(__file__)),
        'data': os.path.dirname(os.path.realpath(__file__)) + '/../POSTPROC/',
        'plastic': os.path.dirname(os.path.realpath(__file__)) + '/../PLASTIC_DATA/',
        'fig': os.path.dirname(os.path.realpath(__file__)) + '/../FIGURES/'}

# FILE HANDLES
fh = {'flux': dirs['data'] + '/terrestrial_flux_' + param['mode'] + '_s' + str(param['us_d']) + '_b' + str(param['ub_d']) + '.nc',
      'source_list': dirs['plastic'] + 'country_list.in',
      'sink_list': dirs['plastic'] + 'sink_list.in',
      'cmap': dirs['fig'] + 'cmap_data.pkl',
      'fig': dirs['fig'] + 'time_series_AldPra_' + param['mode'] + '_s' + str(param['us_d']) + '_b' + str(param['ub_d'])}

##############################################################################
# EXTRACT DATA                                                               #
##############################################################################

fmatrix = xr.open_dataarray(fh['flux'])

##############################################################################
# PLOT SIMPLE TIME SERIES AND SEASONAL CYCLE                                 #
##############################################################################

f, ax = plt.subplots(1, 1, figsize=(20, 6))

# FIRSTLY GENERATE A TIME-SERIES
fmatrix_ts = fmatrix[:, :, :, 24:264].sum(dim=('source_time', 'source'))

# Generate a t-axis from the datetime format
y_min = fmatrix_ts.coords['sink_time'].values[0].astype('datetime64[Y]').astype(int) + 1970
n_year = int(len(fmatrix_ts.coords['sink_time'].values)/12)
t_axis = np.linspace(y_min + 1/24, y_min + n_year - 1/24, num=n_year*12)

ax.set_xlim([y_min, y_min + n_year])

for site_i, site in enumerate(param['sink']):
    ax.plot(t_axis, fmatrix_ts.loc[site, :],
               param['linestyle'][site_i], c='k', label=param['sink_name'][site_i] + ' (' + param['name'] + ')')

legend = ax.legend(loc="upper right", fontsize=20)
legend.get_frame().set_linewidth(0)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.set_ylabel('Monthly beaching debris (kg)', fontsize=22)
ax.set_xlabel('Beaching year', fontsize=22)
ax.set_xticks(np.arange(1997, 2015, 2))
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)
ax.grid(axis='x')
ax.set_yscale('log')
plt.savefig(fh['fig'] + '_full.pdf' , dpi=300, bbox_inches="tight")

# NOW CALCULATE A MONTHLY CLIMATOLOGY FOR  SINK TIME
f, ax = plt.subplots(1, 1, figsize=(20, 6))
fmatrix_sinkt = fmatrix[:, :, :, 24:264].sum(dim=('source_time', 'source'))
abs_max = []

for site_i, site in enumerate(param['sink']):
    norm_ts = fmatrix_sinkt.loc[site, :]/fmatrix_sinkt.loc[site, :].mean()
    monclim = np.zeros((12,))
    monmin = np.zeros((12,))
    monmax = np.zeros((12,))

    for m in range(12):
        monclim[m] = np.mean(norm_ts[m::12])
        monmin[m] = np.min(norm_ts[m::12])
        monmax[m] = np.max(norm_ts[m::12])

    abs_max.append(np.max(monclim))

    ax.plot(np.arange(1, 13), monclim, marker='s', c='k',
                   ls=param['linestyle'][site_i], label=param['sink_name'][site_i] + ' (' + param['name'] + ')')

abs_max_lim = 1.1*np.max(abs_max)

ax.set_xlim([1-0.5, 12+0.5])
legend = ax.legend(loc="upper right", fontsize=20)
legend.get_frame().set_linewidth(0)

ax.fill_between([0, 2.5], [abs_max_lim, abs_max_lim],
               color='none', hatch='/', edgecolor='k', linewidth=0)
ax.fill_between([11.5, 12.5], [abs_max_lim, abs_max_lim],
               color='none', hatch='/', edgecolor='k', linewidth=0)
ax.fill_between([5.5, 8.5], [abs_max_lim, abs_max_lim],
               color='none', hatch='\\', edgecolor='k', linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.set_ylabel('Normalised beaching rate', fontsize=22)
ax.set_xlabel('Beaching month', fontsize=22)
ax.set_xticks(np.arange(1, 13, 1))
ax.set_yticks(np.arange(20))
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)
ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul',
                      'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax.set_ylim([0, abs_max_lim])
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig(fh['fig'] + '_monclim.pdf' , dpi=300, bbox_inches="tight")
