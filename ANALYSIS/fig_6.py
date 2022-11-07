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
from sys import argv

# PARAMETERS
param = {# Analysis parameters
         'us_d': float(argv[1]),    # Sinking timescale (days)
         'ub_d': float(argv[2]),      # Beaching timescale (days)

         # Sink
         'sink': [['Aldabra', 'Assomption', 'Cosmoledo', 'Astove'],
                  ['Mahé', 'Fregate', 'Silhouette', 'Praslin', 'Denis', 'Bird']],
         'sink_name': ['Aldabra Group', 'Seychelles Plateau'],
         'linestyle': ['-', '--'],
         'linecolor': ['firebrick', 'darkblue'],
         'symbol'   : ['o', 's'],

         # Physics
         'mode': argv[3],

         # CMAP
         'cmap': cmr.guppy_r,
         'write_cmap': True, # Whether to write cmap data (good w/ 100/0010)
         'n_source': 10,

         # Export
         'name': argv[4],
         'export': True}

# DIRECTORIES
dirs = {'script': os.path.dirname(os.path.realpath(__file__)),
        'data': os.path.dirname(os.path.realpath(__file__)) + '/../MATRICES/',
        'plastic': os.path.dirname(os.path.realpath(__file__)) + '/../PLASTIC_DATA/',
        'fig': os.path.dirname(os.path.realpath(__file__)) + '/../FIGURES/'}

# FILE HANDLES
fh = {'flux': dirs['data'] + '/terrestrial_flux_' + param['mode'] + '_s' + str(param['us_d']) + '_b' + str(param['ub_d']) + '.nc',
      'source_list': dirs['plastic'] + 'country_list.in',
      'sink_list': dirs['plastic'] + 'sink_list.in',
      'cmap': dirs['fig'] + 'cmap_data.pkl',
      'fig': dirs['fig'] + 'time_series_' + param['mode'] + '_s' + str(param['us_d']) + '_b' + str(param['ub_d'])}

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
    ax.plot(t_axis, fmatrix_ts.loc[site, :].sum('sink'),
               param['linestyle'][site_i], c=param['linecolor'][site_i],
               label=param['sink_name'][site_i] + ' (' + param['name'] + ')')

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
abs_min = []

for site_i, site in enumerate(param['sink']):
    norm_ts = fmatrix_sinkt.loc[site, :].sum('sink')/fmatrix_sinkt.loc[site, :].sum('sink').mean()
    monclim = np.zeros((12,))
    monmin = np.zeros((12,))
    monmax = np.zeros((12,))
    monstd = np.zeros((12,))

    for m in range(12):
        monclim[m] = np.mean(norm_ts[m::12])
        monmin[m] = np.min(norm_ts[m::12])
        monmax[m] = np.max(norm_ts[m::12])
        monstd[m] = np.std(norm_ts[m::12])

    abs_max.append(np.max(monmax))
    abs_min.append(np.min(monmin))

    ax.fill_between(np.arange(1, 13), monmin, monmax,
                    color=param['linecolor'][site_i], alpha=0.1)
    ax.plot(np.arange(1, 13), monclim, marker=param['symbol'][site_i], ms=10,
            c=param['linecolor'][site_i], ls=param['linestyle'][site_i],
            label=param['sink_name'][site_i] + ' (' + param['name'] + ')')


abs_max_lim = 1.05*np.max(abs_max)
abs_min_lim = 0.95*np.min(abs_min)

ax.set_xlim([1, 12])
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
ax.set_yscale('log')
ax.set_xticks(np.arange(1, 13, 1))
# ax.set_yticks(np.arange(20))
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)
ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul',
                      'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax.set_ylim([abs_min_lim, abs_max_lim])
# ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig(fh['fig'] + '_monclim.pdf' , dpi=300, bbox_inches="tight")
