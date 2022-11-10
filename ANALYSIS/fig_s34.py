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
from sys import argv
from tqdm import tqdm

# PARAMETERS
param = {# Analysis parameters
         'us_d': float(argv[1]),    # Sinking timescale (days)
         'ub_d': float(argv[2]),      # Beaching timescale (days)

         # Sink
         'sink': argv[4],
         'linestyle': ['-', '--'],
         'linecolor': ['firebrick', 'darkblue'],
         'symbol'   : ['o', 's'],
         'source'   : ['Tanzania', 'Comoros', 'Seychelles'],

         # Physics
         'mode': argv[3],

         # CMAP
         'cmap': cmr.guppy_r,
         'write_cmap': True, # Whether to write cmap data (good w/ 100/0010)
         'n_source': 10,

         # Export
         'name': argv[5],
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

# Firstly compute the time-mean proportions of major sources
total_mass = fmatrix[:, :, :, 24:264].sum(dim=('source_time', 'sink_time', 'source')).loc[param['sink']]
total_prop = {source: fmatrix[:, :, :, 24:264].sum(dim=('source_time', 'sink_time')).loc[source, param['sink']]/total_mass
              for source in param['source']}

# Now compute the cumulative proportions for many start-times
n_years = 6
cumprop_ts = {}
for source in param['source']:
    cumprop_ts[source] = np.zeros((12*(20-n_years), 12*n_years), dtype=np.float32)

pbar = tqdm(total=len(param['source']*12*(20-n_years)))
for source in param['source']:
    for i in range(12*(20-n_years)):
        # Get the cumulative total and source flux
        cumtot_ts = fmatrix[:, :, :, 24+i:(24+n_years*12)+i].sum(dim=('source_time', 'source')).loc[param['sink']].cumsum(dim='sink_time')
        cumsource_ts = fmatrix[:, :, :, 24+i:(24+n_years*12)+i].sum(dim=('source_time')).loc[source, param['sink']].cumsum(dim='sink_time')
        cumprop_ts[source][i, :] = (cumsource_ts/cumtot_ts).values

        pbar.update(1)

# Now get percentiles
cumpctl_ts = {}

for source in param['source']:
    cumpctl_ts[source] = np.percentile(cumprop_ts[source], [5, 50, 95], axis=0)

# Plot results
f, ax = plt.subplots(1, 1, figsize=(12, 6))
colors = ['firebrick', 'k', 'darkblue']
t_axis = (np.arange(n_years*12)+1)/12

for i, source in enumerate(param['source']):
    ax.fill_between(t_axis, cumpctl_ts[source][0, :],
                    cumpctl_ts[source][2, :], color=colors[i], alpha=0.15)
    ax.plot(t_axis, cumpctl_ts[source][1, :], color=colors[i], label=source)
    ax.plot([t_axis[0], t_axis[-1]], [total_prop[source], total_prop[source]], color=colors[i],
            linestyle='--')

ax.set_xlim([t_axis[0], t_axis[-1]])
ax.set_ylim([0, 1])
ax.set_ylabel('Proportion of debris from source', fontsize=22)
ax.set_xlabel('Accumulation time (years)', fontsize=22)
legend = ax.legend(loc="upper right", fontsize=20)
legend.get_frame().set_linewidth(0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)

ax.set_title(param['name'] + ' debris beaching at ' + param['sink'], fontsize=26)

plt.savefig(fh['fig'] + param['sink'] + '_beaching_accum_' + param['name'] + '.pdf' , dpi=300, bbox_inches="tight")
