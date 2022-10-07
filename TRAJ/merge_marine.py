#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to merge trajectory files
@author: Noam
"""

import xarray as xr
import numpy as np
import os
from glob import glob
from sys import argv
from tqdm import tqdm

##############################################################################
# DEFINE INPUT FILES                                                         #
##############################################################################

# PHYSICS
phys = argv[1]

# DIRECTORIES
dirs = {}
dirs['root'] = os.getcwd() + '/../'
dirs['traj_in'] = dirs['root'] + 'TRAJ/MAR/' + phys + '/'
dirs['traj_out'] = dirs['traj_in'] + 'POSTPROC/'

# PATTERN
pattern = 'FwdMar' + phys + '.nc' 
patern_out = 'Marine_' + phys + '_'
years = np.arange(1993, 2013)

# EVENTS
e_num = 25

##############################################################################
# LOOP THROUGH FILES AND MERGE                                               #
##############################################################################

for year in years:
    # Get list of all files from year
    fh_list_year = sorted(glob(dirs['traj_in'] + str(year) + '_*' + pattern))#               + pattern + str(year) + '_*'))
    months = np.unique([int(item.split('/')[-1].split('_')[-4]) for item in fh_list_year])

    pbar = tqdm(total=len(fh_list_year))

    for month in months:
        # Get list of all files from month
        fh_list_month = sorted(glob(dirs['traj_in'] + str(year) + '_' + str(month) + '_*' + pattern))
        releases = np.unique([int(item.split('/')[-1].split('_')[-3]) for item in fh_list_month])

        for release in releases:
            # Get list of all files from release
            fh_list_release = sorted(glob(dirs['traj_in'] + str(year) + '_' + str(month) + '_' + str(release) + '_*' + pattern))
            part_list = np.unique([int(item.split('/')[-1].split('_')[-2]) for item in fh_list_release])
            assert len(fh_list_release) == part_list[-1]+1 # Check all files are present

            for fhi, fh in enumerate(fh_list_release):
                try:
                    file = xr.open_dataset(fh, mask_and_scale=False)
                except:
                    print('Error opening file: ' + fh)
                    continue

                if fhi == 0:
                    var_dict = {}
                    attr_dict = {}

                if len(file['e_num']) == 0:
                    continue

                for v_name in ['e_num', 'lon0', 'lat0', 'gfw_num']:
                    if fhi == 0:
                        var_dict[v_name] = file[v_name].values
                    else:
                        var_dict[v_name] = np.concatenate((var_dict[v_name], file[v_name].values))

                for attr_name in ['parcels_version']:
                    if fhi == 0:
                        attr_dict[attr_name] = file.attrs[attr_name]
                    else:
                        assert attr_dict[attr_name] == file.attrs[attr_name]

                for var_type in ['e']:
                    for e_idx in range(e_num):
                        v_name = var_type + str(e_idx)

                        if fhi == 0:
                            var_dict[v_name] = file[v_name].values
                        else:
                            var_dict[v_name] = np.concatenate((var_dict[v_name], file[v_name].values))

                pbar.update(1)

            for v_name in var_dict.keys():
                var_dict[v_name] = ('traj', var_dict[v_name].flatten())

            new_file = xr.Dataset(data_vars=var_dict,
                                  coords={'traj': np.arange(len(var_dict['e_num'][1]))},
                                  attrs=attr_dict)

            new_fh = patern_out + str(year) + '_' + str(month) + '_' + str(release) + '.zarr'

            new_file.to_zarr(store=dirs['traj_out'] + new_fh, mode='w',
                             consolidated=True)
