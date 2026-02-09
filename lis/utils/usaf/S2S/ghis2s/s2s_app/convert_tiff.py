#!/usr/bin/env python

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
    Converts TIF formatted files to other formats.
"""
import sys
import argparse
import xarray as xr
import numpy as np
import yaml
from ghis2s.shared import utils
from ghis2s.s2splots import plot_utils

def parse_list(arg):
    """
        Parse argument list.
    """
    try:
        return [int(x) for x in arg.strip('[]').split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError('List must be integers.')

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', required=True, help='config file')
parser.add_argument('-i', '--input_tiff', required=True, help='Input tiff filename with full path')
parser.add_argument('-o', '--output_nc4', required=False, help='Output nc4 filename with full path')
parser.add_argument('-p', '--output_png', required=False, help='[Optional] Output png filename with full path')
parser.add_argument(
    '-r', '--range',
    required=False,
    type=parse_list,
    help='[Optional] Colorbar range as a list of integers, e.g., -r "[-3, 3]"'
)

args = parser.parse_args()
tif_file = args.input_tiff
nc4_file = args.output_nc4

with open(args.config_file, 'r', encoding="utf-8") as file:
    cfg = yaml.safe_load(file)

da = utils.tiff_to_da(tif_file)

if args.output_nc4 is not None:
    data_xr = xr.Dataset()
    data_xr['data'] = (('y', 'x'), np.array(da.values))
    data_xr.coords['y'] = (('y'), da.y.data)
    data_xr.coords['x'] = (('x'), da.x.data)
    encoding = {
        'data': {
            'zlib': True,
            'complevel': 6,
            'shuffle': True,
            'missing_value': -9999.0,
            '_FillValue': -9999.0,
        }
    }
    data_xr.to_netcdf(nc4_file, format="NETCDF4", encoding=encoding)

if args.output_png is None:
    sys.exit()

# plot

domain =  [np.min(da.y.data), np.max(da.y.data), np.min(da.x.data), np.max(da.x.data)]

if args.range is None:
    cartopy_dir = f"{cfg['SETUP']['supplementarydir']}/s2splots/share/cartopy/"
    plot_utils.contours(
        da.x.data,
        da.y.data,
        1,
        1,
        np.expand_dims(da.values, axis=0),
        'clim_reanaly',
        [' '],
        domain,
        args.output_png,
        ['white', 'white'],
        fscale=0.75,
        cartopy_datadir=cartopy_dir
    )
else:
    cartopy_dir = f"{cfg['SETUP']['supplementarydir']}/s2splots/share/cartopy/"
    plot_utils.contours(
        da.x.data,
        da.y.data,
        1,
        1,
        np.expand_dims(da.values, axis=0),
        'clim_reanaly',
        [' '],
        domain,
        args.output_png,
        ['navy', 'black'],
        fscale=0.75,
        cartopy_datadir=cartopy_dir,
        min_val=args.range[0],
        max_val=args.range[1]
    )
