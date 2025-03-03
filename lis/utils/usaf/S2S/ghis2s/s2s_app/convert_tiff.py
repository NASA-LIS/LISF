import sys
import argparse
import xarray as xr
import numpy as np
import yaml
from shared import utils

def parse_list(arg):
    try:
        return [int(x) for x in arg.strip('[]').split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError('List must be integers.')
    
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', required=True, help='config file')
parser.add_argument('-i', '--input_tiff', required=True, help='Input tiff filename with full path')
parser.add_argument('-o', '--output_nc4', required=False, help='Output nc4 filename with full path')
parser.add_argument('-p', '--output_png', required=False, help='[Optional] Output png filename with full path')
parser.add_argument('-r', '--range', required=False, type=parse_list, help='[Optional] Provide the colorbar range as a two-element list of integers in quotes, for e.g. -r  "[-3, 3]"')

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
    data_xr.to_netcdf(nc4_file, format="NETCDF4",
                      encoding = {'data': {"zlib":True, "complevel":6, "shuffle":True, "missing_value": -9999., "_FillValue": -9999.}})

if args.output_png is None:
    sys.exit()

# plot
from s2s_modules.s2splots import plot_utils
domain =  [np.min(da.y.data), np.max(da.y.data), np.min(da.x.data), np.max(da.x.data)]

if args.range is None:
    plot_utils.contours(da.x.data, da.y.data, 1, 1, np.expand_dims(da.values, axis=0),'clim_reanaly',[' '], domain, args.output_png, ['white','white'],
                    fscale=0.75, cartopy_datadir = cfg['SETUP']['supplementarydir'] + '/s2splots/share/cartopy/')
else:
    plot_utils.contours(da.x.data, da.y.data, 1, 1, np.expand_dims(da.values, axis=0),'clim_reanaly',[' '], domain, args.output_png, ['navy','black'],
                        fscale=0.75, cartopy_datadir = cfg['SETUP']['supplementarydir'] + '/s2splots/share/cartopy/', min_val = args.range[0], max_val = args.range[1])
    
                     
