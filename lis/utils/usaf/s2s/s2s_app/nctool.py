#!/usr/bin/env python
'''
A python script to print the summary of a netCDF file OR compare netCDF files in all subdirectories (or just 2 individual netCDF files).
USAGE: (1) print Mean, Maximum (Maxloc) and Minimum (Minloc) of each spatial layer separately for each variable
           nctool.py file_name
       (2) compare 2 netcdf files
           a) check closeness at absolute tolerance of 1.e-3: nctool.py file1 file2
           b) check whether the 2 files are identical: nctool.py file1 file2 zero_diff
       (3) compare all netcdf files with similar names in two directories (including sub directories) 
           a) check closeness at absolute tolerance of 1.e-3: nctool.py directory1 directory2
           b) check whether the 2 are identical: nctool.py directory1 directory2 zero_diff
       (4) This is also good to print variable value[s] at a specific location. First ncdump -h and get exact names of lat, lon and the variable  
           as in the netCDF header. 
           nctool.py lat lon latname longname var_name filename
           for e.g.
           python nctool.py 39.875 25.875 lat lon Evap_tavg /gpfsm/dnb06/projects/p159/GHI_S2S/AFRICOM/E2ES/NRT/lis_darun/output/SURFACEMODEL/202204/LIS_HIST_202204010000.d01.nc
           python nctool.py -0.8486657635296301 32.97690772217679 latitude longitude anom /gpfsm/dnb06/projects/p159/GHI_S2S/AFRICOM/E2ES/NRT/s2smetric/output/202203/DYN_ANOM/AFRICOM/NOAHMP/GNEMO5_Surface-SM_ANOM_init_monthly_03_2022.nc  
           when j, i indexes at the location is known
           python nctool.py j_index i_index lat lon Evap_tavg /gpfsm/dnb06/projects/p159/GHI_S2S/AFRICOM/E2ES/NRT/lis_darun/output/SURFACEMODEL/202204/LIS_HIST_202204010000.d01.nc
'''

import os
import sys
import numpy as np
import xarray as xr
from numpy import unravel_index

zero_diff = False
ATOL = 1.e-03

def file_info(infile):
    
    def print_summary (arr):
        str_out =  'Mean :' + str(np.nanmean(arr)) + ' Max:' + str(np.nanmax(arr)) + str(unravel_index(arr.argmax(), arr.shape)) + \
            ' Min :' + str(np.nanmin(arr)) + str(unravel_index(arr.argmin(), arr.shape))
        return str_out
    
    _a = xr.open_dataset(infile)
    for var in _a.data_vars:
        darr = np.array(_a[var].values)
        print (_a[var].dims)
        if len(_a[var].dims) == 2:
            arr = darr
            print (var, ':=> ', print_summary (arr))
        elif len(_a[var].dims) == 3:
            for j in range (darr.shape[0]):
                arr = darr[j,:,:]
                print (var,j, ':=> ', print_summary (arr))
        elif len(_a[var].dims) == 4:
            for j in range (darr.shape[0]):
                for i in range (darr.shape[1]):                            
                    arr = darr[j,i,:,:]
                    print (var,j,i, ':=> ', print_summary (arr))
        elif len(_a[var].dims) == 5:
            for j in range (darr.shape[0]):
                for i in range (darr.shape[1]):
                    for n in range (darr.shape[2]):
                        arr = darr[j,i,n,:,:]
                        print (var,j,i,n, ':=> ', print_summary (arr))

    sys.exit()

def diff_nc(file1, file2):
    ''' Compare two netCDF files'''
    if file1.endswith('.nc') or file1.endswith('.NC') or file1.endswith('.nc4'):
        _a = xr.open_dataset(file1)
        _b = xr.open_dataset(file2)
        if zero_diff:
            xr.testing.assert_equal(_a, _b)
            print (file1,': is identical.')  
        else:
            xr.testing.assert_allclose(_a, _b, atol=1.e-03)
            print (file1,': is close at absolute tolerance of ' + str(ATOL) + '.') 

def print_var(latp, lonp, latname, lonname, varname, filename):

    def getclosest_ij(lats,lons,latpt,lonpt):
        # find squared distance of every point on grid
        dist_sq = (lats-latpt)**2 + (lons-lonpt)**2
        minindex_flattened = dist_sq.argmin()
        return np.unravel_index(minindex_flattened, lats.shape)
    
    _a = xr.open_dataset(filename)
    lats = np.array(_a[latname].values)
    lons = np.array(_a[lonname].values)
    darr = np.array(_a[varname].values)

    # find indeces closest to latp/lonp
    if latp.find('.') == -1:
        iy, ix = int(latp), int(lonp)
    else:
        if len(lats.shape) == 2:
            iy, ix = getclosest_ij(lats.astype(np.float), lons.astype(np.float), float(latp), float(lonp))
        if len(lats.shape) == 1:
            iy = min(range(len(lats)), key=lambda i: abs(float(lats[i]) - float(latp)))
            ix = min(range(len(lons)), key=lambda i: abs(float(lons[i]) - float(lonp)))

    if len (darr.shape) == 2:
        print (darr [iy,ix])
    if len (darr.shape) == 3:
        print (darr [:,iy,ix])
    if len (darr.shape) == 4:
        for i in range (darr.shape[1]):
            print (i, darr [:,i, iy,ix])
    
if __name__ == "__main__":
    if len(sys.argv) == 2:
        # print file infor
        if sys.argv[1] == 'help':
            print (' ')
            print ("\x1B[4m" + "\033[1m" + 'A python script to: (i) Print summary of a netCDF file; (ii) Compare 2 individual netCDF files; (iii) Compare netCDF files in two directories including all subdirectories; and (iv)  Print variable value[s] at a specific location.' + "\033[0m" + "\x1B[0m")
            print ( ' ')
            print ("\033[1m" +    '(1) print Mean, Maximum (Maxloc) and Minimum (Minloc) of each spatial layer separately for each variable:' +  "\033[0m")
            print ('            nctool.py file_name')
            print (' ')
            print ("\033[1m" +    '(2) Compare 2 netcdf files:' + "\033[0m")
            print ('         a) check closeness at absolute tolerance of 1.e-3:')
            print ('            nctool.py file1 file2')
            print ('         b) check whether the 2 files are identical:')
            print ('            nctool.py file1 file2 zero_diff')
            print (' ')            
            print ("\033[1m" +    '(3) Compare all netcdf files with similar names in two directories (including sub directories):' +  "\033[0m") 
            print ('         a) check closeness at absolute tolerance of 1.e-3:')
            print ('            nctool.py directory1 directory2')
            print ('         b) check whether the 2 are identical:')
            print ('            nctool.py directory1 directory2 zero_diff')
            print (' ')  
            print ("\033[1m" +    '(4) This is also good to print variable value[s] at a specific location specified by lat/lon OR J_INDEX [north-south dimension] and I_INDEX [west-east dimension] as in the netCDF file.' +  "\033[0m")  
            print ('         First ncdump -h and get exact names of lat, lon and the variable from the netCDF header.') 
            print ('         nctool.py lat lon latname longname var_name filename')
            print (' ')       
            print ('         Example 1 (when lat/lon at the place are known - the script will map the place on the grid space and indentify the grid cell that where the place is located.): ')
            print ('         nctool.py -0.8486657635296301 32.97690772217679 latitude longitude anom /gpfsm/dnb06/projects/p159/GHI_S2S/AFRICOM/E2ES/NRT/s2smetric/output/202203/DYN_ANOM/AFRICOM/NOAHMP/GNEMO5_Surface-SM_ANOM_init_monthly_03_2022.nc')
            print (' ')       
            print ('         Example 2 (when J_INDEX and I_INDEX at the location are known):')
            print ('         nctool.py j_index i_index lat lon Evap_tavg /gpfsm/dnb06/projects/p159/GHI_S2S/AFRICOM/E2ES/NRT/lis_darun/output/SURFACEMODEL/202204/LIS_HIST_202204010000.d01.nc')
        else:
            file_info(str(sys.argv[1]))
        
    if len(sys.argv) == 3 or len(sys.argv) == 4:
        # diff 2 dirs or 2 files
        CMDARGS = str(sys.argv)
        DIR1 = str(sys.argv[1])
        DIR2 = str(sys.argv[2])
        if len(sys.argv) == 4:
            zero_diff = True
            
        if os.path.isdir(DIR1):
            for root, dirs, files in os.walk(DIR1, topdown=False):
                for name in files:
                    _f1 = os.path.join(root, name)
                    _f2 = _f1.replace(DIR1, DIR2)
                    diff_nc(_f1, _f2)
        else:
            diff_nc(DIR1, DIR2)
        sys.exit()
        
    if len(sys.argv) == 7:
        # print var_values at lat/lon
        print_var(sys.argv[1], sys.argv[2], str(sys.argv[3]), str(sys.argv[4]), str(sys.argv[5]), str(sys.argv[6]))
        sys.exit()

    
