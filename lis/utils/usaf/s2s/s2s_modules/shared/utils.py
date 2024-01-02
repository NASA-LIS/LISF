#!/usr/bin/env python3

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
#------------------------------------------------------------------------------
#
# SCRIPT: job_script
#
# PURPOSE: writes batch job script
#
# REVISION HISTORY:
# 7 Mar 2022: Sarith Mahanama, first version
#
#------------------------------------------------------------------------------
"""


import glob
import os
import platform
import re
import datetime
import numpy as np
from netCDF4 import Dataset as nc4 #pylint: disable=no-name-in-module
import yaml
#pylint: disable=consider-using-f-string, too-many-statements, too-many-locals, too-many-arguments

def job_script(s2s_configfile, jobfile, job_name, ntasks, hours, cwd, in_command = None,
               command2 = None, command_list = None):
    ''' writes SLURM job script '''
    if in_command is None:
        this_command = 'COMMAND'
    else:
        this_command = in_command

    if command2 is None:
        sec_command = '\n'
    else:
        sec_command = command2

    with open(s2s_configfile, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    sponsor_code = cfg['SETUP']['SPCODE']
    lisf = cfg['SETUP']['LISFDIR']
    lisf_module = cfg['SETUP']['LISFMOD']
    supd = cfg['SETUP']['supplementarydir']

    with open(jobfile, 'w', encoding="utf-8") as _f:

        _f.write('#!/bin/bash' + '\n')
        _f.write('\n')
        _f.write('#######################################################################' + '\n')
        _f.write('#                        Batch Parameters ' + '\n')
        _f.write('#######################################################################' + '\n')
        _f.write('\n')
        _f.write('#SBATCH --account=' + sponsor_code + '\n')
        _f.write('#SBATCH --ntasks=' + ntasks + '\n')
        _f.write('#SBATCH --time=' + hours + ':00:00' + '\n')
        if 'discover' in platform.node() or 'borg' in platform.node():
            _f.write('#SBATCH --constraint=' + cfg['SETUP']['CONSTRAINT'] + '\n')
        else:
#            _f.write('#SBATCH --cluster-constraint=green' + '\n')
            _f.write('#SBATCH --cluster-constraint=' + cfg['SETUP']['CONSTRAINT'] + '\n')
            _f.write('#SBATCH --partition=batch' + '\n')
        _f.write('#SBATCH --job-name=' + job_name + '\n')
        _f.write('#SBATCH --output ' + cwd + '/' + job_name + '%j.out' + '\n')
        _f.write('#SBATCH --error ' + cwd + '/' + job_name + '%j.err' + '\n')
        _f.write('\n')
        _f.write('#######################################################################' + '\n')
        _f.write('#                  Run LIS-Hydro S2S ' + job_name + '\n')
        _f.write('#######################################################################' + '\n')
        _f.write('\n')
        if 'discover' in platform.node() or 'borg' in platform.node():
            _f.write('source /etc/profile.d/modules.sh' + '\n')
            _f.write('module purge' + '\n')
        if os.path.isfile(lisf + '/env/discover/' + lisf_module):
            _f.write('module use -a ' + lisf + '/env/discover/' + '\n')
            _f.write('module --ignore-cache load ' + lisf_module + '\n')
        else:
            _f.write('module use -a ' + supd + '/env/' + '\n')
            _f.write('module load ' + lisf_module + '\n')
        _f.write('ulimit -s unlimited' + '\n')
        _f.write('\n')
        _f.write('cd ' + cwd + '\n')

        if  command_list is None:
            _f.write( this_command + ' || exit 1' + '\n')
            _f.write( sec_command + '\n')
        else:
            for this_command in command_list:
                _f.write( this_command + '\n')
        _f.write('\n')
        _f.write('echo "[INFO] Completed ' + job_name + '!"' + '\n')
        _f.write('\n')
        _f.write('/usr/bin/touch DONE' + '\n')
        _f.write('exit 0' + '\n')
    _f.close()

def update_job_schedule (filename, myid, jobname, afterid):
    ''' writes the SLURM_JOB_SCHEDULE file '''
    with open(filename, "a", encoding="utf-8") as sch_file:
        sch_file.write('{:<10}{:<30}{}\n'.format(myid, jobname, afterid))

def print_status_report (e2es, yyyymm):
    ''' prints status report on the screen '''

    def read_out (this_no, nfiles, ofile, job_file):
        with open(ofile, "r", encoding="utf-8") as file:
            pattern_walltime = "Walltime Used"
            pattern_sbu =  "Estimated SBUs"
            pattern_not1 = "Pct Walltime Used"
            pattern_not2 = "Total CPU-Time Allocated"

            for line in file:
                if re.search(pattern_walltime, line):
                    if (not re.search(pattern_not1, line)) and (not re.search(pattern_not2, line)):
                        _l2 = [int(x) for x in line.split(":")[1:4]]
                if re.search(pattern_sbu, line):
                    sbu = np.float (line.split(":")[1])
        file.close()
        print ('{:>3}/{:>3}  {:<35}{:>2}h {:>2}m {:>2}s'.format
               (this_no,nfiles,job_file,_l2[0],_l2[1],_l2[2]))
        return sbu

    os.chdir(e2es)
    jfiles = glob.glob("scratch/" + yyyymm + "/*/*.j")
    jfiles.sort(key=os.path.getmtime)

    print ("  ")
    print ("#######################################################################")
    print ("                          STATUS OF SLURM JOBS                         ")
    print ("#######################################################################")
    print ("  ")
    print ("            JOB FILE                          WALLTIME ")
    print ("  ")
    total_sbu = 0.

    for file_no, jfile in enumerate(jfiles):
        job_file = jfile.split("/")[3]
        jcut = jfile[:-(len('run.j'))]
        ofile = glob.glob(jcut + "*.out")
        if len(ofile) == 1:
            if file_no == 0:
                efile = ofile[0][:-(len('out'))]+'err'
                time_begin = datetime.datetime.fromtimestamp(os.path.getctime(efile))
            sbu = read_out (file_no + 1, len(jfiles), ofile[0], job_file)
            total_sbu = total_sbu + sbu
            time_end = datetime.datetime.fromtimestamp(os.path.getmtime(ofile[0]))

    print ("  ")
    timedt = time_end - time_begin

    str1 = ' TOTAL SBUs        : {0:.2f}'.format (total_sbu)
    str2 = ' ELAPSED TIME      : {0:.2f} hours'.format (timedt.total_seconds() / 3600.)
    print (str1)
    print (str2)


def job_script_lis(s2s_configfile, jobfile, job_name, cwd, hours=None, in_command=None):
    ''' writes SLURM job scripts for LISF '''
    if in_command is None:
        this_command = 'COMMAND'
    else:
        this_command = in_command
    if hours is None:
        if 'discover' in platform.node() or 'borg' in platform.node():
            thours ='7:15:00'
        else:
            thours ='6:00:00'
    else:
        thours = hours + ':00:00'

    with open(s2s_configfile, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    sponsor_code = cfg['SETUP']['SPCODE']
    lisf = cfg['SETUP']['LISFDIR']
    lisf_module = cfg['SETUP']['LISFMOD']
    supd = cfg['SETUP']['supplementarydir']
    domain=cfg['EXP']['DOMAIN']
    datatype=cfg['SETUP']['DATATYPE']
    numprocx=cfg['FCST']['numprocx']
    numprocy=cfg['FCST']['numprocy']
    ntasks=str(numprocx*numprocy)

    with open(jobfile, 'w', encoding="utf-8") as _f:

        _f.write('#!/bin/bash' + '\n')
        _f.write('\n')
        _f.write('#######################################################################' + '\n')
        _f.write('#                        Batch Parameters ' + '\n')
        _f.write('#######################################################################' + '\n')
        _f.write('\n')
        _f.write('#SBATCH --account=' + sponsor_code + '\n')
        _f.write('#SBATCH --time=' + thours + '\n')
        if 'discover' in platform.node() or 'borg' in platform.node():
            _f.write('#SBATCH --constraint=' + cfg['SETUP']['CONSTRAINT'] + '\n')
        else:
#            _f.write('#SBATCH --cluster-constraint=green' + '\n')
            _f.write('#SBATCH --cluster-constraint=' + cfg['SETUP']['CONSTRAINT'] + '\n')
            _f.write('#SBATCH --partition=batch' + '\n')
        if datatype == 'hindcast':
            _f.write('#SBATCH --ntasks=' + ntasks + '\n')
        else:
            if domain == 'GLOBAL':
                _f.write('#SBATCH  -N 12' + '\n')
                _f.write('#SBATCH --ntasks-per-node=24' + '\n')
            else:
                _f.write('#SBATCH  -N 1' + '\n')
                _f.write('#SBATCH --ntasks-per-node='+ ntasks + '\n')

        _f.write('#SBATCH --job-name=' + job_name + '\n')
        _f.write('#SBATCH --output ' + cwd + '/' + job_name + '%j.out' + '\n')
        _f.write('#SBATCH --error ' + cwd + '/' + job_name + '%j.err' + '\n')
        _f.write('\n')
        _f.write('#######################################################################' + '\n')
        _f.write('#                  Run LIS-Hydro S2S ' + job_name + '\n')
        _f.write('#######################################################################' + '\n')
        _f.write('\n')
        if 'discover' in platform.node() or 'borg' in platform.node():
            _f.write('source /etc/profile.d/modules.sh' + '\n')
            _f.write('module purge' + '\n')
        if os.path.isfile(lisf + '/env/discover/' + lisf_module):
            _f.write('module use -a ' + lisf + '/env/discover/' + '\n')
            _f.write('module --ignore-cache load ' + lisf_module + '\n')
        else:
            _f.write('module use -a ' + supd + '/env/' + '\n')
            _f.write('module load ' + lisf_module + '\n')
        _f.write('ulimit -s unlimited' + '\n')
        _f.write('\n')
        _f.write('cd ' + cwd + '\n')
        _f.write( this_command + ' || exit 1' + '\n')
        _f.write('\n')
        _f.write('echo "[INFO] Completed ' + job_name + '!"' + '\n')
        _f.write('\n')
        _f.write('/usr/bin/touch DONE' + '\n')
        _f.write('exit 0' + '\n')
    _f.close()

def get_domain_info (s2s_configfile, extent=None, coord=None):
    ''' get domain infor from LDTINPUT file'''

    with open(s2s_configfile, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    ldtfile = cfg['SETUP']['supplementarydir'] + '/lis_darun/' + cfg['SETUP']['ldtinputfile']
    ldt = nc4(ldtfile, 'r')

    if extent is not None:
        lon = np.array(ldt['lon'])
        lat = np.array(ldt['lat'])
        return np.int(np.floor(np.min(lat[:,0]))), np.int(np.ceil(np.max(lat[:,0]))), \
            np.int(np.floor(np.min(lon[0,:]))), np.int(np.ceil(np.max(lon[0,:])))

    if coord is not None:
        lon = np.array(ldt['lon'])
        lat = np.array(ldt['lat'])
        return lat[:,0], lon[0,:]
    return None

def tiff_to_da(file):
    import xarray as xr
    import rasterio
    dataset = rasterio.open(file)
    # Read the data from the GeoTIFF using rasterio
    data = dataset.read(1)  # Read the first band, adjust if necessary

    # Extract the metadata
    transform = dataset.transform
    crs = dataset.crs
    x_coords = dataset.bounds.left + transform[0] * np.arange(dataset.width)
    y_coords = dataset.bounds.top + transform[4] * np.arange(dataset.height)
    
    # Create an xarray DataArray
    da = xr.DataArray(data, dims=('y', 'x'), coords={'y': y_coords, 'x': x_coords}, attrs={'crs': crs})
    
    return da

