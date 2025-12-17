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
# SCRIPT: utils.py
#
# PURPOSE: writes batch job script
#
# REVISION HISTORY:
#  7 Mar 2022: Sarith Mahanama, first version
# 22 Oct 2024: K. Arsenault, updated to account for srun submissions on discover
#
#------------------------------------------------------------------------------
"""

import glob
import os
import sys
import platform
import math
import shutil
import psutil
import rasterio
import numpy as np
import xarray as xr
import dask as da
from netCDF4 import Dataset as nc4 #pylint: disable=no-name-in-module
import yaml
#pylint: disable=consider-using-f-string, too-many-statements, too-many-locals, too-many-arguments

def job_script(s2s_configfile, jobfile, job_name, ntasks, hours, cwd,
               parallel_run, in_command = None,
               command2 = None, command_list = None, group_jobs=None):
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
    pythonpath = cfg['SETUP']['LISFDIR'] + 'lis/utils/usaf/S2S/'

    with open(jobfile, 'w', encoding="utf-8") as _f:

        _f.write('#!/bin/bash' + '\n')
        _f.write('\n')
        _f.write('#######################################################################' + '\n')
        _f.write('#                        Batch Parameters ' + '\n')
        _f.write('#######################################################################' + '\n')
        _f.write('\n')
        _f.write('#SBATCH --account=' + sponsor_code + '\n')
        _f.write('#SBATCH --nodes=1' + '\n')
        if parallel_run is not None:
            if parallel_run['TPN'] is None:
                if parallel_run['MP']:
                    _f.write('#SBATCH --ntasks=' + str(ntasks) + '\n')
                else:
                    _f.write('#SBATCH --ntasks-per-node=' + str(ntasks) + '\n')
            else:
                _f.write('#SBATCH --ntasks-per-node=' + str(parallel_run['TPN']) + '\n')
        else:
            if not cfg['SETUP']['CONSTRAINT'] == 'cssrw':
                _f.write('#SBATCH --ntasks-per-node=' + str(ntasks) + '\n')
        if len(hours) > 1:
            if hours in ['12','24']:
                _f.write('#SBATCH --time=' + f'{hours}:00:00' + '\n')
            else:
                _f.write('#SBATCH --time=' + f'00:{hours}:00' + '\n')
        else:
            _f.write('#SBATCH --time=' + hours + ':00:00' + '\n')
        if 'discover' in platform.node() or 'borg' in platform.node():
            _f.write('#SBATCH --constraint=' + cfg['SETUP']['CONSTRAINT'] + '\n')
            if group_jobs:
                mpc = min(math.ceil(240 / ntasks), 100)
                if parallel_run is not None:
                    if parallel_run['MP']:
                        _f.write('#SBATCH --mem=' + parallel_run['MEM'] + '\n')
                    else:
                        if 'mil' in cfg['SETUP']['CONSTRAINT']:
                            _f.write('#SBATCH --partition=packable'  + '\n')
                        mpc = str(math.ceil(240 / parallel_run['TPN'])) + 'GB'
                        _f.write('#SBATCH --mem-per-cpu=' + mpc + '\n')
                    _f.write('#SBATCH --cpus-per-task=' + parallel_run['CPT'] + '\n')
                else:
                    if 'cssrw' in cfg['SETUP']['CONSTRAINT']:
                        _f.write('#SBATCH --partition=datamove'  + '\n')
                    if 'mil' in cfg['SETUP']['CONSTRAINT']:
                        _f.write('#SBATCH --partition=packable'  + '\n')
                        _f.write('#SBATCH --mem-per-cpu=' + str(mpc) + 'GB'  + '\n')
            else:
                _f.write('#SBATCH --mem-per-cpu=40GB'  + '\n')
                if 'mil' in cfg['SETUP']['CONSTRAINT']:
                    _f.write('#SBATCH --partition=packable'  + '\n')

        else:
            _f.write('#SBATCH --cluster-constraint=' + cfg['SETUP']['CONSTRAINT'] + '\n')
            _f.write('#SBATCH --partition=batch' + '\n')
            _f.write('#SBATCH --exclusive' + '\n')
            _f.write('#SBATCH --mem=0' + '\n')

        _f.write('#SBATCH --job-name=' + job_name + '\n')
        _f.write('#SBATCH --output ' + cwd + '/logs/' + job_name + '%j.out' + '\n')
        _f.write('#SBATCH --error ' + cwd + '/logs/' + job_name + '%j.err' + '\n')
        _f.write('\n')
        _f.write('#######################################################################' + '\n')
        _f.write('#                  Run LISF S2S ' + job_name + '\n')
        _f.write('#######################################################################' + '\n')
        _f.write('\n')
        _f.write('export USE_CYLC_ENV=0' + '\n')
        if 'discover' in platform.node() or 'borg' in platform.node():
            _f.write('source /etc/profile.d/modules.sh' + '\n')
            _f.write('module purge' + '\n')
        if os.path.isfile(lisf + '/env/discover/' + lisf_module):
            _f.write('unset LD_LIBRARY_PATH' + '\n')
            _f.write('unset PROJ_DATA' + '\n')
            _f.write('unset PROJ_LIB' + '\n')
            _f.write('module use -a ' + lisf + '/env/discover/' + '\n')
            _f.write('module --ignore-cache load ' + lisf_module + '\n')
        else:
            _f.write('module use -a ' + supd + '/env/' + '\n')
            _f.write('module load ' + lisf_module + '\n')
        _f.write('ulimit -s unlimited' + '\n')
        _f.write('\n')
        _f.write('export PYTHONPATH='+ pythonpath + '\n')
        _f.write('export SCRIPT_NAME='+ job_name + 'run.j' + '\n')
        if parallel_run is not None:
            _f.write('export NUM_WORKERS='+ parallel_run['CPT'] + '\n')

        # To handle MPLCONFIG and Cache home dir write situation,
        #  where home or /tmp are read-only from compute nodes:
        if 'discover' not in platform.node() or 'borg' not in platform.node():
            if 'user_cache_dir' in cfg['SETUP']:
                # print(" -- Using a user-designated cache directory from the s2s_config file --")
                _f.write(f"export USER_CACHE_DIR=\"/{cfg['SETUP']['user_cache_dir']}/ghis2s_python_cache_$$\" \n")
                _f.write('export MPLCONFIGDIR="$USER_CACHE_DIR" \n')
                _f.write('export XDG_CACHE_HOME="$USER_CACHE_DIR" \n')

        _f.write('cd ' + cwd + '\n')
        _f.write("PIDS=()\n")
        if command_list is None and group_jobs is None:
            _f.write(f"{this_command} || exit 1\n")
            _f.write("PIDS+=($!)\n")
            _f.write("\n")
            _f.write(f"{sec_command}\n")
            _f.write("PIDS+=($!)\n")
            _f.write("\n")
        else:
            if group_jobs:
                for cmd in group_jobs:
                    if parallel_run is not None:
                        if 'SKIP_ARG' in parallel_run:
                            _f.write(f"{cmd} &\n")
                        else:
                            _f.write(f"srun --exclusive --cpus-per-task={parallel_run['CPT']}" +
                                     f" --ntasks {parallel_run['NT']} {cmd} &\n")
                        _f.write("PIDS+=($!)\n")
                        _f.write("\n")
                    else:
                        _f.write(f"srun --exclusive --ntasks 1 {cmd} &\n")
                        _f.write("PIDS+=($!)\n")
                        _f.write("\n")

            if command_list:
                for cmd in command_list:
                    _f.write(f"{cmd} \n")
                    _f.write("PIDS+=($!)\n")
                    _f.write("\n")
        _f.write(f"""for pid in "${{PIDS[@]}}"; do
    wait $pid || {{ echo "[ERROR] Process failed. Exiting."; touch logs/{job_name}FAILED; exit 1; }}
done
        """)
        _f.write('\n')
        _f.write('echo "[INFO] Completed ' + job_name + '!"' + '\n')
        _f.write('\n')
        _f.write('sleep 90' + '\n')
        _f.write('\n')
        _f.write('exit 0' + '\n')
    _f.close()

def remove_sbatch_lines(filename):
    ''' the function removes SBATCH lines from the input *.j file'''
    with open(filename, 'r', encoding="utf-8") as file:
        lines = file.readlines()

    filtered_lines = [line for line in lines if not line.strip().startswith('#SBATCH')]
    with open(filename, 'w', encoding="utf-8") as file:
        file.writelines(filtered_lines)

def cylc_job_scripts(job_file, hours, command_list=None, loop_list=None, command2=None):
    ''' writes Cylc specific .sh files without srun'''
    with open(job_file, 'w', encoding="utf-8") as f:
        f.write("#!/bin/bash\n\n")
        f.write("# Run tasks in parallel\n")
        f.write("PIDS=()\n")
        if command2 is not None:
            f.write(f"{command2} &\n")
            f.write("PIDS+=($!)\n")
            f.write("\n")

        if loop_list is None:
            # If loop_list is not provided, loop through command_list
            for cmd in command_list:
                f.write(f"{cmd} &\n")
                f.write("PIDS+=($!)\n")
                f.write("\n")
        else:
            quoted_items = [f"'{item}'" for item in loop_list]
            f.write(f"ITEMS=({' '.join(quoted_items)})\n")
            # First loop
            f.write("for ITEM in \"${ITEMS[@]}\"; do\n")
            f.write(f"    {command_list[0]} &\n")
            f.write("    PIDS+=($!)\n")
            f.write("done\n\n")

            if len(command_list) > 1:
                # Second loop (assuming MODELS is defined elsewhere in your script)
                f.write("for ITEM in \"${ITEMS[@]}\"; do\n")
                f.write(f"    {command_list[1]} &\n")
                f.write("    PIDS+=($!)\n")
                f.write("done\n\n")

        # Set runtime
        f.write("# Set runtime\n")
        f.write("START_TIME=$(date +%s)\n")
        f.write(f"TIME_LIMIT_SECONDS=$(({hours} * 60 * 60))  \n\n")

        # While loop for time limit and process checking
        f.write("""while true; do
sleep 60
CURRENT_TIME=$(date +%s)
ELAPSED_TIME=$((CURRENT_TIME - START_TIME))

if [ $ELAPSED_TIME -ge $TIME_LIMIT_SECONDS ]; then
    echo "[ERROR] Job exceeded time limit ($TIME_LIMIT). Killing processes..."
    for PID in "${PIDS[@]}"; do
        kill $PID 2>/dev/null
        sleep 2
        kill -9 $PID 2>/dev/null
    done
exit 1
fi

ALL_DONE=true
for PID in "${PIDS[@]}"; do
    if kill -0 $PID 2>/dev/null; then
        ALL_DONE=false
        break
    fi
done

if $ALL_DONE; then
    break
fi
done
        """)
        f.write(f"echo [INFO] Completed {job_file} ! \n\n")
        f.write("""
exit 0
        """)

def update_job_schedule (filename, myid, jobname, afterid):
    ''' writes the SLURM_JOB_SCHEDULE file '''
    with open(filename, "a", encoding="utf-8") as sch_file:
        sch_file.write('{:<10}{:<30}{}\n'.format(myid, jobname, afterid))

def job_script_lis(s2s_configfile, jobfile, job_name, cwd, hours=None, in_command=None):
    ''' writes SLURM job scripts for LISF '''
    if in_command is None:
        this_command = 'COMMAND'
    else:
        this_command = in_command
    if hours is None:
        if 'discover' in platform.node() or 'borg' in platform.node():
            thours ='7:00:00'
#            thours ='8:15:00'
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
    cylc_command = f'mpirun -np {ntasks} '
    with open(jobfile, 'w', encoding="utf-8") as _f:

        _f.write('#!/bin/bash' + '\n')
        _f.write('\n')
        _f.write('#######################################################################' + '\n')
        _f.write('#                        Batch Parameters ' + '\n')
        _f.write('#######################################################################' + '\n')
        _f.write('\n')
        _f.write('#SBATCH --account=' + sponsor_code + '\n')
        _f.write('#SBATCH --time=' + thours + '\n')

        # Slurm constaint entry:
        if 'discover' in platform.node() or 'borg' in platform.node():
            _f.write('#SBATCH --constraint=' + cfg['SETUP']['CONSTRAINT'] + '\n')
        else:
            _f.write('#SBATCH --cluster-constraint=' + cfg['SETUP']['CONSTRAINT'] + '\n')
            _f.write('#SBATCH --partition=batch' + '\n')
            _f.write('#SBATCH --exclusive' + '\n')
            _f.write('#SBATCH --mem=0' + '\n')

        # Assign ntasks based on hindcast | forecast LISF step:
        if datatype == 'hindcast':
            if 'mil' in cfg['SETUP']['CONSTRAINT']:
                _f.write('#SBATCH --ntasks=' + ntasks + ' --ntasks-per-socket=48 --ntasks-per-core=1' + '\n')
            else:
                _f.write('#SBATCH --ntasks=' + ntasks + '\n')
        # Forecast runmode (and LIS DA run setup):
        else:
            if domain == 'GLOBAL':
                if 'mil' in cfg['SETUP']['CONSTRAINT']:
                    _f.write('#SBATCH --ntasks=' + ntasks + ' --ntasks-per-socket=48 --ntasks-per-core=1' + '\n')
                    cylc_command = cylc_command + '--map-by socket:PE=48 --bind-to none '
                else:
                    _f.write('#SBATCH  -N 12' + '\n')
                    _f.write('#SBATCH --ntasks-per-node=24' + '\n')
                    cylc_command = cylc_command + '--map-by socket:PE=24 --bind-to none '
            else:
                _f.write('#SBATCH  -N 1' + '\n')
                _f.write('#SBATCH --ntasks-per-node='+ ntasks + '\n')

        _f.write('#SBATCH --job-name=' + job_name + '\n')
        _f.write('#SBATCH --output ' + cwd + '/logs/' + job_name + '%j.out' + '\n')
        _f.write('#SBATCH --error ' + cwd + '/logs/' + job_name + '%j.err' + '\n')
        _f.write('\n')
        _f.write('#######################################################################' + '\n')
        _f.write('#                  Run LISF S2S ' + job_name + '\n')
        _f.write('#######################################################################' + '\n')
        _f.write('\n')

        if 'discover' in platform.node() or 'borg' in platform.node():
            _f.write('source /etc/profile.d/modules.sh' + '\n')
            _f.write('module purge' + '\n')
        if os.path.isfile(lisf + '/env/discover/' + lisf_module):
            _f.write('unset LD_LIBRARY_PATH' + '\n')
            _f.write('module use --append ' + lisf + '/env/discover/' + '\n')
            _f.write('module --ignore-cache load ' + lisf_module + '\n')
        else:
            _f.write('module use -a ' + supd + '/env/' + '\n')
            _f.write('module load ' + lisf_module + '\n')
        _f.write('ulimit -s unlimited' + '\n')
        _f.write('\n')

        if 'mil' in cfg['SETUP']['CONSTRAINT']:
            _f.write('export I_MPI_PMI_LIBRARY=/usr/slurm/lib64/libpmi2.so' + '\n')
            _f.write('export I_MPI_PMI_VALUE_LENGTH_MAX=' + ntasks + '\n')
            _f.write('cd ' + cwd + '\n')
            _f.write('srun --mpi=pmi2 --ntasks=$SLURM_NTASKS \\' + '\n')
            _f.write('     --ntasks-per-socket=$SLURM_NTASKS_PER_SOCKET \\' + '\n')
            _f.write('     --ntasks-per-core=$SLURM_NTASKS_PER_CORE \\' + '\n')
            _f.write('     --cpu-bind="none"  \\' + '\n')
            # Separate out LIS DA run from LIS fcst run:
            if job_name == "lisda_":
                _f.write('     ' + this_command + ' || exit 1' + '\n')
                cylc_command = cylc_command + './LIS'
            else:
                _f.write('     ./LIS -f ' + this_command + ' || exit 1' + '\n')
                cylc_command = cylc_command + './LIS -f ' + this_command.split()[-1]
        else:
            _f.write('cd ' + cwd + '\n')
            _f.write( this_command + ' || exit 1' + '\n')

        _f.write('\n')
        _f.write('echo "[INFO] Completed ' + job_name + '!"' + '\n')
        _f.write('\n')
        _f.write('exit 0' + '\n')
    _f.close()
    shutil.copy(jobfile, job_name + 'run.sh')
    remove_sbatch_lines(job_name + 'run.sh')
    #cylc_job_scripts(job_name + 'run.sh', int(thours.split(':')[0]), command_list=[cylc_command])

def get_domain_info (s2s_configfile, extent=None, coord=None):
    ''' get domain infor from LDTINPUT file'''

    with open(s2s_configfile, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    ldtfile = cfg['SETUP']['supplementarydir'] + '/lis_darun/' + cfg['SETUP']['ldtinputfile']
    ldt = nc4(ldtfile, 'r')

    if extent is not None:
        lon = np.array(ldt['lon'])
        lat = np.array(ldt['lat'])
        return int(np.floor(np.min(lat[:,0]))), int(np.ceil(np.max(lat[:,0]))), \
            int(np.floor(np.min(lon[0,:]))), int(np.ceil(np.max(lon[0,:])))

    if coord is not None:
        lon = np.array(ldt['lon'])
        lat = np.array(ldt['lat'])
        return lat[:,0], lon[0,:]
    return None

def tiff_to_da(file):
    ''' converts TIF files to xarray DataArray'''
    dataset = rasterio.open(file)
    # Read the data from the GeoTIFF using rasterio
    data = dataset.read(1)  # Read the first band, adjust if necessary

    # Extract the metadata
    transform = dataset.transform
    crs = dataset.crs
    x_coords = dataset.bounds.left + transform[0] * np.arange(dataset.width)
    y_coords = dataset.bounds.top + transform[4] * np.arange(dataset.height)

    # Create an xarray DataArray
    tiff2da = xr.DataArray(data, dims=('y', 'x'), coords={'y': y_coords, 'x': x_coords}, attrs={'crs': crs})
    return tiff2da

def load_ncdata(infile, logger,  var_name=None, **kwargs):
    ''' generic function to load letcdf file[s] as a xarray dataset/datarray'''
    try:
        if isinstance(infile, str) and ('*' in infile or '?' in infile):
            matching_files = glob.glob(infile)
            if not matching_files:
                logger[0].error(f"No files found matching pattern: {infile}", subtask=logger[1])
                sys.exit(1)
            elif len(matching_files) == 1:
                infile = matching_files[0]
                multi_file_kwargs = ['combine', 'concat_dim', 'data_vars', 'coords', 'compat', 'join']
                kwargs = {k: v for k, v in kwargs.items() if k not in multi_file_kwargs}
            else:
                infile = matching_files
        if var_name is not None:
            if isinstance(infile, str):
                dataset = xr.open_dataset(infile, **kwargs)
            else:
                dataset = xr.open_mfdataset(infile, **kwargs)
            data = dataset[var_name]
            dataset.close()
            del dataset
            return data
        if isinstance(infile, str):
            return xr.open_dataset(infile, **kwargs)
        return xr.open_mfdataset(infile, **kwargs)

    except Exception as e:
        logger[0].error(f"Couldn't open {infile}", subtask=logger[1])
        logger[0].error(f"xarray error {e}", subtask=logger[1])
        sys.exit(1)

def detect_spatial_dimensions(out_xr):
    """Detect spatial dimension """
    lat_names = ['latitude', 'lat', 'lats', 'y', 'north_south']
    lon_names = ['longitude', 'lon', 'lons', 'x', 'east_west']

    lat_dim = None
    lon_dim = None

    # Find latitude/longiitude dimension names
    for dim in out_xr.dims:
        if dim.lower() in [name.lower() for name in lat_names]:
            lat_dim = dim
            break
    for dim in out_xr.dims:
        if dim.lower() in [name.lower() for name in lon_names]:
            lon_dim = dim
            break

    return lat_dim, lon_dim

def get_chunk_sizes(dataset, dim_in=None):
    ''' returns chuck sizes '''
    if dim_in is None:
        lat_dim, lon_dim = detect_spatial_dimensions(dataset)
        lat_size = dataset.sizes[lat_dim]
        lon_size = dataset.sizes[lon_dim]
    else:
        lat_size, lon_size = dim_in[0], dim_in[1]

    if lat_size == 3600:
        lat_chunk = 600
        lon_chunk = 1200
    elif lat_size == 1800:
        lat_chunk = 450
        lon_chunk = 900
    else:
        lat_chunk = max(100, lat_size // 4)
        lon_chunk = max(100, lon_size // 4)

    return lat_chunk, lon_chunk

def write_ncfile(out_xr, outfile, encoding, logger):
    ''' generic function to write netcdf4 files from xarray datasets'''
    try:
        # Get main data variable (exclude coordinates)
        data_vars = list(out_xr.data_vars.keys())
        if not data_vars:
            logger[0].error("No data variables found in dataset", subtask=logger[1])
            sys.exit(1)

        main_var = data_vars[0]
        is_lazy = hasattr(out_xr[main_var].data, 'chunks')

        if is_lazy:
            logger[0].info(f"Rechunking lazy data for optimal writing: {outfile}", subtask=logger[1])
            lat_dim, lon_dim = detect_spatial_dimensions(out_xr)

            if lat_dim and lon_dim:
                rechunk_dict = {}
                # Handle ensemble/time dimensions
                for dim in out_xr.dims:
                    if dim.lower() in ['ens', 'ensemble', 'member', 'time', 'lead']:
                        rechunk_dict[dim] = 1
                    elif dim == lat_dim:
                        rechunk_dict[dim] = out_xr.sizes[lat_dim]
                    elif dim == lon_dim:
                        rechunk_dict[dim] = out_xr.sizes[lon_dim]

                logger[0].info(f"Rechunking with: {rechunk_dict}", subtask=logger[1])
                updated_encoding = {}
                for var_name in data_vars:
                    if var_name in encoding:
                        updated_encoding[var_name] = encoding[var_name].copy()
                    else:
                        updated_encoding[var_name] = {
                            'dtype': 'float32',
                            'zlib': True,
                            'complevel': 6,
                            'shuffle': True,
                            '_FillValue': -9999.0
                        }

                    updated_encoding[var_name]['chunksizes'] = (
                        1,
                        out_xr.sizes[lat_dim],
                        out_xr.sizes[lon_dim]
                    )

                logger[0].info(f"Updated encoding for variables: {list(updated_encoding.keys())}", subtask=logger[1])

                with da.config.set({'array.chunk-size': '2GB', 'optimization.fuse': {}}):
                    out_xr_rechunked = out_xr.chunk(rechunk_dict)
                out_xr_rechunked.to_netcdf(outfile, format='NETCDF4', encoding=updated_encoding, engine='netcdf4')
            else:
                logger[0].info("Could not detect spatial dimensions, using standard write to NetCDF file", subtask=logger[1])
                out_xr.to_netcdf(outfile, format='NETCDF4', encoding=encoding, engine='netcdf4')

        else:
            logger[0].info("Pre-computed xarray dataset detected, writing directly to NetCDF file", subtask=logger[1])
            out_xr.to_netcdf(outfile, format='NETCDF4', encoding=encoding, engine='netcdf4')

        return
    except Exception as e:
        logger[0].error(f"Error saving file: {e}", subtask=logger[1])
        sys.exit(1)

def write_zarrfile(out_xr, outfile, encoding, logger):
    ''' generic function to write zarr files from xarray datasets'''
    try:
        # Parameters to drop (netCDF4-specific) zarr does not know
        netcdf_only_params = {'zlib', 'complevel', 'shuffle', 'fletcher32', 'contiguous'}

        # zarr encoding dictionary
        zarr_encoding = {}
        for var, enc in encoding.items():
            zarr_encoding[var] = {}
            for key, value in enc.items():
                if key not in netcdf_only_params:
                    if key == 'chunksizes':
                        zarr_encoding[var]['chunks'] = value
                    else:
                        zarr_encoding[var][key] = value
        out_xr.to_zarr(outfile, mode='w', encoding=zarr_encoding)
        return
    except Exception as e:
        logger[0].error(f"Error saving zarr file: {e}", subtask=logger[1])
        sys.exit(1)

def log_memory_usage(message, logger):
    """
    Log current program memory usage with a custom message.
    Args:
        message (str): Custom message describing where/when memory is being checked
        logger (list): List containing (logger_object, subtask_name)
    """
    try:
        process = psutil.Process()
        memory_mb = process.memory_info().rss / 1024 / 1024
        full_message = f"{message} - Memory usage: {memory_mb:.2f} MB"
        logger[0].info(full_message, subtask=logger[1])
    except Exception as e:
        logger[0].error(f"Error logging memory usage: {e}", subtask=logger[1])
