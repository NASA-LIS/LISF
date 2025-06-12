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
#  7 Mar 2022: Sarith Mahanama, first version
# 22 Oct 2024: K. Arsenault, updated to account for srun submissions on discover
#
#------------------------------------------------------------------------------
"""

import glob
import os
import sys
import platform
import re
import datetime
import math
import numpy as np
from netCDF4 import Dataset as nc4 #pylint: disable=no-name-in-module
import yaml
#pylint: disable=consider-using-f-string, too-many-statements, too-many-locals, too-many-arguments

def job_script(s2s_configfile, jobfile, job_name, ntasks, hours, cwd, parallel_run, in_command = None,
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
                _f.write('#SBATCH --ntasks-per-node=' + str(ntasks) + '\n')
            else:
                _f.write('#SBATCH --ntasks-per-node=' + str(parallel_run['TPN']) + '\n')
        else:
            _f.write('#SBATCH --ntasks-per-node=' + str(ntasks) + '\n')
        _f.write('#SBATCH --time=' + hours + ':00:00' + '\n')
        if 'discover' in platform.node() or 'borg' in platform.node():
            _f.write('#SBATCH --constraint=' + cfg['SETUP']['CONSTRAINT'] + '\n')
            if 'mil' in cfg['SETUP']['CONSTRAINT']:
                _f.write('#SBATCH --partition=packable'  + '\n')
            if group_jobs:   
                mpc = min(math.ceil(480 / ntasks), 100)
                if parallel_run is not None:
                    _f.write('#SBATCH --mem-per-cpu=' + parallel_run['MEM'] + '\n')
                    _f.write('#SBATCH --cpus-per-task=' + parallel_run['CPT'] + '\n')
                else:
                    _f.write('#SBATCH --mem-per-cpu=' + str(mpc) + 'GB'  + '\n')
            else:
                _f.write('#SBATCH --mem-per-cpu=40GB'  + '\n')

        else:
            _f.write('#SBATCH --cluster-constraint=' + cfg['SETUP']['CONSTRAINT'] + '\n')
            _f.write('#SBATCH --partition=batch' + '\n')
            _f.write('#SBATCH --exclusive' + '\n')
            _f.write('#SBATCH --mem=0' + '\n')

        _f.write('#SBATCH --job-name=' + job_name + '\n')
        _f.write('#SBATCH --output ' + cwd + '/' + job_name + '%j.out' + '\n')
        _f.write('#SBATCH --error ' + cwd + '/' + job_name + '%j.err' + '\n')
        _f.write('\n')
        _f.write('#######################################################################' + '\n')
        _f.write('#                  Run LISF S2S ' + job_name + '\n')
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
        _f.write('export PYTHONPATH='+ pythonpath + '\n')
        if parallel_run is not None:
            _f.write('export NUM_WORKERS='+ parallel_run['CPT'] + '\n')
        _f.write('cd ' + cwd + '\n')

        if command_list is None and group_jobs is None:
            _f.write(f"{this_command} || exit 1\n")
            _f.write(f"{sec_command}\n")
        else:
            if group_jobs:
                for cmd in group_jobs:
                    if parallel_run is not None:
                        _f.write(f"srun --exclusive --cpus-per-task={parallel_run['CPT']} --ntasks {parallel_run['NT']} {cmd} &\n")
                    else:
                        _f.write(f"srun --exclusive --ntasks 1 {cmd} &\n")
                _f.write("wait\n")
            if command_list:
                for cmd in command_list:
                    _f.write(f"{cmd}\n")
        _f.write('\n')
        _f.write('echo "[INFO] Completed ' + job_name + '!"' + '\n')
        _f.write('\n')
        _f.write('/usr/bin/touch DONE' + '\n')
        _f.write('exit 0' + '\n')
    _f.close()

def cylc_job_scripts(job_file, hours, cwd, command_list=None, loop_list=None, command2=None):
    with open(job_file, 'w') as f:
        f.write("#!/bin/bash\n\n")
        # Set ITEMS to loop_list
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
        f.write("""/usr/bin/touch DONE
exit 0
        """)    

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
                    sbu = float (line.split(":")[1])
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
        _f.write('#SBATCH --output ' + cwd + '/' + job_name + '%j.out' + '\n')
        _f.write('#SBATCH --error ' + cwd + '/' + job_name + '%j.err' + '\n')
        _f.write('\n')
        _f.write('#######################################################################' + '\n')
        _f.write('#                  Run LISF S2S ' + job_name + '\n')
        _f.write('#######################################################################' + '\n')
        _f.write('\n')

        if 'discover' in platform.node() or 'borg' in platform.node():
            _f.write('source /etc/profile.d/modules.sh' + '\n')
            _f.write('module purge' + '\n')
        if os.path.isfile(lisf + '/env/discover/' + lisf_module):
# KRA TESTING:
#            _f.write('module use -a ' + lisf + '/env/discover/' + '\n')
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
        _f.write('/usr/bin/touch DONE' + '\n')
        _f.write('exit 0' + '\n')
    _f.close()
    cylc_job_scripts(job_name + 'run.sh', int(thours.split(':')[0]), cwd, command_list=[cylc_command])

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
    
