#!/usr/bin/env python3
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

import yaml
import argparse

def job_script(s2s_configfile, jobfile, job_name, ntasks, hours, cwd, in_command = None,  command2 = None, command_list = None):
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
    
    with open(jobfile, 'w') as f:
        
        f.write('#!/bin/bash' + '\n')
        f.write('\n')
        f.write('#######################################################################' + '\n')
        f.write('#                        Batch Parameters ' + '\n')
        f.write('#######################################################################' + '\n')
        f.write('\n')
        f.write('#SBATCH --account=' + sponsor_code + '\n')
        f.write('#SBATCH --ntasks=' + ntasks + '\n')
        f.write('#SBATCH --time=' + hours + ':00:00' + '\n')
        f.write('#SBATCH --constraint=cssrw' + '\n')
        f.write('#SBATCH --job-name=' + job_name + '\n')
        f.write('#SBATCH --output ' + cwd + '/' + job_name + '%j.out' + '\n')
        f.write('#SBATCH --error ' + cwd + '/' + job_name + '%j.err' + '\n')
        f.write('\n')
        f.write('#######################################################################' + '\n')
        f.write('#                  Run LIS-Hydro S2S ' + job_name + '\n')
        f.write('#######################################################################' + '\n')
        f.write('\n')
        f.write('source /etc/profile.d/modules.sh' + '\n')
        f.write('module purge' + '\n')
        f.write('module use -a ' + lisf + '/env/discover/' + '\n')
        f.write('module --ignore-cache load ' + lisf_module + '\n')
        f.write('ulimit -s unlimited' + '\n')
        f.write('\n')
        f.write('cd ' + cwd + '\n')

        if  command_list is None:
            f.write( this_command + ' || exit 1' + '\n')
            f.write( sec_command + '\n')
        else:
            for this_command in command_list:
                f.write( this_command + '\n')
        f.write('\n')
        f.write('echo "[INFO] Completed ' + job_name + '!"' + '\n')
        f.write('\n')
        f.write('/usr/bin/touch DONE' + '\n')
        f.write('exit 0' + '\n')
    f.close()

def job_script_long(s2s_configfile, jobfile, job_name, ntasks, cwd, in_command = None):
    if in_command is None:
        this_command = 'COMMAND'
    else:
        this_command = in_command
                
    with open(s2s_configfile, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    sponsor_code = cfg['SETUP']['SPCODE']
    lisf = cfg['SETUP']['LISFDIR']
    lisf_module = cfg['SETUP']['LISFMOD']
    
    with open(jobfile, 'w') as f:
        
        f.write('#!/bin/bash' + '\n')
        f.write('\n')
        f.write('#######################################################################' + '\n')
        f.write('#                        Batch Parameters ' + '\n')
        f.write('#######################################################################' + '\n')
        f.write('\n')
        f.write('#SBATCH --account=' + sponsor_code + '\n')
        f.write('#SBATCH --ntasks=' + ntasks + '\n')
        f.write('#SBATCH --time=23:59:59' + '\n')
        f.write('#SBATCH --qos=long' + '\n')
        f.write('#SBATCH --constraint=cssrw' + '\n')
        f.write('#SBATCH --job-name=' + job_name + '\n')
        f.write('#SBATCH --output ' + cwd + '/' + job_name + '%j.out' + '\n')
        f.write('#SBATCH --error ' + cwd + '/' + job_name + '%j.err' + '\n')
        f.write('\n')
        f.write('#######################################################################' + '\n')
        f.write('#                  Run LIS-Hydro S2S ' + job_name + '\n')
        f.write('#######################################################################' + '\n')
        f.write('\n')
        f.write('source /etc/profile.d/modules.sh' + '\n')
        f.write('module purge' + '\n')
        f.write('module use -a ' + lisf + '/env/discover/' + '\n')
        f.write('ulimit -s unlimited' + '\n')
        f.write('\n')
        f.write('cd ' + cwd + '\n')
        f.write( this_command + ' || exit 1' + '\n')
        f.write('\n')
        f.write('echo "[INFO] Completed ' + job_name + '!"' + '\n')
        f.write('\n')
        f.write('/usr/bin/touch DONE' + '\n')
        f.write('exit 0' + '\n')
    f.close()
    
def update_job_schedule (filename, myid, jobname, afterid):
    with open(filename, "a") as sch_file:
        sch_file.write('{:<10}{:<30}{}\n'.format(myid, jobname, afterid))

def print_status_report (e2es, yyyymm):
    import glob
    import os, time
    import re
    import sys
    import datetime
    import numpy as np

    def read_out (this_no, nfiles, ofile, job_file):
        file = open(ofile, "r")
        pattern_walltime = "Walltime Used"
        pattern_sbu =  "Estimated SBUs"
        pattern_not1 = "Pct Walltime Used"
        pattern_not2 = "Total CPU-Time Allocated"
        
        for line in file:
            if re.search(pattern_walltime, line):
                if (not re.search(pattern_not1, line)) and (not re.search(pattern_not2, line)):
                    l2 = [int(x) for x in line.split(":")[1:4]]
            if re.search(pattern_sbu, line):
                sbu = np.float (line.split(":")[1])
        file.close()
        print ('{:>3}/{:>3}  {:<35}{:>2}h {:>2}m {:>2}s'.format(this_no,  nfiles, job_file, l2[0], l2[1], l2[2]))
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
    if in_command is None:
        this_command = 'COMMAND'
    else:
        this_command = in_command
    if hours is None:
        thours ='12'
    else:
        thours = hours

    with open(s2s_configfile, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    sponsor_code = cfg['SETUP']['SPCODE']
    lisf = cfg['SETUP']['LISFDIR']
    lisf_module = cfg['SETUP']['LISFMOD']
    domain=cfg['EXP']['DOMAIN']
    DATATYPE=cfg['SETUP']['DATATYPE']
    numprocx=cfg['FCST']['numprocx']
    numprocy=cfg['FCST']['numprocy']
    ntasks=str(numprocx*numprocy)
    
    with open(jobfile, 'w') as f:
        
        f.write('#!/bin/bash' + '\n')
        f.write('\n')
        f.write('#######################################################################' + '\n')
        f.write('#                        Batch Parameters ' + '\n')
        f.write('#######################################################################' + '\n')
        f.write('\n')
        f.write('#SBATCH --account=' + sponsor_code + '\n')
        f.write('#SBATCH --constraint=cssrw' + '\n')        
        f.write('#SBATCH --time=' + thours + ':00:00' + '\n')
        if DATATYPE == 'hindcast':
            f.write('#SBATCH --ntasks=' + ntasks + '\n')
        else:
            if domain == 'GLOBAL':
                f.write('#SBATCH  -N 12' + '\n')
                f.write('#SBATCH --ntasks-per-node=24' + '\n')
            else:
                f.write('#SBATCH  -N 1' + '\n')
                f.write('#SBATCH --ntasks-per-node='+ ntasks + '\n')
                
        f.write('#SBATCH --job-name=' + job_name + '\n')
        f.write('#SBATCH --output ' + cwd + '/' + job_name + '%j.out' + '\n')
        f.write('#SBATCH --error ' + cwd + '/' + job_name + '%j.err' + '\n')
        f.write('\n')
        f.write('#######################################################################' + '\n')
        f.write('#                  Run LIS-Hydro S2S ' + job_name + '\n')
        f.write('#######################################################################' + '\n')
        f.write('\n')
        f.write('source /etc/profile.d/modules.sh' + '\n')
        f.write('module purge' + '\n')
        f.write('module use -a ' + lisf + '/env/discover/' + '\n')
        f.write('module --ignore-cache load ' + lisf_module + '\n')
        f.write('ulimit -s unlimited' + '\n')
        f.write('\n')
        f.write('cd ' + cwd + '\n')
        f.write( this_command + ' || exit 1' + '\n')
        f.write('\n')
        f.write('echo "[INFO] Completed ' + job_name + '!"' + '\n')
        f.write('\n')
        f.write('/usr/bin/touch DONE' + '\n')
        f.write('exit 0' + '\n')
    f.close()    

def get_domain_info (s2s_configfile, extent=None, coord=None):
    import numpy as np
    from netCDF4 import Dataset as nc4
    
    with open(s2s_configfile, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    ldtfile = cfg['BCSD']['supplementarydir'] + '/lis_darun/' + cfg['FCST']['ldtinputfile']
    ldt = nc4(ldtfile, 'r')

    if extent is not None:
        lon = np.array(ldt['lon'])
        lat = np.array(ldt['lat'])
        return np.int(np.floor(np.min(lat[:,0]))), np.int(np.ceil(np.max(lat[:,0]))), np.int(np.floor(np.min(lon[0,:]))), np.int(np.ceil(np.max(lon[0,:])))

    if coord is not None:
        lon = np.array(ldt['lon'])
        lat = np.array(ldt['lat'])
        return lat[:,0], lon[0,:]

    

    
    
    

    
    
    
