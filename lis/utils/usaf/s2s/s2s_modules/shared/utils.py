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
    
        
    with open(s2s_configfile, 'r') as file:
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
        f.write('#SBATCH --constraint=sky|cas' + '\n')
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
                
    with open(s2s_configfile, 'r') as file:
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
        f.write('#SBATCH --constraint=sky' + '\n')
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


