import os
import sys
import glob
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import subprocess
import re
import threading
import time
import platform
import shutil
import tempfile
import yaml
import argparse
from ghis2s.s2s_app import s2s_api
from ghis2s.shared import utils
from ghis2s.lis_fcst import generate_lis_config_scriptfiles_fcst
from ghis2s.s2spost import run_s2spost_9months
from ghis2s.s2smetric import postprocess_nmme_job
from ghis2s import bcsd

class DownloadForecasts():
    def __init__(self, year, month, config_file):
        with open(config_file, 'r', encoding="utf-8") as file:
            self.config = yaml.safe_load(file)
        self.config_file = config_file
        self.E2ESDIR = self.config['SETUP']['E2ESDIR']
        self.year = year
        self.month = month
        YYYY = '{:04d}'.format(year)
        MM = '{:02d}'.format(month)
        self.YYYYMM = YYYY + MM
        self.SCRDIR = self.E2ESDIR + 'scratch/' + YYYY + MM + '/'
        if 'discover' in platform.node() or 'borg' in platform.node():
            self.cfsv2datadir = self.config['BCSD']['fcst_download_dir'] + "/Oper_TS/"
        else:
            self.cfsv2datadir = self.config['BCSD']['fcst_download_dir']
        self.patchfile = self.config['SETUP']['supplementarydir'] + "/bcsd_fcst/patch_files/patch_files_list.txt"
        self.patchdir = self.config['SETUP']['supplementarydir'] + "/bcsd_fcst/patch_files/"
        self.cfsv2_log = os.path.join(self.SCRDIR, 'CFSv2_missing_corrupted_files')
        self.srcdir = "https://noaacfs.blob.core.windows.net/cfs"
        self.NMME_RAWDIR =  self.config['BCSD']['nmme_download_dir']
        if os.path.exists(self.cfsv2_log):
            os.remove(self.cfsv2_log)

        def set_month_days(mon):
            if mon == 1:
                print("January ...")
                prevmon = 12
                day1, day2, day3 = 17, 22, 27
            elif mon == 2:
                print("February ...")
                prevmon = 1
                day1, day2, day3 = 21, 26, 31
            elif mon == 3:
                print("March ...")
                prevmon = 2
                day1, day2, day3 = 15, 20, 25
            elif mon == 4:
                print("April ...")
                prevmon = 3
                day1, day2, day3 = 17, 22, 27
            elif mon == 5:
                print("May ...")
                prevmon = 4
                day1, day2, day3 = 16, 21, 26
            elif mon == 6:
                print("June ...")
                prevmon = 5
                day1, day2, day3 = 21, 26, 31
            elif mon == 7:
                print("July ...")
                prevmon = 6
                day1, day2, day3 = 20, 25, 30
            elif mon == 8:
                print("August ...")
                prevmon = 7
                day1, day2, day3 = 20, 25, 30
            elif mon == 9:
                print("September ...")
                prevmon = 8
                day1, day2, day3 = 19, 24, 29
            elif mon == 10:
                print("October ...")
                prevmon = 9
                day1, day2, day3 = 18, 23, 28
            elif mon == 11:
                print("November ...")
                prevmon = 10
                day1, day2, day3 = 18, 23, 28
            elif mon == 12:
                print("December ...")
                prevmon = 11
                day1, day2, day3 = 17, 22, 27
            else:
                print("Invalid month")
                return None

            return '{:02d}'.format(prevmon), day1, day2, day3

        self.prevmon, self.day1, self.day2, self.day3 =  set_month_days(month)

    def NMME_file_checker(self):
        """ Checks NMME forecast files availability """
        Mmm = datetime(int(self.year), int(self.month), 1).strftime("%b")
        prec_files = glob.glob(f"{self.NMME_RAWDIR}/*/*{Mmm}.{self.year}.nc")
        print("[INFO] Current NMME files present: ")
        for file in prec_files:
            print(file)

        with open(self.E2ESDIR + self.config_file, 'r') as f:
            config_lines = f.readlines()
        
        nmme_models = {
            "CanSIPS-IC4": "GNEMO52, CanESM5",
            "COLA-RSMAS-CESM1": "CESM1",
            "GFDL-SPEAR": "GFDL",
            "NASA-GEOSS2S": "GEOSv2",
            "NCEP-CFSv2": "CFSv2"
        }

        have_model = "["

        if len(prec_files) < 5:
            print("[WARN] -- Some NMME model data may be missing or not available ... ")
            command = f"bash {self.E2ESDIR}/s2s_app/download_nmme_precip.sh {self.YYYY} {Mmm} {self.E2ESDIR}/${self.config_file}"
            process = subprocess.run(command, shell=True)
            
            prec_files = glob.glob(f"{self.NMME_RAWDIR}/*/*{Mmm}.{self.year}.nc")
            print(" ... New NMME files downloaded: ")
            for file in prec_files:
                print(file)

            for dir in ["CanSIPS-IC4", "COLA-RSMAS-CESM1", "GFDL-SPEAR", "NASA-GEOSS2S", "NCEP-CFSv2"]:
                fdown = f"{self.NMME_RAWDIR}/{dir}/prec.{dir}.mon_{Mmm}.{self.year}.nc"
                if os.path.isfile(fdown) and os.path.getsize(fdown) > 0:
                    have_model += f"{nmme_models[dir]}, "
                else:
                    if os.path.exists(fdown):
                        os.remove(fdown)

            have_model = have_model.rstrip(', ') + "]"

            if len(have_model.split()) < 7:
                navail = len(have_model.split()) - 1
                print(f"\nPrecipitation forecasts are available for only {navail} NMME models ({have_model}).")
                yesorno = input("Do you want to continue (Y/N)? ")

                if yesorno.lower() == 'y':
                    line2 = next(i for i, line in enumerate(config_lines) if 'NMME_models:' in line)
                    new_line = f"  NMME_models: {have_model}\n"
                    config_lines[line2] = new_line

                    with open(self.E2ESDIR + self.config_file, 'w') as f:
                        f.writelines(config_lines)
                else:
                    exit()

    def CFSv2_download(self):
        """ download CFSv2 forecasts """
        if self.month > 1:
            os.makedirs(self.cfsv2datadir  + '{:04d}'.format(self.year), exist_ok=True)
            os.chdir(self.cfsv2datadir  + '{:04d}'.format(self.year))
            year2 = self.year
        else:
            # - Need to account for Dec/Jan crossover
            os.makedirs(self.cfsv2datadir  + '{:04d}'.format(self.year -1), exist_ok=True)
            os.chdir(self.cfsv2datadir  + '{:04d}'.format(self.year -1))
            year2 = self.year -1

        for prevmondays in [self.day1, self.day2, self.day3]:
            icdate = f"{year2:04d}{self.prevmon}{prevmondays:02d}"
            os.makedirs(icdate, exist_ok=True)
            os.chdir(icdate)
            # Loop over variable type:
            for vartype in ['dlwsfc', 'dswsfc', 'q2m', 'wnd10m', 'prate', 'tmp2m', 'pressfc']:   
                # Forecast cycle (00, 06, 12, 18):
                for cycle in ['00', '06', '12', '18']:
                    file_name = f"{vartype}.01.{icdate}{cycle}.daily.grb2"
                    file = f"{self.srcdir}/cfs.{icdate}/{cycle}/time_grib_01/{vartype}.01.{icdate}{cycle}.daily.grb2"
                    if not os.path.isfile(file_name):
                        command = f"wget {file}"
                        subprocess.run(command)
        
    def CFSv2_file_checker(self):
        def create_cfsv2_log(cfsv2_log):
            # Write the header and instructions to the log file
            with open(cfsv2_log, 'a') as log_file:
                log_file.write(" #####################################################################################\n")
                log_file.write("                                  MISSING/INCOMPLETE CFSV2 FILES                      \n")
                log_file.write(" #####################################################################################\n")
                log_file.write("                         \n")
                log_file.write("  A replacement file is required for each missing or corrupted file. CFSv2 replacement files are saved in:\n")
                log_file.write(f"  {self.patchdir}, \n")
                log_file.write("  and comma-delimited lines in:  \n")
                log_file.write(f"  {self.patchfile} \n")
                log_file.write("  lists the replacement file names for each corrupted file. The table has three columns:  \n")
                log_file.write("  YYYYMMDDHH, bad_file_name, replacement_file_name.     \n")
                log_file.write("                         \n")
                log_file.write(f" (1) cd {self.patchdir}      \n")
                log_file.write(" (2) Each problematic file name in the section below is followed by a list of wget commands to download a suitable replacement file in order of preference.\n")
                log_file.write("     Download the first suggested replacement file and add a new entry to: \n")
                log_file.write(f"     {self.patchfile} \n")
                log_file.write(" (3) Repeat the same procedure to download replacements and update: \n")
                log_file.write(f"     {self.patchfile} \n")
                log_file.write("     for every missing/corrupted file.\n")
                log_file.write(" (4) Relaunch the forecast: s2s_app/s2s_run.sh -y YEAR -m MONTH -c CONFIGFILE\n")
                log_file.write(" (5) If any of the replacement files fail, you will be redirected to this file.\n")
                log_file.write(f"     {self.cfsv2_log}\n")
                log_file.write(" (6) Repeat steps (2) and (3) using a different replacement file for the original bad file.\n")
                log_file.write("                         \n")
                
        def neighb_days(vartype, icdate, cycle, mon):
            count_days = 1
            while count_days <= 4:
                # Calculate the previous date
                nwdate = (datetime.strptime(icdate, '%Y%m%d') - timedelta(days=count_days)).strftime('%Y%m%d')
                with open(self.cfsv2_log, 'a') as f:
                    f.write(f"wget {self.srcdir}/cfs.{nwdate}/{cycle}/time_grib_01/{vartype}.01.{nwdate}{cycle}.daily.grb2\n")
        
                # Calculate the next date
                nwdate = (datetime.strptime(icdate, '%Y%m%d') + timedelta(days=count_days)).strftime('%Y%m%d')
                if datetime.strptime(nwdate, '%Y%m%d').month != int(mon):
                    with open(self.cfsv2_log, 'a') as f:
                        f.write(f"wget {self.srcdir}/cfs.{nwdate}/{cycle}/time_grib_01/{vartype}.01.{nwdate}{cycle}.daily.grb2\n")
        
                count_days += 1
    
            with open(self.cfsv2_log, 'a') as f:
                f.write("   \n")
                
        def print_message(log_file):
            with open(log_file, 'a') as f:
                f.write("Note: If all recommended substitutes are also not available, you could try a different forecast hour from any of above dates.\n")
                f.write("\n")

        create_cfsv2_log(self.cfsv2_log)
        
        print(f"Previous mon, days 1-2-3 :: {self.prevmon}, {self.day1}-{self.day2}-{self.day3}")
        print(" ")
        print("==================================================================================================")
        print(" CFSv2 file checker is running to ensure all forcings files are available and not corrupted.......")
        print("==================================================================================================")
        
        if self.month > 1:
            os.chdir(self.cfsv2datadir  + '{:04d}'.format(self.year))
            year2 = self.year
        else:
            # - Need to account for Dec/Jan crossover
            os.chdir(self.cfsv2datadir  + '{:04d}'.format(self.year -1))
            year2 = self.year -1
            
        # Loop through variables and dates
        # Initial forecast dates:
        ret_code = 0
        for prevmondays in [self.day1, self.day2, self.day3]:
            icdate = f"{year2:04d}{self.prevmon}{prevmondays:02d}"
            os.chdir(icdate)

            # Loop over variable type:
            for vartype in ['dlwsfc', 'dswsfc', 'q2m', 'wnd10m', 'prate', 'tmp2m', 'pressfc']:   
                # Forecast cycle (00, 06, 12, 18):
                for cycle in ['00', '06', '12', '18']:
                    file_name = f"{vartype}.01.{icdate}{cycle}.daily.grb2"
                    # File check 1: missing file
                    if not os.path.isfile(file_name):
                        with open(patchfile, 'r') as f:
                            have_patch = f.read()
                            if file_name not in have_patch:
                                with open(self.cfsv2_log, 'a') as log_file:
                                    log_file.write(f"{file_name}:  MISSING\n")
                                    log_file.write("Possible substitutes in order of preference are:\n")
                                    neighb_days(vartype, icdate, cycle, prevmon)  
                                    ret_code = 1
                                    
                    # File check 2: corrupted file
                    if os.path.isfile(file_name):
                        py_code = s2s_api.cfsv2_file_checker(file_name, self.YYYYMM, py_call=True)
                        if py_code > 0:
                            with open(self.patchfile, 'r') as f:
                                have_patch = f.read()
                                
                            if file_name not in have_patch:
                                with open(self.cfsv2_log, 'a') as log_file:
                                    log_file.write(f"{file_name}: CORRUPTED\n")
                                    log_file.write("Possible substitutes in order of preference are:\n")
                                    neighb_days(vartype, icdate, cycle, prevmon)  
                                ret_code = 1
                                
                            else:
                                supfile = [line.split(',')[2].strip() for line in open(self.patchfile) if file_name in line][0]
                                supfile_path = os.path.join(self.patchdir, supfile)
                                py_code = s2s_api.cfsv2_file_checker(supfile_path, self.YYYYMM, py_call=True)
                        
                                if py_code > 0:
                                    with open(self.cfsv2_log, 'a') as log_file:
                                        log_file.write(f"{file_name}: Replacement {supfile} is also CORRUPTED!\n")
                                        log_file.write(f"Try downloading the next file (DON'T forget to update {patchfile})\n")
                                        neighb_days(vartype, icdate, cycle, prevmon) 
                                        ret_code = 1

            os.chdir('..')

        if ret_code > 0:
            print("*** Missing or Incomplete CFSv2 forcing files were found ***.")
            print("Please follow the instructions in:")
            print(self.cfsv2_log)
            print_message(self.cfsv2_log)
        else:
            with open(self.cfsv2_log, 'a') as log_file:
                log_file.write("**************************************************************\n")
                log_file.write(" SUCCESS ! All CFSv2 forcings files passed the file check.\n")
                log_file.write("**************************************************************\n")

                print("**************************************************************")
                print(" SUCCESS ! All CFSv2 forcings files passed the file check.")
                print("**************************************************************")

        print(" -- Done checking (and/or downloading) CFSv2 Forecast files -- ")

        return ret_code
        
class S2Srun(DownloadForecasts):
    def __init__(self, year, month, config_file):
        super(S2Srun, self).__init__(year, month, config_file)
        with open(config_file, 'r', encoding="utf-8") as file:
            self.config = yaml.safe_load(file)
        self.config_file = config_file
        self.year = year
        self.month = month
        self.YYYY = '{:04d}'.format(year)
        self.MM = '{:02d}'.format(month)
        self.E2ESDIR = self.config['SETUP']['E2ESDIR']
        self.LISFDIR = self.config['SETUP']['LISFDIR']
        self.LISHDIR = self.config['SETUP']['LISFDIR'] + 'lis/utils/usaf/S2S/'
        self.METFORC = self.config['SETUP']['METFORC']
        self.E2ESROOT = self.config['SETUP']['E2ESDIR']
        self.DOMAIN = self.config['EXP']['DOMAIN']
        self.SUPDIR = self.config['SETUP']['supplementarydir']
        self.LDTFILE = self.config['SETUP']['ldtinputfile']
        self.SCRDIR = self.E2ESDIR + 'scratch/' + self.YYYY + self.MM + '/'
        self.MODELS = self.config["EXP"]["NMME_models"]
        self.CONSTRAINT = self.config['SETUP']['CONSTRAINT']
        self.schedule = {}
        
        if not os.path.exists(self.E2ESDIR + 'scratch/'):
            subprocess.run(["setfacl", "-R", "-m", "u::rwx,g::rwx,o::r", self.E2ESDIR], check=True)

        os.makedirs(self.SCRDIR + '/lis_darun', exist_ok=True)
        os.makedirs(self.SCRDIR + '/ldt_ics', exist_ok=True)
        os.makedirs(self.SCRDIR + '/bcsd_fcst', exist_ok=True)
        os.makedirs(self.SCRDIR + '/lis_fcst', exist_ok=True)
        os.makedirs(self.SCRDIR + '/s2spost', exist_ok=True)
        os.makedirs(self.SCRDIR + '/s2smetric', exist_ok=True)
        os.makedirs(self.SCRDIR + '/s2splots', exist_ok=True)

    def create_symlink(self, source, link_name):
        """
        Create a symbolic link and ignore the error if the link already exists.
        :param source: The source file or directory.
        :param link_name: The destination for the symbolic link.
        """
        try:
            os.symlink(source, link_name)
        except FileExistsError:
            print(f"Symbolic link {link_name} already exists, skipping creation.")
        return
            
    def CFSv2_file_checker(self):
        spinner_done = [False]
        def spinner():
            """Display a spinner while waiting for a process to complete."""
            spin_chars = ['\\', '|', '/', '-']
            idx = 0
            while not spinner_done[0]:
                print(f'Please wait... {spin_chars[idx]}', end='\r', flush=True)
                idx = (idx + 1) % len(spin_chars)
                time.sleep(0.1)

        #command = f"bash {self.E2ESDIR}/s2s_app/wget_cfsv2_oper_ts_e2es.sh -y {self.YYYY} -m {self.MM} -c {self.E2ESDIR}/{self.config_file} -d N"
        #process = subprocess.run(command, shell=True)
        #ret_code = process.returncode

        spinner_thread = threading.Thread(target=spinner)
        spinner_thread.start()
        ret_code = super(S2Srun, self).CFSv2_file_checker()
        spinner_done[0] = True
        spinner_thread.join()        
        print('\rDone')
        
        if ret_code > 0:
            print(f"Error return code from the CFSv2 file download checker :: {ret_code}")
            print("> 0 :: Exiting from s2s_run.py --")
            sys.exit(ret_code)
 
        return
    
    def create_dict(self, jobfile, subdir, prev=None):
        self.schedule[jobfile] = {'subdir': subdir, 'jobid': None, 'prev': []}
        if prev is not None:
            if isinstance(prev, list):
                self.schedule[jobfile]['prev'].extend(prev)
            else:
                self.schedule[jobfile]['prev'].append(prev)

    def submit_jobs(self):
        def get_previds(jobfile):
            previd_list = []
            if len(self.schedule[jobfile]['prev']) > 0:
                for file in self.schedule[jobfile]['prev']:
                    previd_list.append(self.schedule[file]['jobid'])
            if len(previd_list) == 0:
                previd_list = None
            return previd_list
        
        def submit_slurm_job(job_script, prev_id=None):
            try:
                sbatch_command = ['/usr/bin/sbatch ']
                # Add dependency if prev_id is provided
                if prev_id is not None:
                    if isinstance(prev_id, list):
                        dependency_ids = ':'.join(map(str, prev_id))
                    else:
                        # If prev_id is a single value, convert it to string
                        dependency_ids = str(prev_id)
                    sbatch_command.extend(['--dependency=afterok:' + dependency_ids])
                    
                sbatch_command.append(' ' + job_script)
                sbatch_command = ''.join(sbatch_command)
                result = subprocess.run(sbatch_command, capture_output=True, text=True, check=True, shell=True)
                match = re.search(r'Submitted batch job (\d+)', result.stdout)
                if match:
                    job_id = match.group(1)
                    print(f"Submitted successfully. Job ID: {job_id}; Job Script: {job_script}")
                    return job_id
                else:
                    print("Failed to extract job ID from sbatch output")
                    return None
            except subprocess.CalledProcessError as e:
                print(f"An error occurred while submitting the job: {e}")
                return None
            
        JOB_SCHEDULE = os.path.join(self.SCRDIR, 'SLURM_JOB_SCHEDULE')
        if os.path.exists(JOB_SCHEDULE):
            os.remove(JOB_SCHEDULE)

        header = [
            "#######################################################################",
            "                         SLURM JOB SCHEDULE                            ",
            "#######################################################################",
            "                         "
        ]
        sub_section = {
            'lis_darun' : ["                         ",
                       "(1) LIS Data Assimilation",
                       "-------------------------",
                       "                         "],
            'ldt_ics' :["              ",
                        "(2) LDT and Initial Conditions",
                        "------------------------------",
                        "              "],

            'bcsd_fcst':["              ",
                         "(3) Bias Correction and Spatial Downscaling",
                         "-------------------------------------------",
                         "                                           "],
            'lis_fcst':["              ",
                        "(4) LIS Forecast Runs                      ",
                        "-------------------------------------------",
                        "                                           "],
            's2spost':["              ",
                       "(5) S2S post-process                       ",
                       "-------------------------------------------",
                       "                                           "],
            's2smetric':["              ",
                         "(6) S2S metric                             ",
                         "-------------------------------------------",
                         "                                           "],
            's2splots':["              ",
                       "(7) S2S plots                              ",
                       "-------------------------------------------",
                       "                                           "]}
        
        with open(JOB_SCHEDULE, "a", encoding="utf-8") as file:
            for line in header:
                file.write(line + "\n")
                
        utils.update_job_schedule(JOB_SCHEDULE, "JOB ID", "JOB SCRIPT", "AFTER")
        dir_list = []
        for jfile in self.schedule.keys():
            subdir = self.schedule[jfile]['subdir']
            if subdir not in dir_list:
                dir_list.append(subdir)
                with open(JOB_SCHEDULE, "a", encoding="utf-8") as file:
                    for line in sub_section.get(subdir):
                        file.write(line + "\n")
            os.chdir(self.SCRDIR + subdir)
            prev_ids = get_previds(jfile)
            job_id = submit_slurm_job(jfile, prev_id=prev_ids)
            self.schedule[jfile]['jobid'] = job_id
            if isinstance(prev_ids, list):
                prev_str = ','.join(map(str, prev_ids))
            else:
                # If prev_id is a single value, convert it to string
                prev_str = str(prev_ids)
            utils.update_job_schedule(JOB_SCHEDULE, job_id, jfile, prev_str)
                                                                            
    def split_list(self, input_list, length_sublist):
        result = []
        for i in range(0, len(input_list), length_sublist):
            sublist = input_list[i:i + length_sublist]
            result.append(sublist)
        return result

    def sublist_to_file(self, sublist, CWD):
        temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt')
        temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt', dir=CWD)
        for item in sublist:
            temp_file.write(f"{item}\n")
            temp_file.flush()
        return temp_file

    def write_cylc_snippet(self):
        """ writes Cylc runtime snippet """
        def extract_slurm_info(slurm_file):
            """ reads *.j SLURM file for directives, prescript and env variables"""
            def convert_time_to_minutes(time_str):
                """Convert time from HH:MM:SS to minutes."""
                parts = time_str.split(':')
                if len(parts) == 3:  
                    hours = int(parts[0])
                    minutes = int(parts[1])
                    return hours * 60 + minutes
                return None
            
            directives = []
            pre_script = []
            environment = []
            
            with open(slurm_file, 'r') as file:
                lines = file.readlines()       
                for line in lines:
                    line = line.strip()
                    # Extract SLURM directives
                    if line.startswith('#SBATCH'):
                        directive = line[8:].strip()
                        if directive.startswith('--time='):
                            # Convert time to minutes
                            time_value = directive.split('=', 1)[1].strip()
                            minutes = convert_time_to_minutes(time_value)
                            if minutes is not None:
                                directives.append(f'--time={minutes}') 
                        else:
                            directives.append(directive)
                            # Stop reading pre-script after 'cd' command
                    elif line.startswith('cd'):
                        break
                    elif line.startswith('#'):
                        continue
                    else:
                        pre_script.append(line)

            # Extract environment variables from the pre-script
            commands_to_remove = []
            for command in pre_script:
                if command.startswith('export'):
                    env_var = command.split('=', 1)
                    if len(env_var) == 2:
                        rm_export = env_var[0].split(' ',1)
                        if len(rm_export) > 1:
                            var_name = rm_export[1].strip()
                            var_value = env_var[1].strip()
                            environment.append(f"{var_name}={var_value}")
                            commands_to_remove.append(command)
                            
            # remove env var from pre script
            for cmd in commands_to_remove:
                if cmd in pre_script:
                    pre_script.remove(cmd)            

            # remove empty lines from pre_script
            pre_script = [cmd for cmd in pre_script if cmd.strip()]           
            return directives, pre_script, environment

        def write_lines(file, jfile, inherit_list, subdir, pre_script, directives, environment):
            file.write(f"    [[{jfile}]]\n")
            if inherit_list is not None:
                if isinstance(inherit_list, list):
                    inherits = ','.join(map(str, inherit_list))
                else:
                    inherits = str(inherit_list)
                file.write(f"        inherit = {inherits}\n")
            sh_script = self.SCRDIR + subdir + '/' + jfile + '.sh'
            file.write(f"        script = {sh_script}\n")
            
            # Write pre-script
            file.write("        pre-script = \n")
            for command in pre_script:
                file.write(f"                     {command}\n")
        
            # Write directives
            file.write("        [[[directives]]]\n")
            for directive in directives:
                file.write(f"            {directive}\n")

            if len(environment) > 0:
                # Write environment variables
                file.write("        [[[environment]]]\n")
                for env in environment:
                    file.write(f"            {env}\n")
        
        ''' the main function '''
        cylc_file = f"{self.SCRDIR}CYLC_workflow.rc"
        with open(cylc_file, 'w') as file:
            file.write("[runtime]\n")
            for jfile in self.schedule.keys():
                subdir = self.schedule[jfile]['subdir']
                directives, pre_script, environment = extract_slurm_info(self.SCRDIR + subdir + '/' + jfile)
                
                inherit_list = []
                if len(self.schedule[jfile]['prev']) > 0:
                    for pfile in self.schedule[jfile]['prev']:
                        inherit_list.append(pfile.removesuffix('.j'))
                if len(inherit_list) == 0:
                    inherit_list = None
                write_lines(file, jfile.removesuffix('.j'), inherit_list, subdir, pre_script, directives, environment)
                
    def lis_darun(self):
        """ LIS DARUN STEP """
        
        os.makedirs(self.E2ESDIR + '/lis_darun/input/', exist_ok=True)
        os.makedirs(self.E2ESDIR + '/lis_darun/output/lis.config_files/', exist_ok=True)
        os.chdir(self.E2ESDIR + '/lis_darun/input/')
        self.create_symlink(self.LISHDIR + '/ghis2s/lis_darun/forcing_variables.txt', 'forcing_variables.txt')
        self.create_symlink(self.LISHDIR + '/ghis2s/lis_darun/noahmp401_parms','noahmp401_parms')
        self.create_symlink(self.LISHDIR + '/ghis2s/lis_darun/template_files','template_files')
        self.create_symlink(self.LISHDIR + '/ghis2s/lis_darun/attribs','attribs')
        self.create_symlink(self.LISHDIR + '/ghis2s/lis_darun/tables','tables')
        self.create_symlink(self.SUPDIR + '/lis_darun/cdf/' +self.DOMAIN,'cdf')
        self.create_symlink(self.SUPDIR + '/lis_darun/RS_DATA','RS_DATA')
        self.create_symlink(self.SUPDIR + '/lis_darun/' + self.LDTFILE, self.LDTFILE)        
        os.chdir(self.E2ESDIR)

        # previous month
        date_obj = datetime.strptime(f"{self.YYYY}-{self.MM}-01", "%Y-%m-%d")
        previous_month_date = date_obj - relativedelta(months=1)
        YYYYMMP = previous_month_date.strftime("%Y%m")
        YYYYP = YYYYMMP[:4]
        MMP = YYYYMMP[4:6]
        monP = int(MMP)    
        PERTMODE = self.config['EXP']['pertmode']
        os.chdir(self.SCRDIR +'/lis_darun')
        
        CWD=self.SCRDIR +'lis_darun'
        self.create_symlink(self.LISFDIR + '/lis/LIS', 'LIS')
        self.create_symlink(self.E2ESDIR + '/lis_darun/input', 'input')
        self.create_symlink(self.E2ESDIR + '/lis_darun/output', 'output')
        self.create_symlink(self.METFORC, self.METFORC.rstrip('/').split('/')[-1])

        # write lisda_run.j
        # -----------------        
        s2s_api.lis_job_file(self.E2ESDIR +'/' + self.config_file, 'lisda_run.j', 'lisda_', CWD, str(5))

        if 'discover' in platform.node() or 'borg' in platform.node():
            if 'mil' in self.CONSTRAINT:
                COMMAND = './LIS'
            else:
                SLURM_NTASKS = os.getenv('SLURM_NTASKS', '1')  
                COMMAND = f'mpirun -np {SLURM_NTASKS} ./LIS'
        else:
            COMMAND = 'srun ./LIS'

        # add LIS command
        # ---------------
        with open('lisda_run.j', 'r') as file:
            filedata = file.read()   
        filedata = filedata.replace('COMMAND', COMMAND)
        with open('lisda_run.j', 'w') as file:
            file.write(filedata)

        self.create_dict('lisda_run.j', 'lis_darun')

        # configure lis.config
        # --------------------
        
        shutil.copy(self.E2ESDIR + '/lis_darun/input/template_files/lis.config_template.' + self.DOMAIN, 'lis.config')
        DAPERTRSTFILE = './output/DAPERT/' + YYYYP + MMP + '/LIS_DAPERT_{}{}010000.d01.bin'.format(YYYYP, MMP)
        NOAHMP401RSTFILE = './output/SURFACEMODEL/' + YYYYP + MMP + '/LIS_RST_NOAHMP401_{}{}010000.d01.nc'.format(YYYYP, MMP)
        HYMAP2RSTFILE = './output/ROUTING/' + YYYYP + MMP + '/LIS_RST_HYMAP2_router_{}{}010000.d01.nc'.format(YYYYP, MMP)
        LSMLISLOGFILE = CWD + '/logs_{}{}/lislog'.format(YYYYP, MMP)

        with open('lis.config', 'r') as file:
            filedata = file.read()

        filedata = filedata.replace('DAPERTRSTFILE', DAPERTRSTFILE)
        filedata = filedata.replace('NOAHMP401RSTFILE', NOAHMP401RSTFILE)
        filedata = filedata.replace('HYMAP2RSTFILE', HYMAP2RSTFILE)
        filedata = filedata.replace('STARTYR', YYYYP)
        filedata = filedata.replace('STARTMO', str(monP))
        filedata = filedata.replace('STARTDA', '1')
        filedata = filedata.replace('FINALYR', self.YYYY)
        filedata = filedata.replace('FINALMO', str(self.month))
        filedata = filedata.replace('FINALDA', '1')
        filedata = filedata.replace('PERTMODE', PERTMODE)
        filedata = filedata.replace('LSMLISLOGFILE', LSMLISLOGFILE)
        
        with open('lis.config', 'w') as file:
            file.write(filedata)

        shutil.copy('lis.config', 'output/lis.config_files/lis.config_darun_{}{}'.format(YYYYP, MMP))
        
        os.chdir(self.E2ESDIR)        
        return

    def ldt_ics(self):
        """ LDT-ICS STEP """
        if 'lisda_run.j' in self.schedule.keys():
            prev = 'lisda_run.j'
        else:
            prev = None
            
        os.makedirs(self.E2ESDIR + '/ldt_ics/input/', exist_ok=True)
        os.chdir(self.E2ESDIR + '/ldt_ics/')
        [os.makedirs(model, exist_ok=True) for model in self.MODELS]
        os.makedirs('ldt.config_files', exist_ok=True)
        os.makedirs('template_files', exist_ok=True)
        shutil.copy(self.LISHDIR +
                    '/ghis2s/ldt_ics/template_files/ldt.config_noahmp401_nmme_TEMPLATE.{}'.format(self.DOMAIN),
                    'template_files/ldt.config_noahmp401_nmme_TEMPLATE')

        os.chdir(self.E2ESDIR + '/ldt_ics/input/')
        self.create_symlink(self.SUPDIR + '/lis_darun/' + self.LDTFILE, self.LDTFILE)
        self.create_symlink(self.SUPDIR + '/LS_PARAMETERS', 'LS_PARAMETERS')
    
        os.chdir(self.SCRDIR + '/ldt_ics')
        CWD=self.SCRDIR + 'ldt_ics'
        self.create_symlink(self.LISFDIR + '/ldt/LDT', 'LDT')
        self.create_symlink(self.E2ESDIR + '/lis_darun/output', 'lisda_output')
        self.create_symlink(self.E2ESDIR + '/ldt_ics/input', 'input')

        for model in self.MODELS:
            self.create_symlink(self.E2ESDIR + '/ldt_ics/' + model, model)

        self.create_symlink(self.E2ESDIR + '/ldt_ics/ldt.config_files', 'ldt.config_files')
        self.create_symlink(self.E2ESDIR + '/ldt_ics/template_files', 'template_files')
    
        # configure batch script
        # ----------------------

        s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, 'ldtics_run.j', 'ldtics_', str(1), str(2), CWD, None) 
        
        COMMAND="python {}/ghis2s/ldt_ics/generate_ldtconfig_files_ensrst_nrt.py -y {} -m {} -i ./lisda_output -w {} -s {}/{}".\
            format(self.LISHDIR, self.YYYY, self.month, CWD, self.E2ESDIR, self.config_file)

        # add command
        # ---------------
        with open('ldtics_run.j', 'r') as file:
            filedata = file.read()   
        filedata = filedata.replace('COMMAND', COMMAND)
        with open('ldtics_run.j', 'w') as file:
            file.write(filedata)

        self.create_dict('ldtics_run.j', 'ldt_ics', prev=prev)

        os.chdir(self.E2ESDIR)        
        return

    def bcsd(self):
        """ BCSD 12 steps """            
        obs_clim_dir='{}/hindcast/bcsd_fcst/CFSv2_25km/raw/Climatology/'.format(self.E2ESROOT)
        nmme_clim_dir='{}/hindcast/bcsd_fcst/NMME/raw/Climatology/'.format(self.E2ESROOT)
        usaf_25km='{}/hindcast/bcsd_fcst/USAF-LIS7.3rc8_25km/raw/Climatology/'.format(self.E2ESROOT)

        os.makedirs(self.E2ESDIR + '/bcsd_fcst', exist_ok=True)
        os.chdir(self.E2ESDIR + '/bcsd_fcst')
        os.makedirs('USAF-LIS7.3rc8_25km/raw', exist_ok=True)
        os.makedirs('CFSv2_25km/raw', exist_ok=True)
        os.makedirs('NMME/raw', exist_ok=True)
    
        # link Climatology directories
        os.chdir(self.E2ESDIR + '/bcsd_fcst/CFSv2_25km/raw')
        self.create_symlink(obs_clim_dir, 'Climatology')

        os.chdir(self.E2ESDIR + '/bcsd_fcst/NMME/raw')
        self.create_symlink(nmme_clim_dir, 'Climatology')

        os.chdir(self.E2ESDIR + '/bcsd_fcst/USAF-LIS7.3rc8_25km/raw')
        self.create_symlink(usaf_25km, 'Climatology')
    
        # manage jobs from SCRATCH
        os.chdir(self.SCRDIR + 'bcsd_fcst')
        self.create_symlink(self.E2ESDIR + 'bcsd_fcst/', 'bcsd_fcst')
        CWD=self.SCRDIR + 'bcsd_fcst'

        date_obj = datetime.strptime(f"{self.YYYY}-{self.MM}-01", "%Y-%m-%d")
        mmm = date_obj.strftime("%b").lower()

        # (1) bcsd01 - regrid CFSv2 files
        # -------------------------------
        jobname='bcsd01_'
        slurm_commands = bcsd.task_01.main(self.E2ESDIR +'/' + self.config_file, self.year, None,
                                           mmm, CWD, jobname, 1, 2, py_call=True)

        # multi tasks per job
        l_sub = 2
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i in range(len(slurm_sub)):
            tfile = self.sublist_to_file(slurm_sub[i], CWD)
            try:
                s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + '{:02d}_run.j'.format(i+1),
                                    jobname+ '{:02d}_'.format(i+1), 1, str(3), CWD, tfile.name)
                self.create_dict(jobname+ '{:02d}_run.j'.format(i+1), 'bcsd_fcst')
            finally:
                tfile.close()
                os.unlink(tfile.name)
                
            utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 3, CWD, command_list=slurm_sub[i])
 
        # (3) bcsd03 regridding NMME
        # --------------------------
        jobname='bcsd03_'
        slurm_commands = bcsd.task_03.main(self.E2ESDIR +'/' + self.config_file, self.year, self.month,
                                           jobname, 1, str(2), CWD, py_call=True)
        tfile = self.sublist_to_file(slurm_commands, CWD)
        try:
            s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(3), CWD, tfile.name)
            self.create_dict(jobname+ 'run.j', 'bcsd_fcst')
        finally:
            tfile.close()
            os.unlink(tfile.name)
            
        utils.cylc_job_scripts(jobname + 'run.sh', 3, CWD, command_list=slurm_commands)

        # (4) bcsd04: Monthly "BC" step applied to CFSv2 (task_04.py, after 1 and 3)
        # --------------------------------------------------------------------------
        jobname='bcsd04_'
        prev = [f"{key}" for key in self.schedule.keys() if 'bcsd01_' in key]
        prev.extend([f"{key}" for key in self.schedule.keys() if 'bcsd03_' in key])
        slurm_commands = bcsd.task_04.main(self.E2ESDIR +'/' + self.config_file, self.year, self.year, mmm, self.MM, jobname,
                              1, 3, CWD, py_call=True)
        # multi tasks per job
        l_sub = 4
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i in range(len(slurm_sub)):
            tfile = self.sublist_to_file(slurm_sub[i], CWD)
            try:
                s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + '{:02d}_run.j'.format(i+1),
                                    jobname + '{:02d}_'.format(i+1), 1, str(4), CWD, tfile.name)
                self.create_dict(jobname + '{:02d}_run.j'.format(i+1), 'bcsd_fcst', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)
                
            utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 4, CWD, command_list=slurm_sub[i])

        # (5) bcsd05: Monthly "BC" step applied to NMME (task_05.py: after 1 and 3)
        # -------------------------------------------------------------------------
        jobname='bcsd05_'
        prev = [f"{key}" for key in self.schedule.keys() if 'bcsd01_' in key]
        prev.extend([f"{key}" for key in self.schedule.keys() if 'bcsd03_' in key])
        slurm_commands = []
        for nmme_model in self.MODELS:
            var1 = bcsd.task_05.main(self.E2ESDIR +'/' + self.config_file, self.year, self.year, mmm, self.MM, jobname,
                                     1, 3, CWD, nmme_model, py_call=True)
            slurm_commands.extend(var1)
        
        # multi tasks per job
        l_sub = 3
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i in range(len(slurm_sub)):
            tfile = self.sublist_to_file(slurm_sub[i], CWD)
            try:
                s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + '{:02d}_run.j'.format(i+1),
                                    jobname + '{:02d}_'.format(i+1), 1, str(4), CWD, tfile.name)
                self.create_dict(jobname + '{:02d}_run.j'.format(i+1), 'bcsd_fcst', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)
                
            utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 4, CWD, command_list=slurm_sub[i])

        # (6) bcsd06: CFSv2 Temporal Disaggregation (task_06.py: after 4 and 5)
        # ---------------------------------------------------------------------
        jobname='bcsd06_'
        prev = [f"{key}" for key in self.schedule.keys() if 'bcsd04_' in key]
        prev.extend([f"{key}" for key in self.schedule.keys() if 'bcsd05_' in key])         
        slurm_commands = bcsd.task_06.main(self.E2ESDIR +'/' + self.config_file, self.year, self.year, mmm, self.MM, jobname,
                                           1, 3, CWD, self.E2ESDIR, py_call=True)

        tfile = self.sublist_to_file(slurm_commands, CWD)
        try:
            s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(5), CWD, tfile.name)
            self.create_dict(jobname + 'run.j', 'bcsd_fcst', prev=prev)
        finally:
            tfile.close()
            os.unlink(tfile.name)
            
        utils.cylc_job_scripts(jobname + 'run.sh', 5, CWD, command_list=slurm_commands)

        # (8) bcsd08: NMME disaagregation
        # -------------------------------
        jobname='bcsd08_'
        slurm_commands = []
        for nmme_model in self.MODELS:
            var1 = bcsd.task_08.main(self.E2ESDIR +'/' + self.config_file, self.year, self.year, mmm, self.MM, jobname,
                                     1, 3, CWD, self.E2ESDIR, nmme_model, py_call=True)
            slurm_commands.extend(var1)
        
        tfile = self.sublist_to_file(slurm_commands, CWD)
        try:
            s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(5), CWD, tfile.name)
            self.create_dict(jobname + 'run.j', 'bcsd_fcst', prev='bcsd06_run.j')
        finally:
            tfile.close()
            os.unlink(tfile.name)

        utils.cylc_job_scripts(jobname + 'run.sh', 5, CWD, command_list=slurm_commands)

        # Task 9: Combine the CFSv2 forcing fields into final format for LIS to read
        # Task 10: Combine the NMME forcing fields into final format for LIS to read
        #          and symbolically link to the reusable CFSv2 met forcings
        # ---------------------------------------------------------------------------
        jobname='bcsd09-10_'
        slurm_9_10, slurm_11_12 = bcsd.task_09.main(self.E2ESDIR +'/' + self.config_file, self.year, self.year, mmm,
                                                    self.MM, jobname,1, 4, CWD, self.E2ESDIR, 'CFSv2', py_call=True)
        # bcsd09-10
        tfile = self.sublist_to_file(slurm_9_10, CWD)
        try:
            s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(6), CWD, tfile.name)
            self.create_dict(jobname + 'run.j', 'bcsd_fcst', prev='bcsd08_run.j')
        finally:
            tfile.close()
            os.unlink(tfile.name)
            
        utils.cylc_job_scripts(jobname + 'run.sh', 6, CWD, command_list=slurm_9_10)

        # bcsd11-12
        jobname='bcsd11-12_'
        tfile = self.sublist_to_file(slurm_11_12, CWD)
        try:
            s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(6), CWD, tfile.name)
            self.create_dict(jobname + 'run.j', 'bcsd_fcst', prev='bcsd09-10_run.j')
        finally:
            tfile.close()
            os.unlink(tfile.name)
            
        utils.cylc_job_scripts(jobname + 'run.sh', 6, CWD, command_list=slurm_11_12)
        
        os.chdir(self.E2ESDIR)
        return

    def lis_fcst(self):
        """ LIS forecast """
        prev = [job for job in ['ldtics_run.j', 'bcsd11-12_run.j'] if job in self.schedule] or None
        jobname='lis_fcst'
        os.makedirs(self.E2ESDIR + 'lis_fcst/', exist_ok=True)
        os.makedirs(self.E2ESDIR + 'lis_fcst/input/LDT_ICs/', exist_ok=True)

        os.chdir(self.E2ESDIR + 'lis_fcst/input/')
        self.create_symlink(self.LISHDIR + '/ghis2s/lis_darun/forcing_variables.txt','forcing_variables.txt')
        self.create_symlink(self.LISHDIR + '/ghis2s/lis_darun/noahmp401_parms','noahmp401_parms')
        self.create_symlink(self.LISHDIR + '/ghis2s/lis_fcst/template_files','template_files')
        self.create_symlink(self.LISHDIR + '/ghis2s/lis_fcst/tables','tables')
        self.create_symlink(self.SUPDIR + '/lis_darun/' + self.LDTFILE, self.LDTFILE)

        os.chdir(self.E2ESDIR + 'lis_fcst/input/LDT_ICs/')
        for model in self.MODELS:
            os.makedirs(self.E2ESDIR + 'lis_fcst/' + self.YYYY + self.MM + '/' + model + '/logs/', exist_ok=True)
            self.create_symlink(self.E2ESDIR + 'ldt_ics/'  + model, model)

        os.chdir(self.SCRDIR + 'lis_fcst/')
        CWD=self.SCRDIR + 'lis_fcst'
        self.create_symlink(self.LISFDIR + 'lis/LIS', 'LIS')
        self.create_symlink(self.E2ESDIR + 'lis_fcst/input', 'input')
        self.create_symlink(self.E2ESDIR + 'lis_fcst/' + self.YYYY + self.MM, self.YYYY + self.MM)
        self.create_symlink(self.E2ESDIR + 'bcsd_fcst', 'bcsd_fcst')

        generate_lis_config_scriptfiles_fcst.main(self.E2ESDIR + self.config_file, self.year, self.month, CWD, jobname)
        for model in self.MODELS:
            job_list = sorted(glob.glob(f"{jobname}_{model}*_run.j"))
            nFiles = len(job_list)
            for FileNo, jfile in enumerate(job_list):
                if nFiles > 1:
                    if FileNo == 0:
                        self.create_dict(jfile, 'lis_fcst', prev=prev)
                    else:
                        self.create_dict(jfile, 'lis_fcst', prev=job_list[FileNo-1])
                else:
                    self.create_dict(jfile, 'lis_fcst', prev=prev)
            
        os.chdir(self.E2ESDIR)
        return
        
    def s2spost(self):
        """ S2SPOST STEP """
        fcst_list = [item for item in self.schedule.keys() if re.match('lis_fcst', item)]
        if (len(fcst_list) > 0):
            prev=[]
            for model in self.MODELS:
                sublist = sorted([item for item in fcst_list if model in item])
                last_item = sublist[-1] if sublist else None
                if last_item is not None:
                    prev.append(last_item)
        else:
            prev = None
        
        [os.makedirs(self.E2ESDIR + 's2spost/' + self.YYYY + self.MM + '/' + model, exist_ok=True) for model in self.MODELS]
        os.chdir(self.SCRDIR + 's2spost')
        self.create_symlink(self.E2ESDIR + 'lis_fcst/' + self.YYYY + self.MM + '/', 'lis_fcst')
        self.create_symlink(self.E2ESDIR + 'lis_fcst/input/', 'input')
        CWD=self.SCRDIR + 's2spost'

        jobname='s2spost_'        
        slurm_commands = []
        for model in self.MODELS:
            self.create_symlink(self.E2ESDIR + 's2spost/' + self.YYYY + self.MM + '/' + model, model)
            var1 = run_s2spost_9months.main(self.E2ESDIR +'/' + self.config_file, self.year, self.month, jobname, 1, str(3), CWD, model, py_call=True)
            slurm_commands.extend(var1)
            
        # multi tasks per job
        l_sub = 27
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i in range(len(slurm_sub)):
            tfile = self.sublist_to_file(slurm_sub[i], CWD)
            try:
                s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + '{:02d}_run.j'.format(i+1),
                                        jobname + '{:02d}_'.format(i+1), 1, str(4), CWD, tfile.name)
                self.create_dict(jobname + '{:02d}_run.j'.format(i+1), 's2spost', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)
                
            utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 4, CWD, command_list=slurm_sub[i])

        os.chdir(self.E2ESDIR)        
        return

    def s2smetric(self):
        """ S2SMETRICS STEP """
        s2spost_list = [item for item in self.schedule.keys() if re.match('s2spost', item)]
        if len(s2spost_list) > 0:
            prev = s2spost_list
        else:
            prev = None
            
        os.makedirs(self.E2ESDIR + 's2smetric/' + self.YYYY + self.MM + '/DYN_ANOM', exist_ok=True)
        os.makedirs(self.E2ESDIR + 's2smetric/' + self.YYYY + self.MM + '/DYN_SANOM', exist_ok=True)
        os.makedirs(self.E2ESDIR + 's2smetric/' + self.YYYY + self.MM + '/metrics_cf', exist_ok=True)

        os.chdir(self.SCRDIR + 's2smetric')
        self.create_symlink(self.E2ESDIR + 's2spost/', 's2spost')
        self.create_symlink(self.E2ESDIR + 's2smetric/', 's2smetric')
        CWD=self.SCRDIR + 's2smetric'

        jobname='s2smetric_'
        slurm_commands = []
        for model in self.MODELS:
            var1 = postprocess_nmme_job.main(self.E2ESDIR +'/' + self.config_file, self.year, self.month, CWD, jobname=jobname, ntasks=1,
                                             hours=str(4), nmme_model= model, py_call=True)
            slurm_commands.extend(var1)

        tfile = self.sublist_to_file(slurm_commands, CWD)
        try:
            s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(4), CWD, tfile.name)
            self.create_dict(jobname + 'run.j', 's2smetric', prev=prev)
        finally:
            tfile.close()
            os.unlink(tfile.name)

        utils.cylc_job_scripts(jobname + 'run.sh', 4, CWD, command_list=slurm_commands)

        # weekly metrics
        jobname='s2smetric_weekly_'
        slurm_commands = []
        weekly_vars = self.config["POST"]["weekly_vars"]
        par_info = {}
        par_info['NPROCS'] = str(len(weekly_vars))
        par_info['MEM']= '8GB'
        for model in self.MODELS:
            var1 = postprocess_nmme_job.main(self.E2ESDIR +'/' + self.config_file, self.year, self.month, CWD, jobname=jobname, ntasks=1,
                                             hours=str(4), nmme_model= model, py_call=True, weekly=True)
            slurm_commands.extend(var1)

        if 'discover' in platform.node() or 'borg' in platform.node():
            tfile = self.sublist_to_file(slurm_commands, CWD)
            try:
                s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + 'run.j',
                                        jobname, 1, str(4), CWD, tfile.name, parallel_run=par_info)
                self.create_dict(jobname + 'run.j', 's2smetric', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)

            utils.cylc_job_scripts(jobname + 'run.sh', 4, CWD, command_list=slurm_commands)
        else:
            # multi tasks per job
            l_sub = 6
            slurm_sub = self.split_list(slurm_commands, l_sub)
            for i in range(len(slurm_sub)):
                tfile = self.sublist_to_file(slurm_sub[i], CWD)
                try:
                    s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + '{:02d}_run.j'.format(i+1),
                                            jobname + '{:02d}_'.format(i+1), 1, str(4), CWD, tfile.name)
                    self.create_dict(jobname + '{:02d}_run.j'.format(i+1), 's2smetric', prev=prev)
                finally:
                    tfile.close()
                    os.unlink(tfile.name)
                
                utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 4, CWD, command_list=slurm_sub[i])

        jobname='s2smetric_tiff_'
        prev = ['s2smetric_run.j']
        prev.extend([f"{key}" for key in self.schedule.keys() if 's2smetric_weekly_' in key])
        
        s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, 's2smetric_tiff_run.j', 's2smetric_tiff_', 1, str(3), CWD, None)
        COMMAND = f"srun --exclusive --ntasks 1 python {self.LISHDIR}/ghis2s/s2smetric/postprocess_nmme_job.py -y {self.YYYY} -m {self.MM} -w {CWD} -c {self.E2ESDIR}{self.config_file}"
        slurm_commands = [f"python {self.LISHDIR}/ghis2s/s2smetric/postprocess_nmme_job.py -y {self.YYYY} -m {self.MM} -w {CWD} -c {self.E2ESDIR}{self.config_file}"]

        # add LIS command
        with open('s2smetric_tiff_run.j', 'r') as file:
            filedata = file.read()   
        filedata = filedata.replace('COMMAND', COMMAND)
        with open('s2smetric_tiff_run.j', 'w') as file:
            file.write(filedata)

        self.create_dict('s2smetric_tiff_run.j', 's2smetric', prev=prev)
        utils.cylc_job_scripts(jobname + 'run.sh', 3, CWD, command_list=slurm_commands)

        os.chdir(self.E2ESDIR)        
        return

    def s2splots(self):
        if 's2smetric_tiff_run.j' in self.schedule:
            prev = 's2smetric_tiff_run.j'
        else:
            prev = None
            
        os.makedirs(self.E2ESDIR + 's2smetric/' + self.YYYY + self.MM, exist_ok=True)
        os.chdir(self.SCRDIR + 's2splots')
        self.create_symlink(self.E2ESDIR + 's2splots/', 's2splots')
        self.create_symlink(self.E2ESDIR + 's2smetric/', 's2smetric')
        CWD=self.SCRDIR + 's2splots'

        jobname='s2splots_01_'
        slurm_commands = []
        slurm_commands.append(f"python {self.LISHDIR}/ghis2s/s2splots/plot_mena.py -y {self.YYYY} -m {self.MM} -w {self.E2ESDIR} -c {self.E2ESDIR}{self.config_file}")
        slurm_commands.append(f"python {self.LISHDIR}/ghis2s/s2splots/plot_anom_verify.py -y {self.YYYY} -m {self.month} -w {self.E2ESDIR} -c {self.E2ESDIR}{self.config_file} -l 1")
        slurm_commands.append(f"python {self.LISHDIR}/ghis2s/s2splots/plot_anom_verify.py -y {self.YYYY} -m {self.month} -w {self.E2ESDIR} -c {self.E2ESDIR}{self.config_file} -l 2")
        slurm_commands.append(f"python {self.LISHDIR}/ghis2s/s2splots/plot_weekly_anom.py -y {self.YYYY} -m {self.month} -w {self.E2ESDIR} -c {self.E2ESDIR}{self.config_file}")

        tfile = self.sublist_to_file(slurm_commands, CWD)
        try:
            s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(2), CWD, tfile.name)
            self.create_dict('s2splots_01_run.j', 's2splots', prev=prev)
        finally:
            tfile.close()
            os.unlink(tfile.name)

        utils.cylc_job_scripts(jobname + 'run.sh', 2, CWD, command_list=slurm_commands)

        # 2nd job
        jobname='s2splots_02_'
        slurm_commands = [f"python {self.LISHDIR}/ghis2s/s2splots/plot_s2smetrics.py -y {self.YYYY} -m {self.MM} -w {self.E2ESDIR} -c {self.E2ESDIR}{self.config_file}"]
        par_info = {}
        par_info['NPROCS'] = str(14)
        par_info['MEM']= '10GB'
        tfile = self.sublist_to_file(slurm_commands, CWD)
        try:
            s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(2), CWD, tfile.name, parallel_run=par_info)
            self.create_dict(jobname + 'run.j', 's2splots', prev=prev)
        finally:
            tfile.close()
            os.unlink(tfile.name)

        utils.cylc_job_scripts(jobname + 'run.sh', 2, CWD, command_list=slurm_commands)

        # 3rd job
        jobname='s2splots_03_'
        slurm_commands = [f"python {self.LISHDIR}/ghis2s/s2splots/plot_hybas.py -y {self.YYYY} -m {self.month} -w {self.E2ESDIR} -c {self.E2ESDIR}{self.config_file}"]
        par_info = {}
        par_info['NPROCS'] = str(5)
        par_info['MEM']= '10GB'
        tfile = self.sublist_to_file(slurm_commands, CWD)
        try:
            s2s_api.python_job_file(self.E2ESDIR +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(2), CWD, tfile.name, parallel_run=par_info)
            self.create_dict(jobname + 'run.j', 's2splots', prev=prev)
        finally:
            tfile.close()
            os.unlink(tfile.name)

        utils.cylc_job_scripts(jobname + 'run.sh', 2, CWD, command_list=slurm_commands)
        os.chdir(self.E2ESDIR)
        
        return
         
    def main(self):
        # (1) Run CFSV2 file checker to ensure downloaded files are not corrupted/
        self.CFSv2_file_checker()
        super(S2Srun, self).NMME_file_checker()
        
        # (2) LISDA run
        self.lis_darun()

        # (3) LDT-ICS
        self.ldt_ics()

        # (4) BCSD
        self.bcsd()

        # (5) LIS FCST
        self.lis_fcst()

        # (6) S2SPOST
        self.s2spost()

        # (7) S2SMETRIC
        self.s2smetric()

        # (8) S2SPLOTS
        self.s2splots()

        # (9) Write CYLC workflow runtime snippet
        # -----------------------------------
        self.write_cylc_snippet() 
       
        return
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, type=str, help='config file')
    parser.add_argument('-y', '--year', required=True, type=int, help='forecast year')
    parser.add_argument('-m', '--month', required=True, type=int, help='forecast month')
    parser.add_argument('-s', '--step', required=False, default=None, type=str, help='S2S step: LISDA, LDTICS, BCSD, FCST, POST, METRICS or PLOTS')
    parser.add_argument('-o', '--one_step', action='store_true', help='Is only one step (default: False)?')
    parser.add_argument('-j', '--submit_job', action='store_true', help='Submit SLURM jobs (default: False)?')
    parser.add_argument('-r', '--report', action='store_true', help='Print report')
    args = parser.parse_args()

    s2s = S2Srun(year=args.year, month=args.month, config_file=args.config_file)

    # Print SLURM job report
    if  args.report:
        command = f"s2s_app/s2s_run.sh -y {args.year} -m {args.month} -c {args.config_file} -r Y"
        process = subprocess.run(command, shell=True)
        sys.exit()
    
    if args.step is not None:
        if args.step == 'LISDA':
            s2s.lis_darun()
        elif args.step == 'LDTICS':
            s2s.ldt_ics()
            if not args.one_step:
                s2s.bcsd()
                s2s.lis_fcst()
                s2s.s2spost()
                s2s.s2smetric()
                s2s.s2splots()
        elif args.step == 'BCSD':
            s2s.bcsd()
            if not args.one_step:
                s2s.lis_fcst()
                s2s.s2spost()
                s2s.s2smetric()
                s2s.s2splots()
        elif args.step == 'FCST':
            s2s.lis_fcst()
            if not args.one_step:
                s2s.s2spost()
                s2s.s2smetric()
                s2s.s2splots()
        elif args.step == 'POST':
            s2s.s2spost()
            if not args.one_step:
                s2s.s2smetric()
                s2s.s2splots()
        elif args.step == 'METRICS':
            s2s.s2smetric()
            if not args.one_step:
                s2s.s2splots()
        elif args.step == 'PLOTS':
            s2s.s2splots()
        
        # Write CYLC workflow runtime snippet
        # -----------------------------------
        s2s.write_cylc_snippet()

    else:
        s2s.main()
        
    # Submit SLURM jobs
    # -----------------
    if  args.submit_job:
        s2s.submit_jobs()
        
