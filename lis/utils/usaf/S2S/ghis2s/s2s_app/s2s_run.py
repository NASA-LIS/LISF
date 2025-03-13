import os
import sys
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import subprocess
import threading
import time
import platform
import shutil
import tempfile
import yaml
import argparse
from ghis2s.s2s_app import s2s_api
from ghis2s.shared import utils
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
        self.BWD = os.getcwd()
        self.SCRDIR = self.E2ESDIR + 'scratch/' + self.YYYY + self.MM + '/'
        self.MODELS = self.config["EXP"]["NMME_models"]
        self.CONSTRAINT = self.config['SETUP']['CONSTRAINT']
        self.schedule = self.job_schedule()
        
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

        #command = f"bash {self.E2ESDIR}/s2s_app/wget_cfsv2_oper_ts_e2es.sh -y {self.YYYY} -m {self.MM} -c {self.BWD}/{self.config_file} -d N"
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

    def job_schedule(self):
        def create_dict(prev):
            return {'jfiles': [], 'jobid': [], 'prev': prev}
            
        schedule = {}
        schedule['lisda'] = create_dict(None)
        schedule['ldtics'] = create_dict(['lisda'])
        schedule['bcsd01'] = create_dict(None)
        schedule['bcsd03'] = create_dict(None)
        schedule['bcsd04'] = create_dict(['bcsd01', 'bcsd03'])
        schedule['bcsd05'] = create_dict(['bcsd01', 'bcsd03'])
        schedule['bcsd06'] = create_dict(['bcsd04', 'bcsd05'])
        schedule['bcsd08'] = create_dict(['bcsd06'])
        schedule['bcsd09-10'] = create_dict(['bcsd08'])
        schedule['bcsd11-12'] = create_dict(['bcsd09-10'])
        schedule['lis_fcst'] = create_dict(['bcsd11-12'])
        schedule['s2spost'] = create_dict(['lis_fcst'])
        schedule['s2smetric'] = create_dict(['s2spost'])
        schedule['s2splots'] = create_dict(['s2smetric'])

        return schedule
    def split_list(self, input_list, length_sublist):
        result = []
        for i in range(0, len(input_list), length_sublist):
            sublist = input_list[i:i + length_sublist]
            result.append(sublist)
        return result
    #def split_list(self, input_list, num_sublists):
    #    """divide a list to sublists"""
    #    sublist_size = (len(input_list) - 1) // (num_sublists - 1)
    #    result = [input_list[i:i + sublist_size] for i in range(0, len(input_list) - 1, sublist_size)]
    #    result.append([input_list[-1]])
    #    return result

    def sublist_to_file(self, sublist, CWD):
        temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt')
        temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt', dir=CWD)
        for item in sublist:
            temp_file.write(f"{item}\n")
            temp_file.flush()
        return temp_file

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
        os.chdir(self.BWD)

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
        s2s_api.lis_job_file(self.BWD +'/' + self.config_file, 'lisda_run.j', 'lisda_', CWD, str(5))

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

        self.schedule['lisda']['jfiles'].append('lisda_run.j')

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
        
        os.chdir(self.BWD)
        
        return self.schedule

    def ldt_ics(self):
        """ LDT-ICS STEP """

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

        s2s_api.python_job_file(self.BWD +'/' + self.config_file, 'ldtics_run.j', 'ldtics_', str(1), str(2), CWD, None) 
        
        COMMAND="python {}/ghis2s/ldt_ics/generate_ldtconfig_files_ensrst_nrt.py -y {} -m {} -i ./lisda_output -w {} -s {}/{}".\
            format(self.LISHDIR, self.YYYY, self.month, CWD, self.BWD, self.config_file)

        # add command
        # ---------------
        with open('ldtics_run.j', 'r') as file:
            filedata = file.read()   
        filedata = filedata.replace('COMMAND', COMMAND)
        with open('ldtics_run.j', 'w') as file:
            file.write(filedata)

        self.schedule['ldtics']['jfiles'].append('ldtics_run.j')

        os.chdir(self.BWD)
        
        return self.schedule

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
        slurm_commands = bcsd.task_01.main(self.BWD +'/' + self.config_file, self.year, None,
                                           mmm, CWD, jobname, 1, 2, py_call=True)

        # multi tasks per job
        l_sub = 2
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i in range(len(slurm_sub)):
            tfile = self.sublist_to_file(slurm_sub[i], CWD)
            try:
                s2s_api.python_job_file(self.BWD +'/' + self.config_file, jobname + '{:02d}_run.j'.format(i+1),
                                    jobname+ '{:02d}_'.format(i+1), 1, str(3), CWD, tfile.name)
                self.schedule['bcsd01']['jfiles'].append(jobname + '{:02d}_run.j'.format(i+1))
            finally:
                tfile.close()
                os.unlink(tfile.name)
                
            utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 3, CWD, command_list=slurm_sub[i])
 
        # (3) bcsd03 regridding NMME
        # --------------------------
        jobname='bcsd03_'
        slurm_commands = \
            bcsd.task_03.main(self.BWD +'/' + self.config_file, self.year, self.MM,
                              jobname, 1, str(2), CWD, py_call=True)
        tfile = self.sublist_to_file(slurm_commands, CWD)
        try:
            s2s_api.python_job_file(self.BWD +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(2), CWD, tfile.name)
            self.schedule['bcsd03']['jfiles'].append(jobname + 'run.j')
        finally:
            tfile.close()
            os.unlink(tfile.name)
            
        utils.cylc_job_scripts(jobname + 'run.sh', 2, CWD, command_list=slurm_commands)

        # (4) bcsd04: Monthly "BC" step applied to CFSv2 (task_04.py, after 1 and 3)
        # --------------------------------------------------------------------------
        jobname='bcsd04_'
        slurm_commands = bcsd.task_04.main(self.BWD +'/' + self.config_file, self.year, self.year, mmm, self.MM, jobname,
                              1, 3, CWD, py_call=True)
        # multi tasks per job
        l_sub = 4
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i in range(len(slurm_sub)):
            tfile = self.sublist_to_file(slurm_sub[i], CWD)
            try:
                s2s_api.python_job_file(self.BWD +'/' + self.config_file, jobname + '{:02d}_run.j'.format(i+1),
                                    jobname + '{:02d}_'.format(i+1), 1, str(3), CWD, tfile.name)
                self.schedule['bcsd04']['jfiles'].append(jobname + '{:02d}_run.j'.format(i+1))
            finally:
                tfile.close()
                os.unlink(tfile.name)
                
            utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 3, CWD, command_list=slurm_sub[i])

        # (5) bcsd05: Monthly "BC" step applied to NMME (task_05.py: after 1 and 3)
        # -------------------------------------------------------------------------
        jobname='bcsd05_'
        slurm_commands = []
        for nmme_model in self.MODELS:
            var1 = bcsd.task_05.main(self.BWD +'/' + self.config_file, self.year, self.year, mmm, self.MM, jobname,
                                     1, 3, CWD, nmme_model, py_call=True)
            slurm_commands.append(var1)
        
        # multi tasks per job
        l_sub = 3
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i in range(len(slurm_sub)):
            tfile = self.sublist_to_file(slurm_sub[i], CWD)
            try:
                s2s_api.python_job_file(self.BWD +'/' + self.config_file, jobname + '{:02d}_run.j'.format(i+1),
                                    jobname + '{:02d}_'.format(i+1), 1, str(3), CWD, tfile.name)
                self.schedule['bcsd05']['jfiles'].append(jobname + '{:02d}_run.j'.format(i+1))
            finally:
                tfile.close()
                os.unlink(tfile.name)
                
            utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 3, CWD, command_list=slurm_sub[i])

        # (6) bcsd06: CFSv2 Temporal Disaggregation (task_06.py: after 4 and 5)
        # ---------------------------------------------------------------------
        jobname='bcsd06_'
        slurm_commands = bcsd.task_06.main(self.BWD +'/' + self.config_file, self.year, self.year, mmm, self.MM, jobname,
                                           1, 3, CWD, self.E2ESDIR, py_call=True)

        tfile = self.sublist_to_file(slurm_commands, CWD)
        try:
            s2s_api.python_job_file(self.BWD +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(5), CWD, tfile.name)
            self.schedule['bcsd06']['jfiles'].append(jobname + 'run.j')
        finally:
            tfile.close()
            os.unlink(tfile.name)
            
        utils.cylc_job_scripts(jobname + 'run.sh', 5, CWD, command_list=slurm_commands)

        # (8) bcsd08: NMME disaagregation
        # -------------------------------
        jobname='bcsd08_'
        slurm_commands = []
        for nmme_model in self.MODELS:
            var1 = bcsd.task_08.main(self.BWD +'/' + self.config_file, self.year, self.year, mmm, self.MM, jobname,
                                     1, 3, CWD, self.E2ESDIR, nmme_model, py_call=True)
            slurm_commands.append(var1)
        
        tfile = self.sublist_to_file(slurm_commands, CWD)
        try:
            s2s_api.python_job_file(self.BWD +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(5), CWD, tfile.name)
            self.schedule['bcsd08']['jfiles'].append(jobname + 'run.j')
        finally:
            tfile.close()
            os.unlink(tfile.name)

        utils.cylc_job_scripts(jobname + 'run.sh', 5, CWD, command_list=slurm_commands)

        # Task 9: Combine the CFSv2 forcing fields into final format for LIS to read
        # Task 10: Combine the NMME forcing fields into final format for LIS to read
        #          and symbolically link to the reusable CFSv2 met forcings
        # ---------------------------------------------------------------------------
        jobname='bcsd09-10_'
        slurm_9_10, slurm_11_12 = bcsd.task_09.main(self.BWD +'/' + self.config_file, self.year, self.year, mmm,
                                                    self.MM, jobname,1, 4, CWD, self.E2ESDIR, 'CFSv2', py_call=True)
        # bcsd09-10
        tfile = self.sublist_to_file(slurm_9_10, CWD)
        try:
            s2s_api.python_job_file(self.BWD +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(6), CWD, tfile.name)
            self.schedule['bcsd09-10']['jfiles'].append(jobname + 'run.j')
        finally:
            tfile.close()
            os.unlink(tfile.name)
            
        utils.cylc_job_scripts(jobname + 'run.sh', 6, CWD, command_list=slurm_9_10)

        # bcsd11-12
        jobname='bcsd11-12_'
        tfile = self.sublist_to_file(slurm_11_12, CWD)
        try:
            s2s_api.python_job_file(self.BWD +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(6), CWD, tfile.name)
            self.schedule['bcsd11-12']['jfiles'].append(jobname + 'run.j')
        finally:
            tfile.close()
            os.unlink(tfile.name)
            
        utils.cylc_job_scripts(jobname + 'run.sh', 6, CWD, command_list=slurm_11_12)

        return self.schedule
             
    def main(self):
        # (1) Run CFSV2 file checker to ensure downloaded files are not corrupted/
        #self.CFSv2_file_checker()

        # (2) LISDA run
        #self.lis_darun()

        # (3) LDT-ICS
        #self.ldt_ics()

        # (4) BCSD
        self.bcsd()

        # (5) LIS FCST
        #self.lis_fcst()
        

        return self.schedule

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, type=str, help='config file')
    parser.add_argument('-y', '--year', required=True, type=int, help='forecast year')
    parser.add_argument('-m', '--month', required=True, type=int, help='forecast month')
    parser.add_argument('-s', '--step', required=False, default=None, type=str, help='S2S step: LISDA, LDTICS, BCSD, FCST, POST, METRICS or PLOTS')
    parser.add_argument('-o', '--one_step', action='store_true', help='Is only one step (default: False)?')
    parser.add_argument('-j', '--submit_job', action='store_true', help='Submit SLURM jobs (default: False)')
    parser.add_argument('-r', '--report', action='store_true', help='Report')
    args = parser.parse_args()

    s2s = S2Srun(year=args.year, month=args.month, config_file=args.config_file)
    
    if args.step is not None:
        if args.step == 'LISDA':
            schedule = s2s.lis_darun()
        elif args.step == 'LDTICS':
            schedule =  s2s.ldt_ics()
        elif args.step == 'BCSD':
            schedule = s2s.bcsd()
            if not args.one_step:
                schedule = s2s.lis_fcst()
                schedule = s2s.s2spost()
                schedule = s2smetrics()
                schedule = s2splots()
        elif args.step == 'FCST':
            schedule = s2s.lis_fcst()
            if not args.one_step:
                schedule = s2s.s2spost()
                schedule = s2smetrics()
                schedule = s2splots()
        elif args.step == 'POST':
            schedule = s2spost()
            if not args.one_step:
                schedule = s2smetrics()
                schedule = s2splots()
        elif args.step == 'METRICS':
            schedule = s2smetrics()
            schedule = s2splots()
    else:
        schedule = s2s.main()

    # Submit SLURM jobs
    # -----------------
    if  args.submit_job:
        print(schedule)
        #s2s.submit_jobs(schedule)
        
