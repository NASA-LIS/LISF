'''
GHI-S2S Main Script
'''

import os
import sys
import glob
import subprocess
import re
import threading
import time
import platform
import shutil
import tempfile
import argparse
from datetime import datetime, timedelta, date
from dateutil.relativedelta import relativedelta
import xarray as xr
import yaml
from ghis2s.s2s_app import s2s_api
from ghis2s.shared import utils, logging_utils
from ghis2s.lis_fcst import generate_lis_config_scriptfiles_fcst
from ghis2s.s2spost import s2spost_driver
from ghis2s.s2smetric import s2smetric_driver
from ghis2s import bcsd
from ghis2s.shared.logging_utils import TaskLogger
from ghis2s.bcsd.bcsd_library.nmme_module import NMMEParams
# pylint: disable=too-many-lines

class DownloadForecasts():
    ''' Contains methods to download CFSv2 and NMME forecasts '''
    def __init__(self, year, month, config_file):
        with open(config_file, 'r', encoding="utf-8") as file:
            self.config = yaml.safe_load(file)
        self.config_file = config_file
        self.e2esdir = self.config['SETUP']['E2ESDIR']
        self.year = year
        self.month = month
        yyyy = f'{year:04d}'
        mm = f'{month:02d}'
        self.yyyymm = yyyy + mm
        self.scrdir = self.e2esdir + 'scratch/' + yyyy + mm + '/'
        if 'discover' in platform.node() or 'borg' in platform.node():
            self.cfsv2datadir = self.config['BCSD']['fcst_download_dir'] + "/Oper_TS/"
        else:
            self.cfsv2datadir = self.config['BCSD']['fcst_download_dir']
        self.patchfile = self.config['SETUP']['supplementarydir'] + \
            "/bcsd_fcst/patch_files/patch_files_list.txt"
        self.patchdir = self.config['SETUP']['supplementarydir'] + "/bcsd_fcst/patch_files/"
        self.cfsv2_log = os.path.join(self.scrdir, 'CFSv2_missing_corrupted_files')
        self.srcdir = "https://noaacfs.blob.core.windows.net/cfs"
        self.nmme_rawdir =  self.config['BCSD']['nmme_download_dir']
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

            return f'{prevmon:02d}', day1, day2, day3

        self.prevmon, self.day1, self.day2, self.day3 =  set_month_days(month)

    def nmme_file_checker(self):
        """ Checks NMME forecast files availability """
        print(" ")
        print("==========================================================================")
        print(" NMME Precip file checker.......")
        print("==========================================================================")
        nmme_path_dict = {
            'CFSv2': ['NCEP-CFSv2','NCEP-CFSv2'],
            'GEOSv2': ['NASA-GEOSS2S','NASA-GEOSS2S'],
            'CCM4': ['CanSIPS-IC3','CanSIPS-IC3'],
            'GNEMO5': ['CanSIPS-IC3','CanSIPS-IC3'],
            'CanESM5': ['CanSIPS-IC4','CanESM5'],
            'GNEMO52': ['CanSIPS-IC4', 'GEM5.2-NEMO'],
            'CCSM4': ['COLA-RSMAS-CCSM4', 'COLA-RSMAS-CCSM4'],
            'CESM1': ['COLA-RSMAS-CESM1', 'COLA-RSMAS-CESM1'],
            'GFDL': ['GFDL-SPEAR', 'GFDL-SPEAR'],
        }

        mon_abbr = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
               'Sep', 'Oct', 'Nov', 'Dec']
        infile_temp = '{}/{}/prec.{}.mon_{}.{:04d}.nc'

        havenot_files = {}

        for model in self.models:
            nmme_path = nmme_path_dict[model]
            infile = infile_temp.format(self.nmme_rawdir, nmme_path[0], nmme_path[1],
                                        mon_abbr[self.month-1], self.year)
            if not os.path.exists(infile):
                havenot_files[model] = infile

        if len(havenot_files) > 0:
            print("Following NMME precip files are missing:")
            for key, value in havenot_files.items():
                print(f"{key}: {value}")

            print("  ")
            print("Please download missing files before launching the forecast.")
            print("Alternatively, if you wish to exclude them from the current forecast:")
            print(f"      Remove {list(havenot_files.keys())} from NMME_models in")
            print(f"      {self.config_file} and launch the forecast.")

            sys.exit()

        print("Found below NMME precipitation files:")
        for model in self.models:
            nmme_path = nmme_path_dict[model]
            infile = infile_temp.format(self.nmme_rawdir, nmme_path[0], nmme_path[1],
                                        mon_abbr[self.month-1], self.year)
            print(f"{model}: {infile}")

    def cfsv2_download(self):
        """ download CFSv2 forecasts """
        if self.month > 1:
            os.makedirs(self.cfsv2datadir  + f'{self.year:04d}', exist_ok=True)
            os.chdir(self.cfsv2datadir  + f'{self.year:04d}')
            year2 = self.year
        else:
            # - Need to account for Dec/Jan crossover
            os.makedirs(self.cfsv2datadir  + f'{self.year -1:04d}', exist_ok=True)
            os.chdir(self.cfsv2datadir  + f'{self.year -1:04d}')
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
                    file = f"{self.srcdir}/cfs.{icdate}/{cycle}/time_grib_01/" + \
                        f"{vartype}.01.{icdate}{cycle}.daily.grb2"
                    if not os.path.isfile(file_name):
                        command = f"wget {file}"
                        subprocess.run(command, check=True)

    def cfsv2_file_checker(self):
        ''' Checks CFSv2 files available and uncorrupted '''
        def create_cfsv2_log(cfsv2_log):
            # Write the header and instructions to the log file
            with open(cfsv2_log, 'a', encoding="utf-8") as log_file:
                log_file.write(" ###############################################################\n")
                log_file.write("                 MISSING/INCOMPLETE CFSV2 FILES                 \n")
                log_file.write(" ###############################################################\n")
                log_file.write("                         \n")
                log_file.write("  A replacement file is required for each missing or corrupted\n")
                log_file.write("     file. CFSv2 replacement files are saved in:\n")
                log_file.write(f"  {self.patchdir}, \n")
                log_file.write("  and comma-delimited lines in:  \n")
                log_file.write(f"  {self.patchfile} \n")
                log_file.write("  lists the replacement file names for each corrupted file.\n")
                log_file.write("  The table has three columns:  \n")
                log_file.write("  YYYYMMDDHH, bad_file_name, replacement_file_name.     \n")
                log_file.write("                         \n")
                log_file.write(f" (1) cd {self.patchdir}      \n")
                log_file.write(" (2) Each problematic file name in the section below is followed\n")
                log_file.write("     by a list of wget commands to download a suitable \n")
                log_file.write("     replacement file in order of preference.\n")
                log_file.write("     Download the first suggested replacement file and add a new\n")
                log_file.write("     entry to: \n")
                log_file.write(f"     {self.patchfile} \n")
                log_file.write(
                    " (3) Repeat the same procedure to download replacements and update:\n")
                log_file.write(f"     {self.patchfile} \n")
                log_file.write("     for every missing/corrupted file.\n")
                log_file.write(
                    " (4) Relaunch forecast: s2s_app/s2s_run.sh -y YEAR -m MONTH -c CONFIGFILE\n")
                log_file.write(" (5) If any replacement file fails, you will be redirected\n")
                log_file.write(f"      to this file. {self.cfsv2_log}\n")
                log_file.write(" (6) Repeat steps (2) and (3) using a different replacement file\n")
                log_file.write("     for the original bad file.\n")
                log_file.write("                         \n")

        def neighb_days(vartype, icdate, cycle, mon):
            count_days = 1
            while count_days <= 4:
                # Calculate the previous date
                nwdate = (datetime.strptime(icdate, '%Y%m%d') - \
                          timedelta(days=count_days)).strftime('%Y%m%d')
                with open(self.cfsv2_log, 'a', encoding="utf-8") as f:
                    f.write(f"wget {self.srcdir}/cfs.{nwdate}/{cycle}/time_grib_01/"
                            f"{vartype}.01.{nwdate}{cycle}.daily.grb2\n")

                # Calculate the next date
                nwdate = (datetime.strptime(icdate, '%Y%m%d') + \
                          timedelta(days=count_days)).strftime('%Y%m%d')
                if datetime.strptime(nwdate, '%Y%m%d').month != int(mon):
                    with open(self.cfsv2_log, 'a', encoding="utf-8") as f:
                        f.write(f"wget {self.srcdir}/cfs.{nwdate}/{cycle}/time_grib_01/"
                                f"{vartype}.01.{nwdate}{cycle}.daily.grb2\n")

                count_days += 1

            with open(self.cfsv2_log, 'a', encoding="utf-8") as f:
                f.write("   \n")

        def print_message(log_file):
            with open(log_file, 'a', encoding="utf-8") as f:
                f.write("Note: If all recommended substitutes are also not available, you could try"
                        "      a different forecast hour from any of above dates.\n")
                f.write("\n")

        create_cfsv2_log(self.cfsv2_log)

        print(f"Previous mon, days 1-2-3 :: {self.prevmon}, {self.day1}-{self.day2}-{self.day3}")
        print(" ")
        print("==========================================================================")
        print(" CFSv2 file checker is running to ensure all forcings files are available ")
        print("       and not corrupted.......")
        print("==========================================================================")

        if self.month > 1:
            os.chdir(self.cfsv2datadir  + f'{self.year:04d}')
            year2 = self.year
        else:
            # - Need to account for Dec/Jan crossover
            os.chdir(self.cfsv2datadir  + f'{self.year -1:04d}')
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
                        with open(self.patchfile, 'r', encoding="utf-8") as f:
                            have_patch = f.read()
                            if file_name not in have_patch:
                                with open(self.cfsv2_log, 'a', encoding="utf-8") as log_file:
                                    log_file.write(f"{file_name}:  MISSING\n")
                                    log_file.write("Possible substitutes in order of preference:\n")
                                    neighb_days(vartype, icdate, cycle, self.prevmon)
                                    ret_code = 1

                    # File check 2: corrupted file
                    if os.path.isfile(file_name):
                        py_code = s2s_api.cfsv2_file_checker(file_name, self.yyyymm, py_call=True)
                        if py_code > 0:
                            with open(self.patchfile, 'r', encoding="utf-8") as f:
                                have_patch = f.read()

                            if file_name not in have_patch:
                                with open(self.cfsv2_log, 'a', encoding="utf-8") as log_file:
                                    log_file.write(f"{file_name}: CORRUPTED\n")
                                    log_file.write("Possible substitutes in order of preference:\n")
                                    neighb_days(vartype, icdate, cycle, self.prevmon)
                                ret_code = 1

                            else:
                                supfile = [line.split(',')[2].strip()
                                           for line in open(self.patchfile, encoding="utf-8")
                                           if file_name in line][0]
                                supfile_path = os.path.join(self.patchdir, supfile)
                                py_code = s2s_api.cfsv2_file_checker(supfile_path, self.yyyymm,
                                                                     py_call=True)

                                if py_code > 0:
                                    with open(self.cfsv2_log, 'a', encoding="utf-8") as log_file:
                                        log_file.write(
                                            f"{file_name}: Replacement {supfile} also CORRUPTED!\n")
                                        log_file.write(
                                            "Try downloading the next file (DON'T forget \n")
                                        log_file.write(f"   to update {self.patchfile})\n")
                                        neighb_days(vartype, icdate, cycle, self.prevmon)
                                        ret_code = 1

            os.chdir('..')

        if ret_code > 0:
            print("*** Missing or Incomplete CFSv2 forcing files were found ***.")
            print("Please follow the instructions in:")
            print(self.cfsv2_log)
            print_message(self.cfsv2_log)
        else:
            with open(self.cfsv2_log, 'a', encoding="utf-8") as log_file:
                log_file.write("**************************************************************\n")
                log_file.write(" SUCCESS ! All CFSv2 forcings files passed the file check.\n")
                log_file.write("**************************************************************\n")

                print("**************************************************************")
                print(" SUCCESS ! All CFSv2 forcings files passed the file check.")
                print("**************************************************************")

        print(" -- Done checking (and/or downloading) CFSv2 Forecast files -- ")

        return ret_code

class S2Srun(DownloadForecasts):
    ''' The main GHI-S2S class that creates monthly forecats '''
    def __init__(self, year, month, config_file, additional_env_vars=None):
        super().__init__(year, month, config_file)
        with open(config_file, 'r', encoding="utf-8") as file:
            self.config = yaml.safe_load(file)
        self.config_file = config_file
        self.year = year
        self.month = month
        self.yyyy = f'{year:04d}'
        self.mm = f'{month:02d}'
        self.e2esdir = self.config['SETUP']['E2ESDIR']
        self.hindcast = False
        if self.config['SETUP']['DATATYPE'] == 'hindcast':
            self.hindcast = True
            self.e2esdir = self.config['SETUP']['E2ESDIR'] + "/hindcast/"
        self.lisfdir = self.config['SETUP']['LISFDIR']
        self.lishdir = self.config['SETUP']['LISFDIR'] + 'lis/utils/usaf/S2S/'
        self.lishmod = self.config['SETUP']['LISFMOD']
        self.metforc = self.config['SETUP']['METFORC']
        self.e2esroot = self.config['SETUP']['E2ESDIR']
        self.domain = self.config['EXP']['DOMAIN']
        self.supdir = self.config['SETUP']['supplementarydir']
        self.ldtfile = self.config['SETUP']['ldtinputfile']
        self.scrdir = self.e2esdir + 'scratch/' + self.yyyy + self.mm + '/'
        self.models = self.config["EXP"]["NMME_models"]
        self.constraint = self.config['SETUP']['CONSTRAINT']
        self.fcst_model = self.config['BCSD']['metforce_source']
        self.schedule = {}
        self.additional_env_vars = additional_env_vars

        if not os.path.exists(self.e2esdir + 'scratch/'):
            subprocess.run(["setfacl", "-R", "-m", "u::rwx,g::rwx,o::r",
                            self.e2esdir], check=True)

        os.makedirs(self.scrdir + '/lis_darun/logs', exist_ok=True)
        os.makedirs(self.scrdir + '/ldt_ics/logs', exist_ok=True)
        os.makedirs(self.scrdir + '/bcsd_fcst/logs', exist_ok=True)
        os.makedirs(self.scrdir + '/lis_fcst/logs', exist_ok=True)
        os.makedirs(self.scrdir + '/s2spost/logs', exist_ok=True)
        os.makedirs(self.scrdir + '/s2smetric/logs', exist_ok=True)
        os.makedirs(self.scrdir + '/s2splots/logs', exist_ok=True)

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

    def delete_forecast(self):
        ''' deletes a forecast from the E2ES directory '''
        e2es, yyyy, yyyymm = self.e2esroot, self.yyyy, self.yyyymm
        lats, _ = utils.get_domain_info(self.config_file, coord=True)
        resol = f'{round((lats[1] - lats[0])*100)}km'
        met_dir = f'{self.fcst_model}_{resol}'
        mm = yyyymm[4:6]
        date_obj = datetime.strptime(f"{yyyy}-{mm}-01", "%Y-%m-%d")
        mmm01 = date_obj.strftime("%b").lower() + '01'
        umm = date_obj.strftime("%b")
        previous_month_date = date_obj - relativedelta(months=1)
        yyyymmp = previous_month_date.strftime("%Y%m")

        dirs = [
            f'{e2es}/lis_darun/output/EnKF/{yyyymmp}',
            f'{e2es}/lis_darun/output/ROUTING/{yyyymmp}/LIS_HIST_*.nc',
            f'{e2es}/lis_darun/output/SURFACEMODEL/{yyyymmp}/LIS_HIST_*.nc',
            f'{e2es}/lis_darun/output/lis.config_files/lis.config_darun_{yyyymmp}']
        dirs.extend([
            f'{e2es}/ldt_ics/{model}/ldtlog_noahmp401_{umm}{yyyy}.0000' for model in self.models])
        dirs.extend([
            f'{e2es}/ldt_ics/{model}/LIS_RST_HYMAP2_router_*{umm}{yyyy}*nc' for model in self.models
        ])
        dirs.extend([
            f'{e2es}/ldt_ics/{model}/LIS_RST_NOAHMP401_*{umm}{yyyy}*nc' for model in self.models])
        dirs.extend([
            f'{e2es}/bcsd_fcst/{met_dir}/raw/6-Hourly/{mmm01}/{yyyy}',
            f'{e2es}/bcsd_fcst/{met_dir}/raw/Monthly/{mmm01}/{yyyy}',
            f'{e2es}/bcsd_fcst/{met_dir}/bcsd/6-Hourly/{mmm01}/{yyyy}',
            f'{e2es}/bcsd_fcst/{met_dir}/bcsd/Monthly/{mmm01}/*{yyyy}_{yyyy}.nc',
            f'{e2es}/bcsd_fcst/{met_dir}/final/6-Hourly/{mmm01}/{yyyy}',
            f'{e2es}/bcsd_fcst/NMME/bcsd/Monthly/{mmm01}/PRECTOT.*{yyyy}_{yyyy}.nc'])
        dirs.extend([
            f'{e2es}/bcsd_fcst/NMME/raw/Monthly/{mmm01}/{model}/{yyyy}/' for model in self.models])
        dirs.extend([
            f'{e2es}/bcsd_fcst/NMME/final/6-Hourly/{model}/{mmm01}/{yyyy}' for model in self.models
        ])
        dirs.extend([f'{e2es}/lis_fcst/{yyyymm}',
                     f'{e2es}/s2spost/{yyyymm}',
                     f'{e2es}/s2smetric/{yyyymm}',
                     f'{e2es}/s2splots/{yyyymm}',
                     f'{e2es}/scratch/{yyyymm}'])

        for path in dirs:
            # Handle wildcards
            if '*' in path:
                matching_items = glob.glob(path)
                if not matching_items:
                    continue

                print(f"\nFound {len(matching_items)} items matching pattern: {path}")
                confirm = input(
                    f"Delete all {len(matching_items)} items above? (y/n): ").strip().lower()

                # delete each matching item
                for item in matching_items:
                    if os.path.exists(item):
                        if confirm.startswith('y'):
                            try:
                                if os.path.isdir(item):
                                    shutil.rmtree(item)
                                    print(f"Directory deleted: {item}")
                                else:
                                    os.remove(item)
                                    print(f"File deleted: {item}")
                            except Exception as e:
                                print(f"Error deleting {item}: {e}")

            # Handle non-wildcard paths
            elif os.path.exists(path):
                confirm = input(f"Delete {path}? (y/n): ").strip().lower()
                if confirm.startswith('y'):
                    try:
                        if os.path.isdir(path):
                            shutil.rmtree(path)
                            print(f"Directory deleted: {path}")
                        else:
                            os.remove(path)
                            print(f"File deleted: {path}")
                    except Exception as e:
                        print(f"Error deleting {path}: {e}")
            else:
                print(f"Not found: {path}")

    def rst_file_checker(self):
        ''' Checks availability of LSM, DA and HyMAP restart files '''
        # previous month
        date_obj = datetime.strptime(f"{self.yyyy}-{self.mm}-01", "%Y-%m-%d")
        previous_month_date = date_obj - relativedelta(months=1)
        yyyymmp = previous_month_date.strftime("%Y%m")
        yyyyp = yyyymmp[:4]
        mmp = yyyymmp[4:6]

        rst_path = self.e2esdir + 'lis_darun/output/'
        dapert = f'DAPERT/{yyyyp}{mmp}/LIS_DAPERT_{yyyyp}{mmp}010000.d01.bin'
        surf_rst = f'SURFACEMODEL/{yyyyp}{mmp}/LIS_RST_NOAHMP401_{yyyyp}{mmp}010000.d01.nc'
        rout_rst = f'ROUTING/{yyyyp}{mmp}/LIS_RST_HYMAP2_router_{yyyyp}{mmp}010000.d01.nc'

        # Check each restart file
        print(" ")
        print("==========================================================================")
        print(" Land Model Restart file checker.......")
        print("==========================================================================")

        restart_files = {
            'DAPERT': rst_path + dapert,
            'SURFACEMODEL': rst_path + surf_rst,
            'ROUTING': rst_path + rout_rst
        }

        missing_files = []
        for file_type, file_path in restart_files.items():
            if not os.path.exists(file_path):
                missing_files.append(f"{file_type}: {file_path}")

        if missing_files:
            print(f"[ERROR] Missing restart files for {yyyyp}-{mmp}:")
            for missing_file in missing_files:
                print(f"[ERROR]   {missing_file}")
            print("[ERROR] Cannot launch forecast without required restart files!")
            sys.exit()
        else:
            print(f"[INFO] All restart files found for {yyyyp}-{mmp}")
            for file_type, file_path in restart_files.items():
                print(f"[INFO]   {file_type}: {os.path.basename(file_path)}")

    def clim_files_checker(self):
        ''' Checks availability of clim files '''
        print(" ")
        print("==========================================================================")
        print(" BCSD Climatological file checker.......")
        print("==========================================================================")

        def check_file(infile, description):
            try:
                dataset = xr.open_dataset(infile)
                print(f"{description} : {infile} is good.")
                dataset.close()
                del dataset
            except Exception as e:
                print(f"PROBLEM openning {infile} {e}")

        fcst_clim_file_template = '{}hindcast/bcsd_fcst/{}_{}/raw/Climatology/{}/{}_fcst_clim.nc'
        date_obj = datetime.strptime(f"{self.yyyy}-{self.mm}-01", "%Y-%m-%d")
        mmm = date_obj.strftime("%b").lower() + '01'
        lats, _ = utils.get_domain_info(self.config_file, coord=True)
        resol = f'{round((lats[1] - lats[0])*100)}km'

        def check_available_file(fcst_var):
            alternate_name = {
                'LWGAB': 'LWS',
                'SWGDN': 'SLRSF',
                'QV2M': 'Q2M',
            }
            if os.path.exists(fcst_clim_file_template.format(
                    self.e2esdir, self.fcst_model, resol, mmm, fcst_var)):
                return fcst_clim_file_template.format(
                    self.e2esdir, self.fcst_model, resol, mmm, fcst_var)

            return fcst_clim_file_template.format(
                self.e2esdir, self.fcst_model, resol, mmm,alternate_name[fcst_var])

        obs_clim_file_template = '{}/{}_obs_clim.nc'

        # Obs clim files
        obs_var_list = ["PRECTOT", "LWGAB", "SWGDN", "PS", "QV2M", "T2M", "U10M"]
        obs_path = self.e2esdir + 'hindcast/bcsd_fcst/USAF-LIS7.3rc8_25km/raw/Climatology/'

        for var in obs_var_list:
            check_file(obs_clim_file_template.format(obs_path, var), 'USAF-LIS7.3rc8_25km Clim')

        # metforce clim
        fcst_var_list = ["PRECTOT", "LWGAB", "SWGDN", "PS", "QV2M", "T2M", "WIND10M"]
        for _fcst_var in fcst_var_list:
            infile = fcst_clim_file_template.format(
                self.e2esdir, self.fcst_model, resol, mmm, _fcst_var)
            if _fcst_var in ['LWGAB', 'SWGDN', 'QV2M']:
                infile = check_available_file(_fcst_var)
            check_file(infile, self.fcst_model+ '_' + resol)

    def cfsv2_file_checker(self):
        spinner_done = [False]
        def spinner():
            """Display a spinner while waiting for a process to complete."""
            spin_chars = ['\\', '|', '/', '-']
            idx = 0
            while not spinner_done[0]:
                print(f'Please wait... {spin_chars[idx]}', end='\r', flush=True)
                idx = (idx + 1) % len(spin_chars)
                time.sleep(0.1)

        spinner_thread = threading.Thread(target=spinner)
        spinner_thread.start()
        ret_code = super().cfsv2_file_checker()
        spinner_done[0] = True
        spinner_thread.join()
        print('\rDone')

        if ret_code > 0:
            print(f"Error return code from the CFSv2 file download checker :: {ret_code}")
            print("> 0 :: Exiting from s2s_run.py --")
            sys.exit(ret_code)

    def create_dict(self, jobfile, subdir, prev=None):
        ''' create dependency dictionary '''
        self.schedule[jobfile] = {'subdir': subdir, 'jobid': None, 'prev': []}
        if prev is not None:
            if isinstance(prev, list):
                self.schedule[jobfile]['prev'].extend(prev)
            else:
                self.schedule[jobfile]['prev'].append(prev)

    def submit_jobs(self):
        ''' submits SLURM jobs '''
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
                result = subprocess.run(sbatch_command, capture_output=True, text=True,
                                        check=True, shell=True)
                match = re.search(r'Submitted batch job (\d+)', result.stdout)
                if match:
                    job_id = match.group(1)
                    print(f"Submitted successfully. Job ID: {job_id}; Job Script: {job_script}")
                    return job_id

                print("Failed to extract job ID from sbatch output")
                return None
            except subprocess.CalledProcessError as e:
                print(f"An error occurred while submitting the job: {e}")
                return None

        job_schedule = os.path.join(self.scrdir, 'SLURM_JOB_SCHEDULE')
        if os.path.exists(job_schedule):
            os.remove(job_schedule)

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

        with open(job_schedule, "a", encoding="utf-8") as file:
            for line in header:
                file.write(line + "\n")

        utils.update_job_schedule(job_schedule, "JOB ID", "JOB SCRIPT", "AFTER")
        dir_list = []
        for jfile in self.schedule.keys():
            subdir = self.schedule[jfile]['subdir']
            if subdir not in dir_list:
                dir_list.append(subdir)
                with open(job_schedule, "a", encoding="utf-8") as file:
                    for line in sub_section.get(subdir):
                        file.write(line + "\n")
            os.chdir(self.scrdir + subdir)
            prev_ids = get_previds(jfile)
            job_id = submit_slurm_job(jfile, prev_id=prev_ids)
            self.schedule[jfile]['jobid'] = job_id
            if isinstance(prev_ids, list):
                prev_str = ','.join(map(str, prev_ids))
            else:
                # If prev_id is a single value, convert it to string
                prev_str = str(prev_ids)
            utils.update_job_schedule(job_schedule, job_id, jfile, prev_str)

    def split_list(self, input_list, length_sublist):
        ''' splits job list '''
        result = []
        for i in range(0, len(input_list), length_sublist):
            sublist = input_list[i:i + length_sublist]
            result.append(sublist)
        return result

    def sublist_to_file(self, sublist, cwd):
        ''' writes sublists to temporary files '''
        temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt')
        temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt', dir=cwd)
        for item in sublist:
            temp_file.write(f"{item}\n")
            temp_file.flush()
        return temp_file

    def write_cylc_snippet(self):
        """ writes Cylc runtime snippet """
        is_nccs = True
        if os.path.isfile(self.lisfdir + 'env/discover/' + self.lishmod):
            modulepath = self.lisfdir + 'env/discover/'
        else:
            modulepath = self.config['SETUP']['supplementarydir'] + '/env/'
            is_nccs = False

        def write_log_monitoring_script():
            """writes the shell script to run log monitoring"""   
            log_script_path = f"{self.scrdir}ghis2s_log.sh"
            with open(log_script_path, 'w', encoding="utf-8") as f:
                f.write("#!/bin/bash\n")
                f.write("set -eu\n")
                f.write("export USE_CYLC_ENV=0\n")
                if is_nccs:
                    f.write("source /etc/profile.d/modules.sh\n")
                f.write(f"module use -a {modulepath}\n")
                f.write(f"module load {self.lishmod}\n")
                f.write("ulimit -s unlimited\n")
                f.write(f"export PYTHONPATH={self.lisfdir}lis/utils/usaf/S2S/\n")
                f.write(f"cd {self.scrdir}\n")
                f.write(f"python {self.lisfdir}lis/utils/usaf/S2S/ghis2s/s2s_app/s2s_run.py "
                        f"-y {self.year} -m {self.month} -c {self.e2esdir}{self.config_file} -l\n")

            os.chmod(log_script_path, 0o755)

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

            with open(slurm_file, 'r', encoding="utf-8") as file:
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

        def write_lines(file, jfile, subdir, directives, environment):
            file.write(f"    [[{jfile}]]\n")

            sh_script = self.scrdir + subdir + '/' + jfile + '.sh'
            file.write(f"        script = {sh_script}\n")

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
        cylc_file = f"{self.scrdir}CYLC_workflow.rc"

        # Build dependency graph based on schedule structure
        dependency_map = {}
        for jfile in self.schedule.keys():
            task_name = jfile.removesuffix('.j')
            prev_tasks = []
            if len(self.schedule[jfile]['prev']) > 0:
                for pfile in self.schedule[jfile]['prev']:
                    prev_task = pfile.removesuffix('.j')
                    prev_tasks.append(prev_task)
                    dependency_map[task_name] = prev_tasks

        # write log monitoring bash script
        write_log_monitoring_script()

        # write flow.cylc
        has_dependencies = any(dependencies for dependencies in dependency_map.values())
        with open(cylc_file, 'w', encoding="utf-8") as file:
            # Write header
            file.write("#!jinja2\n")
            file.write("# S2S Forecast Workflow\n")
            file.write("  \n")
            file.write("[meta]\n")
            file.write("    title = S2S Forecast Workflow\n")
            file.write("    description = S2S Forecast Initialized on YYYY-MM\n")
            file.write("  \n")
            file.write("[scheduler]\n")
            file.write("    UTC mode = True\n")
            file.write("    cycle point format = %Y%m%dT%H%M\n")
            file.write("  \n")

            # Write scheduling section
            file.write("[scheduling]\n")
            file.write("    initial cycle point = now\n")
            file.write("    final cycle point = +P2D\n")
            file.write("    runahead limit = P2D\n")
            file.write("  \n")

            file.write("    [[xtriggers]]\n")
            file.write("        log_check = wall_clock(offset=PT0M)\n")
            file.write("  \n")

            # Update scheduling structure
            file.write("    [[graph]]\n")
            file.write("        R1 = \"\"\"\n")

            # Generate dependency graph from dependency_map
            if has_dependencies:
                for task, dependencies in dependency_map.items():
                    if dependencies:
                        dep_str = " & ".join(dependencies)
                        file.write(f"                {dep_str} => {task}\n")
            else:
                # Single task run: write all tasks as standalone
                for jfile in self.schedule.keys():
                    task_name = jfile.removesuffix('.j')
                    file.write(f"                {task_name}\n")
                    dependency_map[task_name] = task_name
            # Find terminal tasks for final log collection
            all_tasks = set(dependency_map.keys())
            dependency_tasks = set()
            for deps in dependency_map.values():
                dependency_tasks.update(deps)
            terminal_tasks = all_tasks - dependency_tasks

            # Add final log collection dependency
            if terminal_tasks:
                terminal_str = " & ".join(sorted(terminal_tasks))
                file.write(
                    f"                {terminal_str} => final_log_collect => stop_log_monitor\n")
            else:
                all_tasks_str = " & ".join(sorted(all_tasks))
                file.write(
                    f"            {all_tasks_str}:finish-all => final_log_collect => stop_log_monitor\n")
            file.write("                stop_log_monitor => !log_monitor\n")
            file.write("        \"\"\"\n")
            file.write("  \n")
            file.write("        R//PT15M = \"\"\"\n")
            file.write("           @log_check => log_monitor\n")
            file.write("        \"\"\"\n")
            file.write("  \n")

            # Write runtime section
            file.write("[runtime]\n")
            file.write("    [[root]]\n")
            file.write("        platform = slurm-ghi\n")
            file.write("        pre-script = \"\"\"\n")
            if 'discover' in platform.node() or 'borg' in platform.node():
                file.write("            source /etc/profile.d/modules.sh\n")
            file.write(f"            module load {self.lishmod}\n")
            file.write("        \"\"\"\n")
            file.write("        [[[environment]]]\n")
            file.write("            USE_CYLC_ENV = 1\n")
            file.write(f"            MODULEPATH = {modulepath}:$MODULEPATH\n")
            file.write(f"            PYTHONPATH = {self.lisfdir}lis/utils/usaf/S2S/\n")

            # Add additional environment variables if provided
            if self.additional_env_vars:
                for key, value in self.additional_env_vars.items():
                    file.write(f"            {key} = {value}\n")

            file.write("        [[[mail]]]\n")
            file.write("            to = USEREMAIL\n")
            file.write("        [[[events]]]\n")
            file.write("            mail events = failed\n")
            file.write("  \n")

            file.write("    [[log_monitor]]\n")
            file.write(f"        script = {self.scrdir}ghis2s_log.sh\n")
            file.write("        platform = localhost\n")
            file.write("  \n")

            file.write("    [[final_log_collect]]\n")
            file.write(f"        script = {self.scrdir}ghis2s_log.sh\n")
            file.write("        platform = localhost\n")
            file.write("  \n")

            file.write("    [[stop_log_monitor]]\n")
            file.write("        platform = localhost\n")
            file.write("        script = \"\"\"\n")
            file.write("            cylc stop $CYLC_WORKFLOW_ID --now\n")
            file.write("        \"\"\"\n")
            file.write("  \n")

            for jfile in self.schedule.keys():
                subdir = self.schedule[jfile]['subdir']
                directives, pre_script, environment = \
                    extract_slurm_info(self.scrdir + subdir + '/' + jfile)

                #inherit_list = []
                #if len(self.schedule[jfile]['prev']) > 0:
                #    for pfile in self.schedule[jfile]['prev']:
                #        inherit_list.append(pfile.removesuffix('.j'))
                #if len(inherit_list) == 0:
                #    inherit_list = None
                write_lines(file, jfile.removesuffix('.j'), subdir,
                            directives, environment)

    def lis_darun(self):
        """ LIS DARUN STEP """
        os.makedirs(self.e2esdir + '/lis_darun/input/', exist_ok=True)
        os.makedirs(self.e2esdir + '/lis_darun/output/lis.config_files/', exist_ok=True)
        os.chdir(self.e2esdir + '/lis_darun/input/')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/forcing_variables.txt',
                            'forcing_variables.txt')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/noahmp401_parms','noahmp401_parms')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/template_files','template_files')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/attribs','attribs')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/tables','tables')
        self.create_symlink(self.supdir + '/lis_darun/cdf/' +self.domain,'cdf')
        self.create_symlink(self.supdir + '/lis_darun/RS_DATA','RS_DATA')
        self.create_symlink(self.supdir + '/lis_darun/' + self.ldtfile, self.ldtfile)
        os.chdir(self.e2esdir)

        # previous month
        date_obj = datetime.strptime(f"{self.yyyy}-{self.mm}-01", "%Y-%m-%d")
        previous_month_date = date_obj - relativedelta(months=1)
        yyyymmp = previous_month_date.strftime("%Y%m")
        yyyyp = yyyymmp[:4]
        mmp = yyyymmp[4:6]
        monp = int(mmp)
        pertmode = self.config['EXP']['pertmode']
        os.makedirs(self.scrdir + '/lis_darun/input/', exist_ok=True)
        os.chdir(self.scrdir +'/lis_darun/input/')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/forcing_variables.txt',
                            'forcing_variables.txt')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/noahmp401_parms','noahmp401_parms')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/template_files','template_files')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/attribs','attribs')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/tables','tables')
        self.create_symlink(self.supdir + '/lis_darun/cdf/' +self.domain,'cdf')
        self.create_symlink(self.supdir + '/lis_darun/RS_DATA','RS_DATA')
        self.create_symlink(self.supdir + '/lis_darun/' + self.ldtfile, self.ldtfile)

        os.chdir(self.scrdir +'/lis_darun/')
        cwd=self.scrdir +'lis_darun'
        self.create_symlink(self.lisfdir + '/lis/LIS', 'LIS')
        self.create_symlink(self.e2esdir + '/lis_darun/output', 'output')
        self.create_symlink(self.metforc, self.metforc.rstrip('/').split('/')[-1])

        # write lisda_run.j
        # -----------------
        s2s_api.lis_job_file(self.e2esroot +'/' + self.config_file, 'lisda_run.j',
                             'lisda_', cwd, str(5))

        if 'discover' in platform.node() or 'borg' in platform.node():
            if 'mil' in self.constraint:
                command = './LIS'
            else:
                slurm_ntasks = os.getenv('SLURM_NTASKS', '1')
                command = f'mpirun -np {slurm_ntasks} ./LIS'
        else:
            command = 'srun ./LIS'

        # add LIS command
        # ---------------
        with open('lisda_run.j', 'r', encoding="utf-8") as file:
            filedata = file.read()
        filedata = filedata.replace('COMMAND', command)
        with open('lisda_run.j', 'w', encoding="utf-8") as file:
            file.write(filedata)

        shutil.copy('lisda_run.j', 'lisda_run.sh')
        utils.remove_sbatch_lines('lisda_run.sh')
        self.create_dict('lisda_run.j', 'lis_darun')

        # configure lis.config
        # --------------------

        shutil.copy(self.e2esdir + '/lis_darun/input/template_files/lis.config_template.' +\
                    self.domain, 'lis.config')
        dapertrstfile = './output/DAPERT/' + yyyyp + mmp + f'/LIS_DAPERT_{yyyyp}{mmp}010000.d01.bin'
        noahmp401rstfile = './output/SURFACEMODEL/' + yyyyp + mmp +\
            f'/LIS_RST_NOAHMP401_{yyyyp}{mmp}010000.d01.nc'
        hymap2rstfile = './output/ROUTING/' + yyyyp + mmp +\
            f'/LIS_RST_HYMAP2_router_{yyyyp}{mmp}010000.d01.nc'
        lsmlislogfile = cwd + f'/logs_{yyyyp}{mmp}/lislog'

        with open('lis.config', 'r', encoding="utf-8") as file:
            filedata = file.read()

        filedata = filedata.replace('DAPERTRSTFILE', dapertrstfile)
        filedata = filedata.replace('NOAHMP401RSTFILE', noahmp401rstfile)
        filedata = filedata.replace('HYMAP2RSTFILE', hymap2rstfile)
        filedata = filedata.replace('STARTYR', yyyyp)
        filedata = filedata.replace('STARTMO', str(monp))
        filedata = filedata.replace('STARTDA', '1')
        filedata = filedata.replace('FINALYR', self.yyyy)
        filedata = filedata.replace('FINALMO', str(self.month))
        filedata = filedata.replace('FINALDA', '1')
        filedata = filedata.replace('PERTMODE', pertmode)
        filedata = filedata.replace('LSMLISLOGFILE', lsmlislogfile)

        with open('lis.config', 'w', encoding="utf-8") as file:
            file.write(filedata)

        shutil.copy('lis.config', f'output/lis.config_files/lis.config_darun_{yyyyp}{mmp}')
        _ = TaskLogger('lisda_run.j',
                            os.getcwd(),
                            f'LISDA log files are at: \n{cwd}/logs_{yyyymmp[:6]}/')

        os.chdir(self.e2esdir)

    def ldt_ics(self):
        """ LDT-ICS STEP """
        if 'lisda_run.j' in self.schedule.keys():
            prev = 'lisda_run.j'
        else:
            prev = None

        os.makedirs(self.e2esdir + '/ldt_ics/input/', exist_ok=True)
        os.makedirs(self.scrdir + '/ldt_ics/input/', exist_ok=True)
        os.chdir(self.e2esdir + '/ldt_ics/')
        [os.makedirs(model, exist_ok=True) for model in self.models]
        os.makedirs('ldt.config_files', exist_ok=True)
        os.makedirs('template_files', exist_ok=True)
        shutil.copy(
            self.lishdir +
            f'/ghis2s/ldt_ics/template_files/ldt.config_noahmp401_nmme_TEMPLATE.{self.domain}',
            'template_files/ldt.config_noahmp401_nmme_TEMPLATE')

        os.chdir(self.e2esdir + '/ldt_ics/input/')
        self.create_symlink(self.supdir + '/lis_darun/' + self.ldtfile, self.ldtfile)
        self.create_symlink(self.supdir + '/LS_PARAMETERS', 'LS_PARAMETERS')
        os.chdir(self.scrdir + '/ldt_ics/input/')
        self.create_symlink(self.supdir + '/lis_darun/' + self.ldtfile, self.ldtfile)
        self.create_symlink(self.supdir + '/LS_PARAMETERS', 'LS_PARAMETERS')

        os.chdir(self.scrdir + '/ldt_ics')
        cwd=self.scrdir + 'ldt_ics'
        self.create_symlink(self.lisfdir + '/ldt/LDT', 'LDT')
        self.create_symlink(self.e2esdir + '/lis_darun/output', 'lisda_output')

        for model in self.models:
            self.create_symlink(self.e2esdir + '/ldt_ics/' + model, model)

        self.create_symlink(self.e2esdir + '/ldt_ics/ldt.config_files', 'ldt.config_files')
        self.create_symlink(self.e2esdir + '/ldt_ics/template_files', 'template_files')

        # configure batch script
        # ----------------------

        s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, 'ldtics_run.j', 'ldtics_',
                                str(1), str(2), cwd, None)

        command=(f"python {self.lishdir}/ghis2s/ldt_ics/generate_ldtconfig_files_ensrst_nrt.py -y"
                 f" {self.yyyy} -m {self.month} -i ./lisda_output -w {cwd} -s"
                 f" {self.e2esdir}/{self.config_file}")

        # add command
        # ---------------
        with open('ldtics_run.j', 'r', encoding="utf-8") as file:
            filedata = file.read()
        filedata = filedata.replace('COMMAND', command)
        with open('ldtics_run.j', 'w', encoding="utf-8") as file:
            file.write(filedata)

        self.create_dict('ldtics_run.j', 'ldt_ics', prev=prev)
        shutil.copy('ldtics_run.j', 'ldtics_run.sh')
        utils.remove_sbatch_lines('ldtics_run.sh')
        #utils.cylc_job_scripts('ldtics_run.sh', 2, cwd, command_list=[command])
        mon_abbr = date(int(self.yyyy), int(self.mm), 1).strftime('%b')
        _ = TaskLogger('ldtics_run.j',
                            os.getcwd(),
                            (f'LDTICS log files are at: \n{self.e2esdir}/ldt_ics/*/'
                             f'ldtlog_noahmp401_{mon_abbr}{self.yyyy}.0000'))

        os.chdir(self.e2esdir)

    def bcsd(self):
        """ BCSD 12 steps """
        lats, _ = utils.get_domain_info(self.config_file, coord=True)
        resol = f'{round((lats[1] - lats[0])*100)}km'
        fcast_clim_dir=\
            f'{self.e2esroot}/hindcast/bcsd_fcst/{self.fcst_model}_{resol}/raw/Climatology/'
        nmme_clim_dir=f'{self.e2esroot}/hindcast/bcsd_fcst/NMME/raw/Climatology/'
        usaf_25km=f'{self.e2esroot}/hindcast/bcsd_fcst/USAF-LIS7.3rc8_25km/raw/Climatology/'

        os.makedirs(self.e2esdir + '/bcsd_fcst', exist_ok=True)
        os.chdir(self.e2esdir + '/bcsd_fcst')
        os.makedirs('USAF-LIS7.3rc8_25km/raw', exist_ok=True)
        os.makedirs(f'{self.fcst_model}_{resol}/raw', exist_ok=True)
        os.makedirs('NMME/raw', exist_ok=True)

        # link Climatology directories
        os.chdir(self.e2esdir + f'/bcsd_fcst/{self.fcst_model}_{resol}/raw')
        self.create_symlink(fcast_clim_dir, 'Climatology')

        os.chdir(self.e2esdir + '/bcsd_fcst/NMME/raw')
        self.create_symlink(nmme_clim_dir, 'Climatology')

        os.chdir(self.e2esdir + '/bcsd_fcst/USAF-LIS7.3rc8_25km/raw')
        self.create_symlink(usaf_25km, 'Climatology')

        # manage jobs from SCRATCH
        os.chdir(self.scrdir + 'bcsd_fcst')
        self.create_symlink(self.e2esdir + 'bcsd_fcst/', 'bcsd_fcst')
        cwd=self.scrdir + 'bcsd_fcst'

        date_obj = datetime.strptime(f"{self.yyyy}-{self.mm}-01", "%Y-%m-%d")
        mmm = date_obj.strftime("%b").lower()

        # (1) metforce_regridding (formerly bcsd01) - regrid metforce (CFSv2/GEOSv3 files
        # ---------------------------------------------------------------------
        jobname='mf_regrid_'
        resol_info = {
            '25km': {'CPT': str(1), 'MEM':'13GB', 'NT': str(self.config["EXP"]["lead_months"]),
                     'TPN': 2*int(self.config["EXP"]["lead_months"]), 'l_sub': 2,'HOURS': str(1)},
            '10km': {'CPT': str(1), 'MEM':'120GB', 'NT': str(1), 'TPN': 2,
                     'l_sub': 2, 'HOURS': str(3)},
            '5km': {'CPT': str(1), 'MEM':'240GB', 'NT': str(1), 'TPN': 1, 'l_sub': 1,
                    'HOURS': str(self.config["EXP"]["lead_months"])},
        }
        info = resol_info[resol]
        slurm_commands = bcsd.metforce_regridding.main(
            self.e2esroot +'/' + self.config_file, self.year, None,
            mmm, cwd, jobname, 1, 2, py_call=True)

        # multi tasks per job
        l_sub = info['l_sub']
        par_info = {}
        par_info['CPT'] = info['CPT']
        par_info['MEM'] = info['MEM']
        par_info['NT'] = info['NT']
        par_info['TPN'] = info['TPN']
        par_info['MP'] = False
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i, sub_val in enumerate(slurm_sub):
            tfile = self.sublist_to_file(sub_val, cwd)
            try:
                s2s_api.python_job_file(self.e2esroot +'/' + self.config_file,
                                        jobname + f'{i+1:02d}_run.j', jobname+ f'{i+1:02d}_',
                                        info['TPN'], info['HOURS'], cwd,
                                        tfile.name, parallel_run=par_info)
                self.create_dict(jobname+ f'{i+1:02d}_run.j', 'bcsd_fcst')
            finally:
                tfile.close()
                os.unlink(tfile.name)

            shutil.copy(jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_run.sh')
            utils.remove_sbatch_lines(jobname + f'{i+1:02d}_run.sh')
            #utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), info['HOURS'],
            # cwd, command_list=sub_val)

        # (3) precip_regridding (formerly bcsd03) regridding precipitation (NMME)
        # ---------------------------------------
        jobname='pr_regrid_'
        slurm_commands = bcsd.precip_regridding.main(self.e2esroot +'/' + self.config_file,
                                                     self.year, self.month, jobname, 1, str(2),
                                                     cwd, py_call=True)
        tfile = self.sublist_to_file(slurm_commands, cwd)
        try:
            s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(3), cwd, tfile.name)
            self.create_dict(jobname+ 'run.j', 'bcsd_fcst')
        finally:
            tfile.close()
            os.unlink(tfile.name)

        shutil.copy(jobname + 'run.j', jobname + 'run.sh')
        utils.remove_sbatch_lines(jobname + 'run.sh')
        #utils.cylc_job_scripts(jobname + 'run.sh', 3, cwd, command_list=slurm_commands)

        # (4) metforce_biascorrection (formerly bcsd04): Monthly "BC" step applied to CFSv2
        #    (task_04.py, after 1 and 3)
        # --------------------------------------------------------------------------
        jobname='mf_biascorr_'
        par_info = {}
        par_info['CPT'] = '10'
        par_info['NT']= str(1)
        par_info['MEM']= '240GB'
        par_info['TPN'] = None
        par_info['MP'] = True
        prev = [f"{key}" for key in self.schedule.keys() if 'mf_regrid_' in key]
        slurm_commands = bcsd.metforce_biascorrection.main(
            self.e2esroot +'/' + self.config_file, self.year, self.year, mmm,
            self.mm, jobname, 1, 3, cwd, py_call=True)
        # multi tasks per job
        l_sub = 1
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i, sub_val in enumerate(slurm_sub):
            tfile = self.sublist_to_file(sub_val, cwd)
            try:
                s2s_api.python_job_file(self.e2esroot +'/' + self.config_file,
                                        jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_',
                                        1, str(4), cwd, tfile.name, parallel_run=par_info)
                self.create_dict(jobname + f'{i+1:02d}_run.j', 'bcsd_fcst', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)

            shutil.copy(jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_run.sh')
            utils.remove_sbatch_lines(jobname + f'{i+1:02d}_run.sh')
            #utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 4, cwd,
            # command_list=sub_val)

        # (5) precip_biascorrection (bcsd05): Monthly "BC" step applied to NMME
        # -------------------------------------------------------------------------
        jobname='pr_biascorr_'
        prev = [f"{key}" for key in self.schedule.keys() if 'pr_regrid_' in key]
        slurm_commands = []
        for nmme_model in self.models:
            var1 = bcsd.precip_biascorrection.main(self.e2esroot +'/' + self.config_file,
                                                   self.year, self.year, mmm, self.mm, jobname,
                                                   1, 3, cwd, nmme_model, py_call=True)
            slurm_commands.extend(var1)

        # multi tasks per job
        l_sub = 1
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i, sub_val in enumerate(slurm_sub):
            tfile = self.sublist_to_file(sub_val, cwd)
            nmme_model = sub_val[0].split()[7]
            ens_num = NMMEParams(nmme_model).ens_num
            par_info['CPT'] = str(max(10,ens_num))
            try:
                s2s_api.python_job_file(self.e2esroot +'/' + self.config_file,
                                        jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_', 1,
                                        str(4), cwd, tfile.name, parallel_run=par_info)
                self.create_dict(jobname + f'{i+1:02d}_run.j', 'bcsd_fcst', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)

            shutil.copy(jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_run.sh')
            utils.remove_sbatch_lines(jobname + f'{i+1:02d}_run.sh')
            #utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 4, cwd,
            #  command_list=sub_val)

        # (6) metforce_temporal_disaggregation (formerly bcsd06: after 4 and 5)
        # ---------------------------------------------------------------------
        jobname='mf_tempdis_'
        par_info = {}
        par_info['CPT'] = '12'
        par_info['NT']= str(1)
        par_info['MEM']= '240GB'
        par_info['TPN'] = None
        par_info['MP'] = True
        prev = [f"{key}" for key in self.schedule.keys() if 'mf_biascorr_' in key]
        slurm_commands = bcsd.metforce_temporal_disaggregation.main(self.e2esroot +'/' + \
                                                                    self.config_file, self.year,
                                                                    self.year, mmm, self.mm,
                                                                    jobname, 1, 3, cwd,
                                                                    self.e2esdir, py_call=True)

        # multi tasks per job
        l_sub = 3
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i, sub_val in enumerate(slurm_sub):
            tfile = self.sublist_to_file(sub_val, cwd)
            try:
                s2s_api.python_job_file(self.e2esroot +'/' + self.config_file,
                                        jobname + f'{i+1:02d}_run.j',
                                        jobname + f'{i+1:02d}_', 1, str(4), cwd, tfile.name,
                                        parallel_run=par_info)
                self.create_dict(jobname + f'{i+1:02d}_run.j', 'bcsd_fcst', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)

            shutil.copy(jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_run.sh')
            utils.remove_sbatch_lines(jobname + f'{i+1:02d}_run.sh')
            #utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 5, cwd,
            # command_list=sub_val)

        # (8) precip temporal disaggregation (formerly task_08.py)
        # -------------------------------
        jobname='pr_tempdis_'
        par_info = {}
        par_info['CPT'] = '15'
        par_info['NT']= str(1)
        par_info['MEM']= '240GB'
        par_info['TPN'] = None
        par_info['MP'] = True
        slurm_commands = []
        prev = [f"{key}" for key in self.schedule.keys() if 'mf_regrid_' in key]
        prev.extend([f"{key}" for key in self.schedule.keys() if 'pr_biascorr_' in key])
        for nmme_model in self.models:
            var1 = bcsd.precip_temporal_disaggregation.main(self.e2esroot +'/' + self.config_file,
                                                            self.year, self.year, mmm, self.mm,
                                                            jobname, 1, 3, cwd, self.e2esdir,
                                                            nmme_model, py_call=True)
            slurm_commands.extend(var1)

        # multi tasks per job
        l_sub = 2
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i, sub_val in enumerate(slurm_sub):
            tfile = self.sublist_to_file(sub_val, cwd)
            try:
                s2s_api.python_job_file(self.e2esroot +'/' + self.config_file,
                                        jobname + f'{i+1:02d}_run.j',
                                        jobname + f'{i+1:02d}_', 1, str(2), cwd, tfile.name,
                                        parallel_run=par_info)
                self.create_dict(jobname + f'{i+1:02d}_run.j', 'bcsd_fcst', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)

            shutil.copy(jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_run.sh')
            utils.remove_sbatch_lines(jobname + f'{i+1:02d}_run.sh')
            #utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 5, cwd,
            #    command_list=sub_val)

        # Task 9: Combine the CFSv2 forcing fields into final format for LIS to read
        # Task 10: Combine the NMME forcing fields into final format for LIS to read
        #          and symbolically link to the reusable CFSv2 met forcings
        # ---------------------------------------------------------------------------
        prev = [f"{key}" for key in self.schedule.keys() if 'mf_tempdis_' in key]
        slurm_9_10 = bcsd.combine_forcings.main(self.e2esroot +'/' + self.config_file, self.year,
                                                self.year, mmm, self.mm, jobname,1, 4, cwd,
                                                self.e2esdir, self.fcst_model, py_call=True)
        # bcsd09
        jobname='combine_files_'
        par_info = {}
        par_info['CPT'] = str(self.config['BCSD']['nof_raw_ens'])
        par_info['NT'] = str(1)
        par_info['MEM'] = '240GB'
        par_info['TPN'] = None
        par_info['MP'] = True

        tfile = self.sublist_to_file([slurm_9_10[0]], cwd)
        try:
            s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(1), cwd, tfile.name, parallel_run=par_info)
            self.create_dict(jobname + 'run.j', 'bcsd_fcst', prev=prev)
        finally:
            tfile.close()
            os.unlink(tfile.name)

        shutil.copy(jobname + 'run.j', jobname + 'run.sh')
        utils.remove_sbatch_lines(jobname + 'run.sh')
        #utils.cylc_job_scripts(jobname + 'run.sh', 6, cwd, command_list=slurm_9_10[0])

        os.chdir(self.e2esdir)

    def lis_fcst(self):
        """ LIS forecast """
        def check_recommend_nseg(lead_mons, nseg):
            """Check validity and print recommendation if invalid"""
            def check_jobseg(lead_mons, nseg):
                """Check if nseg is valid """
                if lead_mons % nseg == 0:
                    return True

                base_length = lead_mons // nseg + 1
                remainder = base_length * nseg - lead_mons
                return remainder < base_length

            def try_these_nseg(lead_mons):
                """Find all valid nseg values"""
                good_nseg = []
                for i in range(1, lead_mons + 1):
                    if check_jobseg(lead_mons, i):
                        good_nseg.append(i)
                return good_nseg

            is_valid = check_jobseg(lead_mons, nseg)
            good_nseg = []
            if not is_valid:
                good_nseg = try_these_nseg(lead_mons)
            return is_valid, good_nseg

        # checks validity of JOB_SEGMENTS
        for model in self.models:
            is_valid, good_nseg = check_recommend_nseg(
                self.config["EXP"]["lead_months"],
                self.config['FCST']['JOB_SEGMENTS'][0].get(model))
            if not is_valid:
                print(f"[ERROR] {model} Unsupported nof job segments: "
                      f"{self.config['FCST']['JOB_SEGMENTS'][0].get(model)}")
                print(f'Try these instead: {good_nseg}')
                sys.exit()

        prev = [job for job in ['ldtics_run.j', 'combine_files_run.j']
                if job in self.schedule] or None

        if prev is not None:
            prev.extend([f"{key}" for key in self.schedule.keys() if 'pr_tempdis_' in key])
        else:
            temp_keys = [f"{key}" for key in self.schedule.keys() if 'pr_tempdis_' in key]
            prev = temp_keys or None

        jobname='lis_fcst'
        os.makedirs(self.e2esdir + 'lis_fcst/', exist_ok=True)
        os.makedirs(self.e2esdir + 'lis_fcst/input/LDT_ICs/', exist_ok=True)

        os.chdir(self.e2esdir + 'lis_fcst/input/')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/forcing_variables.txt',
                            'forcing_variables.txt')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/noahmp401_parms','noahmp401_parms')
        self.create_symlink(self.lishdir + '/ghis2s/lis_fcst/template_files','template_files')
        self.create_symlink(self.lishdir + '/ghis2s/lis_fcst/tables','tables')
        self.create_symlink(self.supdir + '/lis_darun/' + self.ldtfile, self.ldtfile)

        os.chdir(self.e2esdir + 'lis_fcst/input/LDT_ICs/')
        for model in self.models:
            os.makedirs(self.e2esdir + 'lis_fcst/' + self.yyyy + self.mm + '/' + model + '/logs/',
                        exist_ok=True)
            self.create_symlink(self.e2esdir + 'ldt_ics/'  + model, model)

        os.makedirs(self.scrdir + 'lis_fcst/input', exist_ok=True)
        os.chdir(self.scrdir + 'lis_fcst/input/')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/forcing_variables.txt',
                            'forcing_variables.txt')
        self.create_symlink(self.lishdir + '/ghis2s/lis_darun/noahmp401_parms','noahmp401_parms')
        self.create_symlink(self.lishdir + '/ghis2s/lis_fcst/template_files','template_files')
        self.create_symlink(self.lishdir + '/ghis2s/lis_fcst/tables','tables')
        self.create_symlink(self.e2esdir + '/ldt_ics','LDT_ICs')
        self.create_symlink(self.supdir + '/lis_darun/' + self.ldtfile, self.ldtfile)

        os.chdir(self.scrdir + 'lis_fcst/')
        cwd=self.scrdir + 'lis_fcst'
        self.create_symlink(self.lisfdir + 'lis/LIS', 'LIS')
        self.create_symlink(self.e2esdir + 'lis_fcst/input', 'input')
        self.create_symlink(self.e2esdir + 'lis_fcst/' + self.yyyy + self.mm, self.yyyy + self.mm)
        self.create_symlink(self.e2esdir + 'bcsd_fcst', 'bcsd_fcst')

        generate_lis_config_scriptfiles_fcst.main(self.e2esroot + self.config_file, self.year,
                                                  self.month, cwd, jobname)

        for model in self.models:
            lindex = ''
            job_list = sorted(glob.glob(f"{jobname}_{model}*_run.j"))
            n_files = len(job_list)
            for file_no, jfile in enumerate(job_list):
                if n_files > 1:
                    if file_no == 0:
                        self.create_dict(jfile, 'lis_fcst', prev=prev)
                    else:
                        self.create_dict(jfile, 'lis_fcst', prev=job_list[file_no-1])
                        lindex = f'_{file_no:02d}'
                else:
                    self.create_dict(jfile, 'lis_fcst', prev=prev)
                _ = TaskLogger(jfile,
                                    os.getcwd(),
                                    f'{model} log files are at: \n{self.e2esdir}/lis_fcst/'
                                    f'{self.yyyy}{self.mm}/{model}/logs/lislog{lindex}.')

        os.chdir(self.e2esdir)

    def s2spost(self):
        """ S2SPOST STEP """
        fcst_list = [item for item in self.schedule.keys() if re.match('lis_fcst', item)]
        if len(fcst_list) > 0:
            prev=[]
            for model in self.models:
                sublist = sorted([item for item in fcst_list if model in item])
                last_item = sublist[-1] if sublist else None
                if last_item is not None:
                    prev.append(last_item)
        else:
            prev = None

        if self.hindcast:
            [os.makedirs(self.e2esdir + 's2spost/' + self.mm + '/' + self.yyyy + self.mm + '/' + \
                         model, exist_ok=True) for model in self.models]
        else:
            [os.makedirs(self.e2esdir + 's2spost/' + self.yyyy + self.mm + '/' + model,
                         exist_ok=True) for model in self.models]
        os.chdir(self.scrdir + 's2spost')
        self.create_symlink(self.e2esdir + 'lis_fcst/' + self.yyyy + self.mm + '/', 'lis_fcst')
        self.create_symlink(self.e2esdir + 'lis_fcst/input/', 'input')
        cwd=self.scrdir + 's2spost'

        jobname='s2spost_'
        slurm_commands = []
        monthly_commands = []
        weekly_commands = []
        for model in self.models:
            if self.hindcast:
                self.create_symlink(self.e2esdir + 's2spost/' + self.mm + '/' + self.yyyy + \
                                    self.mm + '/' + model, model)
            else:
                self.create_symlink(self.e2esdir + 's2spost/' + self.yyyy + self.mm + '/' + model,
                                    model)
            var1, var2, var3 = s2spost_driver.main(self.e2esroot +'/' + self.config_file, self.year,
                                                   self.month, jobname, 1, str(3), cwd, model,
                                                   py_call=True)
            slurm_commands.extend(var1)
            monthly_commands.extend(var2)
            weekly_commands.extend(var3)

        # post process hindcast
        if self.hindcast:
            jobname='s2spost_weekly_'
            l_sub = 6
            slurm_sub = self.split_list(weekly_commands, l_sub)
            for i, sub_val in enumerate(slurm_sub):
                tfile = self.sublist_to_file(sub_val, cwd)
                try:
                    s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, jobname + \
                                            f'{i+1:02d}_run.j',
                                            jobname + f'{i+1:02d}_', 1, str(30), cwd, tfile.name)
                    self.create_dict(jobname + f'{i+1:02d}_run.j', 's2spost', prev=prev)
                finally:
                    tfile.close()
                    os.unlink(tfile.name)

                shutil.copy(jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_run.sh')
                utils.remove_sbatch_lines(jobname + f'{i+1:02d}_run.sh')
                #utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 1, cwd,
                #  command_list=sub_val)
            os.chdir(self.e2esdir)
            return

        # processing dailies multi tasks per job
        l_sub = 27
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i, sub_val in enumerate(slurm_sub):
            tfile = self.sublist_to_file(sub_val, cwd)
            try:
                s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, jobname + \
                                        f'{i+1:02d}_run.j',
                                        jobname + f'{i+1:02d}_', 1, str(4), cwd, tfile.name)
                self.create_dict(jobname + f'{i+1:02d}_run.j', 's2spost', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)

            shutil.copy(jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_run.sh')
            utils.remove_sbatch_lines(jobname + f'{i+1:02d}_run.sh')
            #utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 4, cwd,
            #   command_list=sub_val)

        # processing monthlies multi tasks per job
        jobname='s2spost_mon_'
        prev = sorted(glob.glob("s2spost_0*_run.j"))
        l_sub = 9
        slurm_sub = self.split_list(monthly_commands, l_sub)
        for i, sub_val in enumerate(slurm_sub):
            tfile = self.sublist_to_file(sub_val, cwd)
            try:
                s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, jobname + \
                                        f'{i+1:02d}_run.j',
                                        jobname + f'{i+1:02d}_', 1, str(3), cwd, tfile.name)
                self.create_dict(jobname + f'{i+1:02d}_run.j', 's2spost', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)

            shutil.copy(jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_run.sh')
            utils.remove_sbatch_lines(jobname + f'{i+1:02d}_run.sh')
            #utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 1, cwd,
            #  command_list=sub_val)

        # processing dailies multi tasks per job
        jobname='s2spost_weekly_'
        l_sub = 18
        slurm_sub = self.split_list(weekly_commands, l_sub)
        for i, sub_val in enumerate(slurm_sub):
            tfile = self.sublist_to_file(sub_val, cwd)
            try:
                s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, jobname + \
                                        f'{i+1:02d}_run.j',
                                        jobname + f'{i+1:02d}_', 1, str(1), cwd, tfile.name)
                self.create_dict(jobname + f'{i+1:02d}_run.j', 's2spost', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)

            shutil.copy(jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_run.sh')
            utils.remove_sbatch_lines(jobname + f'{i+1:02d}_run.sh')
            #utils.cylc_job_scripts(jobname + '{:02d}_run.sh'.format(i+1), 1, cwd,
            #  command_list=sub_val)
        os.chdir(self.e2esdir)
        return

    def s2smetric(self):
        """ S2SMETRICS STEP """
        s2spost_list = [item for item in self.schedule.keys() if re.match('s2spost', item)]
        if len(s2spost_list) > 0:
            prev = s2spost_list
        else:
            prev = None

        os.makedirs(self.e2esdir + 's2smetric/' + self.yyyy + self.mm, exist_ok=True)
        os.chdir(self.scrdir + 's2smetric')
        cwd=self.scrdir + 's2smetric'
        self.create_symlink(self.e2esdir + 's2spost/', 's2spost')

        # post process hindcast
        if self.hindcast:
            jobname='s2smetric_weekly_'
            par_info = {}
            par_info['CPT'] = str(len(self.models))
            par_info['MEM']= '240GB'
            par_info['NT']= str(1)
            par_info['TPN'] = None
            par_info['MP'] = True
            slurm_commands = []
            var1 = s2smetric_driver.main(self.e2esroot +'/' + self.config_file, self.year,
                                         self.month, cwd, jobname=jobname, ntasks=1,
                                         hours=str(1), nmme_model= 'all_models', py_call=True,
                                         weekly=True)
            slurm_commands.extend(var1)
            tfile = self.sublist_to_file(slurm_commands, cwd)
            try:
                s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, jobname + 'run.j',
                                        jobname, 1, str(6), cwd, tfile.name, parallel_run=par_info)
                self.create_dict(jobname+ 'run.j', 's2smetric')
            finally:
                tfile.close()
                os.unlink(tfile.name)
            os.chdir(self.e2esdir)
            return

        self.create_symlink(self.e2esdir + 's2smetric/', 's2smetric')
        if len(s2spost_list) > 0:
            prev = [item for item in self.schedule.keys() if re.match('s2spost_mon', item)]
        jobname='s2smetric_'
        metric_vars = self.config["POST"]["metric_vars"]
        par_info = {}
        par_info['CPT'] = str(2*len(metric_vars))
        par_info['MEM']= '240GB'
        par_info['NT']= str(1)
        par_info['TPN'] = None
        par_info['MP'] = True
        for i, model in enumerate(self.models):
            slurm_commands = s2smetric_driver.main(self.e2esroot +'/' + self.config_file, self.year,
                                                   self.month, cwd, jobname=jobname, ntasks=1,
                                             hours=str(1), nmme_model= model, py_call=True)

            tfile = self.sublist_to_file(slurm_commands, cwd)
            try:
                s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, jobname + \
                                        f'{i+1:02d}_run.j',
                                        jobname + f'{i+1:02d}_', 1, str(1), cwd, tfile.name,
                                        parallel_run=par_info)
                self.create_dict(jobname + f'{i+1:02d}_run.j', 's2smetric', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)

            shutil.copy(jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_run.sh')
            utils.remove_sbatch_lines(jobname + f'{i+1:02d}_run.sh')
            #utils.cylc_job_scripts(jobname + 'run.sh', 4, cwd, command_list=slurm_commands)

        # weekly metrics
        jobname='s2smetric_weekly_'
        if len(s2spost_list) > 0:
            prev = [item for item in self.schedule.keys() if re.match('s2spost_weekly', item)]
        slurm_commands = []
        weekly_vars = self.config["POST"]["weekly_vars"]
        par_info = {}
        par_info['CPT'] = str(2*len(weekly_vars))
        par_info['MEM']= '240GB'
        par_info['NT']= str(1)
        par_info['TPN'] = None
        par_info['MP'] = True
        for i, model in enumerate(self.models):
            slurm_commands = s2smetric_driver.main(self.e2esroot +'/' + self.config_file, self.year,
                                                   self.month, cwd, jobname=jobname, ntasks=1,
                                                   hours=str(4), nmme_model= model, py_call=True,
                                                   weekly=True)

            tfile = self.sublist_to_file(slurm_commands, cwd)
            try:
                s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, jobname + \
                                        f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_', 1, str(1),
                                        cwd, tfile.name, parallel_run=par_info)
                self.create_dict(jobname + f'{i+1:02d}_run.j', 's2smetric', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)

            shutil.copy(jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_run.sh')
            utils.remove_sbatch_lines(jobname + f'{i+1:02d}_run.sh')
            #utils.cylc_job_scripts(jobname + 'run.sh', 4, cwd, command_list=slurm_commands)

        # write tiff files
        jobname='s2smetric_tiff_'
        prev = [f"{key}" for key in self.schedule.keys() if 's2smetric_0' in key]
        par_info = {}
        par_info['CPT'] = str(2*len(metric_vars))
        par_info['MEM']= '240GB'
        par_info['NT']= str(1)
        par_info['TPN'] = None
        par_info['MP'] = True

        command = [f"python {self.lishdir}/ghis2s/s2smetric/s2smetric_driver.py -y {self.yyyy}"
                   f" -m {self.mm} -w {cwd} -c {self.e2esroot}{self.config_file}"]
        tfile = self.sublist_to_file(command, cwd)
        try:
            s2s_api.python_job_file(
                self.e2esroot +'/' + self.config_file, 's2smetric_tiff_run.j',
                's2smetric_tiff_', 1, str(1), cwd, tfile.name, parallel_run=par_info)
            self.create_dict('s2smetric_tiff_run.j', 's2smetric', prev=prev)
        finally:
            tfile.close()
            os.unlink(tfile.name)

        shutil.copy(jobname + 'run.j', jobname + 'run.sh')
        utils.remove_sbatch_lines(jobname + 'run.sh')
        #utils.cylc_job_scripts(jobname + 'run.sh', 3, cwd, command_list=command)

        # write tiff files
        jobname='s2smetric_weekly_tiff_'
        prev = [f"{key}" for key in self.schedule.keys() if 's2smetric_weekly_0' in key]
        par_info = {}
        par_info['CPT'] = str(2*len(weekly_vars))
        par_info['MEM']= '240GB'
        par_info['NT']= str(1)
        par_info['TPN'] = None
        par_info['MP'] = True

        command = [f"python {self.lishdir}/ghis2s/s2smetric/s2smetric_driver.py -y {self.yyyy}"
                   f" -m {self.mm} -w {cwd} -c {self.e2esroot}{self.config_file} -W"]
        tfile = self.sublist_to_file(command, cwd)
        try:
            s2s_api.python_job_file(self.e2esroot +'/' + self.config_file,
                                    's2smetric_weekly_tiff_run.j', 's2smetric_weekly_tiff_', 1,
                                    str(1), cwd, tfile.name, parallel_run=par_info)
            self.create_dict('s2smetric_weekly_tiff_run.j', 's2smetric', prev=prev)
        finally:
            tfile.close()
            os.unlink(tfile.name)

        shutil.copy(jobname + 'run.j', jobname + 'run.sh')
        utils.remove_sbatch_lines(jobname + 'run.sh')
        #utils.cylc_job_scripts(jobname + 'run.sh', 3, cwd, command_list=command)

        os.chdir(self.e2esdir)
        return

    def s2splots(self):
        ''' S2SPLOTS method '''
        if 's2smetric_tiff_run.j' in self.schedule:
            prev = ['s2smetric_tiff_run.j', 's2smetric_weekly_tiff_run.j']
        else:
            prev = None

        os.makedirs(self.e2esdir + 's2splots/' + self.yyyy + self.mm, exist_ok=True)
        os.chdir(self.scrdir + 's2splots')
        self.create_symlink(self.e2esdir + 's2splots/', 's2splots')
        self.create_symlink(self.e2esdir + 's2smetric/', 's2smetric')
        cwd=self.scrdir + 's2splots'

        jobname='s2splots_01_'
        slurm_commands = []
        slurm_commands.append(
            f"python {self.lishdir}/ghis2s/s2splots/plot_mena.py -y {self.yyyy}"
            f" -m {self.mm} -w {self.e2esdir} -c {self.e2esroot}{self.config_file}")
        slurm_commands.append(
            f"python {self.lishdir}/ghis2s/s2splots/plot_anom_verify.py -y {self.yyyy}"
            f" -m {self.month} -w {self.e2esdir} -c {self.e2esroot}{self.config_file} -l 1")
        slurm_commands.append(
            f"python {self.lishdir}/ghis2s/s2splots/plot_anom_verify.py"
            f" -y {self.yyyy} -m {self.month} -w {self.e2esdir}"
            f" -c {self.e2esroot}{self.config_file} -l 2")
        slurm_commands.append(
            f"python {self.lishdir}/ghis2s/s2splots/plot_weekly_anom.py -y {self.yyyy}"
            f" -m {self.month} -w {self.e2esdir} -c {self.e2esroot}{self.config_file}")

        tfile = self.sublist_to_file(slurm_commands, cwd)
        try:
            s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(2), cwd, tfile.name)
            self.create_dict('s2splots_01_run.j', 's2splots', prev=prev)
        finally:
            tfile.close()
            os.unlink(tfile.name)

        shutil.copy(jobname + 'run.j', jobname + 'run.sh')
        utils.remove_sbatch_lines(jobname + 'run.sh')
        #utils.cylc_job_scripts(jobname + 'run.sh', 2, cwd, command_list=slurm_commands)

        # 2nd job
        jobname='s2splots_02_'
        slurm_commands = [
            f"python {self.lishdir}/ghis2s/s2splots/plot_s2smetrics.py -y {self.yyyy}"
            f" -m {self.mm} -w {self.e2esdir} -c {self.e2esroot}{self.config_file} -M ANOM"]
        par_info = {}
        par_info['CPT'] = str(7)
        par_info['MEM']= '240GB'
        par_info['NT']= str(1)
        par_info['TPN'] = None
        par_info['MP'] = True
        tfile = self.sublist_to_file(slurm_commands, cwd)
        try:
            s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(2), cwd, tfile.name, parallel_run=par_info)
            self.create_dict(jobname + 'run.j', 's2splots', prev=prev)
        finally:
            tfile.close()
            os.unlink(tfile.name)

        shutil.copy(jobname + 'run.j', jobname + 'run.sh')
        utils.remove_sbatch_lines(jobname + 'run.sh')
        #utils.cylc_job_scripts(jobname + 'run.sh', 2, cwd, command_list=slurm_commands)

       # 3rd job
        jobname='s2splots_03_'
        slurm_commands = [
            f"python {self.lishdir}/ghis2s/s2splots/plot_s2smetrics.py -y {self.yyyy}"
            f" -m {self.mm} -w {self.e2esdir} -c {self.e2esroot}{self.config_file} -M SANOM"]
        par_info = {}
        par_info['CPT'] = str(7)
        par_info['MEM']= '240GB'
        par_info['NT']= str(1)
        par_info['TPN'] = None
        par_info['MP'] = True
        tfile = self.sublist_to_file(slurm_commands, cwd)
        try:
            s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(2), cwd, tfile.name, parallel_run=par_info)
            self.create_dict(jobname + 'run.j', 's2splots', prev=prev)
        finally:
            tfile.close()
            os.unlink(tfile.name)

        shutil.copy(jobname + 'run.j', jobname + 'run.sh')
        utils.remove_sbatch_lines(jobname + 'run.sh')
        #utils.cylc_job_scripts(jobname + 'run.sh', 2, cwd, command_list=slurm_commands)

        # 4th job
        jobname='s2splots_04_'
        slurm_commands = [
            f"python {self.lishdir}/ghis2s/s2splots/plot_hybas.py -y {self.yyyy}"
            f" -m {self.month} -w {self.e2esdir} -c {self.e2esroot}{self.config_file}"]
        par_info = {}
        par_info['CPT'] = str(5)
        par_info['MEM']= '240GB'
        par_info['NT']= str(1)
        par_info['TPN'] = None
        par_info['MP'] = True
        tfile = self.sublist_to_file(slurm_commands, cwd)
        try:
            s2s_api.python_job_file(self.e2esroot +'/' + self.config_file, jobname + 'run.j',
                                    jobname, 1, str(2), cwd, tfile.name, parallel_run=par_info)
            self.create_dict(jobname + 'run.j', 's2splots', prev=prev)
        finally:
            tfile.close()
            os.unlink(tfile.name)

        shutil.copy(jobname + 'run.j', jobname + 'run.sh')
        utils.remove_sbatch_lines(jobname + 'run.sh')
        #utils.cylc_job_scripts(jobname + 'run.sh', 2, cwd, command_list=slurm_commands)
        os.chdir(self.e2esdir)

    def main(self):
        ''' S2SRun driver '''
        # (1) File checkers to ensure RST files are available & downloaded files are not corrupted.
        self.rst_file_checker()
        super().nmme_file_checker()
        self.clim_files_checker()
        self.cfsv2_file_checker()

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, type=str, help='config file')
    parser.add_argument('-y', '--year', required=True, type=int, help='forecast year')
    parser.add_argument('-m', '--month', required=True, type=int, help='forecast month')
    parser.add_argument('-s', '--step', required=False, default=None, type=str,
                        help='S2S step: LISDA, LDTICS, BCSD, FCST, POST, METRICS or PLOTS')
    parser.add_argument('-o', '--one_step', action='store_true',
                        help='Is only one step (default: False)?')
    parser.add_argument('-j', '--submit_job', action='store_true',
                        help='Submit SLURM jobs (default: False)?')
    parser.add_argument('-r', '--report', action='store_true',
                        help='Print report')
    parser.add_argument('-l', '--logging', action='store_true',
                        help='Write centralized log file')
    parser.add_argument('-d', '--delete_forecast', action='store_true',
                        help='delete forecast')

    args = parser.parse_args()

    s2s = S2Srun(year=args.year, month=args.month, config_file=args.config_file)

    # Print SLURM job report
    if  args.report:
        JOB_SCHEDULE = os.path.join(s2s.scrdir, 'SLURM_JOB_SCHEDULE')
        if os.path.exists(JOB_SCHEDULE):
            CMD = f"s2s_app/s2s_run.sh -y {args.year} -m {args.month} -c {args.config_file} -r Y"
            process = subprocess.run(CMD, shell=True, check=True)
        else:
            if args.step is None:
                CYLC_WORKFLOW = f'cylc_e2e_{s2s.yyyy}{s2s.mm}'
            else:
                CYLC_WORKFLOW = f'cylc_{args.step.lower()}_{s2s.yyyy}{s2s.mm}'
            CMD = f"sh s2s_app/cylc_walltime.sh {CYLC_WORKFLOW}"
            process = subprocess.run(CMD, shell=True, check=True)
        sys.exit()

    # Write LOG file
    if args.logging:
        logging_utils.write_centralized_logging(s2s.scrdir)
        sys.exit()

    # Delete forecast
    if args.delete_forecast:
        s2s.delete_forecast()
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

    # save S2S schedule dictionary
    logging_utils.save_schedule(s2s.scrdir, s2s.schedule)

    # Submit SLURM jobs
    # -----------------
    if  args.submit_job:
        s2s.submit_jobs()
