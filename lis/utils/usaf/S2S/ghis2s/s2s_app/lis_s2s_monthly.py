#!/usr/bin/env python3
"""
##########################################################################################
#
# NAME: lis_s2s_monthly.py
#
# PURPOSE
# =======
# Prepares the environment and executes the monthly Subseasonal-to-seasonal (S2S) component of LIS
#
#
# USAGE
# =====
# lis_s2s_monthly.py [-m yyyymm] [-r] [-e] [-d]
#
#     where:
#
#         yyyymm (optional) is the month for which S2S will be run (default: current month)
#         report (optional) will simply provide a status update on the S2S run
#         email  (optional) will send an email with the results of the report
#         delete (optional) will delete old S2S data
#
#     Notes: -r can be used without -e (will print/log a report but not send an email)
#            If -e is used, script will force -r to be True
#
#
# Updates:
# ========
# 20230906 Initial version ..................................... Eddy Hildebrand/SAIC/TADS
#
##########################################################################################
"""

#-----------------------------------------------------------------------------------------
# Import section
#-----------------------------------------------------------------------------------------
import copy
import glob
import logging
import os
import sys
import shutil
import subprocess
from datetime import datetime as dt
from datetime import timedelta
import ymdate
import specm.restart_program
from util import copy_and_flag
sys.path.append("/ccs/proj/nwp303/green/lis_nasa/lis_nasa-7.5.13/lis/utils/usaf/S2S/")
from ghis2s.s2s_app import s2s_api 

#-----------------------------------------------------------------------------------------
# Set up the configuration definition
#-----------------------------------------------------------------------------------------
_CONFIG_DEFINITION = {
  "default": os.path.join(os.path.dirname(__file__), "../cfg/lis_s2s_monthly.cfg"),
  "pos_args": [
  ],
  "opt_args": [
    ["m", "month", "S2S month to run", str],
    ["r", "report", "Status check for S2S run", bool],
    ["e", "email", "Email the status check results", bool],
    ["d", "delete", "Delete old month of data", bool],
    ],
  "config": {
    "batch_options": {
      "batch_system": "slurm",
      "program": "lis_s2s_monthly",
      "project": "nwp303",
      "queue": "development",
      "cluster_constraint": "green",
      "time_limit": 120,
      "nodes": 1,
    },
    "dirs": {
      "BIN": "/ccs/proj/nwp303/lis-local/lis-local-master/bin",
      "CFG": "/ccs/proj/nwp303/lis-local/lis-local-master/cfg",
      "DATA": "/lustre/active/nwp303/proj-shared/data/lis-local/7.5",
      "WORK": "/lustre/active/nwp303/proj-shared/data/lis-local/7.5/work/s2s_monthly",
      "CTRL": "/lustre/active/nwp303/proj-shared/data/lis-local/7.5/ctrl_times",
    },
    "misc": {
      "category": "CRITICAL",
      "service": "LIS",
      "contacts": "/ccs/proj/nwp303/lis-local/lis-local-master/cfg/specm_email.cfg",
      "log_dir": "/lustre/active/nwp303/proj-shared/data/lis-local/7.5/log",
      "log_level": "INFO",
      "splunk_cp": "1",
      "splunk_dir": "/lustre/active/nwp602/pdata/557WW_Outbound/test/dropbox",
    },
    "lis": {
      "lis_nasa": "/ccs/proj/nwp303/green/lis_nasa/lis_nasa-7.5.13",
      "s2s_app": "/ccs/proj/nwp303/green/lis_nasa/lis_nasa-7.5.13/lis/utils/usaf/S2S/s2s_app",
      "hindcast": "/lustre/active/nwp303/world-shared/LIS75/AdHOC/E2ES_DIR/hindcast/",
      "SMAP_L2_SM_NRT": "/lustre/active/nwp303/proj-shared/data/lis-local/7.5/SMAP",
      "usaf_lis75s2s_gfs2galwem": "/lustre/active/nwp303/proj-shared/data/lis-local/7.5/S2S",
      "retention_time": "180",
      "WS_copy": "1",
      "WS_dir": "/lustre/active/nwp303/world-shared/LIS75/S2S",
      "Out557_copy": "1",
      "Out557_dir": "/lustre/active/nwp303/proj-shared/outbound/inbound",
    }
  }
}

class LisS2sMonthly(specm.restart_program.RestartProgram):
    """Run lis_s2s_monthly"""

    #---------------------------------------------------------------------------------------
    # init method - set up class variables
    #---------------------------------------------------------------------------------------
    def __init__(self, name, config_definition = None):
        """Initialize the process"""
        if config_definition is None:
            config_definition = copy.deepcopy(_CONFIG_DEFINITION)
        super().__init__(name, config_definition)

        self.config.log()
        self.cfgdir  = self.config.get("dirs", "CFG")
        self.logdir  = self.config.get("misc", "log_dir")
        self.workdir = self.config.get("dirs", "WORK")
        self.ctrldir = self.config.get("dirs", "CTRL")

    #---------------------------------------------------------------------------------------
    # run method - runs the program
    #---------------------------------------------------------------------------------------
    def _run(self):

        #------------------------------------------------------------------------
        # Determine if this is a real S2S run or just a status report or deletion
        #------------------------------------------------------------------------
        report = self.config.get("user", "report")
        email = self.config.get("user", "email")
        delete = self.config.get("user", "delete")
        if email and not report:
            logging.info("Email option was specified.  Setting report to True")
            report = True

        #--------------------------------------------------------
        # Check beginning date
        # Value defined by the caller will be used if it is valid
        # Value not defined by caller will be set automatically
        #--------------------------------------------------------
        yyyymm = self.check_begdtg(self.config.get("user", "month"), delete)

        #------------------------------------
        # Log the arguments that will be used
        #------------------------------------
        logging.info("Calling Arguments")
        logging.info("=================")
        logging.info("Month: %s", yyyymm)
        logging.info("Report: %s", report)
        logging.info("Email: %s", email)
        logging.info("Delete: %s", delete)

        #--------------------------------
        # Define log and work directories
        #--------------------------------
        os.makedirs(self.logdir, 0o755, exist_ok=True)
        os.makedirs(os.path.join(self.workdir,'E2ES'), 0o755, exist_ok=True)
        os.chdir(self.workdir)

        #------------------------------
        # If we are doing a real run...
        #------------------------------
        if not report and not delete:

            #----------------
            # Set up symlinks
            #----------------
            self.link_setup()

            #-----------------------------
            # Copy files to work directory
            #-----------------------------
            logging.info('Copying LIS NASA files to working directory')
            shutil.copytree(self.config.get("lis", "lis_nasa"),
                            os.path.join(self.workdir,'lis_nasa'), dirs_exist_ok=True)

            #------------------------
            # Get the S2S config file
            #------------------------
            self.cfg_setup()

        #--------
        # Run S2S
        #--------
        output = self.run_s2s(yyyymm, report, delete)

        #------------------------------------------------
        # Copy output to world-shared and 557 WW Outbound
        #------------------------------------------------
        if report:

            #----------------
            # Generate report
            #----------------
            self.generate_report(yyyymm, output, email)

            #--------------------------------------------------------
            # Get control times file to see if we already copied this
            # month's output to world-shared and 557WW_Outbound
            #--------------------------------------------------------
            with open(os.path.join(self.ctrldir, 'lis.control.times.S2S'), "r",
                encoding="utf-8") as fil:
                lines = fil.readlines()
            last_copy = None
            for line in lines:
                line = line.strip().split()
                if "last_s2s_output_copy" in line:
                    last_copy = line[0]
                    break
            logging.info('last_s2s_output_copy from lis.control.times.S2S file is %s', last_copy)

            #----------------------
            # If S2S run is done...
            #----------------------
            if 'ELAPSED TIME' in output:

                #-------------------------------------------
                # ...and we haven't copied the output yet...
                #-------------------------------------------
                if int(last_copy) < int(yyyymm):

                    #-----------------------------------------------------------
                    # ...then copy the output to world-shared and 557WW_Outbound
                    # Note: Plots copy intentionally commented out
                    #-----------------------------------------------------------
                    logging.info('S2S run appears to be complete.  Calling copy_output')
                    self.copy_output(yyyymm, 'metric')   # about 9 GB per month
                    self.copy_output(yyyymm, 'post')     # north of 700 GB per month!
                    # self.copy_output(yyyymm, 'plots')  # about 200 MB per month

                    #-----------------------------------------------------------------------------
                    # Update control times file to indicate we finished copying this month's files
                    #-----------------------------------------------------------------------------
                    logging.info('Updating lis.control.times.S2S')
                    with open(os.path.join(self.ctrldir, 'lis.control.times.S2S'), "w",
                        encoding="utf-8") as fil:
                        fil.write(f'{yyyymm} last_s2s_output_copy\n')
                else:
                    logging.info('%s %s %s %s %s',
                                 'S2S output files for this month have already been copied to',
                                 'world-shared and 557WW Outbound.  Not recopying as it is a',
                                 'large amount of data.  If you would like to recopy the output,',
                                 'manually set lis.control.times.S2S to the previous month and',
                                 f'rerun: lis_s2s_monthly.py -m {yyyymm} -r')
            else:
                logging.info('S2S run not yet complete')

        #---------------------------------
        # Copy log files to Splunk dropbox
        #---------------------------------
        self.splunk_copy()

        #----------
        # Finish up
        #----------
        self.logger.finish()


    #------------------
    # link_setup method
    #------------------
    def link_setup(self):
        """
        Create symlinks
        """
        os.chdir(os.path.join(self.workdir,'E2ES'))
        for link_name in ["s2s_app", "hindcast"]:
            try:
                os.unlink(link_name)
            except OSError:
                pass
            os.symlink(self.config.get("lis", link_name), link_name)

        os.makedirs(os.path.join(self.workdir, 'USAF_FORCING'), mode=0o775, exist_ok=True)
        os.chdir(os.path.join(self.workdir, 'USAF_FORCING'))
        for link_name in ["SMAP_L2_SM_NRT", "usaf_lis75s2s_gfs2galwem"]:
            try:
                os.unlink(link_name)
            except OSError:
                pass
            os.symlink(self.config.get("lis", link_name), link_name)

    #-----------------
    # cfg_setup method
    #-----------------
    def cfg_setup(self):
        """
        Copy the shell config to the work directory
        """
        os.chdir(os.path.join(self.workdir,'E2ES'))
        shell_cfg = os.path.join(self.cfgdir,"lis.config_SHELL_s2s_monthly")
        logging.info('Copying LIS shell config to work directory')
        shutil.copy(shell_cfg, "s2s_config_global_fcast")

    #---------------
    # run_s2s method
    #---------------
    def run_s2s(self, yyyymm, report, delete):
        """
        Call NASA's shell script to initiate an S2S run
        Example: ./s2s_app/s2s_run.sh -y 2023 -m 5 -c s2s_config_global_fcast [-r Y] [-d Y]
            -y: 4 digit year
            -m: 1 or 2 digit month
            -c: Name of S2S configuration file
            -r: Report (Y/N) This optional argument shows a progress report on a currently
                running or completed S2S run.  Do not use the -r option when initially
                starting a monthly S2S run.
            -d: Delete an old month of data
        """
        os.chdir(os.path.join(self.workdir, 'E2ES'))

        # If month < 10, we have to use a single digit month
        month = yyyymm[-1] if yyyymm[-2] == '0' else yyyymm[-2:]

        # initialize s2s_run instance
        s2s_run = s2s_api.S2SRun(config_file=args.config_file, model_type=args.model_type)
        s2s_run.main()
        
        # Build the command
        cmd = f'./s2s_app/s2s_run.sh -y {yyyymm[0:4]} -m {month} -c s2s_config_global_fcast'
        if report:
            cmd += ' -r Y'
        if delete:
            cmd += ' -d Y'

        # Call the shell script
        logging.info('Calling S2S run shell script for %s', yyyymm)
        logging.info('Command is: %s', cmd)
        output = subprocess.run(cmd, shell=True, check=False, capture_output=True, text=True).stdout

        return output

    #-----------------------
    # generate_report method
    #-----------------------
    def generate_report(self, yyyymm, output, email):
        """ Generate the report """

        # Check SLURM_JOB_SCHEDULE
        failed_jobs = self.check_slurm_jobs(yyyymm)

        # If we found any failed jobs, force email and abort to true
        abort = False
        if len(failed_jobs) > 0:
            email = True
            abort = True
            workdirs = []
            for failed_job in failed_jobs:
                with subprocess.Popen([f'sacct -n --format WorkDir%150 -j {failed_job}'],
                    shell=True, stdout=subprocess.PIPE) as proc:
                    retcode = proc.communicate()
                workdirs.append(retcode[0].decode('UTF-8').strip())

        # Print and email handling
        if not email:
            logging.info('Status report:\n%s', output)
            logging.info('Found 0 failed Slurm jobs')
        else:
            email_body = "\n"
            email_body += "".join(output)
            if not abort:
                logging.info("Emailing the status report")
                email_body += "\nFound 0 failed jobs"
                self.warn(f'LIS S2S Status Report for {yyyymm} Monthly Run', email_body)
            else:
                email_body += f"ERROR: FAILED SLURM JOBS: {failed_jobs}"
                email_body += f"\nWork Dirs for failed jobs: {workdirs}"
                logging.error("Status report: %s", email_body)
                self.errhnd.abort(f'LIS S2S Status Report for {yyyymm} Monthly Run')

    #------------------------
    # check_slurm_jobs method
    #------------------------
    def check_slurm_jobs(self, yyyymm):
        """Check Slurm job schedule for failed jobs"""
        with open(os.path.join(self.workdir,"E2ES","scratch",yyyymm,"SLURM_JOB_SCHEDULE"), "r",
            encoding="utf-8") as fil:
            lines = fil.readlines()
        failed_jobs = []
        for line in lines:
            if line.strip() == '':
                continue # skip blank lines
            jobid = line.strip().split()[0]
            if jobid.isnumeric() is False:
                continue # skip header lines
            with subprocess.Popen([f'sacct -n --format state -j {jobid}'], shell=True,
                stdout=subprocess.PIPE) as proc:
                retcode = proc.communicate()
            status = retcode[0].decode('UTF-8').strip().split('\n',maxsplit=1)[0].split()[-1]
            if status == 'FAILED':
                failed_jobs.append(jobid)
        return failed_jobs

    #-------------------
    # copy_output method
    #-------------------
    def copy_output(self, yyyymm, s2stype):
        """Copy S2S output to world-shared and/or 557WW Outbound"""

        # Get a list of files
        if s2stype == 'metric':
            files = sorted(glob.glob(os.path.join(self.workdir, 'E2ES', 's2smetric', yyyymm,
                                                  'metrics_cf', 'NOAHMP', 'PS.557WW*.*')))
        elif s2stype == 'post':
            files = sorted(glob.glob(os.path.join(self.workdir, 'E2ES', 's2spost', yyyymm,
                                                  '*', 'PS.557WW*.NC')))
        elif s2stype == 'plots':
            files = sorted(glob.glob(os.path.join(self.workdir, 'E2ES', 's2splots', yyyymm,
                                                  'NOAHMP', '*.png')))
        else:
            logging.info('%s is an unrecognized data type.  Not copying any files.', s2stype)
            return

        if len(files) == 0:
            logging.info('Could not find any S2S %s output files to copy', s2stype)
            return

        # Copy to world-shared
        if int(self.config.get("lis", "WS_copy")) == 1:
            wshared_dir = os.path.join(self.config.get("lis", "WS_dir"), yyyymm)
            os.makedirs(wshared_dir, mode=0o775, exist_ok=True)

            logging.info('Copying S2S %s files to %s', s2stype, wshared_dir)
            for file in files:
                logging.info('Copying %s to %s', file, wshared_dir)
                try:
                    shutil.copy(file, wshared_dir)
                except OSError as err:
                    raise OSError('ERROR: Could not copy S2S file to world-shared') from err

        # Copy to 557WW_Outbound
        if int(self.config.get("lis", "Out557_copy")) == 1:
            out557_dir = self.config.get("lis", "Out557_dir")

            logging.info('Copying S2S %s files to %s', s2stype, out557_dir)
            for file in files:
                logging.info('Copying %s to %s', file, out557_dir)
                try:
                    copy_and_flag(file, out557_dir)
                except OSError as err:
                    raise OSError('ERROR: Could not copy S2S file to 557WW_Outbound') from err

    #-------------------
    # splunk_copy method
    #-------------------
    def splunk_copy(self):
        """Copy log and err files to Splunk dropbox"""
        if int(self.config.get("misc", "splunk_cp")) == 1:

            os.chdir(self.logdir)
            splunk_dir = self.config.get("misc", "splunk_dir")
            logging.info('Done with lis_s2s_monthly.  Copying log files to Splunk dropbox')

            # Get the log file from lis_s2s_monthly.py
            logfile = sorted(glob.glob(os.path.join('lis_s2s_monthly*'+str(os.getpid())+'.log')))
            if len(logfile) > 0:
                copy_and_flag(logfile[0], splunk_dir)
            else:
                logging.info('Could not find %s.  Skipping copy to Splunk dropbox', logfile)

            # Get the err file from lis_s2s_monthly.py
            errfile = sorted(glob.glob(os.path.join('lis_s2s_monthly*'+str(os.getpid())+'.err')))
            if len(errfile) > 0:
                copy_and_flag(errfile[0], splunk_dir)
            else:
                logging.info('Could not find %s.  Skipping copy to Splunk dropbox', errfile)
        else:
            logging.info('Copying log files to Splunk dropbox is disabled')


    #--------------------
    # check_begdtg method
    #--------------------
    def check_begdtg(self, beg_dtg, delete):
        """
        If a beg date (yyyymm) was specified as a calling argument, make sure it is valid
        """
        if beg_dtg is not None:
            logging.info('Beg date provided via calling argument. Determining if it is valid')
            if ymdate.is_valid(beg_dtg+'01') and len(beg_dtg) == 6:
                logging.info('yyyymm appears to be valid')
            else:
                raise ValueError('ERROR: Invalid beg date - must be in YYYYMM format')
        else:
            if delete:
                backup = self.config.get("lis", "retention_time")
                logging.info('Delete option specified. Backing up %s days', backup)
                beg_dtg = (dt.now()-timedelta(days=int(backup))).strftime('%Y%m')
            else:
                logging.info('No yyyymm provided via calling arg.  Using current month.')
                beg_dtg = (dt.now()-timedelta(days=1)).strftime('%Y%m')
        logging.info('yyyymm is %s', beg_dtg)
        return beg_dtg


#-----------------------------------------------------------------------------------------
# main method
#-----------------------------------------------------------------------------------------
def main():
    """Define the application and run it"""
    app = LisS2sMonthly(__file__, _CONFIG_DEFINITION)

    try:
        app.run()
    except(OSError, ValueError) as err:
        app.error(err)

if __name__ == "__main__":
    main()
