''' 
Python main script for HydroSFS
'''
import os
import sys
import glob
import subprocess
import shutil
import argparse
from datetime import datetime
from dateutil.relativedelta import relativedelta
from ghis2s.main.hcast_module import S2SHcast
from ghis2s.main.experiment_setup import S2Srun
from ghis2s.main import s2s_api, walltime
from ghis2s.shared import utils, logging_utils
from ghis2s.shared.logging_utils import TaskLogger

class HydroSFS(S2Srun):
    ''' HydroSFS class'''
    def __init__(self, month, config_file, pre_process=False, year=2025):
        super().__init__(year, month, config_file)
        mm = f'{month:02d}'
        if pre_process:
            self.scrdir = self.e2esdir + 'scratch/' + mm + '/'
            os.makedirs(self.scrdir + '/bcsd_fcst/logs', exist_ok=True)
        else:
            self.scrdir = self.e2esdir + 'scratch/' + self.yyyy + self.mm + '/'
            os.makedirs(self.scrdir + '/bcsd_fcst/logs', exist_ok=True)

    def monthly_hydroscs(self):
        ''' compute HydroSCS monthlies from hourly data '''
        # manage jobs from SCRATCH
        outdir='/discover/nobackup/projects/ghilis/smahanam/SERVIR/HydroSCS_monthly/'

        clim_syr = self.config['BCSD']['clim_start_year']
        clim_eyr = self.config['BCSD']['clim_end_year']

        os.chdir(self.scrdir + 'bcsd_fcst')
        self.create_symlink(self.e2esdir + 'bcsd_fcst/', 'bcsd_fcst')
        cwd=self.scrdir + 'bcsd_fcst'

        info = {}
        info['CPT'] = str(1)
        info['MEM']= '480GB'
        info['NT']= str(1)
        info['TPN'] = str(4)
        info['MP'] = True
        info['SKIP_ARG'] = True
        jobname = 'monthly_hydroscs_'
        cmdlist = []
        for year in range(clim_syr, clim_eyr + 1):
            for var in [1,2,3,4]:
                cmdlist.append(f"python {self.lisfdir}/lis/utils/usaf/S2S/ghis2s/bcsd/"
                               f"bcsd_library/calc_and_write_hydroscs_climatology.py {var} "
                               f"{self.e2esroot}/{self.config_file} {outdir} {year}")

        # write job file
        slurm_sub = self.split_list(cmdlist, 4)
        for i, sub_val in enumerate(slurm_sub):
            tfile = self.sublist_to_file(sub_val, cwd)
            try:
                s2s_api.python_job_file(f'{self.e2esroot}/{self.config_file}', f'{jobname}{i+1:02d}_run.j',
                                        jobname + f'{i+1:02d}_', str(4), str(6), cwd, tfile.name,
                                        parallel_run=info)
                self.create_dict(f'{jobname}{i+1:02d}_run.j', 'bcsd_fcst')
            finally:
                tfile.close()
                os.unlink(tfile.name)

    def clim_hydroscs(self):
        ''' Create HydoSCS Climatology '''
        lats, _ = utils.get_domain_info(self.e2esroot +'/' + self.config_file, coord=True)
        resol = f'{round((lats[1] - lats[0])*100)}km'

        # manage jobs from SCRATCH
        os.chdir(self.scrdir + 'bcsd_fcst')
        self.create_symlink(self.e2esdir + 'bcsd_fcst/', 'bcsd_fcst')
        cwd=self.scrdir + 'bcsd_fcst'
        outdir=f'{self.e2esdir}/bcsd_fcst/{self.obs_model}_{resol}/raw/Climatology/'
        os.makedirs(outdir, exist_ok=True)

        info = {}
        info['CPT'] = str(1)
        info['MEM']= '480GB'
        info['NT']= str(1)
        info['TPN'] = None
        info['MP'] = True

        jobname = 'clim_hydroscs_'
        cmdlist = []
        for var in ['LWdown', 'Rainf', 'Psurf', 'Qair',
                    'SWdown', 'Tair', 'Wind']:
            cmdlist.append(f"python {self.lisfdir}/lis/utils/usaf/S2S/ghis2s/bcsd/"
                           f"bcsd_library/calc_and_write_hydroscs_climatology.py {var} "
                           f"{self.e2esroot}/{self.config_file} {outdir}")

        # write job file
        slurm_sub = self.split_list(cmdlist, 1)
        for i, sub_val in enumerate(slurm_sub):
            tfile = self.sublist_to_file(sub_val, cwd)
            try:
                s2s_api.python_job_file(f'{self.e2esroot}/{self.config_file}', f'{jobname}{i+1:02d}_run.j',
                                        jobname + f'{i+1:02d}_', str(1), str(6), cwd, tfile.name,
                                        parallel_run=info)
                self.create_dict(f'{jobname}{i+1:02d}_run.j', 'bcsd_fcst')
            finally:
                tfile.close()
                os.unlink(tfile.name)

    def write_hydrosfs_jobs(self):
        """ writes job files for write_hydrosfs_files """
        jobname='hydrosfs_'
        par_info = {}
        par_info['CPT'] = '20'
        par_info['NT']= str(1)
        par_info['MEM']= '480GB'
        par_info['TPN'] = None
        par_info['MP'] = True
        par_info['SKIP_ARG'] = True
        prev = [f"{key}" for key in self.schedule.keys() if 'mf_tempdis_' in key]

        lats, _ = utils.get_domain_info(self.e2esroot +'/' + self.config_file, coord=True)
        resol = f'{round((lats[1] - lats[0])*100)}km'

        indir_template = '{}/{:04d}/ens{:01d}'
        outdir_template = '{}/{:04d}{:02d}/ens{:01d}'
        forcedir = f"{self.e2esdir}/bcsd_fcst/{self.fcst_model}_{resol}"
        subdaily_raw_fcst_dir = f"{forcedir}/bcsd/6-Hourly/{self._mmm}01"

        slurm_commands = []
        for ens in range(self.config['BCSD']['nof_raw_ens']):
            indir = indir_template.format(subdaily_raw_fcst_dir, self.year, ens+1)
            outdir = outdir_template.format(self.config['SETUP']['HYDROSFSDIR'], self.year, self.month, ens+1)
            os.makedirs(outdir, exist_ok=True)
            slurm_commands.append(f"python {self.e2esroot}/s2s_hydrosfs.py -m {self.month} -y {self.year} "
                                  f"-c {self.e2esroot}/{self.config_file} "
                                  f'-w "ens{ens+1:02d} {indir}/ {outdir}/"')

        # multi tasks per job
        l_sub = 5
        slurm_sub = self.split_list(slurm_commands, l_sub)
        for i, sub_val in enumerate(slurm_sub):
            tfile = self.sublist_to_file(sub_val, self._cwd)
            try:
                s2s_api.python_job_file(self.e2esroot +'/' + self.config_file,
                                        jobname + f'{i+1:02d}_run.j',
                                        jobname + f'{i+1:02d}_', 1, str(6), self._cwd, tfile.name,
                                        parallel_run=par_info)
                self.create_dict(jobname + f'{i+1:02d}_run.j', 'bcsd_fcst', prev=prev)
            finally:
                tfile.close()
                os.unlink(tfile.name)

            shutil.copy(jobname + f'{i+1:02d}_run.j', jobname + f'{i+1:02d}_run.sh')
            utils.remove_sbatch_lines(jobname + f'{i+1:02d}_run.sh')

    def write_hydrosfs_files(self, cmd):
        """ move BCSD output to HydroSFS directory """
        in_var_list = ["PRECTOT", "LWGAB", "SWGDN", "PS", "QV2M", "T2M", "U10M"]
        compress_encoding = {'dtype': 'int16', "zlib": True, "complevel": 6, "shuffle": True, "missing_value": -32767}

        # Split the command string into 3 variables
        try:
            parts = cmd.split()
            if len(parts) != 3:
                raise ValueError(f"Expected 3 space-delimited arguments, got {len(parts)}")
            ens, indir, outdir = parts
        except ValueError as e:
            print(f"Error parsing command: {e}")
            sys.exit(1)

        task_name = os.environ.get('SCRIPT_NAME')
        logger = TaskLogger(task_name,
                    os.getcwd(),
                    f's2s_app/s2s_hydrosfs.py copying ens{ens}')
        for in_var in in_var_list:
            file_pattern = f'{indir}/{in_var}.*.nc4'
            file_list = sorted(glob.glob(file_pattern))
            filenames = [os.path.basename(f) for f in file_list]
            subtask = f'ens{ens}_{in_var}'
            for file in filenames:
                ds = utils.load_ncdata(indir + file, [logger, subtask],
                                       chunks={'time': 'auto', 'lat': 'auto', 'lon': 'auto'})
                # Drop float specific attributes
                ds[in_var].attrs.pop('_FillValue', None)
                ds[in_var].attrs.pop('least_significant_digit', None)
                vname = in_var
                outfile = file
                if in_var == "U10M":
                    vname = "WIND10M"
                    ds = ds.rename({'U10M': 'WIND10M'})
                    date_str = file.split('.')[1]
                    outfile = f'{vname}.{date_str}.nc4'

                packed_ds, packing_params = utils.pack_dataset_to_int16(ds, [vname], logger=[logger, subtask])
                utils.write_ncfile(packed_ds, outdir+outfile, {vname: compress_encoding}, [logger, subtask])
                utils.add_packing_attributes(outdir+outfile, packing_params, [logger, subtask])

            # create the link for the next month
            date_str = outfile.split('.')[1]
            current_date = datetime.strptime(date_str, '%Y%m')
            next_month = current_date + relativedelta(months=1)
            next_mon_str = next_month.strftime('%Y%m')
            nextfile = f"{vname}.{next_mon_str}.nc4"
            cwd = os.getcwd()
            os.chdir(outdir)
            cmd = f"ln -sfn  {outfile} {nextfile}"
            rc = subprocess.call(cmd, shell=True)
            os.chdir(cwd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, type=str, help='config file')
    parser.add_argument('-y', '--year', required=False, type=int, help='forecast year')
    parser.add_argument('-m', '--month', required=True, type=int, help='forecast month')
    parser.add_argument('-p', '--preprocess', required=False, default=None, type=str,
                        help='preprocess options: regrid, clim, monthly')
    parser.add_argument('-w', '--write_hydrosfs', required=False, type=str, default=None,
                        help='write hydrosfs files command (3 space-delimited ens indir outdir)')
    parser.add_argument('-j', '--submit_job', action='store_true',
                        help='Submit SLURM jobs (default: False)?')
    parser.add_argument('-r', '--report', action='store_true',
                        help='Print report')
    parser.add_argument('-l', '--logging', action='store_true',
                        help='Write centralized log file')

    args = parser.parse_args()

    if  args.preprocess is not None:
        s2s = HydroSFS(month=args.month, config_file=args.config_file, pre_process=True)
        hcst = S2SHcast(month=args.month, config_file=args.config_file, pre_process=True)
    else:
        s2s = HydroSFS(month=args.month, config_file=args.config_file, year=args.year)

    # Print job report
    if  args.report:
        walltime.slurm(s2s.scrdir)
        sys.exit()

    # Write LOG file
    if args.logging:
        logging_utils.write_centralized_logging(s2s.scrdir)
        sys.exit()

    # Write HydroSFS files standalone
    if args.write_hydrosfs is not None:
        s2s.write_hydrosfs_files(args.write_hydrosfs)
        sys.exit()

    if  args.preprocess is not None:
        if args.preprocess.lower() == 'regrid':
            s2s.bcsd(clim=True)
        elif args.preprocess.lower() == 'monthly':
            s2s.monthly_hydroscs()
        elif args.preprocess.lower() == 'clim':
            hcst.clim_metforce()
            if not glob.glob(f'{s2s.e2esroot}/hindcast/bcsd_fcst/{s2s.obs_model}_*/raw/Climatology/PRECTOT_obs_clim.nc'):
                s2s.clim_hydroscs()
            s2s.schedule.update(hcst.schedule)
        else:
            print(f"Invalid preprocessor {args.preprocess}")
            sys.exit()
    else:
        s2s.bcsd()
        s2s.write_hydrosfs_jobs()

    # save S2S schedule dictionary
    logging_utils.save_schedule(s2s.scrdir, s2s.schedule)

    # Submit SLURM jobs
    # -----------------
    if  args.submit_job:
        s2s.write_log_monitoring_script()
        s2s.submit_jobs()
