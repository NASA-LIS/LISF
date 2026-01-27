''' 
Python main script for HydroSFS
'''
import os
import sys
import argparse
from ghis2s.main.experiment_setup import S2Srun
from ghis2s.main import s2s_api, walltime
from ghis2s.shared import utils, logging_utils

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
        outdir=f'{self.e2esdir}/bcsd_fcst/hydroscs_{resol}/raw/Climatology/'
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, type=str, help='config file')
    parser.add_argument('-y', '--year', required=False, type=int, help='forecast year')
    parser.add_argument('-m', '--month', required=True, type=int, help='forecast month')
    parser.add_argument('-p', '--preprocess', required=False, default=None, type=str,
                        help='preprocess options: regrid, clim, monthly')
    parser.add_argument('-j', '--submit_job', action='store_true',
                        help='Submit SLURM jobs (default: False)?')
    parser.add_argument('-r', '--report', action='store_true',
                        help='Print report')
    parser.add_argument('-l', '--logging', action='store_true',
                        help='Write centralized log file')

    args = parser.parse_args()

    if  args.preprocess is not None:
        s2s = HydroSFS(month=args.month, config_file=args.config_file, pre_process=True)
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

    if  args.preprocess is not None:
        if args.preprocess.lower() == 'regrid':
            s2s.bcsd(clim=True)
        elif args.preprocess.lower() == 'monthly':
            s2s.monthly_hydroscs()
        elif args.preprocess.lower() == 'clim':
            #hcast.clim_metforce()
            s2s.clim_hydroscs()
        else:
            print(f"Invalid preprocessor {args.preprocess}")
            sys.exit()
        s2s.schedule.update(hcast.schedule)
    else:
        s2s.bcsd()

    # save S2S schedule dictionary
    logging_utils.save_schedule(s2s.scrdir, s2s.schedule)

    # Submit SLURM jobs
    # -----------------
    if  args.submit_job:
        s2s.submit_jobs()
