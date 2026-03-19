''' 
Python main script for hindcast processing
'''
import os
import sys
import argparse
from datetime import datetime
from ghis2s.s2s_app.s2s_run import S2Srun
from ghis2s.s2s_app import s2s_api, walltime
from ghis2s.shared import utils, logging_utils

class S2SHcast(S2Srun):
    ''' Hind cast processing class'''
    def __init__(self, month, config_file, pre_process=False, year=2025):
        super().__init__(year, month, config_file)
        mm = f'{month:02d}'
        if pre_process:
            self.scrdir = self.e2esdir + 'scratch/' + mm + '/'
            os.makedirs(self.scrdir + '/bcsd_fcst/logs', exist_ok=True)
        else:
            self.scrdir = self.e2esdir + 'scratch/' + self.yyyy + self.mm + '/'
            os.makedirs(self.scrdir + '/bcsd_fcst/logs', exist_ok=True)
            os.makedirs(self.scrdir + '/lis_fcst/logs', exist_ok=True)
            os.makedirs(self.scrdir + '/s2spost/logs', exist_ok=True)

    def regrid_cfsv2(self):
        '''CFSv2 regridding'''
        lats, _ = utils.get_domain_info(self.e2esroot +'/' + self.config_file, coord=True)
        resol = f'{round((lats[1] - lats[0])*100)}km'

        fcast_clim_dir=f'{self.e2esdir}/bcsd_fcst/{self.fcst_model}_{resol}/raw/Climatology/'
        os.makedirs(fcast_clim_dir, exist_ok=True)

        # manage jobs from SCRATCH
        os.chdir(self.scrdir + 'bcsd_fcst')
        self.create_symlink(self.e2esdir + 'bcsd_fcst/', 'bcsd_fcst')
        cwd=self.scrdir + 'bcsd_fcst'
        date_obj = datetime.strptime(f"{self.yyyy}-{self.mm}-01", "%Y-%m-%d")
        mmm = date_obj.strftime("%b").lower()

        clim_syr = self.config['BCSD']['clim_start_year']
        clim_eyr = self.config['BCSD']['clim_end_year']
        self.bcsd_mf_regrid(cwd, mmm, resol, clim_syr, clim_eyr=clim_eyr)

    def clim_cfsv2(self):
        ''' create CFSv2 climatology '''
        lats, _ = utils.get_domain_info(self.e2esroot +'/' + self.config_file, coord=True)
        resol = f'{round((lats[1] - lats[0])*100)}km'

        fcast_clim_dir=f'{self.e2esdir}/bcsd_fcst/{self.fcst_model}_{resol}/raw/Climatology/'
        os.makedirs(fcast_clim_dir, exist_ok=True)

        # manage jobs from SCRATCH
        os.chdir(self.scrdir + 'bcsd_fcst')
        self.create_symlink(self.e2esdir + 'bcsd_fcst/', 'bcsd_fcst')
        cwd=self.scrdir + 'bcsd_fcst'
        date_obj = datetime.strptime(f"{self.yyyy}-{self.mm}-01", "%Y-%m-%d")
        mmm = date_obj.strftime("%b").lower()

        fcst_indir=f'{self.e2esdir}/bcsd_fcst/{self.fcst_model}_{resol}/raw/Monthly/'
        outdir=f'{self.e2esdir}/bcsd_fcst/{self.fcst_model}_{resol}/raw/Climatology/{mmm}01'
        os.makedirs(outdir, exist_ok=True)

        jobname = f'clim_{self.fcst_model}_'
        cmdlist = []
        info = {}
        info['CPT'] = str(1)
        info['MEM']= '480GB'
        info['NT']= str(1)
        info['TPN'] = None
        info['MP'] = True
        for var in ['PRECTOT', 'LWGAB', 'PS', 'QV2M', 'SWGDN', 'T2M', 'WIND10M']:
            cmdlist.append(f"python {self.lisfdir}/lis/utils/usaf/S2S/ghis2s/bcsd/"
                           f"bcsd_library/calc_and_write_forecast_climatology.py {var} "
                           f"{self.mm} {self.e2esroot}/{self.config_file} {fcst_indir} {outdir}")

        # write job file
        slurm_sub = self.split_list(cmdlist, 1)
        for i, sub_val in enumerate(slurm_sub):
            tfile = self.sublist_to_file(sub_val, cwd)
            try:
                s2s_api.python_job_file(f'{self.e2esroot}/{self.config_file}',
                                        f'{jobname}{i+1:02d}_run.j',
                                        jobname + f'{i+1:02d}_', str(1), str(3), cwd,
                                        tfile.name, parallel_run=info)
                self.create_dict(f'{jobname}{i+1:02d}_run.j', 'bcsd_fcst')
            finally:
                tfile.close()
                os.unlink(tfile.name)

    def regrid_nmme(self):
        '''NMME regridding'''
        nmme_output_dir=f'{self.e2esdir}/bcsd_fcst/NMME/raw/Monthly/'
        os.makedirs(nmme_output_dir, exist_ok=True)
        os.chdir(self.scrdir + 'bcsd_fcst')
        self.create_symlink(self.e2esdir + 'bcsd_fcst/', 'bcsd_fcst')
        cwd=self.scrdir + 'bcsd_fcst'
        self.bcsd_pr_regrid(cwd)

    def clim_nmme(self):
        ''' NMME climatology '''
        # manage jobs from SCRATCH
        os.chdir(self.scrdir + 'bcsd_fcst')
        self.create_symlink(self.e2esdir + 'bcsd_fcst/', 'bcsd_fcst')
        cwd=self.scrdir + 'bcsd_fcst'
        date_obj = datetime.strptime(f"{self.yyyy}-{self.mm}-01", "%Y-%m-%d")
        mmm = date_obj.strftime("%b").lower()

        climdir=f'{self.e2esdir}/bcsd_fcst/NMME/raw/Climatology/{mmm}01'
        [os.makedirs(f'{climdir}/{model}', exist_ok=True) for model in self.models]

        jobname = 'clim_nmme_'
        cmdlist = []
        info = {}
        info['CPT'] = str(1)
        info['MEM']= '480GB'
        info['NT']= str(1)
        info['TPN'] = None
        info['MP'] = True
        for model in self.models:
            nmme_indir=f'{self.e2esdir}/bcsd_fcst/NMME/raw/Monthly/{mmm}01/{model}/'
            outdir = f'{climdir}/{model}'
            cmdlist.append(f"python {self.lisfdir}/lis/utils/usaf/S2S/ghis2s/bcsd/"
                           f"bcsd_library/calc_and_write_nmme_forecast_climatology.py PRECTOT "
                           f"{self.mm} {model} {self.e2esroot}/{self.config_file} "
                           f"{nmme_indir} {outdir}")

        # write job file
        slurm_sub = self.split_list(cmdlist, 1)
        for i, sub_val in enumerate(slurm_sub):
            tfile = self.sublist_to_file(sub_val, cwd)
            try:
                s2s_api.python_job_file(f'{self.e2esroot}/{self.config_file}',
                                        f'{jobname}{i+1:02d}_run.j',
                                        jobname + f'{i+1:02d}_', str(1), str(3), cwd,
                                        tfile.name, parallel_run=info)
                self.create_dict(f'{jobname}{i+1:02d}_run.j', 'bcsd_fcst')
            finally:
                tfile.close()
                os.unlink(tfile.name)

    def clim_nafpa(self):
        ''' Create NAFPA Climatologies '''
        lats, _ = utils.get_domain_info(self.e2esroot +'/' + self.config_file, coord=True)
        resol = f'{round((lats[1] - lats[0])*100)}km'

        # manage jobs from SCRATCH
        os.chdir(self.scrdir + 'bcsd_fcst')
        self.create_symlink(self.e2esdir + 'bcsd_fcst/', 'bcsd_fcst')
        cwd=self.scrdir + 'bcsd_fcst'
        outdir=f'{self.e2esdir}/bcsd_fcst/USAF-LIS7.3rc8_{resol}/raw/Climatology/'
        os.makedirs(outdir, exist_ok=True)

        jobname = 'clim_nafpa_'
        if resol != '25km':
            ''' in the absense of a LIS simulation temporarily regridding '''
            info = {}
            info['CPT'] = str(1)
            info['MEM']= '480GB'
            info['NT']= str(1)
            info['TPN'] = None
            info['MP'] = True
            cmd = [f"python {self.lisfdir}/lis/utils/usaf/S2S/ghis2s/bcsd/"
                   f"regrid2fine.py {resol}"]
            tfile = self.sublist_to_file(cmd, cwd)
            try:
                s2s_api.python_job_file(f'{self.e2esroot}/{self.config_file}', f'{jobname}run.j',
                                        jobname, str(1), str(6), cwd, tfile.name, parallel_run=info)
                self.create_dict(jobname + 'run.j', 'bcsd_fcst')
            finally:
                tfile.close()
                os.unlink(tfile.name)
            return

        if os.path.exists(outdir + 'PRECTOT_obs_clim.nc'):
            print ("NAFPA Climatology files exist !")
            return

        cmdlist = []
        for var in ['LWdown_f_tavg', 'Rainf_f_tavg', 'Psurf_f_tavg', 'Qair_f_tavg',
                    'SWdown_f_tavg', 'Tair_f_tavg', 'Wind_f_tavg']:
            cmdlist.append(f"python {self.lisfdir}/lis/utils/usaf/S2S/ghis2s/bcsd/"
                           f"bcsd_library/calc_and_write_observational_climatology.py {var} "
                           f"{self.e2esroot}/{self.config_file} {outdir}")

        # write job file
        tfile = self.sublist_to_file(cmdlist, cwd)
        try:
            s2s_api.python_job_file(f'{self.e2esroot}/{self.config_file}', f'{jobname}run.j',
                                    jobname, str(1), str(3), cwd, tfile.name)
            self.create_dict(jobname + 'run.j', 'bcsd_fcst')
        finally:
            tfile.close()
            os.unlink(tfile.name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, type=str, help='config file')
    parser.add_argument('-y', '--year', required=False, type=int, help='forecast year')
    parser.add_argument('-m', '--month', required=True, type=int, help='forecast month')
    parser.add_argument('-p', '--preprocess', required=False, default=None, type=str,
                        help='preprocess options: regrid, clim, stats')
    parser.add_argument('-j', '--submit_job', action='store_true',
                        help='Submit SLURM jobs (default: False)?')
    parser.add_argument('-r', '--report', action='store_true',
                        help='Print report')
    parser.add_argument('-l', '--logging', action='store_true',
                        help='Write centralized log file')

    args = parser.parse_args()
    if  args.preprocess is not None:
        s2s = S2SHcast(month=args.month, config_file=args.config_file, pre_process=True)
    else:
        s2s = S2SHcast(month=args.month, config_file=args.config_file, year=args.year)

    # Print job report
    if  args.report:
        JOB_SCHEDULE = os.path.join(s2s.scrdir, 'SLURM_JOB_SCHEDULE')
        if os.path.exists(JOB_SCHEDULE):
            walltime.slurm(s2s.scrdir)
        else:
            if args.step is None:
                CYLC_WORKFLOW = f'cylc_e2e_{s2s.yyyy}{s2s.mm}'
            else:
                CYLC_WORKFLOW = f'cylc_{args.step.lower()}_{s2s.yyyy}{s2s.mm}'
            walltime.cylc(CYLC_WORKFLOW)
        sys.exit()

    # Write LOG file
    if args.logging:
        logging_utils.write_centralized_logging(s2s.scrdir)
        sys.exit()

    if  args.preprocess is not None:
        if args.preprocess.lower() == 'regrid':
            s2s.regrid_cfsv2()
            s2s.regrid_nmme()
        elif args.preprocess.lower() == 'clim':
            s2s.clim_cfsv2()
            s2s.clim_nmme()
            s2s.clim_nafpa()
            s2s.write_cylc_snippet()
        elif args.preprocess.lower() == 'stats':
            s2s.s2smetric()
        else:
            print(f"Invalid preprocessor {args.preprocess}")
            sys.exit()
    else:
        s2s.bcsd()
        s2s.lis_fcst()
        s2s.s2spost()
        s2s.write_cylc_snippet()

    # save S2S schedule dictionary
    logging_utils.save_schedule(s2s.scrdir, s2s.schedule)

    # Submit SLURM jobs
    # -----------------
    if  args.submit_job:
        s2s.submit_jobs()
