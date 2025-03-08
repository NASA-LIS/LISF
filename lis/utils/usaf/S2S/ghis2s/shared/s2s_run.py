import os
import sys
from datetime import datetime
from dateutil.relativedelta import relativedelta
import subprocess
import platform
import shutil
import yaml
from ghis2s.s2s_app import s2s_api

class S2Srun():
    def __init__(self, year, month, config_file):
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
        self.NODE_NAME = platform.node()
        self.CONSTRAINT = self.config['SETUP']['CONSTRAINT']
        
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
        command = f"bash {self.E2ESDIR}/s2s_app/wget_cfsv2_oper_ts_e2es.sh -y {self.YYYY} -m {self.MM} -c {self.BWD}/{self.config_file} -d N"
        process = subprocess.run(command, shell=True)
        ret_code = process.returncode

        if ret_code > 0:
            print(f"Error return code from the CFSv2 file download checker :: {ret_code}")
            print("> 0 :: Exiting from s2s_run.py --")
            sys.exit(ret_code)
        else:
            print('SUCCESS')

        return

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
        
        # submit Job
        # ----------

        os.chdir(self.BWD)
        
        return

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

        # submit job
        # ----------

        os.chdir(self.BWD)
        
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
        
    
    def main(self):
        # (1) Run CFSV2 file checker to ensure downloaded files are not corrupted/
        #self.CFSv2_file_checker()

        # (2) LISDA run
        #self.lis_darun()

        # (3) LDT-ICS
        #self.ldt_ics()

        # (4) BCSD
        self.bcsd()

        return
