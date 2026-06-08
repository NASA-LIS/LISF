'''
HYMAP dobdomain decomposition
'''
import glob
import os
import re
import platform
import shutil
from datetime import datetime, date
from dateutil.relativedelta import relativedelta
import numpy as np
import xarray as xr
import dask
from dask.diagnostics import ProgressBar
from ghis2s.main import s2s_api
from ghis2s.shared import utils
from ghis2s.lis_fcst import generate_lis_config_scriptfiles_fcst
from ghis2s.shared.logging_utils import TaskLogger

class HyMAPDomains():
    ''' Contains methods to decompose and stitch HyMAP domains '''
    def __init__(self, subdomain_path):
        self.domains = [
            'africa_1',
            'africa_2',
            'africa_3',
            'africa_4',
            'asia_1',
            'asia_4',
            'asia_5',
            'australia',
            'central_america',
            'europe_1',
            'europe_2',
            'greenland',
            'hawaii',
            'north_america_1',
            'north_america_2',
            'north_america_3',
            'russia_tip',
            'siberia_1',
            'siberia_2',
            'siberia_3',
            'south_america_1',
            'south_america_2',
            'south_asia']

        self.subdomain_path = subdomain_path
        self.domain_masks = self.read_file(self.subdomain_path + 'HyMAP_sub-domains.nc',
                                           var_name='HyMAP_sub-domains')

    def read_file(self, file, var_name=None):
        ''' reads mask file '''
        mask_xr = xr.open_dataset(file)
        domain_masks = mask_xr.rename({'north_south': 'lat', 'east_west': 'lon'})
        domain_masks = domain_masks.assign_coords(lat=mask_xr['lat'].values[:,0],
                                                  lon=mask_xr['lon'].values[0,:])
        if var_name is not None:
            domain_masks = domain_masks.sel(time=0)['HyMAP_sub-domains']
        return domain_masks

    def lis_darun(self, setup, prev=None):
        '''lis_darun on subdomains'''
        # previous month
        date_obj = datetime.strptime(f"{setup.yyyy}-{setup.mm}-01", "%Y-%m-%d")
        previous_month_date = date_obj - relativedelta(months=1)
        yyyymmp = previous_month_date.strftime("%Y%m")
        yyyyp = yyyymmp[:4]
        mmp = yyyymmp[4:6]
        monp = int(mmp)
        pertmode = setup.config['EXP']['pertmode']

        # create subdomain job irectories
        for domain in self.domains:
            ldtfile = f'{self.subdomain_path}lis_input_{domain}.nc'
            os.makedirs(f'{setup.e2esdir}/lis_darun/output/{domain}/lis.config_files/',
                        exist_ok=True)
            cwd = f'{setup.scrdir}/lis_darun/{domain}/'
            os.makedirs(cwd + '/input/', exist_ok=True)
            os.chdir(cwd +'/input/')

            setup.create_symlink(setup.lishdir + '/ghis2s/lis_darun/forcing_variables.txt',
                                 'forcing_variables.txt')
            setup.create_symlink(setup.lishdir + '/ghis2s/lis_darun/noahmp401_parms',
                                 'noahmp401_parms')
            setup.create_symlink(setup.lishdir + '/ghis2s/lis_darun/template_files',
                                 'template_files')
            setup.create_symlink(setup.lishdir + '/ghis2s/lis_darun/attribs','attribs')
            setup.create_symlink(setup.lishdir + '/ghis2s/lis_darun/tables','tables')
            setup.create_symlink(setup.supdir + '/lis_darun/cdf/' +setup.domain,'cdf')
            setup.create_symlink(setup.supdir + '/lis_darun/RS_DATA','RS_DATA')
            setup.create_symlink(ldtfile, f'lis_input_{domain}.nc')

            os.chdir(cwd)
            setup.create_symlink(setup.lisfdir + '/lis/LIS', 'LIS')
            setup.create_symlink(setup.e2esdir + '/lis_darun/output/' + domain, 'output')
            setup.create_symlink(setup.metforc, setup.metforc.rstrip('/').split('/')[-1])
            setup.create_symlink(f'{setup.scrdir}/lis_darun/logs/', 'logs')

            # write lisda_run.j
            # -----------------
            s2s_api.lis_job_file(setup.e2esroot +'/' + setup.config_file, f'lisda_{domain}_run.j',
                                 f'lisda_{domain}_', cwd, str(5))

            if 'discover' in platform.node() or 'borg' in platform.node():
                if 'mil' in setup.constraint:
                    command = './LIS'
                else:
                    slurm_ntasks = os.getenv('SLURM_NTASKS', '1')
                    command = f'mpirun -np {slurm_ntasks} ./LIS'
            else:
                command = 'srun ./LIS'

            # add LIS command
            # ---------------
            with open(f'lisda_{domain}_run.j', 'r', encoding="utf-8") as file:
                filedata = file.read()
            filedata = filedata.replace('COMMAND', command)
            if domain == 'russia_tip':
                filedata = filedata.replace(str(setup.config['FCST']['numprocy']), '60')
            with open(f'lisda_{domain}_run.j', 'w', encoding="utf-8") as file:
                file.write(filedata)

            shutil.copy(f'lisda_{domain}_run.j', f'lisda_{domain}_run.sh')
            utils.remove_sbatch_lines(f'lisda_{domain}_run.sh')
            setup.create_dict(f'{domain}/lisda_{domain}_run.j', 'lis_darun')

            # configure lis.config
            # --------------------
            shutil.copy(setup.e2esdir + '/lis_darun/input/template_files/lis.config_template.' +\
                        setup.domain + '_MERRA2', 'lis.config')
            dapertrstfile = './output/DAPERT/' + yyyyp + \
                mmp + f'/LIS_DAPERT_{yyyyp}{mmp}010000.d01.bin'
            noahmp401rstfile = './output/SURFACEMODEL/' + yyyyp + mmp +\
                f'/LIS_RST_NOAHMP401_{yyyyp}{mmp}010000.d01.nc'
            hymap2rstfile = './output/ROUTING/' + yyyyp + mmp +\
                f'/LIS_RST_HYMAP2_router_{yyyyp}{mmp}010000.d01.nc'
            lsmlislogfile = cwd + f'/logs_{yyyyp}{mmp}/lislog'

            with open('lis.config', 'r', encoding="utf-8") as file:
                filedata = file.read()

            filedata = filedata.replace('LDTINPUTFILE', f'lis_input_{domain}.nc')
            filedata = filedata.replace('DAPERTRSTFILE', dapertrstfile)
            filedata = filedata.replace('NOAHMP401RSTFILE', noahmp401rstfile)
            filedata = filedata.replace('HYMAP2RSTFILE', hymap2rstfile)
            filedata = filedata.replace('STARTYR', yyyyp)
            filedata = filedata.replace('STARTMO', str(monp))
            filedata = filedata.replace('STARTDA', '1')
            filedata = filedata.replace('FINALYR', setup.yyyy)
            filedata = filedata.replace('FINALMO', str(setup.month))
            filedata = filedata.replace('FINALDA', '1')
            filedata = filedata.replace('PERTMODE', pertmode)
            filedata = filedata.replace('LSMLISLOGFILE', lsmlislogfile)
            filedata = filedata.replace('NUMPROCX', str(setup.config['FCST']['numprocx']))
            if domain == 'russia_tip':
                filedata = filedata.replace('NUMPROCY', '60')
            else:
                filedata = filedata.replace('NUMPROCY', str(setup.config['FCST']['numprocy']))

            with open('lis.config', 'w', encoding="utf-8") as file:
                file.write(filedata)

            shutil.copy('lis.config', f'output/lis.config_files/lis.config_darun_{yyyyp}{mmp}')
            _ = TaskLogger(f'lisda_{domain}_run.j',
                           f'{setup.scrdir}/lis_darun/',
                           f'LISDA log files are at: \n{cwd}/logs_{yyyymmp[:6]}/')

        # job file to stitch up LIS_HIST
        jobname = 'stitch_lisda_'
        par_info = {}
        par_info['CPT'] = '1'
        par_info['NT']= str(1)
        par_info['MEM']= '120GB'
        par_info['TPN'] = None
        par_info['MP'] = True
        par_info['SKIP_ARG'] = True
        cwd = f'{setup.scrdir}/lis_darun/'
        os.chdir(cwd)
        command=[(f"python {setup.lishdir}/ghis2s/main/hymap_basins.py -c "
                  f"{setup.e2esroot}/{setup.config_file} "
                  f"-m {yyyyp}{mmp} -s LISDA")]
        tfile = setup.sublist_to_file(command, cwd)
        s2s_api.python_job_file(setup.e2esroot +'/' + setup.config_file, jobname + 'run.j', jobname,
                                str(1), str(2), cwd,  tfile.name, parallel_run=par_info)
        tfile.close()
        os.unlink(tfile.name)
        shutil.copy(f'{jobname}run.j', f'{jobname}run.sh')
        utils.remove_sbatch_lines(f'{jobname}run.sh')
        prev = [f"{key}" for key in setup.schedule.keys() if 'lisda_' in key]
        setup.create_dict(f'{jobname}run.j', 'lis_darun', prev=prev)
        os.chdir(setup.e2esdir)
        return setup.schedule

    def ldt_ics(self, setup, prev):
        ''' runs LDTICS in subdomains '''
        # create subdomain job irectories
        for domain in self.domains:
            ldtfile = f'{self.subdomain_path}lis_input_{domain}.nc'
            ldt_xr = xr.open_dataset(ldtfile)
            latmin = np.round(ldt_xr['lat'].min().values, 2)
            latmax = np.round(ldt_xr['lat'].max().values, 2)
            lonmin = np.round(ldt_xr['lon'].min().values, 2)
            lonmax = np.round(ldt_xr['lon'].max().values, 2)
            cwd = f'{setup.scrdir}/ldt_ics/{domain}/'
            os.makedirs(setup.e2esdir + f'/ldt_ics/{domain}/input/', exist_ok=True)
            os.chdir(setup.e2esdir + f'/ldt_ics/{domain}/')
            [os.makedirs(model, exist_ok=True) for model in setup.models]
            os.makedirs('ldt.config_files', exist_ok=True)
            os.makedirs('template_files', exist_ok=True)
            shutil.copy(
                setup.lishdir +
                f'/ghis2s/ldt_ics/template_files/ldt.config_noahmp401_nmme_TEMPLATE.{setup.resol}',
                'template_files/ldt.config_noahmp401_nmme_TEMPLATE')

            with open('template_files/ldt.config_noahmp401_nmme_TEMPLATE',
                      'r', encoding="utf-8") as file:
                filedata = file.read()
            filedata = filedata.replace('LATMIN', str(latmin))
            filedata = filedata.replace('LATMAX', str(latmax))
            filedata = filedata.replace('LONMIN', str(lonmin))
            filedata = filedata.replace('LONMAX', str(lonmax))
            with open('template_files/ldt.config_noahmp401_nmme_TEMPLATE',
                      'w', encoding="utf-8") as file:
                file.write(filedata)

            os.chdir(setup.e2esdir + f'/ldt_ics/{domain}/input/')
            setup.create_symlink(ldtfile, f"{setup.config['SETUP']['ldtinputfile']}")
            setup.create_symlink(setup.supdir + '/LS_PARAMETERS', 'LS_PARAMETERS')
            os.makedirs(cwd + '/input/', exist_ok=True)
            os.chdir(cwd +'/input/')
            setup.create_symlink(ldtfile, f"{setup.config['SETUP']['ldtinputfile']}")
            setup.create_symlink(setup.supdir + '/LS_PARAMETERS', 'LS_PARAMETERS')
            os.chdir(cwd)
            setup.create_symlink(setup.lisfdir + '/ldt/LDT', 'LDT')
            setup.create_symlink(setup.e2esdir + f'/lis_darun/output/{domain}/', 'lisda_output')
            setup.create_symlink(f'{setup.scrdir}/ldt_ics/logs/', 'logs')

            for model in setup.models:
                setup.create_symlink(setup.e2esdir + f'/ldt_ics/{domain}/{model}', model)

            setup.create_symlink(setup.e2esdir + f'/ldt_ics/{domain}/ldt.config_files',
                                 'ldt.config_files')
            setup.create_symlink(setup.e2esdir + f'/ldt_ics/{domain}/template_files',
                                 'template_files')

            # configure batch script
            # ----------------------
            s2s_api.python_job_file(setup.e2esroot +'/' + setup.config_file,
                                    f'ldtics_{domain}_run.j',
                                    f'ldtics_{domain}_', str(1), str(2), cwd, None)

            command=(f"python {setup.lishdir}/ghis2s/ldt_ics/generate_ldtconfig_files_ensrst_nrt.py -y"
                     f" {setup.yyyy} -m {setup.month} -i ./lisda_output -w {cwd} -s"
                     f" {setup.e2esdir}/{setup.config_file}")

            # add command
            # ---------------
            with open(f'ldtics_{domain}_run.j', 'r', encoding="utf-8") as file:
                filedata = file.read()
                filedata = filedata.replace('COMMAND', command)
            with open(f'ldtics_{domain}_run.j', 'w', encoding="utf-8") as file:
                file.write(filedata)

            setup.create_dict(f'{domain}/ldtics_{domain}_run.j', 'ldt_ics', prev=prev)
            shutil.copy(f'ldtics_{domain}_run.j', f'ldtics_{domain}_run.sh')
            utils.remove_sbatch_lines(f'ldtics_{domain}_run.sh')
            #utils.cylc_job_scripts(f'ldtics_{domain}_run.sh', 2, cwd, command_list=[command])
            mon_abbr = date(int(setup.yyyy), int(setup.mm), 1).strftime('%b')
            _ = TaskLogger(f'ldtics_{domain}_run.j',
                           f'{setup.scrdir}/ldt_ics/',
                           (f'LDTICS log files are at: \n{setup.e2esdir}/ldt_ics/*/'
                            f'ldtlog_noahmp401_{mon_abbr}{setup.yyyy}.0000'))
        os.chdir(setup.e2esdir)
        return setup.schedule

    def lis_fcst(self, setup, prev):
        ''' lis_fcst in subdomains '''
        # create subdomain job irectories
        for domain in self.domains:
            jobname= f'lis_fcst_{domain}'
            ldtfile = f'{self.subdomain_path}lis_input_{domain}.nc'
            os.makedirs(setup.scrdir + f'lis_fcst/{domain}/input/LDT_ICs/', exist_ok=True)
            os.chdir(setup.scrdir + f'lis_fcst/{domain}/input/')
            setup.create_symlink(setup.lishdir + '/ghis2s/lis_darun/forcing_variables.txt',
                                 'forcing_variables.txt')
            setup.create_symlink(setup.lishdir + '/ghis2s/lis_darun/noahmp401_parms',
                                 'noahmp401_parms')
            setup.create_symlink(setup.lishdir + '/ghis2s/lis_fcst/template_files','template_files')
            setup.create_symlink(setup.lishdir + '/ghis2s/lis_fcst/tables','tables')
            setup.create_symlink(ldtfile, f"{setup.config['SETUP']['ldtinputfile']}")

            os.chdir(setup.scrdir + f'lis_fcst/{domain}/input/LDT_ICs/')
            for model in setup.models:
                os.makedirs(setup.scrdir + \
                            f'lis_fcst/{domain}/{setup.yyyy}{setup.mm}/{model}/logs/',
                            exist_ok=True)
                setup.create_symlink(f'{setup.e2esdir}/ldt_ics/{domain}/{model}', model)

            os.chdir(setup.scrdir + f'lis_fcst/{domain}/input/')
            setup.create_symlink(setup.lishdir + '/ghis2s/lis_darun/forcing_variables.txt',
                                 'forcing_variables.txt')
            setup.create_symlink(setup.lishdir + '/ghis2s/lis_darun/noahmp401_parms',
                                 'noahmp401_parms')
            setup.create_symlink(setup.lishdir + '/ghis2s/lis_fcst/template_files','template_files')
            setup.create_symlink(setup.lishdir + '/ghis2s/lis_fcst/tables','tables')
            setup.create_symlink(ldtfile, f"{setup.config['SETUP']['ldtinputfile']}")

            os.chdir(setup.scrdir + f'lis_fcst/{domain}/')
            cwd = setup.scrdir + f'lis_fcst/{domain}/'
            setup.create_symlink(setup.lisfdir + 'lis/LIS', 'LIS')
            setup.create_symlink(setup.e2esdir + 'bcsd_fcst', 'bcsd_fcst')

            generate_lis_config_scriptfiles_fcst.main(setup.e2esroot + setup.config_file,
                                                      setup.year, setup.month, cwd, jobname)
            for model in setup.models:
                lindex = ''
                job_list = sorted(glob.glob(f"{jobname}_{model}*_run.j"))
                n_files = len(job_list)
                for file_no, jfile in enumerate(job_list):
                    if n_files > 1:
                        if file_no == 0:
                            setup.create_dict(f'{domain}/{jfile}', 'lis_fcst', prev=prev)
                        else:
                            setup.create_dict(f'{domain}/{jfile}', 'lis_fcst',
                                              prev=f'{domain}/{job_list[file_no-1]}')
                            lindex = f'_{file_no:02d}'
                    else:
                        setup.create_dict(f'{domain}/{jfile}', 'lis_fcst', prev=prev)
                        _ = TaskLogger(jfile,
                                       os.getcwd(),
                                       f'{model} log at: \n{setup.scrdir}/lis_fcst/{domain}/'
                                       f'{setup.yyyy}{setup.mm}/{model}/logs/lislog{lindex}.')


        # job file to stitch up LIS_HIST
        jobname = 'stitch_lisfcst_'
        par_info = {}
        par_info['CPT'] = '1'
        par_info['NT']= str(1)
        par_info['MEM']= '120GB'
        par_info['TPN'] = None
        par_info['MP'] = True
        par_info['SKIP_ARG'] = True
        cwd = f'{setup.scrdir}/lis_fcst/'
        os.chdir(cwd)

        # lead months string
        start_date = datetime(setup.year, setup.month, 1)
        month_strings = [(start_date + relativedelta(months=i)).strftime("%Y%m")
                         for i in range(setup.config["EXP"]["lead_months"])]
        lead_months = ','.join(month_strings)

        command=[(f"python {setup.lishdir}/ghis2s/main/hymap_basins.py "
                  f"-c {setup.e2esroot}/{setup.config_file} "
                  f"-m {lead_months} -s FCST")]
        tfile = setup.sublist_to_file(command, cwd)
        s2s_api.python_job_file(setup.e2esroot +'/' + setup.config_file, jobname + 'run.j', jobname,
                                str(1), str(3), cwd,  tfile.name, parallel_run=par_info)
        tfile.close()
        os.unlink(tfile.name)
        shutil.copy(f'{jobname}run.j', f'{jobname}run.sh')
        utils.remove_sbatch_lines(f'{jobname}run.sh')

        prev = None
        fcst_list = [job for job in setup.schedule
                     if 'lis_fcst' in job] or None

        if len(fcst_list) > 0:
            prev=[]
            for domain in self.domains:
                sublist1 = sorted([item for item in fcst_list if domain in item])
                for model in setup.models:
                    sublist = sorted([item for item in sublist1 if model in item])
                    last_item = sublist[-1] if sublist else None
                    if last_item is not None:
                        prev.append(last_item)

        setup.create_dict(f'{jobname}run.j', 'lis_fcst', prev=prev)
        os.chdir(setup.e2esdir)
        return setup.schedule

    def stitch_lis(self, input_path, model, filename, output_file, logger,
                   fill_value=-9999.0, new_dims=None):
        '''
        Stitch LIS output netCDF files (LSM or Routing) from subdomains to global grid.
        Handles both:
          - LSM files:     variables with dims (north_south, east_west) and
                           (profile_dim, north_south, east_west)
          - Routing files: variables with dims (north_south, east_west) only

        Parameters
        ----------
        input_path : str  config["SETUP"]["E2ESDIR"]}/lis_darun/output/
        model : str /SURFACEMODEL/202512/  or  /ROUTING/202512/
        filename : str LIS_HIST_202512020000.d01.nc
        output_file : str Full path for the output global netCDF file
        fill_value : float Fill/missing value used in LIS files (default: -9999.0)
        new_dims : dict, optional Dimension rename map, e.g. {'north_south': 'lat', 
                   'east_west': 'lon'}
        '''

        def get_subdomain_path(dom_name):
            return os.path.join(input_path, f'{dom_name}{model}{filename}')

        # ----------------------------------------------------------------
        # 1. Global grid from domain_masks
        # ----------------------------------------------------------------
        global_lat = self.domain_masks['lat']
        global_lon = self.domain_masks['lon']
        n_lat      = len(global_lat)
        n_lon      = len(global_lon)

        global_lat_vals = np.round(global_lat.values, 2)
        global_lon_vals = np.round(global_lon.values, 2)
        lat_to_idx = {float(v): i for i, v in enumerate(global_lat_vals)}
        lon_to_idx = {float(v): i for i, v in enumerate(global_lon_vals)}

        # ----------------------------------------------------------------
        # 2. Locate first existing subdomain file
        # ----------------------------------------------------------------
        first_file = None
        for dom_name in self.domains:
            candidate = get_subdomain_path(dom_name)
            if os.path.exists(candidate):
                first_file = candidate
                #print(f"\nReading structure from : {candidate}")
                break

        if first_file is None:
            raise FileNotFoundError(
                f"No subdomain files found.\n"
                f"  input_path : {input_path}\n"
                f"  model      : {model}\n"
                f"  filename   : {filename}"
            )

        # ----------------------------------------------------------------
        # 3. Classify variables, cxclude 'lat' and 'lon'  coordinate variables
        # ----------------------------------------------------------------
        lis_coord_vars = {'lat', 'lon'}

        with xr.open_dataset(first_file) as first_ds:

            vars_to_stitch = []
            vars_to_keep   = []
            all_dims       = {}
            var_attrs      = {}

            for var_name, var in first_ds.data_vars.items():
                if var_name in lis_coord_vars:
                    continue

                all_dims[var_name]  = var.dims
                var_attrs[var_name] = dict(var.attrs)

                if 'north_south' in var.dims and 'east_west' in var.dims:
                    vars_to_stitch.append(var_name)
                else:
                    vars_to_keep.append(var_name)

            # Collect profile/extra dimension sizes SoilMoist_profiles
            extra_dim_sizes  = {}
            extra_dim_coords = {}

            for var_name in vars_to_stitch:
                for dim in all_dims[var_name]:
                    if dim in ('north_south', 'east_west'):
                        continue
                    if dim not in extra_dim_sizes:
                        extra_dim_sizes[dim] = first_ds.sizes[dim]
                        if dim in first_ds.coords:
                            extra_dim_coords[dim] = first_ds[dim].values

            kept_vars_data  = {}
            kept_vars_attrs = {}
            kept_vars_dims  = {}
            for var_name in vars_to_keep:
                kept_vars_data[var_name]  = first_ds[var_name].values.copy()
                kept_vars_attrs[var_name] = dict(first_ds[var_name].attrs)
                kept_vars_dims[var_name]  = first_ds[var_name].dims

            # Save lat/lon attrs to reproduce on output
            if 'ROUTING' in model:
                surf_model_path = model.replace('ROUTING', 'SURFACEMODEL')
                surf_first_file = first_file.replace(model, surf_model_path)
                with xr.open_dataset(surf_first_file) as surf_ds:
                    lat_attrs = dict(surf_ds['lat'].attrs) if 'lat' in surf_ds else {}
                    lon_attrs = dict(surf_ds['lon'].attrs) if 'lon' in surf_ds else {}
            else:
                lat_attrs         = dict(first_ds['lat'].attrs) if 'lat' in first_ds else {}
                lon_attrs         = dict(first_ds['lon'].attrs) if 'lon' in first_ds else {}

            global_file_attrs = dict(first_ds.attrs)

        # ----------------------------------------------------------------
        # 4. Initialise global arrays
        # ----------------------------------------------------------------
        global_vars = {}

        for var_name in vars_to_stitch:
            var_dims = all_dims[var_name]
            shape = []
            for dim in var_dims:
                if dim == 'north_south':   shape.append(n_lat)
                elif dim == 'east_west':   shape.append(n_lon)
                else:                      shape.append(extra_dim_sizes[dim])

            global_vars[var_name] = {
                'array': np.full(shape, fill_value, dtype=np.float32),
                'dims' : list(var_dims)
            }

        # ----------------------------------------------------------------
        # 5. Stitch loop
        #   5.1. Get domain bounding box in global index space from mask
        #   5.2. Map subdomain lat/lon coords → global indices via lookup
        #   5.3. Use subdomain's ACTUAL coordinate coverage as the crop
        #        extent  handles cases where subdomain file is slightly
        #        smaller than the mask bounding box (e.g. siberia_2 where
        #        subdomain starts at lon=58.05 but mask bbox starts at
        #        lon=57.95)
        #   5.4. Build local_mask from subdomain's actual coverage
        #   5.5  Apply mask + NaN filter → write into global array
        # ----------------------------------------------------------------

        n_missing = 0
        n_skipped = 0
        files_to_delete = []

        for i, dom_name in enumerate(self.domains):
            subdomain_file = get_subdomain_path(dom_name)

            if not os.path.exists(subdomain_file):
                logger[0].warning(f"  not found → {subdomain_file}", subtask=logger[1])
                n_missing += 1
                continue

            domain_mask = self.domain_masks == (i + 1)
            lat_indices, _ = domain_mask.values.nonzero()

            if len(lat_indices) == 0:
                logger[0].warning("  No pixels in mask, skipping...", subtask=logger[1])
                n_skipped += 1
                continue

            files_to_delete.append(subdomain_file)

            with xr.open_dataset(subdomain_file) as sub_ds:
                # ----------------------------------------------------------
                # Map every subdomain lat/lon value → global index.
                # ----------------------------------------------------------
                if 'ROUTING' in model:
                    surf_model_path = model.replace('ROUTING', 'SURFACEMODEL')
                    surf_subdomain_file = subdomain_file.replace(model, surf_model_path)

                    if not os.path.exists(surf_subdomain_file):
                        logger[0].warning(f"  SURFACEMODEL file missing → {surf_subdomain_file}",
                                          subtask=logger[1])
                        n_skipped += 1
                        continue

                    with xr.open_dataset(surf_subdomain_file) as surf_ds:
                        sub_lat_r = np.round(surf_ds['lat'].values, 2)
                        sub_lon_r = np.round(surf_ds['lon'].values, 2)
                else:
                    sub_lat_r = np.round(sub_ds['lat'].values, 2)
                    sub_lon_r = np.round(sub_ds['lon'].values, 2)

                sub_lat_global = np.array(
                    [lat_to_idx.get(float(v), -1) for v in sub_lat_r]
                )
                sub_lon_global = np.array(
                    [lon_to_idx.get(float(v), -1) for v in sub_lon_r]
                )

                # Subdomain rows/cols that successfully map to global grid
                valid_sub_rows = np.where(sub_lat_global >= 0)[0]
                valid_sub_cols = np.where(sub_lon_global >= 0)[0]

                if len(valid_sub_rows) == 0 or len(valid_sub_cols) == 0:
                    logger[0].warning(f"  skipping: no coordinate overlap for {dom_name}",
                                      subtask=logger[1])
                    n_skipped += 1
                    continue

                # Contiguous subdomain row/col range covering the global grid
                sub_r0 = int(valid_sub_rows[0])
                sub_r1 = int(valid_sub_rows[-1])
                sub_c0 = int(valid_sub_cols[0])
                sub_c1 = int(valid_sub_cols[-1])

                # Corresponding global row/col range
                gbl_r0 = int(sub_lat_global[sub_r0])
                gbl_r1 = int(sub_lat_global[sub_r1])
                gbl_c0 = int(sub_lon_global[sub_c0])
                gbl_c1 = int(sub_lon_global[sub_c1])

                sub_lat_slice = slice(sub_r0, sub_r1 + 1)
                sub_lon_slice = slice(sub_c0, sub_c1 + 1)
                local_mask = domain_mask.isel(
                    lat=slice(gbl_r0, gbl_r1 + 1),
                    lon=slice(gbl_c0, gbl_c1 + 1)
                ).values

                # ------------------------------------------------------
                # Stitch each variable
                # ------------------------------------------------------
                for var_name in vars_to_stitch:
                    if var_name not in sub_ds.data_vars:
                        logger[0].warning(f"  {var_name} missing from {dom_name}",
                                          subtask=logger[1])
                        continue

                    var_data = sub_ds[var_name].values
                    var_dims = global_vars[var_name]['dims']
                    ns_idx   = var_dims.index('north_south')
                    ew_idx   = var_dims.index('east_west')

                    # Crop subdomain data to its valid extent
                    sub_slices         = [slice(None)] * var_data.ndim
                    sub_slices[ns_idx] = sub_lat_slice
                    sub_slices[ew_idx] = sub_lon_slice
                    var_data_crop      = var_data[tuple(sub_slices)]

                    # Global slice  use actual subdomain coverage extent
                    slices_global         = [slice(None)] * len(var_dims)
                    slices_global[ns_idx] = slice(gbl_r0, gbl_r1 + 1)
                    slices_global[ew_idx] = slice(gbl_c0, gbl_c1 + 1)

                    # Expand local_mask to match variable dimensions
                    # 2D vars: (crop_ns, crop_ew)          → no expansion
                    # 3D vars: (n_profiles, crop_ns, crop_ew) → prepend axis
                    expanded_mask = local_mask
                    for dim_idx, dim in enumerate(var_dims):
                        if dim not in ('north_south', 'east_west'):
                            if dim_idx < ns_idx:
                                expanded_mask = np.expand_dims(
                                    expanded_mask, axis=0)
                            else:
                                expanded_mask = np.expand_dims(
                                    expanded_mask, axis=-1)

                    expanded_mask = np.broadcast_to(
                        expanded_mask, var_data_crop.shape)

                    # Write where domain mask is True and data is not NaN
                    valid_mask = expanded_mask & ~np.isnan(var_data_crop)

                    global_array = global_vars[var_name]['array']
                    global_array[tuple(slices_global)][valid_mask] = \
                        var_data_crop[valid_mask]

        # ----------------------------------------------------------------
        # 6. Build global xarray Dataset
        # ----------------------------------------------------------------
        # lat/lon as auxiliary coordinate variables on north_south/east_west
        global_coords = {
            'lat': xr.Variable(
                'north_south', global_lat.values,
                attrs=lat_attrs if lat_attrs else {
                    'units'         : 'degree_north',
                    'standard_name' : 'latitude',
                    'long_name'     : 'latitude',
                    'missing_value' : fill_value,
                    '_FillValue'    : fill_value,
                }
            ),
            'lon': xr.Variable(
                'east_west', global_lon.values,
                attrs=lon_attrs if lon_attrs else {
                    'units'         : 'degree_east',
                    'standard_name' : 'longitude',
                    'long_name'     : 'longitude',
                    'missing_value' : fill_value,
                    '_FillValue'    : fill_value,
                }
            ),
        }

        for dim_name, coord_vals in extra_dim_coords.items():
            global_coords[dim_name] = coord_vals

        global_data_vars = {}

        for var_name in vars_to_stitch:
            var_info = global_vars[var_name]
            global_data_vars[var_name] = xr.Variable(
                dims  = var_info['dims'],
                data  = var_info['array'],
                attrs = var_attrs[var_name]
            )

        for var_name in vars_to_keep:
            if var_name in kept_vars_data:
                global_data_vars[var_name] = xr.Variable(
                    dims  = kept_vars_dims[var_name],
                    data  = kept_vars_data[var_name],
                    attrs = kept_vars_attrs[var_name]
                )

        # ----------------------------------------------------------------
        # 7. Update global attributes for the stitched domain
        # ----------------------------------------------------------------
        global_file_attrs['SOUTH_WEST_CORNER_LAT'] = float(global_lat.values[-1])
        global_file_attrs['SOUTH_WEST_CORNER_LON'] = float(global_lon.values[0])

        global_ds = xr.Dataset(
            data_vars = global_data_vars,
            coords    = global_coords,
            attrs     = global_file_attrs
        )

        if new_dims is not None:
            global_ds = global_ds.rename(new_dims)

        # ----------------------------------------------------------------
        # 8. Write to netCDF
        # ----------------------------------------------------------------
        encoding = {
            var_name: {
                'zlib'      : True,
                'complevel' : 6,
                '_FillValue': fill_value,
                'dtype'     : 'float32',
            }
            for var_name in vars_to_stitch
        }

        global_ds.to_netcdf(output_file, encoding=encoding)
        global_ds.close()

        # delete subdomain files
        [os.remove(f) for f in files_to_delete if os.path.exists(f)]

    def decompose(self, indata, output_path):
        ''' decompose global domain to subdomains - parallel processing '''

        os.makedirs(output_path, exist_ok=True)
        # Ensure data is chunked
        if not indata.chunks:
            print("Chunking indata dataset...")
            indata = indata.chunk({'lat': 100, 'lon': 100})

        spatial_dims = {'lat', 'lon'}
        tasks = []
        for i, dom_name in enumerate(self.domains):
            print(f"Preparing domain {i+1}/{len(self.domains)}: {dom_name}")

            domain_mask = self.domain_masks == (i + 1)
            mask_computed = domain_mask.compute()
            lat_indices, lon_indices = mask_computed.values.nonzero()

            if len(lat_indices) == 0:
                print(f"  WARNING: No pixels for {dom_name}, skipping...")
                continue

            lat_min, lat_max = lat_indices.min(), lat_indices.max()
            lon_min, lon_max = lon_indices.min(), lon_indices.max()

            subset_mask = domain_mask.isel(
                lat=slice(lat_min, lat_max + 1),
                lon=slice(lon_min, lon_max + 1)
            )

            subset_vars = {}
            subset_coords = {}
            for var_name in indata.data_vars:
                var = indata[var_name]
                var_dims = set(var.dims)

                if spatial_dims.issubset(var_dims):
                    subset_var = var.isel(
                        lat=slice(lat_min, lat_max + 1),
                        lon=slice(lon_min, lon_max + 1)
                    )
                    subset_vars[var_name] = subset_var.where(subset_mask)
                elif not var_dims.intersection(spatial_dims):
                    subset_vars[var_name] = var
                else:
                    print(f"  WARNING: {var_name} has unusual dimensions {var.dims}, keeping as-is")
                    subset_vars[var_name] = var

            for coord_name in indata.coords:
                coord = indata.coords[coord_name]
                coord_dims = set(coord.dims) if coord.dims else set()
                if coord_name == 'lat':
                    subset_coords['lat'] = indata['lat'].isel(lat=slice(lat_min, lat_max + 1))
                elif coord_name == 'lon':
                    subset_coords['lon'] = indata['lon'].isel(lon=slice(lon_min, lon_max + 1))
                elif spatial_dims.issubset(coord_dims):
                    subset_coords[coord_name] = coord.isel(
                        lat=slice(lat_min, lat_max + 1),
                        lon=slice(lon_min, lon_max + 1)
                    )
                elif not coord_dims.intersection(spatial_dims):
                    subset_coords[coord_name] = coord

            subset_ds = xr.Dataset(subset_vars, coords=subset_coords, attrs=indata.attrs.copy())
            output_file = os.path.join(output_path, f'{dom_name}.nc')
            encoding = {var: {'zlib': True, 'complevel': 4}
                       for var in subset_ds.data_vars}

            # Create delayed task
            task = dask.delayed(subset_ds.to_netcdf)(output_file, encoding=encoding, compute=True)
            tasks.append(task)

        print(f"\nWriting {len(tasks)} subdomain files in parallel...")
        with ProgressBar():
            dask.compute(*tasks)

        print(f"\nDecomposition complete! {len(tasks)} files written to {output_path}")

if __name__ == "__main__":
    import argparse
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, type=str, help='config file')
    parser.add_argument('-s', '--stitch', required=False, default=None, type=str,
                        help='stitch options: LISDA, FCST')
    parser.add_argument('-m', '--yyyymm', required=False, default=None, type=str,
                    help='month(s) YYYYMM comma-separated (e.g. 202603 or 202603,202604,202605)')
    args = parser.parse_args()

    with open(args.config_file, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    hymap = HyMAPDomains(config['SETUP']['supplementarydir'] + '/hymap_domains/10km/')

    if args.stitch is not None:
        task_name = os.environ.get('SCRIPT_NAME')
        LOG_MSG = 'main/hymap_basins.py is stitching subdomains to global'
        logger_ = TaskLogger(task_name, os.getcwd(), LOG_MSG)
        if args.yyyymm:
            yyyymm_list = [ym.strip() for ym in args.yyyymm.split(',')]
            for ym in yyyymm_list:
                if not re.match(r'^\d{6}$', ym):
                    logger_.error(
                        f"Invalid YYYYMM format: '{ym}'. Expected 6-digit string like 202603.")
        else:
            logger_.error("YYYYMM list must be provided.")

        # loop through forecast months
        for yyyymm in yyyymm_list:
            for model_ in ['ROUTING', 'SURFACEMODEL']:
                if args.stitch.lower() == 'lisda':
                    INDIR = f'{config["SETUP"]["E2ESDIR"]}/lis_darun/output/'
                    OUTDIR = f'{config["SETUP"]["E2ESDIR"]}/lis_darun/output/{model_}/{yyyymm}/'
                    os.makedirs(OUTDIR, exist_ok=True)
                    file_list = sorted(glob.glob(
                        f'{INDIR}/{hymap.domains[0]}/{model_}/{yyyymm}/LIS_HIST_*nc'))
                    basenames = [os.path.basename(f) for f in file_list]
                    for _file in basenames:
                        OUTFILE = f'{OUTDIR}/{_file}'
                        logger_.info(OUTFILE, subtask=f'{model_}-{yyyymm}')
                        hymap.stitch_lis(INDIR, f'/{model_}/{yyyymm}/', _file, OUTFILE,
                                         [logger_, f'{model_}-{yyyymm}'])

                if args.stitch.lower() == 'fcst':
                    INDIR_SCR = f'{config["SETUP"]["E2ESDIR"]}/scratch/{yyyymm_list[0]}/lis_fcst/'
                    for nmme_model in config["EXP"]["NMME_models"]:
                        INDIR = f'{INDIR_SCR}/{hymap.domains[0]}/{yyyymm_list[0]}/{nmme_model}/'
                        OUTDIR = (f'{config["SETUP"]["E2ESDIR"]}/lis_fcst/{yyyymm_list[0]}/'
                                  f'{nmme_model}/{model_}/{yyyymm}/')
                        os.makedirs(OUTDIR, exist_ok=True)
                        file_list = sorted(glob.glob(
                            f'{INDIR}/{model_}/{yyyymm}/LIS_HIST_*nc'))
                        basenames = [os.path.basename(f) for f in file_list]
                        for _file in basenames:
                            OUTFILE = f'{OUTDIR}/{_file}'
                            logger_.info(OUTFILE, subtask=f'{nmme_model}-{model_}-{yyyymm}')
                            hymap.stitch_lis(INDIR, f'/{model_}/{yyyymm}/', _file, OUTFILE,
                                             [logger_, f'{model_}-{yyyymm}'])
