# ghis2s Python Package for Cylc Implementation

# ACRONYMS

**E2ESDIR**: The GHI-S2S forecast directory where S2S forecasts reside  
**GHIS2S**: GHI-S2S software system developed by the LIS team  
**GHIREPOS**: Operational software developed by 16WS  
**LISFDIR**: The path to LISF installation  

# (1) The ghis2s Python Package

We present `ghis2s` as a Python package that can be efficiently coupled with the GHIREPOS operational software at 16WS. The core component is the [*s2s_run.py*](https://github.com/smahanam/LISF-1/blob/parallelizing/lis/utils/usaf/S2S/ghis2s/s2s_app/s2s_run.py) Python module -- the master script containing `S2Srun` class.

## Configuration

The `S2Srun` class requires only one input file: a YAML [*configuration file*](https://github.com/smahanam/LISF-1/blob/parallelizing/lis/utils/usaf/S2S/ghis2s/s2s_app/s2s_config_global_fcast) containing system and experiment-related parameters/paths for the forecast. This configuration file must be placed in the E2ESDIR directory.

## Key Capabilities

The `S2Srun` class provides the following functionality:

1. **Experiment Setup**: Creates S2S experiments for the specific month, including:
   - Run directories
   - Symbolic links
   - Job scripts for forecast-related tasks

2. **Workflow Generation**: Writes the `flow.cylc` file for Cylc workflow management based on ([for example](https://github.com/smahanam/FileSharing/blob/main/E2ESDIR/scratch/202401/CYLC-202401/flow.cylc)):
   - Selected NMME models 
   - Requirements specified in the forecast configuration file 
   **Note:** The flow.cylc file uses [this global.cylc](https://github.com/smahanam/LISF-1/blob/parallelizing/lis/utils/usaf/S2S/ghis2s/cylc_script/global.cylc) configuration file

## Methods

- **Main method**: Sets up end-to-end S2S forecasts for a particular month
- **Individual methods**: Allows calling specific forecast process steps:
  - `lisda_run`
  - `ldt_ics`
  - `bcsd`
  - `lis_fcst`
  - `post`
  - `s2smetrics`
  - `s2splots`

These methods correspond to the 7 main steps of the end-to-end forecast process.

## Integration

The `ghis2s` Python package includes a supplementary program, `ghis2s_program.py` designed to couple GHI-S2S with operational GHIREPOS systems.

## Command Line Usage

In addition, the `s2s_run.py` script offers two useful command-line features:

### 1. Standard execution (setting up the experiment):
```
s2s_run.py -y YYYY -m M -c CONFIG_FILE
```
### 2. SLURM job submission
```
s2s_run.py -y YYYY -m M -c CONFIG_FILE -j
```
# (2) GHIREPOS-GHIS2S coupler: ghis2s_program.py 

## Program input environment variables

- **"PYTHONPATH"**: (str, LISFDIR/lis/utils/usaf/S2S/)
- **"E2ESDIR"**: (str, E2ESDIR)
- **"CONFIG_FILE"**: (str, Config file name – must be located in E2ESDIR)
- **"FORECAST_YEAR"**: (int, year)
- **"FORECAST_MONTH"**: (int, month)
- **"USER_EMAIL"**: (str, email_address)
- **"S2S_STEP"**: (str, "E2E")

### Acceptable keys for S2S_STEP:
- `E2E`: end-to-end S2S forecast
- `LISDA`
- `LDTICS` 
- `BCSD`
- `FCST`
- `POST`
- `METRICS`
- `PLOTS`

### Optional environment variables:
The default setting of the below two environment variables are `False`:

- **"ONE_STEP"**: `True` allows to run only the above specified S2S_STEP (bool, False)
- **"SUBMIT_JOB"**: `True` submits the job to the SLURM job management system instead of Cylc (bool, False)

## Passing additional environment variables to flow.cylc file
The `additional_env_vars = {}` dictionary in `ghis2s_program.py` allows the user to pass any additional environment variables as key-value pairs to the flow.cylc file. The `s2s_run.py` script will write those variables in flow.cylc.

# (3) Installing Cylc Workflow
```bash
WORKFLOW_NAME="CYLC-${FORECAST_YEAR}${MM} 
LOGDIR=${E2ESDIR}/scratch/${FORECAST_YEAR}${MM}/${WORKFLOW_NAME} 
cylc install --symlink-dirs=run=$LOGDIR
```
This redirects the Cylc workflow logs to the E2ESDIR scratch directory rather than the default /home/$USER/cylc-run location. 

# (4) Centralized logging

# (5) Operational Notes and Cylc Design Rationale  

**a) Why should ghis2s’s ghis2s_program.py be executed every month?**  
  
Each month requires customized LIS input/configuration files, job script arguments, and month-specific symbolic links, all of which must be placed under scratch/YYYYMM.   Therefore, the script, **ghis2s_program.py**, must be executed monthly to install and configure the Cylc {WORKFLOW_NAME} accordingly.  
  
**b) Can the same flow.cylc file be reused each month?**  
  
No. As stated above, monthly differences in input files and configurations require that **ghis2s_program.py** be executed each time. The **ghis2s** package programmatically generates the 800+ line flow.cylc file to minimize human error and improve efficiency.  
  
**c) Why is [[dependencies]] → [[[R1]]] necessary?**  
  
Launching an operational S2S forecast requires human oversight each month due to the following reasons:  
i) CFSv2 data latency is typically a few days  
ii) NMME precipitation data are delivered by the 8 to 10th day of each month    
iii) Occasionally, a particular NMME model may be unavailable  
iv) **ghis2s** performs checks for CFSv2 and NMME file availability before launching the forecast.  
  
Therefore, users must manually verify the availability of all required meteorological forcings before initiating the forecast.  
The {WORKFLOW_NAME} workflow runs only once, covering the duration specified by the lead_time parameter (in months) in the s2s_config_global_fcast file.  
Note that the **ghis2s Cylc** implementation does **NOT** utilize Cylc’s built-in trigger, clock, or timer features to manage the workflow.  
  
**d) Is ghis2s's Cylc implementation 100% system-agnostic?**  
  
Yes and no.  
  
While the shell scripts (*.sh) avoid hardcoded SLURM directives, certain tasks benefit significantly from using SLURM’s **srun** for resource allocation. For example, Python tasks that use multiprocessing perform better when invoked as:  
```srun --exclusive --cpus-per-task=5 --ntasks 1 python script.py```  

This approach has been more effective than using Cylc’s native mechanisms in such cases.  
That said, ghis2s includes a feature to generate fully system-agnostic shell scripts (i.e., no **srun**), although this feature is currently disabled for performance reasons.
An example of a **SLURM's srun-free** shell script can be found here:   
*https://github.com/smahanam/FileSharing/blob/main/E2ESDIR/scratch/202502/s2splots/s2splots_01_run.sh*  

**e) How does ghis2s differ from other GHI subsystems (GHI-NRT, GHI-MR)?**  
  
Although the GHI-S2S workflow includes over 150 tasks and is more complex than other subsystems, **the master script of the ghis2s software tool, [*s2s_run.py*](https://github.com/smahanam/LISF-1/blob/support/lisf_557ww_7.7_s2srf/lis/utils/usaf/S2S/ghis2s/s2s_app/s2s_run.py)**, simplifies execution by consolidating all tasks into a single command driven by a unified configuration file.  
The script automates the execution of all tasks based on their dependencies, effectively eliminating the need for manual intervention.
  








