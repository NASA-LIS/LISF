# ghis2s Python Package for Cylc Implementation

# ACRONYMS

**E2ESDIR**: The GHI-S2S forecast directory where S2S forecasts reside  
**GHIS2S**: GHI-S2S software system developed by the LIS team  
**GHIREPOS**: Operational software developed by 16WS  
**LISFDIR**: The path to LISF installation  

# The ghis2s Python Package

We present `ghis2s` as a Python package that can be efficiently coupled with the GHIREPOS operational software at 16WS. The core component is the `s2s_run.py` Python module -- the master script containing `S2Srun` class.

## Configuration

The `S2Srun` class requires only one input file: a YAML configuration file containing system and experiment-related parameters/paths for the forecast. This configuration file must be placed in the E2ESDIR directory.

## Key Capabilities

The `S2Srun` class provides the following functionality:

1. **Experiment Setup**: Creates S2S experiments for the specific month, including:
   - Run directories
   - Symbolic links
   - Job scripts for forecast-related tasks

2. **Workflow Generation**: Writes the `flow.cylc` file for Cylc workflow management based on:
   - Selected NMME models
   - Requirements specified in the forecast configuration file

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
# GHIREPOS-GHIS2S coupler: ghis2s_program.py 

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

## Installing Cylc Workflow
```bash
WORKFLOW_NAME="CYLC-$${FORECAST_YEAR}$${MM}" 
LOGDIR="$${E2ESDIR}/scratch/$${FORECAST_YEAR}$${MM}/$${WORKFLOW_NAME}" 
cylc install --symlink-dirs=run=$LOGDIR

This redirects the Cylc workflow logs to the E2ESDIR scratch directory rather than the default /home/$USER/cylc-run location.






This repository contains this README file and a mock-up forecast directory (**"E2ESDIR"**) of the GHI-S2S forecasting system.
GHI-S2S consists of approximately 150 tasks that follow a predefined schedule. These tasks have been grouped into 50+ job files to optimize computer resources.  
A comprehensive description of GHI-S2S can be found at:  
*https://github.com/smahanam/LISF-1/blob/support/lisf_557ww_7.7_s2srf/lis/utils/usaf/S2S/README_GHI-S2S_LIS7.7*  
  
The master script of ghis2s ([*s2s_run.py*](https://github.com/smahanam/LISF-1/blob/support/parallelizing/lis/utils/usaf/S2S/ghis2s/s2s_app/s2s_run.py)) creates run-directories, establishes necessary links, generates bash script files, and sets up the complete S2S forecast experiment each month.

We present **ghis2s** as a Python package that can be called either from within an external Python program or run from the command line. The **ghis2s** package provides a Python script, [*ghis2s_program.py*](https://github.com/smahanam/LISF-1/blob/parallelizing/lis/utils/usaf/S2S/ghis2s/cylc_script/ghis2s_program.py) , for importing **ghis2s**, and installing **Cylc** workflow. Optionally, **ghis2s_program.py** can also execute monthly forecast runs in the SLURM system. Users are welcome to copy **[ghis2s_program.py](https://github.com/smahanam/LISF-1/blob/parallelizing/lis/utils/usaf/S2S/ghis2s/cylc_script/ghis2s_program.py)** from the LISF repository into their working directories and modify it as needed.  

The package requires only one main input file: a YAML configuration file containing system and experiment-related parameters/paths. Example configuration:   
*https://github.com/smahanam/FileSharing/blob/main/E2ESDIR/s2s_config_global_fcast*  
  
Below, the two paths are specified among SETUP parameters in the configuration file, **s2s_config_global_fcast**:  
**E2ESDIR:** The GHI-S2S forecast directory where the S2S forecast resides that must contain the above configuration file.  
**LISFDIR:** The path to LISF installation. 

*The following description uses the S2S forecast initialized on January 1, 2025 as an example.*  
  
## 1) Setting up the Forecast Run Directory
### a) Check the main E2ESDIR Directory
i) Ensure **E2ESDIR** and **LISFDIR** in the s2s_config_global_fcast configuration file are correct.  
ii) Ensure the **hindcast** directory and the land initial conditions ("lis_darun/output model restart files) are available and linked.

### b) Set up the Working directory and install the Cylc Workflow

i) **Cylc** creates the /home/$USER/cylc-run directory during the monthly forecast installation. Therefore, /home/$USER/cylc-run  is *NOT* a recommended name for the user's working directory.  
ii) Copy **ghis2s_program.py** to your working directory and edit *E2ESDIR* parameter to specify user's E2ESDIR.  
iii) Load the LISF Python module and set the ENVIRONMENT variable PYTHONPATH  
  
```
on Discover e.g.
module use -a {LISFDIR}/env/discover/
module --ignore-cache load lisf_7.5_intel_2023.2.1_s2s
     
export PYTHONPATH={LISFDIR}/lis/utils/usaf/S2S/  
OR  
setenv PYTHONPATH {LISFDIR}/lis/utils/usaf/S2S/
```
### c) Test the Script
The **ghis2s_program.py** script imports the **“S2Srun”** class from the **“s2s_run”** module in the **“ghis2s”** package, and instantiates as **s2s**.
The **s2s** instance has the following methods: s2s.main() [end-to-end 7-steps]; s2s.lis_darun(); s2s.ldt_ics(); s2s.bcsd(); s2s.lis_fcst(); s2s.s2spost(); s2s.s2smetric(); s2s.s2splots()
s2s.write_cylc_snippet() [writes CYLC_workflow.rc in scratch/YYYYMM]; and s2s.submit_jobs()  [submits to the SLURM queue].
   
To display options, run the help option:
``` python ghis2s_program.py -h ```
This will print:
```
usage: ghis2s_program.py [-h] -c CONFIG_FILE -y YEAR -m MONTH -e EMAIL [-s STEP] [-o] [-j]

options:
  -h, --help            show this help message and exit
  -c CONFIG_FILE, --config_file CONFIG_FILE
                        config file
  -y YEAR, --year YEAR  forecast year
  -m MONTH, --month MONTH
                        forecast month
  -e EMAIL, --email EMAIL
                        user email
  -s STEP, --step STEP  S2S step: LISDA, LDTICS, BCSD, FCST, POST, METRICS or PLOTS
  -o, --one_step        Is only one step (default: False)?
  -j, --submit_job      Submit SLURM jobs (default: False -> CYLC)?
  
Note: CONFIG_FILE, YEAR, MONTH and EMAIL are mandatory input arguments.  
```  

## 2) Creating E2ES Scratch (Run) Directories, and Job Files for the Forecast Month
Run the following command to initialize the forecast environment:  
  
```python ghis2s_program.py -y 2025 -m 1 -c s2s_config_global_fcast -e USER_EMAIL```  

This command performs the following:

**(1) Create the directories and links:**   
This step will create the main **"scratch"** run-directory of the target initial forecast month **(scratch/202501)**, and 7 run-directories for the seven steps: lis_darun, ldt_ics, bcsd_fcst, lis_fcst, s2spost, s2smetric, and s2splots under the main **"scratch"** run-directory. For example:  
*https://github.com/smahanam/FileSharing/tree/main/E2ESDIR/scratch/202501*  
  
**Note:** The temporary **scratch** directory was introduced to store miscellaneous files, including standard *.err and *.out logs, in a location that can be safely deleted. This helps keep the main E2ESDIR and forecast output directories clean and organized.
  
**(2) Bash job scripts:**   
Each run-directory will be populated with two separate sets of bash job scripts for each task (approximately 50 jobs per forecast). The distribution of jobs per each step is as follows:  
lis_darun: 1; ldt_ics: 1; bcsd_fcst: 15; lis_fcst: 25; s2spost: 2; s2smetric: 3; and s2splots: 3.

The two different job scripts that are created for each S2S task are:
  
**\*.j files:** These are used with the SLURM job management system.    
**\*.sh files:** These contain NO SLURM directives, and are designed for **Cylc**. For example:  
*https://github.com/smahanam/FileSharing/blob/main/E2ESDIR/scratch/202501/s2spost/s2spost_01_run.sh*  
  
**(3) Cylc implementation:** This step will also generate the **CYLC_workflow.rc** file for **Cylc** implementation which defines directives, environmental variables, task dependencies, and the order of execution of job files. For example:  
*https://github.com/smahanam/FileSharing/blob/main/E2ESDIR/scratch/202501/CYLC_workflow.rc*  

### Optional Features for ghis2s_program.py:
**i) -s STEP**  
The STEP option allows the user to resume the process from the last completed step of the seven E2ES steps.  
-s STEP directs **ghis2s_program.py** to start from a specific STEP (valid inputs: LISDA, LDTICS, BCSD, FCST, POST, METRICS or PLOTS).  
**ii) -s STEP -o**  
*-o* flag used to run only the above -s STEP, the process will exit upon completion of above STEP.  
**iii) -j** for SLURM submission  
*-j* flag is used to submit jobs to the SLURM system.

## 3) Running Monthly Forecasts using *ghis2s_program.py* 
### a) In a Cylc Environment  
  
```
python ghis2s_program.py -y 2025 -m 1 -c s2s_config_global_fcast -e USER_EMAIL
```

In addition to creating monthly forecast-specific directories, links, and files, the command above also customizes and installs **Cylc-related** files and directories for the specified forecast month.  
**ghis2s_program.py** uses the WORKFLOW_NAME (e.g. **S2S-202501**) variable to name and organize these **Cylc-specific** resources, which include:  

i) A {WORKFLOW_NAME} directory under the user’s working directory (where ghis2s_program.py is executed).  
ii) A customized **Cylc flow.cylc** file containing the user’s email address, placed inside the {WORKFLOW_NAME} directory.  
iii) A **Cylc log/run** directory under the forecast month’s run directory (e.g., **scratch/202501**). Example:    
*https://github.com/smahanam/FileSharing/tree/main/E2ESDIR/scratch/202501/cylc-run/S2S-202501/run1*

Useful **Cylc** commands to run, monitor, inspect, and stop the workflow:    
```
Run {WORKFLOW_NAME}: cylc play {WORKFLOW_NAME}  
Monitor {WORKFLOW_NAME}: cylc tui {WORKFLOW_NAME}  
Show status {WORKFLOW_NAME}: cylc show {WORKFLOW_NAME}  
Stop {WORKFLOW_NAME}: cylc stop --now {WORKFLOW_NAME}  
Cat log: cylc cat-log {WORKFLOW_NAME}  
```
### End-to-end S2S Flow Chart
![End-to-end S2S Flow Chart](https://github.com/smahanam/FileSharing/blob/main/E2ESDIR/scratch/202402/CYLC-202402/ghis2s_workflow.png)

### b) In the SLURM System
  
```
python ghis2s_program.py -y 2025 -m 1 -c s2s_config_global_fcast -e USER_EMAIL -j
```  
Optionally, this command will submit all generated job scripts (~50 \*.j files) to the SLURM system. The jobs will be executed according to their predefined workflow dependencies.  
For example:  
*https://github.com/smahanam/FileSharing/blob/main/E2ESDIR/scratch/202501/SLURM_JOB_SCHEDULE*

## 4) Operational Notes and Cylc Design Rationale  

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
  








