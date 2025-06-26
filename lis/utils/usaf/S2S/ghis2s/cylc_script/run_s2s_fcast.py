import os
import sys
import shutil
import subprocess
import argparse
from ghis2s.s2s_app.s2s_run import S2Srun

CYLCHOME = '/home/user/cylc8-workflows/'
E2ESDIR = '/discover/nobackup/projects/ghilis/S2S/GLOBAL/Cylc_test2/'
WORKFLOW_NAME = 'S2S'
 
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', required=True, type=str, help='config file')
parser.add_argument('-y', '--year', required=True, type=int, help='forecast year')
parser.add_argument('-m', '--month', required=True, type=int, help='forecast month')
parser.add_argument('-e', '--email', required=True, type=str, help='user email')
parser.add_argument('-s', '--step', required=False, default=None, type=str, help='S2S step: LISDA, LDTICS, BCSD, FCST, POST, METRICS or PLOTS')
parser.add_argument('-o', '--one_step', action='store_true', help='Is only one step (default: False)?')
parser.add_argument('-j', '--submit_job', action='store_true', help='Submit SLURM jobs (default: False -> CYLC)?')

args = parser.parse_args()

# go to E2ESDIR
# -------------
LOGDIR = E2ESDIR + f'scratch/{args.year:04d}{args.month:02d}/' 
print(f'Creating working directories, and job files for the S2S forecast initialized on {args.year:04d}-{args.month:02d}-01')
os.chdir(E2ESDIR)
s2s = S2Srun(year=args.year, month=args.month, config_file=args.config_file)

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
else:
    s2s.main()

# Write CYLC workflow runtime snippet
# -----------------------------------
s2s.write_cylc_snippet()

# Submit SLURM jobs
# -----------------
if  args.submit_job:
    s2s.submit_jobs()
    sys.exit()

# Return to CYLCHOME
# ------------------
os.chdir(CYLCHOME)
os.makedirs(WORKFLOW_NAME, exist_ok=True)
os.chdir(WORKFLOW_NAME)
shutil.copy(LOGDIR + 'CYLC_workflow.rc', 'flow.cylc')
result = subprocess.run(f'chmod 755 {LOGDIR}/*/*sh', shell=True, 
                        capture_output=True, text=True)

# Update email
# ------------
with open('flow.cylc', 'r') as file:
    filedata = file.read()
filedata = filedata.replace('USEREMAIL', args.email)
filedata = filedata.replace('YYYY-MM', f'{args.year:04d}-{args.month:02d}')
filedata = filedata.replace('--output ', '--output=')
filedata = filedata.replace('--error ', '--error=')

with open('flow.cylc', 'w') as file:
    file.write(filedata)

# Install Cylc WORKFLOW_NAME
# --------------------------
command = f"cylc install --symlink-dirs=run={LOGDIR}"
process = subprocess.run(command, shell=True)

# Run Cylc
# --------
print(f'Useful CYLC commands from {os.getcwd()}')
print("   ")
print(f'Run {WORKFLOW_NAME}: cylc play {WORKFLOW_NAME}')
print(f'Monitor {WORKFLOW_NAME}: cylc tui {WORKFLOW_NAME}')
print(f'Show status {WORKFLOW_NAME}: cylc show {WORKFLOW_NAME}')
print(f'Stop {WORKFLOW_NAME}: cylc stop --now {WORKFLOW_NAME}')
print(f'Cat log: cylc cat-log {WORKFLOW_NAME}')






