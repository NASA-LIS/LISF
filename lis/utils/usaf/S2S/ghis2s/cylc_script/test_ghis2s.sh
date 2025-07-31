#!/bin/bash

# Set required environment variables
export CONFIG_FILE="s2s_config_global_fcast"
export FORECAST_YEAR=2025
export FORECAST_MONTH=3
export USER_EMAIL="kristi.r.arsenault@nasa.gov"
#export USER_EMAIL="sarith.p.mahanama@nasa.gov"
export E2ESDIR="/discover/nobackup/projects/ghilis/S2S/GLOBAL/cylc_test1/"
export OUTPUT_ROOT="/discover/nobackup/projects/ghilis/S2S/GHI-repos/s2s_workflows/"
export LISFDIR=`grep LISFDIR $E2ESDIR/$CONFIG_FILE | cut -d':' -f2 | tr -d "[:space:]"`
export LISFMOD=`grep LISFMOD $E2ESDIR/$CONFIG_FILE | cut -d':' -f2 | tr -d "[:space:]"`
export PYTHONPATH="${LISFDIR}/lis/utils/usaf/S2S/"
mkdir -p -m 775 $OUTPUT_ROOT

# Optional variables
export S2S_STEP="POST"
export ONE_STEP=false
export SUBMIT_JOB=false

source /etc/profile.d/modules.sh
module purge
module use -a "${LISFDIR}/env/discover/"
module --ignore-cache load $LISFMOD

# Run S2S Python program
cd /discover/nobackup/projects/ghilis/S2S/GHI-repos/ghi-apps/bin/
python ghis2s_program.py

# if SUBMIT JOBS to SLURM exit
[ "$SUBMIT_JOB" = "true" ] && exit

# Run Cylc
module purge
module use -a "${LISFDIR}/env/discover/"
export USE_CYLC_ENV=1
module --ignore-cache load $LISFMOD

if [ $FORECAST_MONTH -lt 10 ]; then
    MM="0${FORECAST_MONTH}"
else
    MM="${FORECAST_MONTH}"
fi

WORKFLOW_NAME="S2S-${FORECAST_YEAR}${MM}"
WORKFLOW_DIR="${OUTPUT_ROOT}/${WORKFLOW_NAME}"
LOGDIR="${E2ESDIR}/scratch/${FORECAST_YEAR}${MM}"

# Install CYLC workflow
cd $WORKFLOW_DIR
CWD=`pwd`
cylc install --symlink-dirs=run=$LOGDIR

echo

echo "======================================================================================================="
echo "Useful CYLC commands from ${CWD}"
echo "======================================================================================================="
echo
echo "First load Cylc module"
echo
echo "setenv USE_CYLC_ENV 1"
echo "module use -a ${LISFDIR}/env/discover/"
echo "module --ignore-cache load $LISFMOD"
echo
echo "Run ${WORKFLOW_NAME}: cylc play ${WORKFLOW_NAME}"
echo "Monitor ${WORKFLOW_NAME}: cylc tui ${WORKFLOW_NAME}"
echo "Show status ${WORKFLOW_NAME}: cylc show ${WORKFLOW_NAME}"
echo "Stop ${WORKFLOW_NAME}: cylc stop --now ${WORKFLOW_NAME}"
echo "Cat log: cylc cat-log ${WORKFLOW_NAME}"
echo
echo "You might want to module purge and unsetenv USE_CYLC_ENV after running above Cylc commands." 


