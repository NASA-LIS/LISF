#!/bin/bash

# Set required environment variables
export CONFIG_FILE="s2s_config_global_fcast"
export FORECAST_YEAR=2025
export FORECAST_MONTH=1
export USER_EMAIL="sarith.p.mahanama@nasa.gov"
export E2ESDIR="/discover/nobackup/projects/ghilis/smahanam/ghi-coupling/"
export OUTPUT_ROOT="/discover/nobackup/projects/ghilis/smahanam/GHI-repos/s2s_workflows/"
export LISFDIR=`grep LISFDIR $E2ESDIR/$CONFIG_FILE | cut -d':' -f2 | tr -d "[:space:]"`
export LISFMOD=`grep LISFMOD $E2ESDIR/$CONFIG_FILE | cut -d':' -f2 | tr -d "[:space:]"`
export PYTHONPATH="${LISFDIR}/lis/utils/usaf/S2S/"
mkdir -p -m 775 $OUTPUT_ROOT

# Optional variables
export S2S_STEP="E2E"
export ONE_STEP=false
export SUBMIT_JOB=false

source /etc/profile.d/modules.sh
module purge
module use -a "${LISFDIR}/env/discover/"
module --ignore-cache load $LISFMOD

# Run S2S Python program
cd /discover/nobackup/projects/ghilis/smahanam/GHI-repos/ghi-apps/bin/
python ghis2s_program.py

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
CWD=`pwd`

# Install CYLC workflow
cd $WORKFLOW_DIR
cylc install --symlink-dirs=run=$LOGDIR

echo "Useful CYLC commands from ${CWD}"
echo
echo "Run ${WORKFLOW_NAME}: cylc play ${WORKFLOW_NAME}"
echo "Monitor ${WORKFLOW_NAME}: cylc tui ${WORKFLOW_NAME}"
echo "Show status ${WORKFLOW_NAME}: cylc show ${WORKFLOW_NAME}"
echo "Stop ${WORKFLOW_NAME}: cylc stop --now ${WORKFLOW_NAME}"
echo "Cat log: cylc cat-log ${WORKFLOW_NAME}"



