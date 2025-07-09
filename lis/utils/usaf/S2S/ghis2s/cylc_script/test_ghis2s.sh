# Set required environment variables
export CONFIG_FILE="s2s_config_global_fcast"
export FORECAST_YEAR=2025
export FORECAST_MONTH=1
export USER_EMAIL="sarith.p.mahanama@nasa.gov"
export PYTHONPATH="/discover/nobackup/projects/usaf_lis/smahanam/S2S/LISF-1/lis/utils/usaf/S2S/"
export E2ESDIR="/discover/nobackup/projects/ghilis/smahanam/ghi-coupling/"

# Optional variables
export S2S_STEP="E2E"
export ONE_STEP=False
export SUBMIT_JOB=False

source /etc/profile.d/modules.sh
module purge
source /home/$USER/venvs/cylc-env/bin/activate
module use -a /discover/nobackup/projects/usaf_lis/smahanam/S2S/LISF-1//env/discover/
module --ignore-cache load lisf_7.5_intel_2023.2.1_s2s

# Run the program
cd /discover/nobackup/projects/ghilis/smahanam/GHI-repos/ghi-apps/bin/
python ghis2s_program.py
