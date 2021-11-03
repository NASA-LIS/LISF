#!/bin/sh
#SBATCH --job-name=lisdarun
#SBATCH --time=2:00:00
#SBATCH --ntasks=28
#SBATCH --account s1189
#SBATCH --output lisdarun.slurm.out
#SBATCH --mail-user=kristi.r.arsenault@nasa.gov
#SBATCH --mail-type=ALL
#------------------------------------------------------------------------------
#
# SCRIPT:  run_lis_darun_nrt_slurm.sh
#
# DESCRIPTION: Batch script for running LIS data assimilation (DA) simulation  
#  to generate the NoahMP401 and HYMAP initial conditions S2S LIS forecast runs. 
#  LIS is run for 4-days, but in operations, the suggestion is to run the 
#  LIS DA simulation for the entire past month, when 00Z of the first of each 
#  month is reached (e.g., setup in that way via cron).
#
# USAGE: sbatch run_lis_darun_nrt_slurm.sh $STYYYYMMDD  $EDYYYYMMDD
#
#     where $STYYYYMMDD is the start date of the LIS run,
#       and $EDYYYYMMDD is the end date of the LIS run.
#
# REVISION HISTORY:
# 22 Sep 2021: Eric Kemp (SSAI), first version.
# 26 Oct 2021: Kristi Arsenault (SAIC), modified for LIS DA run.
# 03 Nov 2021: Kristi Arsenault (SAIC), updated with python config file entries
#------------------------------------------------------------------------------

ulimit -s unlimited

# When a batch script is started, it starts in user's home directory.
# Change to the directory where job was submitted.
if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# Enviroment
module purge
unset LD_LIBRARY_PATH
module use --append /home/karsenau/privatemodules
module load lisf_7_intel_19_1_3_304

# Paths on local system
SCRIPTDIR=/discover/nobackup/projects/E2ES_Test/lis_darun
LISCFGDIR=$SCRIPTDIR/lis.config_files/
PYCFGFILE=lis_darun.cfg
OUTDIR=./output


# Get the shell command line arguments:
if [ -z "$1" ] ; then
    echo "[ERR] Missing the start date (YYYYMMDD) of the LIS run!"
    exit 1
fi
STYYYYMMDD=$1

if [ -z "$2" ] ; then
    echo "[ERR] Missing the end date (YYYYMMDD) of the LIS run!"
    exit 1
fi
EDYYYYMMDD=$2


# Customize lis.config file, with starting and end dates:
if [ ! -e $SCRIPTDIR/generate_lisda_config_nrt.py ] ; then
    echo "[ERR], $SCRIPTDIR/generate_lisda_config_nrt.py does not exist!" && exit 1
fi
echo "[INFO] Customizing lis.config..."
$SCRIPTDIR/generate_lisda_config_nrt.py $PYCFGFILE $STYYYYMMDD $EDYYYYMMDD || exit 1


# Run LIS
if [ ! -e ./LIS ] ; then
    echo "[ERR] ./LIS executable does not exist!" && exit 1
fi
echo "[INFO] Running LIS..."

cp $LISCFGDIR/lis.config_darun_${STYYYYMMDD} $SCRIPTDIR/lis.config
mpirun -np $SLURM_NTASKS ./LIS || exit 1

echo "[INFO] Completed LIS run!"
exit 0
