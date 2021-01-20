#!/bin/sh
#SBATCH --job-name=autotune
#SBATCH --time=1:00:00
#SBATCH --account s1189
#SBATCH --output autotune.slurm.out
#Adjust node, core, and hardware constraints here
#SBATCH --ntasks=4 --mem-per-cpu=4G --constraint=hasw
#Substitute your e-mail here
#SBATCH --mail-user=eric.kemp@nasa.gov
#SBATCH --mail-type=ALL
#Set quality of service, if needed.
#SBATCH --qos=debug
#------------------------------------------------------------------------------
#
# SCRIPT: run_autotune_discover.sh
#
# Batch script for running autotune software to update error covariance
# settings for LIS Air Force Bratseth scheme.  Customized for NASA Discover
# supercomputer running SLURM batch queueing system.
#
# USAGE:  sbatch run_autotune_discover.sh $YYYYMMDDHH $DD
#           where $YYYYMMDDHH is the end date/time of the OBA training period
#           and $DD is the length (in days) of the OBA training period.
#
#
# REVISION HISTORY:
# 17 Dec 2020:  Eric Kemp.  Initial specification.
#
#------------------------------------------------------------------------------

ulimit -s unlimited

# When a batch script is started, it starts in the user's home directory.
# Change to the directory where job was submitted.
if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# Environment
module purge
unset LD_LIBRARY_PATH
module use --append /home/emkemp/privatemodules
module load lisf_7_intel_19_1_3_304_traceback-work-around

# Define the variables to be processed
NWPVARS=(gage rh2m spd10m t2m)
SATVARS=(cmorph geoprecip imerg ssmi)

# Paths on local system
SCRIPTDIR=/discover/nobackup/emkemp/AFWA/autoretune/LISF/lis/utils/afwa/retune_bratseth/scripts
CFGDIR=/discover/nobackup/emkemp/AFWA/autoretune/LISF/lis/utils/afwa/retune_bratseth/cfgs
BINDIR=/discover/nobackup/emkemp/AFWA/autoretune/LISF/lis/utils/afwa/retune_bratseth/src

# Get the command line arguments to specify the training period.
if [ -z "$1" ] ; then
    echo "Missing end date time for autotune software!"
    exit 1
fi
if [ -z "$2" ] ; then
    echo "Missing training day range for autotune software!"
    exit 1
fi
enddt=$1
dayrange=$2

# Customize config files for procOBA_NWP, including blacklist creation.
# These will be run in the background to parallelize the work.
if [ ! -e $SCRIPTDIR/customize_procoba_nwp.py ] ; then
    echo "ERROR, $SCRIPTDIR/customize_procoba_nwp.py does not exist!" && exit 1
fi
if [ ! -e $CFGDIR/autotune.cfg ] ; then
    echo "ERROR, $CFGDIR/autotune.cfg does not exist!" && exit 1
fi
i=0
for varname in "${NWPVARS[@]}" ; do
    echo "Task $i:  Calling customize_procoba_nwp.py for $varname at `date`"
    srun --ntasks=1 \
         $SCRIPTDIR/customize_procoba_nwp.py $CFGDIR/autotune.cfg \
         $enddt $dayrange $varname &
    PIDS+=($!)
    ((i+=1))
    sleep 1
done
i=0
for pid in "${PIDS[@]}"; do
    wait ${pid}
    st=($?)
    if [[ ${st} -ne 0 ]]; then
        echo "ERROR, Task $i failed, ABORTING..."
        exit 1
    else
        echo "Task $i finished with no reported error"
    fi
    ((i+=1))
done
unset PIDS

# Construct empirical semivariograms.
# These will be run in the background to parallelize the work.
if [ ! -e $BINDIR/procOBA_NWP ] ; then
    echo "ERROR, $BINDIR/procOBA_NWP does not exist!" && exit 1
fi
i=0
for varname in "${NWPVARS[@]}" ; do
    echo "Task $i:  Calling procOBA_NWP for $varname at `date`"
    if [ ! -e procOBA_NWP.$varname.config ] ; then
        echo "ERROR, procOBA_NWP.$varname.config does not exist!" && exit 1
    fi
    srun --ntasks=1 \
         $BINDIR/procOBA_NWP procOBA_NWP.$varname.config \
         procOBA_NWP.$varname.log &
        pid=$!
    PIDS+=($!)
    ((i+=1))
    sleep 1
done
i=0
for pid in "${PIDS[@]}"; do
    wait ${pid}
    st=($?)
    if [[ ${st} -ne 0 ]]; then
        echo "ERROR, Task $i failed, ABORTING..."
        exit 1
    else
        echo "Task $i finished with no reported error"
    fi
    ((i+=1))
done
unset PIDS

# Create best-fits to the empirical semivariograms.
# These will be run in the background to parallelize the work.
if [ ! -e $SCRIPTDIR/fit_semivariogram.py ] ; then
    echo "ERROR, $SCRIPTDIR/fit_semivariogram.py does not exist!" && exit 1
fi
i=0
for varname in "${NWPVARS[@]}" ; do
    echo "Task $i:  Calling fit_semivariogram.py for $varname at `date`"
    if [ $varname = gage ] ; then
        if [ ! -e $CFGDIR/${varname}_nwp.cfg ] ; then
            echo "ERROR, $CFGDIR/${varname}_nwp.cfg does not exist!" && exit 1
        fi
        srun --ntasks=1 \
             $SCRIPTDIR/fit_semivariogram.py $CFGDIR/${varname}_nwp.cfg \
            ${varname}_nwp.param &
    else
        if [ ! -e $CFGDIR/${varname}.cfg ] ; then
            echo "ERROR, $CFGDIR/${varname}.cfg does not exist!" && exit 1
        fi
        srun --ntasks=1 \
             $SCRIPTDIR/fit_semivariogram.py $CFGDIR/$varname.cfg \
            $varname.param &
    fi
    PIDS+=($!)
    ((i+=1))
    sleep 1
done
i=0
for pid in "${PIDS[@]}"; do
    wait ${pid}
    st=($?)
    if [[ ${st} -ne 0 ]]; then
        echo "ERROR, Task $i failed, ABORTING..."
        exit 1
    else
        echo "Task $i finished with no reported error"
    fi
    ((i+=1))
done
unset PIDS

# Customize config files for procOBA_Sat.
# These will be run in the background to parallelize the work.
if [ ! -e $SCRIPTDIR/customize_procoba_sat.py ] ; then
    echo "ERROR, $SCRIPTDIR/customize_procoba_sat.py does not exist!" && exit 1
fi
i=0
for varname in "${SATVARS[@]}" ; do
    echo "Task $i:  Calling customize_procoba_sat.py for $varname at `date`"
    srun --ntasks=1 \
         $SCRIPTDIR/customize_procoba_sat.py $CFGDIR/autotune.cfg \
         $enddt $dayrange $varname &
    PIDS+=($!)
    ((i+=1))
    sleep 1
done
i=0
for pid in "${PIDS[@]}"; do
    wait ${pid}
    st=($?)
    if [[ ${st} -ne 0 ]]; then
        echo "ERROR, Task $i failed, ABORTING..."
        exit 1
    else
        echo "Task $i finished with no reported error"
    fi
    ((i+=1))
done
unset PIDS


# Construct empirical semivariograms for the satellite data w/r/t gages.
# These will be run in the background to parallelize the work.
if [ ! -e $BINDIR/procOBA_Sat ] ; then
    echo "ERROR, $BINDIR/procOBA_Sat does not exist!" && exit 1
fi
i=0
for varname in "${SATVARS[@]}" ; do
    echo "Task $i:  Calling procOBA_Sat for $varname at `date`"
    if [ ! -e procOBA_Sat.$varname.config ] ; then
        echo "ERROR, procOBA_Sat.$varname.config does not exist!" && exit 1
    fi
    srun --ntasks=1 \
         $BINDIR/procOBA_Sat procOBA_Sat.$varname.config \
         procOBA_Sat.$varname.log &
    PIDS+=($!)
    ((i+=1))
    sleep 1
done
i=0
for pid in "${PIDS[@]}"; do
    wait ${pid}
    st=($?)
    if [[ ${st} -ne 0 ]]; then
        echo "ERROR, Task $i failed, ABORTING..."
        exit 1
    else
        echo "Task $i finished with no reported error"
    fi
    ((i+=1))
done
unset PIDS

# Create best fits to the satellite semivariograms w/r/t gages.
# These will be run in the background to parallelize the work.
i=0
for varname in "${SATVARS[@]}" ; do
    echo "Task $i:  Calling fit_semivariogram.py for $varname at `date`"
    if [ ! -e $CFGDIR/gage_$varname.cfg ] ; then
        echo "ERROR, $CFGDIR/gage_$varname.cfg does not exist!" && exit 1
    fi
    srun --ntasks=1 \
         $SCRIPTDIR/fit_semivariogram.py $CFGDIR/gage_$varname.cfg \
         gage_$varname.param &
    PIDS+=($!)
    ((i+=1))
    sleep 1
done
i=0
for pid in "${PIDS[@]}"; do
    wait ${pid}
    st=($?)
    if [[ ${st} -ne 0 ]]; then
        echo "ERROR, Task $i failed, ABORTING..."
        exit 1
    else
        echo "Task $i finished with no reported error"
    fi
    ((i+=1))
done
unset PIDS

# Rescale the satellite error variances to be w/r/t NWP.
# These will be run in the background to parallelize the work.
if [ ! -e $SCRIPTDIR/rescale_sat_sigma2.py ] ; then
    echo "ERROR, $SCRIPTDIR/rescale_sat_sigma2.py does not exist!" && exit 1
fi
i=0
for varname in "${SATVARS[@]}" ; do
    echo "Task $i:  Calling rescale_sat_sigma2.py for $varname at `date`"
    srun --ntasks=1 \
         $SCRIPTDIR/rescale_sat_sigma2.py $varname &
    PIDS+=($!)
    ((i+=1))
    sleep 1
done
i=0
for pid in "${PIDS[@]}"; do
    wait ${pid}
    st=($?)
    if [[ ${st} -ne 0 ]]; then
        echo "ERROR, Task $i failed, ABORTING..."
        exit 1
    else
        echo "Task $i finished with no reported error"
    fi
    ((i+=1))
done
unset PIDS

# Update the lis.config error settings.
# Single task, no need to parallelize.
echo "Calling customize_lis_config.py at `date`"
if [ ! -e $SCRIPTDIR/customize_lis_config.py ] ; then
    echo "ERROR, $SCRIPTDIR/customize_lis_config.py does not exist!" && exit 1
fi
srun --ntasks=1 --kill-on-bad-exit=1 \
     $SCRIPTDIR/customize_lis_config.py $CFGDIR/autotune.cfg \
     $enddt $dayrange || exit 1

# The end
echo "INFO, Completed autotuning!"
exit 0
