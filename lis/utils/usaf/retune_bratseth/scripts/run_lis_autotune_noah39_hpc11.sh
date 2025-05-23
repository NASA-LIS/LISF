#!/bin/sh
#SBATCH --job-name=autotune
#SBATCH --time=1:00:00
#SBATCH --account=NWP601
#SBATCH --output autotune.slurm.out
#SBATCH --ntasks=5 --ntasks-per-node=1
#SBATCH --cluster-constraint=blue
#SBATCH --exclusive
#SBATCH --mem=0
#------------------------------------------------------------------------------
#
# SCRIPT: run_lis_autotune_noah39_hpc11.sh
#
# Batch script for running autotune software to update error covariance
# settings for LIS Air Force Bratseth scheme.  Customized for USAF HPC11
# supercomputer running SLURM batch queueing system. NOTE: Paths must
# be customized before running script.
#
# USAGE:  sbatch run_lis_autotune_noah39_hpc11.sh $YYYYMMDDHH $DD
#           where $YYYYMMDDHH is the end date/time of the OBA training period
#           and $DD is the length (in days) of the OBA training period.
#
#
# REVISION HISTORY:
# 05 Nov 2024: Eric Kemp. Initial revision.
#
#------------------------------------------------------------------------------

ulimit -s unlimited

# When a batch script is started, it starts in the user's home directory.
# Change to the directory where job was submitted.
if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# Environment
module use --append /ccs/home/emkemp/hpc11/privatemodules
module load lisf_7.6_prgenv_cray_8.5.0_cpe_23.12
module load afw-python/3.11-202406

# Define the variables to be processed
NWPVARS=(gage rh2m spd10m t2m)
SATVARS=(imerg)

# Paths on local system. Customize before running this script.
SCRIPTDIR=/lustre/storm/nwp601/proj-shared/emkemp/LISFV7.6/lis/noah39_autotune/scripts
CFGDIR=/lustre/storm/nwp601/proj-shared/emkemp/LISFV7.6/lis/noah39_autotune/cfgs
BINDIR=/lustre/storm/nwp601/proj-shared/emkemp/LISFV7.6/lis/noah39_autotune/bin

# Get the command line arguments to specify the training period.
if [ -z "$1" ] ; then
    echo "ERROR, Missing end date time for autotune software!"
    exit 1
fi
if [ -z "$2" ] ; then
    echo "ERROR, Missing training day range for autotune software!"
    exit 1
fi
enddt=$1
dayrange=$2

echo "INFO, Started autotuning at `date`"

# Customize config files for procOBA_NWP and procOBA_Sat, including blacklist
# creation. These will be run in the background to parallelize the work.
echo "---Customizing procOBA config files, and creating blacklists---"
if [ ! -e $SCRIPTDIR/customize_procoba_nwp.py ] ; then
    echo "ERROR, $SCRIPTDIR/customize_procoba_nwp.py does not exist!" && exit 1
fi
if [ ! -e $SCRIPTDIR/customize_procoba_sat.py ] ; then
    echo "ERROR, $SCRIPTDIR/customize_procoba_sat.py does not exist!" && exit 1
fi
if [ ! -e $CFGDIR/autotune.cfg ] ; then
    echo "ERROR, $CFGDIR/autotune.cfg does not exist!" && exit 1
fi
i=0
# We submit the NWP jobs first since they generate blacklists, which is
# time consuming.
for varname in "${NWPVARS[@]}" ; do
    echo "INFO, Task $i:  Calling customize_procoba_nwp.py for $varname at `date`"
    srun --ntasks=1 --nodes=1 --exclusive \
         $SCRIPTDIR/customize_procoba_nwp.py $CFGDIR/autotune.cfg \
         $enddt $dayrange $varname &
    PIDS+=($!)
    actives+=(1)
    ((i+=1))
    sleep 1
done
for varname in "${SATVARS[@]}" ; do
    echo "INFO, Task $i:  Calling customize_procoba_sat.py for $varname at `date`"
    srun --ntasks=1 --nodes=1 --exclusive \
         $SCRIPTDIR/customize_procoba_sat.py $CFGDIR/autotune.cfg \
         $enddt $dayrange $varname &
    PIDS+=($!)
    actives+=(1)
    ((i+=1))
    sleep 1
done

count=$i
echo "INFO, Waiting for $count task(s) to complete..."
while true; do
    for i in "${!PIDS[@]}"; do
        if [ "${actives[$i]}" -eq 0 ]; then
            continue
        fi
        pid="${PIDS[$i]}"
        # See if pid is still running
        ps --pid "$pid" > /dev/null
        if [ "$?" -ne 0 ]; then # ps doesn't see it, so it terminated
            wait "$pid"
            return_code="$?"
            actives[$i]=0
            count=$((count-1))
            if [ "${return_code}" -ne 0 ]; then
                echo "ERROR, Task $i failed, ABORTING at `date`"
                exit 1
            else
                echo "INFO, Task $i finished with no reported error at `date`"
            fi
            if [ "$count" -gt 0 ] ; then
                echo "INFO, Waiting for $count task(s) to complete..."
            fi
        fi
    done
    alldone=1
    for i in "${!actives[@]}" ; do
        if [ "${actives[$i]}" -eq 1 ] ; then
            alldone=0
        fi
    done
    if [ "${alldone}" -eq 1 ] ; then
        break
    fi
done
unset PIDS
unset actives

# Construct empirical semivariograms.
# These will be run in the background to parallelize the work.
echo "---Creating empirical semivariograms---"
if [ ! -e $BINDIR/procOBA_Sat ] ; then
    echo "ERROR, $BINDIR/procOBA_Sat does not exist!" && exit 1
fi
if [ ! -e $BINDIR/procOBA_NWP ] ; then
    echo "ERROR, $BINDIR/procOBA_NWP does not exist!" && exit 1
fi
i=0
for varname in "${SATVARS[@]}" ; do
    echo "INFO, Task $i:  Calling procOBA_Sat for $varname at `date`"
    if [ ! -e procOBA_Sat.$varname.config ] ; then
        echo "ERROR, procOBA_Sat.$varname.config does not exist!" && exit 1
    fi
    srun --ntasks=1 --nodes=1 --exclusive \
         $BINDIR/procOBA_Sat procOBA_Sat.$varname.config \
         procOBA_Sat.$varname.log &
    PIDS+=($!)
    actives+=(1)
    ((i+=1))
    sleep 1
done
for varname in "${NWPVARS[@]}" ; do
    echo "INFO, Task $i:  Calling procOBA_NWP for $varname at `date`"
    if [ ! -e procOBA_NWP.$varname.config ] ; then
        echo "ERROR, procOBA_NWP.$varname.config does not exist!" && exit 1
    fi
    srun --ntasks=1 --nodes=1 --exclusive \
         $BINDIR/procOBA_NWP procOBA_NWP.$varname.config \
         procOBA_NWP.$varname.log &
        pid=$!
    PIDS+=($!)
    actives+=(1)
    ((i+=1))
    sleep 1
done

count=$i
echo "INFO, Waiting for $count task(s) to complete..."
while true; do
    for i in "${!PIDS[@]}"; do
        if [ "${actives[$i]}" -eq 0 ]; then
            continue
        fi
        pid="${PIDS[$i]}"
        # See if pid is still running
        ps --pid "$pid" > /dev/null
        if [ "$?" -ne 0 ]; then # ps doesn't see it, so it terminated
            wait "$pid"
            return_code="$?"
            actives[$i]=0
            count=$((count-1))
            if [ "${return_code}" -ne 0 ]; then
                echo "ERROR, Task $i failed, ABORTING at `date`"
                exit 1
            else
                echo "INFO, Task $i finished with no reported error at `date`"
            fi
            if [ "$count" -gt 0 ] ; then
                echo "INFO, Waiting for $count task(s) to complete..."
            fi
        fi
    done
    alldone=1
    for i in "${!actives[@]}" ; do
        if [ "${actives[$i]}" -eq 1 ] ; then
            alldone=0
        fi
    done
    if [ "${alldone}" -eq 1 ] ; then
        break
    fi
done
unset PIDS
unset actives

# Create best-fits to the empirical semivariograms.
# These will be run in the background to parallelize the work.
echo "---Estimating best-fit covariance parameters---"
if [ ! -e $SCRIPTDIR/fit_semivariogram.py ] ; then
    echo "ERROR, $SCRIPTDIR/fit_semivariogram.py does not exist!" && exit 1
fi
i=0
for varname in "${SATVARS[@]}" ; do
    echo "INFO, Task $i:  Calling fit_semivariogram.py for $varname at `date`"
    if [ ! -e $CFGDIR/gage_$varname.cfg ] ; then
        echo "ERROR, $CFGDIR/gage_$varname.cfg does not exist!" && exit 1
    fi
    srun --ntasks=1 --nodes=1 --exclusive \
         $SCRIPTDIR/fit_semivariogram.py $CFGDIR/gage_$varname.cfg \
         gage_$varname.param &
    PIDS+=($!)
    actives+=(1)
    ((i+=1))
    sleep 1
done
for varname in "${NWPVARS[@]}" ; do
    echo "INFO, Task $i:  Calling fit_semivariogram.py for $varname at `date`"
    if [ $varname = gage ] ; then
        if [ ! -e $CFGDIR/${varname}_nwp.cfg ] ; then
            echo "ERROR, $CFGDIR/${varname}_nwp.cfg does not exist!" && exit 1
        fi
        srun --ntasks=1 --nodes=1 --exclusive \
             $SCRIPTDIR/fit_semivariogram.py $CFGDIR/${varname}_nwp.cfg \
            ${varname}_nwp.param &
    else
        if [ ! -e $CFGDIR/${varname}.cfg ] ; then
            echo "ERROR, $CFGDIR/${varname}.cfg does not exist!" && exit 1
        fi
        srun --ntasks=1 --nodes=1 --exclusive \
             $SCRIPTDIR/fit_semivariogram.py $CFGDIR/$varname.cfg \
            $varname.param &
    fi
    PIDS+=($!)
    actives+=(1)
    ((i+=1))
    sleep 1
done

count=$i
echo "INFO, Waiting for $count task(s) to complete..."
while true; do
    for i in "${!PIDS[@]}"; do
        if [ "${actives[$i]}" -eq 0 ]; then
            continue
        fi
        pid="${PIDS[$i]}"
        # See if pid is still running
        ps --pid "$pid" > /dev/null
        if [ "$?" -ne 0 ]; then # ps doesn't see it, so it terminated
            wait "$pid"
            return_code="$?"
            actives[$i]=0
            count=$((count-1))
            if [ "${return_code}" -ne 0 ]; then
                echo "ERROR, Task $i failed, ABORTING at `date`"
                exit 1
            else
                echo "INFO, Task $i finished with no reported error at `date`"
            fi
            if [ "$count" -gt 0 ] ; then
                echo "INFO, Waiting for $count task(s) to complete..."
            fi
        fi
    done
    alldone=1
    for i in "${!actives[@]}" ; do
        if [ "${actives[$i]}" -eq 1 ] ; then
            alldone=0
        fi
    done
    if [ "${alldone}" -eq 1 ] ; then
        break
    fi
done
unset PIDS
unset actives

# Rescale the satellite error variances to be w/r/t NWP.
# These will be run in the background to parallelize the work.
echo "---Rescaling satellite error variances---"
if [ ! -e $SCRIPTDIR/rescale_sat_sigma2.py ] ; then
    echo "ERROR, $SCRIPTDIR/rescale_sat_sigma2.py does not exist!" && exit 1
fi
i=0
for varname in "${SATVARS[@]}" ; do
    echo "INFO, Task $i:  Calling rescale_sat_sigma2.py for $varname at `date`"
    srun --ntasks=1 --nodes=1 --exclusive \
         $SCRIPTDIR/rescale_sat_sigma2.py $varname &
    PIDS+=($!)
    actives+=(1)
    ((i+=1))
    sleep 1
done

count=$i
echo "INFO, Waiting for $count task(s) to complete..."
while true; do
    for i in "${!PIDS[@]}"; do
        if [ "${actives[$i]}" -eq 0 ]; then
            continue
        fi
        pid="${PIDS[$i]}"
        # See if pid is still running
        ps --pid "$pid" > /dev/null
        if [ "$?" -ne 0 ]; then # ps doesn't see it, so it terminated
            wait "$pid"
            return_code="$?"
            actives[$i]=0
            count=$((count-1))
            if [ "${return_code}" -ne 0 ]; then
                echo "ERROR, Task $i failed, ABORTING at `date`"
                exit 1
            else
                echo "INFO, Task $i finished with no reported error at `date`"
            fi
            if [ "$count" -gt 0 ] ; then
                echo "INFO, Waiting for $count task(s) to complete..."
            fi
        fi
    done
    alldone=1
    for i in "${!actives[@]}" ; do
        if [ "${actives[$i]}" -eq 1 ] ; then
            alldone=0
        fi
    done
    if [ "${alldone}" -eq 1 ] ; then
        break
    fi
done
unset PIDS
unset actives

# Update the lis.config error settings.
# Single task, no need to parallelize.
echo "---Customizing lis.config file with new error covariances---"
echo "INFO, Calling customize_lis_config.py at `date`"
if [ ! -e $SCRIPTDIR/customize_lis_config.py ] ; then
    echo "ERROR, $SCRIPTDIR/customize_lis_config.py does not exist!" && exit 1
fi
srun --ntasks=1 --nodes=1 --exclusive --kill-on-bad-exit=1 \
     $SCRIPTDIR/customize_lis_config.py $CFGDIR/autotune.cfg \
     $enddt $dayrange || exit 1

# The end
echo "INFO, Completed Noah39 autotuning at `date`"
touch noah39.autotune.job.done
exit 0
