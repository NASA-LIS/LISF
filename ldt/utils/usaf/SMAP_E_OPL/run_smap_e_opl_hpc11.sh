#!/bin/sh
#SBATCH --job-name=smapeopl
#SBATCH --time=0:30:00
#SBATCH --account=NWP601
#SBATCH --output smapeopl.slurm.out
#SBATCH --ntasks=13 --ntasks-per-node=13
#SBATCH --cluster-constraint=blue
#SBATCH --exclusive
#SBATCH --mem=0
#--------------------------------------------------------------------------
#
# SCRIPT: run_smap_e_opl_hpc11.sh
#
# Batch script for running LDT to generate SMAP_E_OPL retrievals with
# LDT.  Customized for USAF HPC11 supercomputer running SLURM batch
# queueing system.
#
# USAGE: run_smap_e_opl_hpc11.sh $startdate $starthour $enddate \
#          $endhour $lsm
#          where $startdate and $enddate specify the inclusive UTC
#          date range to run SMAP_E_OPL retrievals (formatted YYYY-MM-DD
#          or YYYYMMDD); $starthour and $endhour specify the hours of the
#          respective dates (formatted HH); and $lsm is the LSM name, used
#          to select the ldt.config template file.
#
# Based on Korn shell script provided by Pang-Wei Liu.
#
# REVISION HISTORY:
# 29 Mar 2024: Eric Kemp.  Initial specification, based on Discover script.
#
#--------------------------------------------------------------------------

ulimit -s unlimited

# When a batch script is started, it starts in the user's home directory.
# Change to the directory where job was submitted.
if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# Environment
module use --append /ccs/home/emkemp/hpc11/privatemodules
module load lisf_7.6_prgenv_cray_8.5.0_cpe_23.12_draft
module load afw-python/3.10-202312

# Paths on local system
SCRIPTDIR=/lustre/storm/nwp601/proj-shared/emkemp/lis76_smap_e_opl/scripts
BINDIR=/lustre/storm/nwp601/proj-shared/emkemp/lis76_smap_e_opl/bin
TMPLDIR=/lustre/storm/nwp601/proj-shared/emkemp/lis76_smap_e_opl/tmpl

# Get the command line arguments to specify the training period
if [ -z "$1" ] ; then
    echo "ERROR, Missing start date for SMAP_E_OPL retrievals!"
    exit 1
fi
if [ -z "$2" ] ; then
    echo "ERROR, Missing start hour for SMAP_E_OPL retrievals!"
    exit 1
fi
if [ -z "$3" ] ; then
    echo "ERROR, Missing end date for SMAP_E_OPL retrievals!"
    exit 1
fi
if [ -z "$4" ] ; then
    echo "ERROR, Missing end hour for SMAP_E_OPL retrievals!"
    exit 1
fi
if [ -z "$5" ] ; then
    echo "ERROR, missing LSM option!" && exit 1
fi

# Use the command line arguments to set start and end datetime limits.
# Sanity check and report errors if found.

istart="$1 $2"
iend="$3 $4"
start_dt=$(date "+%F %H" -d "$istart")
if [ "$?" -ne 0 ] ; then
    echo "ERR, Bad startdate and hour" && exit 1
fi
end_dt=$(date "+%F %H" -d "$iend")
if [ "$?" -ne 0 ] ; then
    echo "ERR, Bad enddate and hour" && exit 1
fi
d="$start_dt"
cur_yyyymmddhh=$(date -d "$start_dt" +%Y%m%d%H)
end_yyyymmddhh=$(date -d "$end_dt" +%Y%m%d%H)

# Use the LSM command line argument to select the appropriate ldt.config
# template.
lsm="$5"
if [ "$lsm" = "noah" ] ; then
    tmplfile="$TMPLDIR/ldt.config.smapeopl.noah.tmpl"
elif [ "$lsm" = "noahmp" ] ; then
    tmplfile="$TMPLDIR/ldt.config.smapeopl.noahmp.tmpl"
elif [ "$lsm" = "jules" ] ; then
    tmplfile="$TMPLDIR/ldt.config.smapeopl.jules.tmpl"
else
    echo "ERR, invalid LSM option, expected noah, noahmp, or jules; got $lsm"
    exit 1
fi
if [ ! -e $tmplfile ] ; then
    echo "ERR, $tmplfile not found!" && exit 1
fi

# Other sanity checks
if [ ! -e $BINDIR/LDT ] ; then
    echo "ERR, $BINDIR/LDT does not exist!" && exit 1
fi
if [ ! -e $SCRIPTDIR/customize_ldt_smapeopl_config.py ] ; then
    echo "ERR, $SCRIPTDIR/customize_ldt_smapeopl_config.py not found!"
    exit 1
fi

# Loop through all times.  This requires parallel invocations using
# srun.
echo "INFO, started SMAP_E_OPL retrievals at `date`"
i=0
while [ "$cur_yyyymmddhh" -le "$end_yyyymmddhh" ] ; do

    echo "INFO, working on $cur_yyyymmddhh at `date`"

    # Update ldt.config file
    $SCRIPTDIR/customize_ldt_smapeopl_config.py "$cur_yyyymmddhh" \
                 "$tmplfile" $i || exit 1

    # Run LDT
    srun --ntasks=1 --cpus-per-task=1 \
         $BINDIR/LDT ldt.config.smapeopl."$cur_yyyymmddhh" &

    PIDS+=($!)
    actives+=(1)
    yyyymmddhh+=("$cur_yyyymmddhh")
    ((i+=1))
    sleep 1 # Don't overwhelm SLURM

    # Update valid time
    d=$(date -d "$d + 1 hours")
    cur_yyyymmddhh=$(date -d "$d" +%Y%m%d%H)

done

# Wait until all the jobs complete
count=$i
echo "INFO, waiting for $count task(s) to complete..."
while true; do
    for i in "${!PIDS[@]}"; do
        if [ "${actives[$i]}" -eq 0 ] ; then
            continue
        fi
        pid="${PIDS[$i]}"
        # See if pid is still running
        echo `ps --pid "$pid"`

        ps --pid "$pid" > /dev/null
        if [ "$?" -ne 0 ] ; then # ps doesn't see it, so it terminated
           wait "$pid"
           return_code="$?"
           actives[$i]=0
           count=$((count-1))
           if [ "${return_code}" -ne 0 ] ; then
               echo "WARN, SMAP_E_OPL ${yyyymmddhh[$i]} failed around `date`"
           else
               echo "INFO, SMAP_E_OPL ${yyyymmddhh[$i]} finished with no reported error at `date`"
           fi
           if [ "$count" -gt 0 ] ; then
               echo "INFO, Waiting for $count task(s) to complete..."
           fi
        fi
    done
    # See if we are all done.  Assume we are, and correct if we are not.
    alldone=1
    for i in "${!actives[@]}"; do
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
unset yyyymmddhh

# The end
touch smapeopl.job.done
exit 0
