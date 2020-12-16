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

NWPVARS=(gage rh2m spd10m t2m)
SATVARS=(cmorph geoprecip imerg ssmi)

SCRIPTDIR=/discover/nobackup/emkemp/AFWA/autoretune/LISF/lis/utils/afwa/retune_bratseth/scripts
CFGDIR=/discover/nobackup/emkemp/AFWA/autoretune/LISF/lis/utils/afwa/retune_bratseth/cfgs
BINDIR=/discover/nobackup/emkemp/AFWA/autoretune/LISF/lis/utils/afwa/retune_bratseth/src

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

# Customize config files for procOBA_NWP, including blacklist creation
if [ ! -e $SCRIPTDIR/customize_procoba_nwp.py ] ; then
    echo "ERROR, $SCRIPTDIR/customize_procoba_nwp.py does not exist!" && exit 1
fi
if [ ! -e $CFGDIR/autotune.cfg ] ; then
    echo "ERROR, $CFGDIR/autotune.cfg does not exist!" && exit 1
fi
for varname in "${NWPVARS[@]}" ; do
    echo `date`
    echo "Calling customize_procoba_nwp.py for $varname"
    srun --label --ntasks=1 --kill-on-bad-exit=1 \
         $SCRIPTDIR/customize_procoba_nwp.py $CFGDIR/autotune.cfg \
         $enddt $dayrange $varname &
    sleep 1
done
wait

# Next, construct empirical semivariograms
if [ ! -e $BINDIR/procOBA_NWP ] ; then
    echo "ERROR, $BINDIR/procOBA_NWP does not exist!" && exit 1
fi
for varname in "${NWPVARS[@]}" ; do
    echo `date`
    echo "Calling procOBA_NWP for $varname"
    if [ ! -e procOBA_NWP.$varname.config ] ; then
        echo "ERROR, procOBA_NWP.$varname.config does not exist!" && exit 1
    fi
    srun --label --ntasks=1 --kill-on-bad-exit=1 \
         $BINDIR/procOBA_NWP procOBA_NWP.$varname.config \
         procOBA_NWP.$varname.log &
    sleep 1
done
wait

# Now fit the semivariogram functions
if [ ! -e $SCRIPTDIR/fit_semivariogram.py ] ; then
    echo "ERROR, $SCRIPTDIR/fit_semivariogram.py does not exist!" && exit 1
fi
for varname in "${NWPVARS[@]}" ; do
    echo `date`
    echo "Calling fit_semivariogram.py for $varname"
    if [ $varname = gage ] ; then
        if [ ! -e $CFGDIR/${varname}_nwp.cfg ] ; then
            echo "ERROR, $CFGDIR/${varname}_nwp.cfg does not exist!" && exit 1
        fi
        srun --label --ntasks=1 --kill-on-bad-exit=1 \
             $SCRIPTDIR/fit_semivariogram.py $CFGDIR/${varname}_nwp.cfg \
            ${varname}_nwp.param &
    else
        if [ ! -e $CFGDIR/${varname}.cfg ] ; then
            echo "ERROR, $CFGDIR/${varname}.cfg does not exist!" && exit 1
        fi
        srun --label --ntasks=1 --kill-on-bad-exit=1 \
             $SCRIPTDIR/fit_semivariogram.py $CFGDIR/$varname.cfg \
            $varname.param &
    fi
    sleep 1
done
wait

# Next, customize config files for procOBA_Sat
if [ ! -e $SCRIPTDIR/customize_procoba_sat.py ] ; then
    echo "ERROR, $SCRIPTDIR/customize_procoba_sat.py does not exist!" && exit 1
fi
for varname in "${SATVARS[@]}" ; do
    echo `date`
    echo "Calling customize_procoba_sat.py for $varname"
    srun --label --ntasks=1 --kill-on-bad-exit=1 \
         $SCRIPTDIR/customize_procoba_sat.py $CFGDIR/autotune.cfg \
         $enddt $dayrange $varname &
    sleep 1
done
wait

# Next, construct empirical semivariograms
if [ ! -e $BINDIR/procOBA_Sat ] ; then
    echo "ERROR, $BINDIR/procOBA_Sat does not exist!" && exit 1
fi
for varname in "${SATVARS[@]}" ; do
    echo `date`
    echo "Calling procOBA_Sat for $varname"
    if [ ! -e procOBA_Sat.$varname.config ] ; then
        echo "ERROR, procOBA_Sat.$varname.config does not exist!" && exit 1
    fi
    srun --label --ntasks=1 --kill-on-bad-exit=1 \
         $BINDIR/procOBA_Sat procOBA_Sat.$varname.config \
         procOBA_Sat.$varname.log &
    sleep 1
done
wait

# Next, fit semivariograms for satellite data
for varname in "${SATVARS[@]}" ; do
    echo `date`
    echo "Calling fit_semivariogram.py for $varname"
    if [ ! -e $CFGDIR/gage_$varname.cfg ] ; then
        echo "ERROR, $CFGDIR/gage_$varname.cfg does not exist!" && exit 1
    fi
    srun --label --ntasks=1 --kill-on-bad-exit=1 \
         $SCRIPTDIR/fit_semivariogram.py $CFGDIR/gage_$varname.cfg \
         gage_$varname.param &
    sleep 1
done
wait

# Next, rescale the satellite error variance to compare with NWP.
if [ ! -e $SCRIPTDIR/rescale_sat_sigma2.py ] ; then
    echo "ERROR, $SCRIPTDIR/rescale_sat_sigma2.py does not exist!" && exit 1
fi
for varname in "${SATVARS[@]}" ; do
    echo `date`
    echo "Calling rescale_sat_sigma2.py for $varname"
    srun --label --ntasks=1 --kill-on-bad-exit=1 \
         $SCRIPTDIR/rescale_sat_sigma2.py $varname &
    sleep 1
done
wait

# Update the lis.config error settings
echo `date`
echo "Calling customize_lis_config.py"
if [ ! -e $SCRIPTDIR/customize_lis_config.py ] ; then
    echo "ERROR, $SCRIPTDIR/customize_lis_config.py does not exist!" && exit 1
fi
$SCRIPTDIR/customize_lis_config.py $CFGDIR/autotune.cfg \
                                   $enddt $dayrange || exit 1

# The end
echo "INFO, Completed autotuning!"
exit 0
