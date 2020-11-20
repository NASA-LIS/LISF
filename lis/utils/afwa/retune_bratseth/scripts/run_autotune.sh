#!/bin/sh
#SBATCH --job-name=autotune
#SBATCH --time=1:00:00
##SBATCH --time=3:00:00
#SBATCH --account s1189
#SBATCH --output autotune.slurm.out
#Adjust node, core, and hardware constraints here
#SBATCH --ntasks=1 --constraint=hasw
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

# First, handle NWP semivariographic analyses
for varname in gage rh2m spd10m t2m ; do

    echo "Calling customize_procoba_nwp.py for $varname"
    $SCRIPTDIR/customize_procoba_nwp.py $CFGDIR/autotune.cfg \
                                        $enddt $dayrange $varname || exit 1

    echo "Calling procOBA_NWP for $varname"
    mpirun -np $SLURM_NTASKS $BINDIR/procOBA_NWP \
           procOBA_NWP.$varname.config || exit 1

    echo "Calling fit_semivariogram.py for $varname"
    if [ $varname = gage ] ; then
        $SCRIPTDIR/fit_semivariogram.py $CFGDIR/${varname}_nwp.cfg \
                                        ${varname}_nwp.param || exit 1
    else
        $SCRIPTDIR/fit_semivariogram.py $CFGDIR/$varname.cfg \
                                        $varname.param || exit 1
    fi

done


# Next, handle satellite variographic analyses
for varname in cmorph geoprecip imerg ssmi ; do

    echo "Calling customize_procoba_sat.py for $varname"
    $SCRIPTDIR/customize_procoba_sat.py $CFGDIR/autotune.cfg \
                                        $enddt $dayrange $varname || exit 1

    echo "Calling procOBA_Sat for $varname"
    mpirun -np $SLURM_NTASKS $BINDIR/procOBA_Sat \
           procOBA_Sat.$varname.config || exit 1

    echo "Calling fit_semivariogram.py for $varname"
    $SCRIPTDIR/fit_semivariogram.py $CFGDIR/gage_$varname.cfg \
                                    gage_$varname.param || exit 1

    # Rescale the satellite error variance to compare with NWP.
    echo "Calling rescale_sat_sigma2.py for $varname"
    $SCRIPTDIR/rescale_sat_sigma2.py $varname || exit 1
done

# Update the lis.config error settings
echo "Calling customize_lis_config.py"
$SCRIPTDIR/customize_lis_config.py $CFGDIR/autotune.cfg \
    $enddt $dayrange || exit 1

# The end
exit 0
