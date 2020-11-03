#!/bin/sh
#SBATCH --job-name=autotune
#SBATCH --time=1:00:00
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
module load lisf_7_intel_19_1_0_166

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

# First, handle non-precipitation
# for varname in rh2m spd10m t2m ; do
#
#     $SCRIPTDIR/customize_procoba_nwp.py $CFGDIR/autotune.cfg \
#                                         $varname $enddt $dayrange || exit 1
#
#     mpirun -np $SLURM_NTASKS $BINDIR/procOBA_NWP \
#            procOBA_NWP.$varname.config || exit 1
#
#     $SCRIPTDIR/fit_semivariogram.py $CFGDIR/$varname.cfg \
#                                     $varname.param || exit 1
#
# done

# # Next, handle gages versus NWP
#
# $SCRIPTDIR/customize_procoba_nwp.py $CFGDIR/autotune.cfg \
#                                     gage $enddt $dayrange || exit 1
#
# mpirun -np $SLURM_NTASKS $BINDIR/procOBA_NWP \
#        procOBA_NWP.gage.config || exit 1
#
# $SCRIPTDIR/fit_semivariogram.py $CFGDIR/gage_nwp.cfg \
#                                 gage_nwp.param || exit 1


# Next, handle gages versus satellite data
#for varname in cmorph geoprecip imerg ssmi ; do
for varname in imerg ; do

    $SCRIPTDIR/customize_procoba_sat.py $CFGDIR/autotune.cfg \
                                        $varname $enddt $dayrange || exit 1

    mpirun -np $SLURM_NTASKS $BINDIR/procOBA_Sat \
           procOBA_Sat.$varname.config || exit 1

    $SCRIPTDIR/fit_semivariogram.py $CFGDIR/gage_$varname.cfg \
                                    gage_$varname.param || exit 1

    # Rescale the satellite error variance to compare with NWP.
    $SCRIPTDIR/rescale_sat_sigma2.py $varname || exit 1
done

# Update the lis.config error settings

$SCRIPTDIR/customize_lis_config.py $CFGDIR/autotune.cfg \
    imerg $enddt $dayrange || exit 1


exit 0
