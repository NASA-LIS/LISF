#!/bin/bash

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# -------------------------END NOTICE -- DO NOT EDIT-----------------------

# source and load shared functions from s2s_run.sh
source s2s_app/s2s_run.sh --source-only
# LIS-Hydro_S2S hindcast pre-processor

######################################################################
#                     PROCESS COMMAND LINE ARGUMENTS
######################################################################

REPORT=
STEP=

while getopts ":m:c:s:r:" opt; do
    case $opt in
    m) mon="$OPTARG"
       if [ $mon -lt 10 ]; then
           MM="0"$mon
       else
           MM=$mon
       fi
       ;;
    c) CFILE="$OPTARG"
       ;;
    r) REPORT="$OPTARG"
       ;;
    s) STEP="$OPTARG"
       ;;
    *) echo "     "
       echo "USAGE: s2s_app/s2s_hcst_preprocess.sh -m MONTH -c CONFIG_FILE -r REPORT -s STEP"
       echo "     "
       echo "where MANDATORY input parameters:"
       echo "---------------------------------"
       echo "  MONTH:       start month"
       echo "  CONFIG_FILE: config file (for hindcast or forecast)"
       echo "  Thus, s2s_app/s2s_hcst_preprocess.sh -m MONTH -c CONFIG_FILE is good to run the complete hindcast preprocess for MONTH."
       echo "     "       
       echo "with OPTIONAL flags:"
       echo "--------------------"
       echo "  REPORT:   Once the E2ES process has begun (jobs have been submitted), the REPORT flag can be used to check the progress of SLURM jobs (valid inputs: Y or N)"
       echo "  STEP:     The E2ES process includes three steps that are run sequentially DOWNLOAD (or DOWNLOADNMME), REORG, CLIM"
       echo "            However, the STEP option allows the user to kick start the process from the last completed step."
       echo "            -s STEP directs s2s_run.sh to start from a specific STEP (valid inputs: DOWNLOAD, DOWNLOADNMME, REORG, CLIM)."
       echo "     "                       
       exit 1 
       ;;    
    esac
done

# Exit if any mandatory argument is missing.

shift "$(( OPTIND - 1 ))"

if [ -z "$mon" ] || [ -z "$CFILE" ]; then
    echo "Missing mandatory arguments. "
    echo "Please enter s2s_app/s2s_hcst_preprocess.sh -h for usage instructions."
    exit 1
fi

#######################################################################
#    System Settings and Architecture Specific Environment Variables
#######################################################################

umask 022
ulimit -s unlimited
export ARCH=`uname`

export LISFDIR=`grep LISFDIR $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
export LISHDIR=${LISFDIR}/lis/utils/usaf/s2s/
export METFORC=`grep METFORC $CFILE | cut -d':' -f2 | tr -d "[:space:]"`    
export LISFMOD=`grep LISFMOD $CFILE | cut -d':' -f2 | tr -d "[:space:]"`    
export SPCODE=`grep SPCODE  $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
export E2ESDIR=`grep E2ESDIR $CFILE | cut -d':' -f2 | tr -d "[:space:]"`"/hindcast/"    
MODELS=`grep NMME_models $CFILE | cut -d'[' -f2 | cut -d']' -f1 | sed 's/,//g'`

export clim_syr=`grep clim_start_year  $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
export clim_eyr=`grep clim_end_year  $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
export NODE_NAME=`uname -n`
mmm=`date -d "1991-$MM-01" +%b | tr '[:upper:]' '[:lower:]'`

unset LD_LIBRARY_PATH
source /etc/profile.d/modules.sh
module use -a $LISFDIR/env/discover/
module --ignore-cache load $LISFMOD

BWD=`pwd`

#**********************************************************************
#                           FUNCTIONS
#**********************************************************************

download_cfsv2(){
    
    #######################################################################
    #                     Download CFSv2 forecasts
    #######################################################################

    # CFSv2 forecast
    cfsv2datadir=`grep cfsv2datadir $CFILE | cut -d':' -f2 | tr -d "[:space:]"`

    for ((YEAR=clim_syr; YEAR<=clim_eyr; YEAR++)); do
        sh s2s_app/wget_cfsv2_oper_ts_e2es.sh -y ${YEAR} -m ${MM} -c ${cfsv2datadir}
        ret_code=$?
        if [ $ret_code -gt 0 ]; then
            exit
        fi
    done
}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

download_nmme(){
    
    #######################################################################
    #                     Download NMME forecasts
    #######################################################################
            
    # NMME Precipitation
    Mmm=`date -d "2001-${MM}-01" +%b`
        
    # download NMME precip forecasts
    s2s_app/download_nmme_hindcasts.sh ${Mmm} ${CFILE}    
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

reorg_cfsv2(){
    
    #######################################################################
    # Reorganize CFSv2 Data
    #######################################################################

    clim_mid=$((clim_syr+clim_eyr))
    clim_mid=$((clim_mid / 2))
    
    cd ${SCRDIR}/reorg/
    CWD=`pwd`
    /bin/ln -s ${E2ESDIR}/bcsd_fcst/
    
    jobname=reorg_cfsv2
    python $LISHDIR/s2s_modules/bcsd_fcst/forecast_task_01.py -s $clim_syr -e $clim_mid -m $mmm -c $BWD/$CFILE -w ${CWD} -t 1 -H 10 -j ${jobname}_set1
    ((clim_mid++))
    python $LISHDIR/s2s_modules/bcsd_fcst/forecast_task_01.py -s $clim_mid -e $clim_eyr -m $mmm -c $BWD/$CFILE -w ${CWD} -t 1 -H 10 -j ${jobname}_set2
    
    job_list="$jobname*.j"
    for jfile in $job_list
    do
	thisID=$(submit_job "" "${jfile}")
	reorg_cfsv2_ID=`echo $reorg_cfsv2_ID`' '$thisID
    done
    reorg_cfsv2_ID=`echo $reorg_cfsv2_ID | sed "s| |:|g"`	
    
    cd ${BWD}
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

reorg_nmme(){
    
    #######################################################################
    # Reorganize NMME Data
    #######################################################################

    cd ${SCRDIR}/reorg/
    CWD=`pwd`
    /bin/ln -s ${E2ESDIR}/bcsd_fcst/
    
    jobname=reorg_nmme
    nmme_output_dir=${E2ESDIR}/bcsd_fcst/NMME/raw/Monthly/
    mkdir -p -m 775 $nmme_output_dir
        
    python $LISHDIR/s2s_app/s2s_api.py -c $BWD/$CFILE -f ${jobname}_run.j -t 1 -H 4 -j ${jobname}_ -w ${CWD}
    this_model=0
    for nmme_model in $MODELS; do	
	COMMAND="python $LISHDIR/s2s_modules/bcsd_fcst/bcsd_library/nmme_reorg_h.py $MM $nmme_output_dir $nmme_model $BWD/$CFILE"
	if [ $this_model == 0 ]; then
	    sed -i "s|COMMAND|${COMMAND}|g" ${jobname}_run.j
	else
	    LLINE=`grep -n nmme_reorg_h.py ${jobname}_run.j | tail -1 | cut -d':' -f1`
	    ((LLINE++))
	    sed -i "${LLINE}i ${COMMAND}" ${jobname}_run.j
	fi
	((this_model++))
    done;

    reorg_nmme_ID=$(submit_job "" "${jobname}_run.j")
    cd ${BWD}
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clim_nafpa(){
    
    #######################################################################
    # Create NAFPA Climatologies
    #######################################################################

    cd ${SCRDIR}/clim/
    CWD=`pwd`
    /bin/ln -s ${E2ESDIR}/bcsd_fcst/
    
    jobname=clim_nafpa
    outdir=${E2ESDIR}/bcsd_fcst/USAF-LIS7.3rc8_25km/raw/Climatology/
    mkdir -p -m 775 $outdir
     
    python $LISHDIR/s2s_app/s2s_api.py -c $BWD/$CFILE -f ${jobname}_run.j -t 1 -H 4 -j ${jobname}_ -w ${CWD}    
    this_var=0
    for var in LWdown_f_tavg Rainf_f_tavg Psurf_f_tavg  Qair_f_tavg SWdown_f_tavg Tair_f_tavg Wind_f_tavg; do
        COMMAND="python $LISHDIR/s2s_modules/bcsd_fcst/bcsd_library/calc_and_write_observational_climatology.py $var $BWD/$CFILE $outdir"
	if [ $this_var == 0 ]; then
	    sed -i "s|COMMAND|${COMMAND}|g" ${jobname}_run.j
	else
	    LLINE=`grep -n calc_and_write_observational_climatology.py ${jobname}_run.j | tail -1 | cut -d':' -f1`
	    ((LLINE++))
	    sed -i "${LLINE}i ${COMMAND}" ${jobname}_run.j
	fi	
	((this_var++))
    done;

    clim_nafpa_ID=$(submit_job "" "${jobname}_run.j")
    cd ${BWD}

}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clim_cfsv2(){
    
    #######################################################################
    # Create CFSv2 Climatologies
    #######################################################################

    cd ${SCRDIR}/clim/
    CWD=`pwd`
    /bin/ln -s ${E2ESDIR}/bcsd_fcst/
    
    jobname=clim_cfsv2
    fcst_indir=${E2ESDIR}/bcsd_fcst/CFSv2_25km/raw/Monthly/
    outdir=${E2ESDIR}/bcsd_fcst/CFSv2_25km/raw/Climatology/${mmm}01/
    mkdir -p -m 775 ${outdir}
 
    python $LISHDIR/s2s_app/s2s_api.py -c $BWD/$CFILE -f ${jobname}_run.j -t 1 -H 5 -j ${jobname}_ -w ${CWD}
    this_var=0
    for var in PRECTOT  LWS  PS  Q2M  SLRSF  T2M  WIND10M; do
        COMMAND="python $LISHDIR/s2s_modules/bcsd_fcst/bcsd_library/calc_and_write_forecast_climatology.py $var $MM $BWD/$CFILE $fcst_indir $outdir"
	if [ $this_var == 0 ]; then
	    sed -i "s|COMMAND|${COMMAND}|g" ${jobname}_run.j
	else
	    LLINE=`grep -n calc_and_write_forecast_climatology.py ${jobname}_run.j | tail -1 | cut -d':' -f1`
	    ((LLINE++))
	    sed -i "${LLINE}i ${COMMAND}" ${jobname}_run.j
	fi	
	((this_var++))
    done;
    
    clim_cfsv2_ID=$(submit_job $reorg_cfsv2_ID "${jobname}_run.j")
    cd ${BWD}
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clim_nmme(){
    
    #######################################################################
    # Create NMME Climatologies
    #######################################################################

    cd ${SCRDIR}/clim/
    CWD=`pwd`
    /bin/ln -s ${E2ESDIR}/bcsd_fcst/
    
    jobname=clim_nmme
    outdir=${E2ESDIR}/bcsd_fcst/NMME/raw/Climatology/${mmm}01/
    mkdir -p -m 775 ${outdir}
    cd ${outdir}
    mkdir -p -m 775 $MODELS
    cd ${SCRDIR}/clim/
 
    python $LISHDIR/s2s_app/s2s_api.py -c $BWD/$CFILE -f ${jobname}_run.j -t 1 -H 4 -j ${jobname}_ -w ${CWD}
    this_model=0
    for model in $MODELS; do

        nmme_indir=${E2ESDIR}/bcsd_fcst/NMME/raw/Monthly/${mmm}01/${model}/
        outdir=${E2ESDIR}/bcsd_fcst/NMME/raw/Climatology/${mmm}01/${model}/
	
	COMMAND="python $LISHDIR/s2s_modules/bcsd_fcst/bcsd_library/calc_and_write_nmme_forecast_climatology.py 'PRECTOT'  $MM ${model} $BWD/$CFILE $nmme_indir $outdir"
		
	if [ $this_model == 0 ]; then
	    sed -i "s|COMMAND|${COMMAND}|g" ${jobname}_run.j
	else
	    LLINE=`grep -n calc_and_write_nmme_forecast_climatology.py ${jobname}_run.j | tail -1 | cut -d':' -f1`
	    ((LLINE++))
	    sed -i "${LLINE}i ${COMMAND}" ${jobname}_run.j
	fi
	((this_model++))
    done;

    clim_nmme_ID=$(submit_job $reorg_nmme_ID "${jobname}_run.j")
    
    cd ${BWD}
}

#######################################################################
#**********************************************************************
#                           MAIN SCRIPT
#**********************************************************************

SCRDIR=${E2ESDIR}/scratch/${MM}/
if [ "$REPORT" = 'Y' ] || [ "$REPORT" = 'y' ]; then
    print_walltimes
    exit
fi

mkdir -p -m 775 $SCRDIR/reorg
mkdir -p -m 775 $SCRDIR/clim
mkdir -p -m 775 ${E2ESDIR}/bcsd_fcst

JOB_SCHEDULE=${SCRDIR}/SLURM_JOB_SCHEDULE
/bin/rm -f $JOB_SCHEDULE

echo "#######################################################################" >> $JOB_SCHEDULE
echo "                         SLURM JOB SCHEDULE                            " >> $JOB_SCHEDULE
echo "#######################################################################" >> $JOB_SCHEDULE
echo "                         " >> $JOB_SCHEDULE
python $LISHDIR/s2s_app/s2s_api.py -s $JOB_SCHEDULE -m "JOB ID" -f "JOB SCRIPT" -a "AFTER" -c $BWD/$CFILE

#######################################################################
#                               Submit jobs
#######################################################################

reorg_cfsv2_ID=
reorg_nmme_ID=
clim_nafpa_ID=

case $STEP in
    DOWNLOAD)
        download_cfsv2
	download_nmme
        ;;
    DOWNLOADNMME)
	download_nmme
        ;;
    REORG)
	echo "              " >> $JOB_SCHEDULE
	echo "(1) Reorganize forcings" >> $JOB_SCHEDULE
	echo "-----------------------" >> $JOB_SCHEDULE
	echo "              " >> $JOB_SCHEDULE
        reorg_cfsv2
	reorg_nmme
        ;;                
    CLIM)
	echo "              " >> $JOB_SCHEDULE
	echo "(2) Write climgatological files" >> $JOB_SCHEDULE
	echo "------------------------------" >> $JOB_SCHEDULE
	echo "              " >> $JOB_SCHEDULE	
        clim_nafpa
	clim_cfsv2
	clim_nmme
        ;;            
    *)
	echo "              " >> $JOB_SCHEDULE
	echo "(1) Reorganize forcings" >> $JOB_SCHEDULE
	echo "-----------------------" >> $JOB_SCHEDULE
	echo "              " >> $JOB_SCHEDULE	
        reorg_cfsv2
        reorg_nmme
	echo "              " >> $JOB_SCHEDULE
	echo "(2) Write climatological files" >> $JOB_SCHEDULE
	echo "------------------------------" >> $JOB_SCHEDULE
	echo "              " >> $JOB_SCHEDULE		
        clim_nafpa
        clim_cfsv2
        clim_nmme
        ;;   
esac                                            
