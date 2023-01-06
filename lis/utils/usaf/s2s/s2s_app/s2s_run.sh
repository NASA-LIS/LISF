#!/bin/bash

# LIS-Hydro_S2S subsystem run script 

######################################################################
#                     PROCESS COMMAND LINE ARGUMENTS
######################################################################

ONE="N"
REPORT=
DELETE=
STEP=

while getopts ":y:m:c:d:r:s:o:" opt; do
    case $opt in
	y) YYYY="$OPTARG"
	   ;;
	m) mon="$OPTARG"
	   if [ $mon -lt 10 ]; then
	       MM="0"$mon
	   else
	       MM=$mon
	   fi
	   ;;
	c) CFILE="$OPTARG"
	   ;;
	d) DELETE="$OPTARG"
	   ;;
	r) REPORT="$OPTARG"
	   ;;
	s) STEP="$OPTARG"
	   ;;
	o) ONE="$OPTARG"
	   ;;	
	*) echo "     "
	   echo "USAGE: s2s_app/s2s_run.sh -y YEAR -m MONTH -c CONFIG_FILE -d DELETE -r REPORT -s STEP -o ONE_STEP"
	   echo "     "
	   echo "where MANDATORY input parameters:"
	   echo "---------------------------------"
	   echo "  YEAR:        start year"
	   echo "  MONTH:       start month [1 to 12]"
	   echo "  CONFIG_FILE: config file (for hindcast or forecast)"
	   echo "  Thus, s2s_app/s2s_run.sh -y YEAR -m MONTH -c CONFIG_FILE is good to run the complete E2ES process for YEAR/MONTH."
	   echo "     "
	   echo "with OPTIONAL flags:"
	   echo "--------------------"
	   echo "  DELETE:   delete YEAR/MONTH (valid inputs: Y or N)"
	   echo "  REPORT:   Once the E2ES process has begun (jobs have been submitted), the REPORT flag is good to check the progress of SLURM jobs (valid inputs: Y or N)"
	   echo "  STEP:     The E2ES process includes seven steps that are run sequentially LISDA, LDTICS, BCSD, FCST, POST, METRICS, PLOTS."
	   echo "            However, the STEP option allows user to kick start the process from the last completed step."
	   echo "            -s STEP is good to ask s2s_run.sh to start from a specific STEP (valid inputs: LISDA, LDTICS, BCSD, FCST, POST, METRICS or PLOTS)."
           echo "  ONE_STEP: flag is good to run only the above -s STEP (valid inputs: Y or N). If ONE_STEP is set to Y, the process will exit upon completion of above STEP"  	      	   	      	   
	   exit 1 
	   ;;	 
    esac
done

# Exit if any mandatory argument is missing.

shift "$(( OPTIND - 1 ))"

if [ -z "$YYYY" ] || [ -z "$mon" ] || [ -z "$CFILE" ]; then
    echo "Missing mandatory arguments. "
    echo "Please enter s2s_app/s2s_run.sh -h for usage instructions."
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
export AF10KM=`grep AF10KM  $CFILE | cut -d':' -f2 | tr -d "[:space:]"`    
export LISFMOD=`grep LISFMOD $CFILE | cut -d':' -f2 | tr -d "[:space:]"`    
export SPCODE=`grep SPCODE  $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
export S2STOOL=$LISHDIR/
export DATATYPE=`grep DATATYPE  $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
export E2ESROOT=`grep E2ESDIR $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
export DOMAIN=`grep DOMAIN $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
if [ $DATATYPE == "hindcast" ]; then
    export E2ESDIR=`grep E2ESDIR $CFILE | cut -d':' -f2 | tr -d "[:space:]"`"/hindcast/"    
else
    export E2ESDIR=`grep E2ESDIR $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
fi
export SUPDIR=`grep supplementarydir $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
export LDTFILE=`grep ldtinputfile $CFILE | cut -d':' -f2 | tr -d "[:space:]"`

unset LD_LIBRARY_PATH
source /etc/profile.d/modules.sh
module use -a $LISFDIR/env/discover/
module --ignore-cache load $LISFMOD

BWD=`pwd`

#**********************************************************************
#                           FUNCTIONS
#**********************************************************************

submit_job(){
    if [[ $1 == "" ]] || [[ $1 == "," ]]; then
	submit_ID="`sbatch $2 |  cut -d' ' -f4`"
	python $LISHDIR/s2s_app/write_to_file.py -s $JOB_SCHEDULE -m $submit_ID -f $2
    else
	submit_ID="`sbatch --dependency=afterok:$1 $2 |  cut -d' ' -f4`"
	python $LISHDIR/s2s_app/write_to_file.py -s $JOB_SCHEDULE -m $submit_ID -f $2 -a `echo $1`
    fi
    echo $submit_ID
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

set_permission(){

    cd ${E2ESDIR}
    /bin/rm -f set_permission.j
    cat << EOF > ${E2ESDIR}/set_permission.j
#!/bin/bash

#######################################################################
#                        Set Read/Write permission 
#######################################################################

#SBATCH --account=${SPCODE}
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
#SBATCH --job-name=set_permission_
#SBATCH --output ${SCRDIR}/set_permission_%j.out
#SBATCH --error ${SCRDIR}/set_permission_%j.err

cd ${E2ESDIR}

find . -type d -exec chmod 0775 {} \;
find . -name "*.nc" -exec chmod 0644 {} \;
find . -name "*.NC" -exec chmod 0644 {} \;
find . -name "*.NC4" -exec chmod 0644 {} \;
find . -name "*.nc4" -exec chmod 0644 {} \;
find . -name "*.TIF" -exec chmod 0644 {} \;
find . -name "*.png" -exec chmod 0644 {} \;

EOF
   perm_ID=$(submit_job $1 "set_permission.j") 
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

delete_forecast(){

    delete_files(){
	echo 'deleting ....  ' $1
	chmod 755 -R $1
	/bin/rm -rf $1
	
    }

    cd ${E2ESDIR}/
    YYYY=$1
    MM=$2
    mon=`echo $MM |bc`
    YYYYMMP=`date -d "$YYYY-$MM-01 -1 month" +%Y%m`
    YYYYP=`echo $YYYYMMP | cut -c1-4`
    MMP=`echo $YYYYMMP | cut -c5-6`    
    mon_names=('jan' 'feb' 'mar' 'apr' 'may' 'jun' 'jul' 'aug' 'sep' 'oct' 'nov' 'dec')
    Mon=`echo ${mon_names[$mon-1]}`
    
    # delete scratch
    /bin/rm -rf ${E2ESDIR}/scratch/${YYYY}${MM}/ldt.config_noahmp401_nmme_\*_${YYYY}${MM}
    
    # delete LISDA
    delete_files ${E2ESDIR}/lis_darun/output/ROUTING/${YYYY}${MM}/
    delete_files ${E2ESDIR}/lis_darun/output/SURFACEMODEL/${YYYY}${MM}/
    delete_files ${E2ESDIR}/lis_darun/output/ROUTING/${YYYYP}${MMP}/LIS_HIST_\*.nc
    delete_files ${E2ESDIR}/lis_darun/output/SURFACEMODEL/${YYYYP}${MMP}/LIS_HIST_\*.nc
    delete_files ${E2ESDIR}/lis_darun/output/lis.config_files/lis.config_darun_${YYYY}${MM}
    
    # delete LDTICS
    delete_files ${E2ESDIR}/ldt_ics/ldt.config_files/ldt.config_noahmp401_nmme_\*_${YYYY}${MM}
    delete_files ${E2ESDIR}/ldt_ics/\*/\*${Mon^}${YYYY}*\
    
    # delete BCSD
    delete_files ${E2ESDIR}/bcsd_fcst/CFSv2_25km/bcsd/6-Hourly/${Mon}01/${YYYY}
    delete_files ${E2ESDIR}/bcsd_fcst/CFSv2_25km/final/6-Hourly/${Mon}01/${YYYY}
    delete_files ${E2ESDIR}/bcsd_fcst/CFSv2_25km/raw/6-Hourly/${Mon}01/${YYYY}
    delete_files ${E2ESDIR}/bcsd_fcst/CFSv2_25km/bcsd/Monthly/${Mon}01/\*_${YYYY}_${YYYY}.nc
    delete_files ${E2ESDIR}/bcsd_fcst/CFSv2_25km/raw/Monthly/${Mon}01/${YYYY}

    delete_files ${E2ESDIR}/bcsd_fcst/NMME/bcsd/6-Hourly/${Mon}01/\*/${YYYY}
    delete_files ${E2ESDIR}/bcsd_fcst/NMME/final/6-Hourly/\*/${Mon}01/${YYYY}
    delete_files ${E2ESDIR}/bcsd_fcst/NMME/bcsd/Monthly/${Mon}01/\*_${YYYY}_${YYYY}.nc
    delete_files ${E2ESDIR}/bcsd_fcst/NMME/raw/Monthly/${Mon}01/\*/${YYYY}
    delete_files ${E2ESDIR}/bcsd_fcst/NMME/linked_cfsv2_precip_files/${Mon}01/${YYYY}
        
    # delete FCST
    delete_files ${E2ESDIR}/lis_fcst/${YYYY}${MM}/
    delete_files ${E2ESDIR}/lis_fcst/input/\*/\*/lis.config.s2sglobal.noahmp401.hymap2.da_ics_forecast_\*_${YYYY}${MM}
    
    # delete POST
    delete_files ${E2ESDIR}/s2spost/${YYYY}${MM}/
    # delete METRICS
    delete_files ${E2ESDIR}/s2smetric/${YYYY}${MM}/
    # delete PLOTS
    delete_files ${E2ESDIR}/s2splots/${YYYY}${MM}/
    exit    
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

download_forecasts(){
    
    #######################################################################
    #                     Download CFSv2 and NMME forecasts
    #######################################################################

    # CFSv2 forecast
    cfsv2datadir=`grep fcst_download_dir $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
    sh s2s_app/wget_cfsv2_oper_ts_e2es.sh -y ${YYYY} -m ${MM} -c ${cfsv2datadir}
    ret_code=$?
    if [ $ret_code -gt 0 ]; then
	exit
    fi
    
    # NMME Precipitation
    Mmm=`date -d "${YYYY}-${MM}-01" +%b`
    NMME_RAWDIR=`grep nmme_download_dir $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
    prec_files=`ls $NMME_RAWDIR/*/*${Mmm}.${YYYY}.nc`
    all_models=`grep -n NMME_ALL: $CFILE | cut -d':' -f3`
    LINE2=`grep -n NMME_models: $CFILE | cut -d':' -f1`
    new_line="  NMME_models: "`echo ${all_models}`
    
    sed -i "${LINE2}s/.*/${new_line}/1" $CFILE
    
    declare -A nmme_models=( ["CanSIPS-IC3"]="GNEMO5, CCM4" ["COLA-RSMAS-CCSM4"]="CCSM4" ["GFDL-SPEAR"]="GFDL" ["NASA-GEOSS2S"]="GEOSv2" ["NCEP-CFSv2"]="CFSv2" )
    have_model="["
    
    if [ `echo $prec_files | wc -w` -lt 5 ]; then
	
	# download NMME precip forecasts
	s2s_app/download_nmme_precip.sh ${YYYY} ${Mmm} ${CFILE}
	prec_files=`ls $NMME_RAWDIR/*/*${Mmm}.${YYYY}.nc`
	
	for dir in CanSIPS-IC3 COLA-RSMAS-CCSM4 GFDL-SPEAR NASA-GEOSS2S NCEP-CFSv2
	do
	    fdown=$NMME_RAWDIR/$dir/prec.$dir.mon_${Mmm}.${YYYY}.nc
	    if [ `file $fdown | rev | cut -d' ' -f1 | rev` == "data" ]; then
		have_model=${have_model}${nmme_models[$dir]}', '
	    else
		/bin/rm ${fdown}
	    fi
	done
	have_model=${have_model}"]"
	if [ `echo ${have_model} | wc -w` -lt 7 ]; then
	    NAVAIL=`echo ${have_model} | wc -w`
	    ((NAVAIL--))
	    echo 
	    read -p "Prpecipitation forecasts are available for only ${NAVAIL} NMME models (${have_model}). Do you want to continue (Y/N)?" YESORNO
	    
	    if [ "$YESORNO" = 'Y' ] || [ "$YESORNO" = 'y' ]; then
		LINE2=`grep -n NMME_models: $CFILE | cut -d':' -f1`
		new_line="  NMME_models: "`echo ${have_model}`
		sed -i "${LINE2}s/.*/${new_line}/1" $CFILE
	    else
		exit
	    fi    
	fi
    fi
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

lis_darun(){

    #######################################################################
    # (1) LIS DA run starts on the 1st of the previous month to generate
    #     initial conditions
    #######################################################################
    
    echo "                         " >> $JOB_SCHEDULE
    echo "(1) LIS Data Assimilation" >> $JOB_SCHEDULE
    echo "-------------------------" >> $JOB_SCHEDULE
    echo "                         " >> $JOB_SCHEDULE

    # set up input directory
    mkdir -p -m 775 ${E2ESDIR}/lis_darun/input/
    cd ${E2ESDIR}/lis_darun/input/
    /bin/ln -s ${LISHDIR}/s2s_modules/lis_darun/forcing_variables.txt
    /bin/ln -s ${LISHDIR}/s2s_modules/lis_darun/noahmp401_parms
    /bin/ln -s ${LISHDIR}/s2s_modules/lis_darun/template_files
    /bin/ln -s ${LISHDIR}/s2s_modules/lis_darun/attribs
    /bin/ln -s ${LISHDIR}/s2s_modules/lis_darun/tables
    /bin/ln -s ${SUPDIR}/lis_darun/cdf/${DOMAIN} cdf
    /bin/ln -s ${SUPDIR}/lis_darun/RS_DATA
    /bin/ln -s ${SUPDIR}/lis_darun/${LDTFILE}
    cd ${BWD}
    
    # previous month
    YYYYMMP=`date -d "$YYYY-$MM-01 -1 month" +%Y%m`
    YYYYP=`echo $YYYYMMP | cut -c1-4`
    MMP=`echo $YYYYMMP | cut -c5-6`
    monP=`echo $MMP |bc`
    
    PERTMODE=`grep pertmode $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
    cd ${SCRDIR}/lis_darun
    CWD=`pwd`
    /bin/ln -s ${LISFDIR}/lis/LIS
    /bin/ln -s ${E2ESDIR}/lis_darun/input
    /bin/ln -s ${E2ESDIR}/lis_darun/output
    mkdir -p -m 775 output/lis.config_files/
    mkdir -p -m 775 ${CWD}/logs_${YYYYP}${MMP}
    
    # configure batch script
    # ----------------------
    
    python $LISHDIR/s2s_app/write_to_file.py -c ${BWD}/${CFILE} -f lisda_run.j -H 2 -j lisda_ -w ${CWD} -L Y
    COMMAND='mpirun -np $SLURM_NTASKS ./LIS'
    sed -i "s|COMMAND|${COMMAND}|g" lisda_run.j
    
    # configure lis.config
    # --------------------
    /bin/cp ${E2ESDIR}/lis_darun/input/template_files/lis.config_template.${DOMAIN} lis.config
    DAPERTRSTFILE=./output/DAPERT/${YYYYP}${MMP}/LIS_DAPERT_${YYYYP}${MMP}010000.d01.bin
    NOAHMP401RSTFILE=./output/SURFACEMODEL/${YYYYP}${MMP}/LIS_RST_NOAHMP401_${YYYYP}${MMP}010000.d01.nc
    HYMAP2RSTFILE=./output/ROUTING/${YYYYP}${MMP}/LIS_RST_HYMAP2_router_${YYYYP}${MMP}010000.d01.nc
    LSMLISLOGFILE=${CWD}/logs_${YYYYP}${MMP}'/lislog'
    
    #sed -i "s|DAPERTRSTFILE|${DAPERTRSTFILE}|g" lis.config
    sed -i "s|NOAHMP401RSTFILE|${NOAHMP401RSTFILE}|g" lis.config
    sed -i "s|HYMAP2RSTFILE|${HYMAP2RSTFILE}|g" lis.config
    sed -i "s|STARTYR|${YYYYP}|g" lis.config
    sed -i "s|STARTMO|${monP}|g" lis.config
    sed -i "s|STARTDA|1|g" lis.config
    sed -i "s|FINALYR|${YYYY}|g" lis.config
    sed -i "s|FINALMO|${mon}|g" lis.config
    sed -i "s|FINALDA|1|g" lis.config
    sed -i "s|PERTMODE|${PERTMODE}|g" lis.config
    sed -i "s|LSMLISLOGFILE|${LSMLISLOGFILE}|g" lis.config
    
    /bin/cp lis.config output/lis.config_files/lis.config_darun_${YYYYP}${MMP}
    
    # submit job
    # ----------
    lisda_ID=$(submit_job "" "lisda_run.j")
    cd ${BWD}
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ldt_ics(){
    
    #######################################################################
    # (2) LDT ICS run to LIS input files
    #######################################################################

    echo "              " >> $JOB_SCHEDULE
    echo "(2) LDT and Initial Conditions" >> $JOB_SCHEDULE
    echo "------------------------------" >> $JOB_SCHEDULE
    echo "              " >> $JOB_SCHEDULE

    mkdir -p ${E2ESDIR}/ldt_ics/input
    cd ${E2ESDIR}/ldt_ics
    mkdir -p -m 775 $MODELS
    mkdir -p -m 775 ldt.config_files
    mkdir -p -m 775 template_files    
    /bin/cp -p ${LISHDIR}/s2s_modules/ldt_ics/template_files/ldt.config_noahmp401_nmme_TEMPLATE.${DOMAIN} template_files/ldt.config_noahmp401_nmme_TEMPLATE
    cd ${E2ESDIR}/ldt_ics/input
    /bin/ln -s ${SUPDIR}/lis_darun/${LDTFILE}
    /bin/ln -s ${SUPDIR}/LS_PARAMETERS
    
    cd ${SCRDIR}/ldt_ics
    CWD=`pwd`
    /bin/ln -s ${LISFDIR}/ldt/LDT
    /bin/ln -s ${E2ESDIR}/lis_darun/output lisda_output
    /bin/ln -s ${E2ESDIR}/ldt_ics/input
    
    for model in $MODELS
    do
	/bin/ln -s ${E2ESDIR}/ldt_ics/$model
    done
    /bin/ln -s ${E2ESDIR}/ldt_ics/ldt.config_files
    /bin/ln -s ${E2ESDIR}/ldt_ics/template_files
    
    # configure batch script
    # ----------------------
    
    python $LISHDIR/s2s_app/write_to_file.py -c $BWD/$CFILE -f ldtics_run.j -t 1 -H 2 -j ldtics_ -w ${CWD}
    COMMAND="python ${LISHDIR}/s2s_modules/ldt_ics/generate_ldtconfig_files_ensrst_nrt.py -y ${YYYY} -m ${mon} -i ./lisda_output -w $CWD -s $BWD/$CFILE"
    sed -i "s|COMMAND|${COMMAND}|g" ldtics_run.j
    
    # submit job
    # ----------
    ldtics_ID=$(submit_job "$lisda_ID" "ldtics_run.j")	    
    cd ${BWD}
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bcsd_fcst(){
    
    #######################################################################
    # (3) BCSD 
    #######################################################################
    
    echo "              " >> $JOB_SCHEDULE
    echo "(3) Bias Correction and Spatial Downscaling" >> $JOB_SCHEDULE
    echo "-------------------------------------------" >> $JOB_SCHEDULE
    echo "                                           " >> $JOB_SCHEDULE
    
    obs_clim_dir=${E2ESROOT}/hindcast/bcsd_fcst/CFSv2_25km/raw/Climatology/
    nmme_clim_dir=${E2ESROOT}/hindcast/bcsd_fcst/NMME/raw/Climatology/
    usaf_25km=${E2ESROOT}/hindcast/bcsd_fcst/USAF-LIS7.3rc8_25km/raw/Climatology/
    
    mkdir -p -m 775 ${E2ESDIR}/bcsd_fcst
    cd ${E2ESDIR}/bcsd_fcst
    
    mkdir -p -m 775 USAF-LIS7.3rc8_25km/raw
    mkdir -p -m 775 CFSv2_25km/raw
    mkdir -p -m 775 NMME/raw
    
    # link Climatology directories
    cd ${E2ESDIR}/bcsd_fcst/CFSv2_25km/raw
    /bin/ln -s $obs_clim_dir Climatology
    
    cd ${E2ESDIR}/bcsd_fcst/NMME/raw
    /bin/ln -s $nmme_clim_dir Climatology
    
    cd ${E2ESDIR}/bcsd_fcst/USAF-LIS7.3rc8_25km/raw
    /bin/ln -s $usaf_25km Climatology
    
    # manage jobs from SCRATCH
    cd ${SCRDIR}/bcsd_fcst
    /bin/ln -s ${E2ESDIR}/bcsd_fcst/
    CWD=`pwd`
    
    mmm=`date -d "$YYYY-$MM-01" +%b | tr '[:upper:]' '[:lower:]'`
    
    # Task 1: Generate and rescale 6-hourly files to 25 KM (forecast_task_01.py)
    # --------------------------------------------------------------------------
    jobname=bcsd01
    python $LISHDIR/s2s_modules/bcsd_fcst/forecast_task_01.py -s $YYYY -e $YYYY -m $mmm -c $BWD/$CFILE -w ${CWD} -t 1 -H 10 -j $jobname
    bcsd01_ID=$(submit_job "" "${jobname}_run.j")
    
    # Task 3: Rescale and reorganize NMME Data (forecast_task_03.py)
    # --------------------------------------------------------------
    jobname=bcsd03
    python $LISHDIR/s2s_modules/bcsd_fcst/forecast_task_03.py -s $YYYY -m $MM -c $BWD/$CFILE -w ${CWD} -t 1 -H 2 -j $jobname

    job_list="$jobname*.j"
    bcsd03_ID=
    for jfile in $job_list
    do
	thisID=$(submit_job "" "${jfile}")
	bcsd03_ID=`echo $bcsd03_ID`' '$thisID
    done
    bcsd03_ID=`echo $bcsd03_ID | sed "s| |,|g"`
    
    # Task 4: Monthly "BC" step applied to CFSv2 (forecast_task_04.py, after 1 and 3)
    # -------------------------------------------------------------------------------
    jobname=bcsd04
    python $LISHDIR/s2s_modules/bcsd_fcst/forecast_task_04.py -s $YYYY -e $YYYY -m $mmm -n $MM -c $BWD/$CFILE -w ${CWD} -t 1 -H 4 -j $jobname
    
    unset job_list
    job_list="$jobname*.j"
    bcsd04_ID=
    for jfile in $job_list
    do
	thisID=$(submit_job "$bcsd01_ID,$bcsd03_ID" "${jfile}")
	bcsd04_ID=`echo $bcsd04_ID`' '$thisID
    done
    bcsd04_ID=`echo $bcsd04_ID | sed "s| |,|g"`
    
    # Task 5: Monthly "BC" step applied to NMME (forecast_task_05.py: after 1 and 3)
    # ------------------------------------------------------------------------------
    jobname=bcsd05
    for model in $MODELS
    do
	python $LISHDIR/s2s_modules/bcsd_fcst/forecast_task_05.py -s $YYYY -e $YYYY -m $mmm -n $MM -c $BWD/$CFILE -w ${CWD} -t 1 -H 4 -M $model -j $jobname    
    done
    
    unset job_list
    job_list="$jobname*.j"
    bcsd05_ID=
    for jfile in $job_list
    do
	thisID=$(submit_job "$bcsd01_ID,$bcsd03_ID" "${jfile}")
	bcsd05_ID=`echo $bcsd05_ID`' '$thisID
    done
    bcsd05_ID=`echo $bcsd05_ID | sed "s| |,|g"`
    
    # Task 6: CFSv2 Temporal Disaggregation (forecast_task_06.py: after 4 and 5)
    # --------------------------------------------------------------------------
    jobname=bcsd06
    python $LISHDIR/s2s_modules/bcsd_fcst/forecast_task_06.py -s $YYYY -e $YYYY -m $mmm -n $MM -c $BWD/$CFILE -w ${CWD} -t 1 -H 4 -j $jobname
    
    unset job_list
    job_list=`ls $jobname*.j`
    bcsd06_ID=
    for jfile in $job_list
    do
	thisID=$(submit_job "$bcsd04_ID,$bcsd05_ID" "${jfile}")
	bcsd06_ID=`echo $bcsd06_ID`' '$thisID
    done
    bcsd06_ID=`echo $bcsd06_ID | sed "s| |,|g"` 
    
    # Task 7: Generate symbolic links to sub-daily CFSv2 BC forecasts for NMME
    # temporal disaggregation due to an uneven number of ensembles between the datasets
    #  (forecast_task_07.py: after 4 and 5)
    # --------------------------------------------------------------------------
    jobname=bcsd07
    # python $LISHDIR/s2s_modules/bcsd_fcst/forecast_task_07.py -s $YYYY -m $mmm -c -w ${CWD}
    # NOTE: This is run inside a Task 6 job
    
    # Task 8: NMME Temporal Disaggregation (forecast_task_08.py: after 6, 7)
    # ----------------------------------------------------------------------------
    jobname=bcsd08
    for model in $MODELS
    do
	python $LISHDIR/s2s_modules/bcsd_fcst/forecast_task_08.py -s $YYYY -e $YYYY -m $mmm -n $MM -c $BWD/$CFILE -w ${CWD} -t 1 -H 2 -M $model -j $jobname    
    done
    
    unset job_list
    job_list=`ls $jobname*.j`
    bcsd08_ID=
    for jfile in $job_list
    do
	thisID=$(submit_job "$bcsd06_ID" "${jfile}")
	bcsd08_ID=`echo $bcsd08_ID`' '$thisID
    done
    bcsd08_ID=`echo $bcsd08_ID | sed "s| |,|g"`
    
    # Task 9: Combine the CFSv2 forcing fields into final format for LIS to read
    #         (forecast_task_09.py: after 8)
    # ---------------------------------------------------------------------------
    jobname=bcsd09
    python $LISHDIR/s2s_modules/bcsd_fcst/forecast_task_09.py -s $YYYY -e $YYYY -m $mmm -n $MM -M CFSv2 -c $BWD/$CFILE -w ${CWD} -j $jobname -t 1 -H 4

    bcsd09_ID=$(submit_job "$bcsd08_ID" "${jobname}_run.j")
    
    # Task 10: Combine the NMME forcing fields into final format for LIS to read
    #          and symbolically link to the reusable CFSv2 met forcings
    #         (forecast_task_10.py: after 8)
    # ---------------------------------------------------------------------------
    jobname=bcsd10
    # NOTE : Task 10  Job scripts are written by forecast_task_09.py to execute: 
    # python $LISHDIR/s2s_modules/bcsd_fcst/forecast_task_10.py -s $YYYY -m $mmm -n $MM -w ${CWD} -M NMME_MODEL 
    
    unset job_list
    job_list=`ls $jobname*.j`
    bcsd10_ID=
    for jfile in $job_list
    do
	thisID=$(submit_job "$bcsd08_ID" "${jfile}")
	bcsd10_ID=`echo $bcsd10_ID`' '$thisID
    done
    bcsd10_ID=`echo $bcsd10_ID | sed "s| |,|g"`
    
    # Task 11: Copy 9th forecast lead file as 10th forecast lead for LIS runs
    #         (forecast_task_11.py: after 9 and 10)
    # ---------------------------------------------------------------------------
    jobname=bcsd11
    # NOTE : Task 11  Job scripts are written by forecast_task_09.py to execute: 
    # python $LISHDIR/s2s_modules/bcsd_fcst/forecast_task_11.py -s $YYYY -m $mmm -n $MM -c $BWD/$CFILE -w ${CWD}
    bcsd11_ID=$(submit_job "$bcsd09_ID,$bcsd10_ID" "${jobname}_run.j")
    
    # Task 12:  Task to introduce an all-zero variable V10M due to the way wind
    #           is handled in the USAF forcing
    #         (forecast_task_12.py: after 9 and 10)
    # ---------------------------------------------------------------------------
    jobname=bcsd12
    # NOTE : Task 12  Job scripts are written by forecast_task_09.py to execute: 
    # python $LISHDIR/s2s_modules/bcsd_fcst/forecast_task_12.py -s $YYYY -m $mmm -n $MM -c $BWD/$CFILE -w ${CWD}
    bcsd12_ID=$(submit_job "$bcsd09_ID,$bcsd10_ID" "${jobname}_run.j")
    
    cd ${BWD}
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

lis_fcst(){
    #######################################################################
    # (4) LIS FORECAST
    #######################################################################
    
    jobname=lis_fcst
    
    echo "              " >> $JOB_SCHEDULE
    echo "(4) LIS Forecast                           " >> $JOB_SCHEDULE
    echo "-------------------------------------------" >> $JOB_SCHEDULE
    echo "                                           " >> $JOB_SCHEDULE
    
    Mmm1=`date -d "$YYYY-$MM-01" +%b`1

    mkdir -p  -m 775 ${E2ESDIR}/lis_fcst
    cd ${E2ESDIR}/lis_fcst
    mkdir -p -m 775 input/LDT_ICs/
    cd ${E2ESDIR}/lis_fcst/input/
    /bin/ln -s ${LISHDIR}/s2s_modules/lis_darun/forcing_variables.txt
    /bin/ln -s ${LISHDIR}/s2s_modules/lis_darun/noahmp401_parms
    /bin/ln -s ${LISHDIR}/s2s_modules/lis_fcst/template_files
    /bin/ln -s ${LISHDIR}/s2s_modules/lis_fcst/tables
    /bin/ln -s ${SUPDIR}/lis_darun/${LDTFILE}
    
    cd ${E2ESDIR}/lis_fcst/input/LDT_ICs/    
    for model in $MODELS
    do
	mkdir -m 775 -p ${E2ESDIR}/lis_fcst/${YYYY}${MM}/${model}/logs/
	/bin/ln -s ${E2ESDIR}/ldt_ics/$model
    done
    
    cd ${SCRDIR}/lis_fcst
    CWD=`pwd`
    /bin/ln -s ${LISFDIR}/lis/LIS
    /bin/ln -s ${E2ESDIR}/lis_fcst/input
    /bin/ln -s ${E2ESDIR}/lis_fcst/${YYYY}${MM}
    /bin/ln -s ${E2ESDIR}/bcsd_fcst
    
    # write SLURM job scripts
    python $LISHDIR/s2s_modules/lis_fcst/generate_lis_config_scriptfiles_fcst.py -c $BWD/$CFILE -y $YYYY -m $MM -w $CWD -j $jobname
    
    lisfcst_ID=
    
    for model in $MODELS
    do
	unset job_list
	unset nFiles
	unset FileNo
	job_list=`ls ${jobname}_${model}*_run.j`
	nFiles=`echo $job_list | wc -w`
	FileNo=1
	for jfile in $job_list
	do
	    if [ $nFiles -gt 1 ]; then
		if [ $FileNo  -eq 1 ]; then
		    thisID=$(submit_job "$bcsd11_ID,$bcsd12_ID" "$jfile")
		    lisfcst_ID=`echo $lisfcst_ID`' '$thisID
		    prevID=$thisID
		else
		    thisID=$(submit_job "$prevID" "$jfile")
		    lisfcst_ID=`echo $lisfcst_ID`' '$thisID
		    prevID=$thisID		
		fi
	    else
		thisID=$(submit_job "$bcsd11_ID,$bcsd12_ID" "$jfile")
		lisfcst_ID=`echo $lisfcst_ID`' '$thisID
	    fi	
	    ((FileNo++))
	done
    done
    lisfcst_ID=`echo $lisfcst_ID | sed "s| |,|g"`

    cd ${BWD}
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

s2spost(){
    
    #######################################################################
    # (5) S2S Post-process
    #######################################################################
    
    jobname=s2spost
    
    echo "              " >> $JOB_SCHEDULE
    echo "(5) S2S post-process                       " >> $JOB_SCHEDULE
    echo "-------------------------------------------" >> $JOB_SCHEDULE
    echo "                                           " >> $JOB_SCHEDULE
    
    for model in $MODELS
    do
	if [ $DATATYPE == "hindcast" ]; then
	    mkdir -p -m 775 ${E2ESDIR}/s2spost/${MM}/${YYYY}${MM}/$model/
	else
	    mkdir -p -m 775 ${E2ESDIR}/s2spost/${YYYY}${MM}/$model/
	fi
    done
    
    cd ${SCRDIR}/s2spost
    
    /bin/ln -s ${E2ESDIR}/lis_fcst/${YYYY}${MM}/ lis_fcst
    /bin/ln -s ${E2ESDIR}/lis_fcst/input
    
    CWD=`pwd`
    for model in $MODELS
    do
	if [ $DATATYPE == "hindcast" ]; then
	    /bin/ln -s ${E2ESDIR}/s2spost/${MM}/${YYYY}${MM}/$model
	else
	    /bin/ln -s ${E2ESDIR}/s2spost/${YYYY}${MM}/$model
	fi
	python $LISHDIR/s2s_modules/s2spost/run_s2spost_9months.py -y ${YYYY} -m ${MM} -w ${CWD} -c $BWD/$CFILE -j $jobname -t 1 -H 4 -M $model
    done
    
    job_list=`ls $jobname*.j`
    s2spost_ID=
    for jfile in $job_list
    do
	thisID=$(submit_job "$lisfcst_ID" "$jfile")
	s2spost_ID=`echo $s2spost_ID`' '$thisID
    done
    s2spost_ID=`echo $s2spost_ID | sed "s| |,|g"`    
    cd ${BWD}
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

s2smetrics(){
    
    #######################################################################
    # (6) S2S metrics
    #######################################################################
    
    jobname=s2smetric
    
    echo "              " >> $JOB_SCHEDULE
    echo "(6) S2S metric                             " >> $JOB_SCHEDULE
    echo "-------------------------------------------" >> $JOB_SCHEDULE
    echo "                                           " >> $JOB_SCHEDULE
    
    mkdir -p -m 775 ${E2ESDIR}/s2smetric/${YYYY}${MM}/DYN_ANOM/
    mkdir -p -m 775 ${E2ESDIR}/s2smetric/${YYYY}${MM}/DYN_SANOM/
    mkdir -p -m 775 ${E2ESDIR}/s2smetric/${YYYY}${MM}/metrics_cf/
    
    cd ${SCRDIR}/s2smetric
    /bin/ln -s ${E2ESDIR}/s2spost
    /bin/ln -s ${E2ESDIR}/s2smetric
    
    CWD=`pwd`
    for model in $MODELS
    do
	python $LISHDIR/s2s_modules/s2smetric/postprocess_nmme_job.py -y ${YYYY} -m ${MM} -w ${CWD} -c $BWD/$CFILE -j $jobname -t 1 -H 12 -M $model
    done
    
    job_list=`ls $jobname*.j`
    s2smetric_ID=
    for jfile in $job_list
    do
	thisID=$(submit_job "$s2spost_ID" "$jfile")
	s2smetric_ID=`echo $s2smetric_ID`' '$thisID
    done
    s2smetric_ID=`echo $s2smetric_ID | sed "s| |,|g"`
    
    # write tiff file
    python $LISHDIR/s2s_app/write_to_file.py -c $BWD/$CFILE -f ${jobname}_tiff_run.j -t 1 -H 2 -j ${jobname}_tiff_ -w ${CWD}
    COMMAND="python $LISHDIR/s2s_modules/s2smetric/postprocess_nmme_job.py -y ${YYYY} -m ${MM} -w ${CWD} -c $BWD/$CFILE"
    sed -i "s|COMMAND|${COMMAND}|g" ${jobname}_tiff_run.j

    s2smetric_tiff_ID=$(submit_job "$s2smetric_ID" "${jobname}_tiff_run.j")
    cd ${BWD}
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

s2splots(){
    
    #######################################################################
    # (7) S2S plots
    #######################################################################
    
    jobname=s2splots
    
    echo "              " >> $JOB_SCHEDULE
    echo "(7) S2S plots                              " >> $JOB_SCHEDULE
    echo "-------------------------------------------" >> $JOB_SCHEDULE
    echo "                                           " >> $JOB_SCHEDULE
    
    mkdir -p -m 775 ${E2ESDIR}/s2splots/${YYYY}${MM}/
    cd ${SCRDIR}/s2splots
    CWD=`pwd`
    /bin/ln -s ${E2ESDIR}/s2splots/
    /bin/ln -s ${E2ESDIR}/s2smetric/ 
    
    python $LISHDIR/s2s_app/write_to_file.py -c $BWD/$CFILE -f ${jobname}_run.j -t 1 -H 12 -j ${jobname}_ -w ${CWD}
    COMMAND="python ${LISHDIR}/s2s_modules/s2splots/plot_s2smetrics.py -y ${YYYY} -m ${MM} -w ${CWD} -c $BWD/$CFILE"
    sed -i "s|COMMAND|${COMMAND}|g" s2splots_run.j
    
    PLINE=`grep -n plot_s2smetrics.py s2splots_run.j | cut -d':' -f1`
    ((PLINE++))
    SEC_COMMAND="python ${LISHDIR}/s2s_modules/s2splots/plot_streamflow_anom.py -y ${YYYY} -m ${mon} -w ${CWD} -c $BWD/$CFILE"
    sed -i "${PLINE}i ${SEC_COMMAND}" s2splots_run.j
    ((PLINE++))
    THIRD_COMMAND="python ${LISHDIR}/s2s_modules/s2splots/plot_hybas.py -y ${YYYY} -m ${MM} -w ${CWD} -c $BWD/$CFILE"
    sed -i "${PLINE}i ${THIRD_COMMAND}" s2splots_run.j
    #((PLINE++))
    #FOURTH_COMMAND="python ${LISHDIR}/s2s_modules/s2splots/plot_precip.py -y ${YYYY} -m ${mon} -w ${CWD} -c $BWD/$CFILE"
    #sed -i "${PLINE}i ${FOURTH_COMMAND}" s2splots_run.j
    #((PLINE++))
    #FIFTH_COMMAND="python ${LISHDIR}/s2s_modules/s2splots/plot_ccdi.py -y ${YYYY} -m ${mon} -w ${CWD} -c $BWD/$CFILE"
    #sed -i "${PLINE}i ${FIFTH_COMMAND}" s2splots_run.j

    s2splots_ID=$(submit_job "$s2smetric_tiff_ID" "${jobname}_run.j")
}

#######################################################################
#**********************************************************************
#                           MAIN SCRIPT
#**********************************************************************

if [ "$REPORT" = 'Y' ] || [ "$REPORT" = 'y' ]; then
    # Print status report
    python $LISHDIR/s2s_app/write_to_file.py -r Y -w ${E2ESDIR} -d ${YYYY}${MM}
    exit
fi

if [ "$DELETE" = 'Y' ] || [ "$DELETE" = 'y' ]; then
    # delete month
    read -p "Are you sure you want to delete ${YYYY}${MM} forecast files entirely (Y/N)?" YESORNO
    if [ "$YESORNO" = 'Y' ] || [ "$YESORNO" = 'y' ]; then
	read -p "I want to double check that I heard you correctly. We want to delete ${YYYY}${MM} forecast files (Y/N)?" YESORNO
	if [ "$YESORNO" = 'Y' ] || [ "$YESORNO" = 'y' ]; then
	    delete_forecast ${YYYY} ${MM}
	fi
    fi
    exit
fi

#######################################################################
#                        Set up scratch directory
#######################################################################
   
SCRDIR=${E2ESDIR}/scratch/${YYYY}${MM}/
mkdir -p -m 775 ${SCRDIR}/global_usaf_forc
mkdir -p -m 775 ${SCRDIR}/lis_darun
mkdir -p -m 775 ${SCRDIR}/ldt_ics
mkdir -p -m 775 ${SCRDIR}/bcsd_fcst
mkdir -p -m 775 ${SCRDIR}/lis_fcst
mkdir -p -m 775 ${SCRDIR}/s2spost

if [ $DATATYPE  == "forecast" ]; then
    mkdir -p -m 775 ${SCRDIR}/s2smetric
    mkdir -p -m 775 ${SCRDIR}/s2splots
    download_forecasts
fi
MODELS=`grep NMME_models $CFILE | cut -d'[' -f2 | cut -d']' -f1 | sed 's/,//g'`

cd ${SCRDIR}/global_usaf_forc
/bin/ln -s $AF10KM usaf_lis73rc8_10km
cd ${SCRDIR}
/bin/ln -s $METFORC
cd ${BWD}
JOB_SCHEDULE=${SCRDIR}/SLURM_JOB_SCHEDULE
/bin/rm -f $JOB_SCHEDULE

echo "#######################################################################" >> $JOB_SCHEDULE
echo "                         SLURM JOB SCHEDULE                            " >> $JOB_SCHEDULE
echo "#######################################################################" >> $JOB_SCHEDULE
echo "                         " >> $JOB_SCHEDULE
python $LISHDIR/s2s_app/write_to_file.py -s $JOB_SCHEDULE -m "JOB ID" -f "JOB SCRIPT" -a "AFTER"

#######################################################################
#                               Submit jobs
#######################################################################

lisda_ID=
ldtics_ID=
bcsd11_ID=
bcsd12_ID=
lisfcst_ID=
s2spost_ID=
s2smetric_ID=
s2splots_ID=

case $STEP in
    LISDA)
	lis_darun
	if [ $ONE == "N" ] || [ $ONE == "n" ]; then
	    ldt_ics
	    bcsd_fcst
	    lis_fcst
	    s2spost
	    if [ $DATATYPE == "forecast" ]; then
		s2smetrics
		s2splots
		set_permission $s2splots_ID
		exit
	    fi
	    set_permission $s2spost_ID
	    exit
	fi
	set_permission $lisda_ID
    ;;
    LDTICS)
	ldt_ics
	if [ $ONE == "N" ] || [ $ONE == "n" ]; then
	    bcsd_fcst
	    lis_fcst
	    s2spost
	    if [ $DATATYPE == "forecast" ]; then
		s2smetrics
		s2splots
		set_permission $s2splots_ID
		exit
	    fi
	    set_permission $s2spost_ID
	    exit
	fi
	set_permission $ldtics_ID
    ;;    
    BCSD)
	bcsd_fcst
	if [ $ONE == "N" ] || [ $ONE == "n" ]; then
	    lis_fcst
	    s2spost
	    if [ $DATATYPE == "forecast" ]; then
		s2smetrics
		s2splots
		set_permission $s2splots_ID
		exit
	    fi
	    set_permission $s2spost_ID
	    exit
	fi
	set_permission $bcsd12_ID
    ;;
    FCST)
	lis_fcst
	if [ $ONE == "N" ] || [ $ONE == "n" ]; then
	    s2spost
	    if [ $DATATYPE == "forecast" ]; then
		s2smetrics
		s2splots
		set_permission $s2splots_ID
		exit
	    fi
	    set_permission $s2spost_ID
	    exit
	fi
	set_permission $lisfcst_ID
    ;;
    POST)
	s2spost
	if [ $ONE == "N" ] || [ $ONE == "n" ]; then
	    if [ $DATATYPE == "forecast" ]; then
		s2smetrics
		s2splots
		set_permission $s2splots_ID
		exit
	    fi
	    set_permission $s2spost_ID
	    exit
	fi
	set_permission $s2spost_ID
    ;;
    METRICS)
	s2smetrics
	if [ $ONE == "N" ] || [ $ONE == "n" ]; then
	    s2splots
	    set_permission $s2splots_ID
	    exit
	fi
	set_permission $s2smetric_ID
    ;;
    PLOTS)
	s2splots
	set_permission $s2splots_ID
	exit
    ;;
    *)
	lis_darun
	ldt_ics
	bcsd_fcst
	lis_fcst
	s2spost
	if [ $DATATYPE == "forecast" ]; then
	    s2smetrics
	    s2splots
	    set_permission $s2splots_ID
	    exit
	fi
	set_permission $s2spost_ID
    ;;
esac
    

