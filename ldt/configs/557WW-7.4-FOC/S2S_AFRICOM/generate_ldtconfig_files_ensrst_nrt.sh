#!/bin/sh
#
#  Script Description: Generates ldt.config files for 
#   adjusting NoahMP401 ensemble restart files for the six  
#   different NMME models, and copies+renames restart files
#   for the 1-member HyMAP2 restart files for use as ICs
#   for the forecast runs.
#
#  Author:  Kristi Arsenault (Sept-2021)
# ________________________________________________________

  WORKDIR="./"
  echo " Working directory:  "${WORKDIR}
  cd ${WORKDIR}

# USER LDT CONFIG INPUTS:

# Name of LDT executable
  LDTexec=LDT

# LIS DA Run - Input directory
  inputdir="../lis_darun/output_target/"

# LSM and Routing models:
  lsmname="noahmp401"
  routingname="hymap2"

# LDT Input parameter file:
  ldtinputfile="./input/lis_input.s2s_africom.${lsmname}_${routingname}.25km.nc"

# Final Config script file directory:
  configs_outputdir="./ldt.config_files"
  mkdir -p ${configs_outputdir}

# LIS - Input NoahMP401 557 DA Run 12-member restart files:
  lsm_rstdir="./${inputdir}/SURFACEMODEL/${rst_prevyr}${rst_prevmon}/"
#  HyMAP2 == 1 member
  hymap_rstdir="./${inputdir}/ROUTING/${rst_prevyr}${rst_prevmon}/"

# Template and file naming conventions:
  template_dir="./template_files"
  ldtconfig_lsm_template=${template_dir}"/ldt.config_${lsmname}_nmme_TEMPLATE"


## SCRIPT AUTOMATION SECTION ##

 # Config file start / end dates:
  currentmon=`date +%m`
  currentyear=`date +%Y`
  echo " == Current Date: "$currentmon", "$currentyear

  # Initial month
  input_startmon=`date -d "${currentmon}/01/${currentyear}" +%m`
  echo " Input Startmon = "${input_startmon}
  # Previous month date for restart file:
  rst_prevyr=`date -d "${currentmon}/01/${currentyear} - 1 day" +%Y`
  echo " Previous Start Year: "${rst_prevyr}
  rst_prevmonname=`date -d "${currentmon}/01/${currentyear} - 1 day" +%b`
  echo " Rst Start Month Name: "${rst_prevmonname}
  rst_prevmon=`date -d "${currentmon}/01/${currentyear} - 1 day" +%m`
  echo " Rst Start Month: "${rst_prevmon}
  rst_prevday=`date -d "${currentmon}/01/${currentyear} - 1 day" +%d`
  echo " Rst Start Day: "${rst_prevday}

  start_year=${rst_prevyr}
  final_year=${currentyear}
  echo " Start and final years: "${start_year}", "${final_year}

### ________________________________________________________

# Number of forecast months to run over
  input_numfcstmons=9     

  if [ ${lsmname} == "noahmp401" ]; then
    lsm="NOAHMP401"
  fi

# Generate each set of LDT config files and run LDT to produce final forecast
#  Ensemble Restart Files:

# Loop over each NMME forecast model and set the number of ensemble members:
  for nmmemodel in CCM4  CCSM4  CFSv2  GEOSv2  GFDL  GNEMO; do
   
  # Number of NMME forecast-based ensemble members to expand to:
    case "$nmmemodel" in
      CCM4) num_ensmems="10" ;;
      CCSM4) num_ensmems="10" ;;
      GNEMO) num_ensmems="10" ;;
      GEOSv2) num_ensmems="10" ;;
      CFSv2) num_ensmems="24" ;;
      GFDL) num_ensmems="30" ;;
    esac
    echo "Number of members :: "${num_ensmems}", for model == "${nmmemodel}

  # LDT upscaling or downscaling of ensemble restart file for ICs:
    case "$nmmemodel" in
      CCM4) ldtrstgen="downscale" ;;
      CCSM4) ldtrstgen="downscale" ;;
      GNEMO) ldtrstgen="downscale" ;;
      GEOSv2) ldtrstgen="downscale" ;;
      CFSv2) ldtrstgen="upscale" ;;
      GFDL) ldtrstgen="upscale" ;;
    esac
    echo "LDT ensemble restart generation:: "${ldtrstgen}", for model == "${nmmemodel}

  # Lower-case the names of the 6-NMME models
    lcnmmemodel="${nmmemodel,,}"

  # Generate LDT config files and slurm script naming convention used:
    ldtconfig_nameconv_lsm="ldt.config_${lsmname}_nmme_"${lcnmmemodel}

    echo " Processing LDT Ensemble Restart file for :: "${lsmname}

  # Create the *target* ldt.config file name:
    ldtconfig_lsm_target=${configs_outputdir}"/"${ldtconfig_nameconv_lsm}"_"${currentyear}${currentmon}
    cp ${ldtconfig_lsm_template} ${ldtconfig_lsm_target}

  # Replace TEMPLATE String placeholders with automated info:

  # Input restart filename to be read in by LDT:
    inputfname=${lsm_rstdir}"/LIS_RST_NOAHMP401_"${rst_prevyr}${rst_prevmon}${rst_prevday}".d01.nc"
    echo $inputfname

  # Output restart filename to be written out by LDT:
    rstdate=${rst_prevyr}${rst_prevmon}${rst_prevday}
    mkdir -p ${nmmemodel}
    outputfname=${nmmemodel}"/LIS_RST_NOAHMP401_"${rstdate}"2345.ICS_"${rst_prevmonname}${rst_prevyr}".ens"${num_ensmems}".nc"

  # Log file and output directory entries:
    lsmlogfile=${nmmemodel}"/ldtlog_"${lsmname}"_"${rst_prevmonname}${rst_prevyr}
    maskparmlogfile=${nmmemodel}"/MaskParamFill.log"

  # Replace main strings
    find ${ldtconfig_lsm_target} -type f | xargs perl -pi -e 's|LDTINPUTFILE|'${ldtinputfile}'|g'
    find ${ldtconfig_lsm_target} -type f | xargs perl -pi -e 's|LDTRSTGENOPT|'${ldtrstgen}'|g'
    find ${ldtconfig_lsm_target} -type f | xargs perl -pi -e 's|INPUTRSTFILE|'${inputfname}'|g'
    find ${ldtconfig_lsm_target} -type f | xargs perl -pi -e 's|OUTPUTRSTFILE|'./${outputfname}'|g'
    find ${ldtconfig_lsm_target} -type f | xargs perl -pi -e 's|OUTENSMEMS|'${num_ensmems}'|g'

    find ${ldtconfig_lsm_target} -type f | xargs perl -pi -e 's|LSMLDTLOGFILE|'./${lsmlogfile}'|g'
    find ${ldtconfig_lsm_target} -type f | xargs perl -pi -e 's|PARAMLOGFILE|'./${maskparmlogfile}'|g'
    find ${ldtconfig_lsm_target} -type f | xargs perl -pi -e 's|MODELDIR|'./${nmmemodel}/'|g'

  # Run the LDT executable with the newly generated ldt.config files:
   ./${LDTexec} ${ldtconfig_lsm_target}


  # -- HyMAP2 -- #

  # Copy HyMAP2 restart file and rename, placing 1-member restart file per 
  #  NMME model directory:
    hymap_inrstfile=${hymap_rstdir}"/LIS_RST_HYMAP2_router_"${rstdate}"2345.d01.nc"
    hymap_outrstfile="./"${nmmemodel}"/LIS_RST_HYMAP2_router_"${rstdate}"2345.ICS_"${rst_prevmonname}${rst_prevyr}".ens1.nc"

    cp ${hymap_inrstfile} ${hymap_outrstfile}

  done

 ##### 

  echo " Done generating LDT config files and output restart files. "

