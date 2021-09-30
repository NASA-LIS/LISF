#!/bin/sh

# currentmon=`date +%m`
 currentmon=10

 echo ${currentmon}

# currentyear=`date +%Y`
currentyear=2020
 echo ${currentyear}

 CSYR=2008
 CEYR=2019

 year=`date -d "${currentmon}/01/${currentyear} - 1 day" +%Y`
 echo ${year}
 
 prev_mon=`date -d "${currentmon}/01/${currentyear} - 1 day" +%m`
 echo ${prev_mon}

 prev_yr=`date -d "${currentmon}/01/${currentyear} - 1 day" +%Y`
 echo ${prev_yr}


 rundir='/discover/nobackup/projects/nca/karsenau/S2S/AFRICOM/share/final_code/lvt/utils/usaf/s2smetric/' 

 echo $rundir
 cd $rundir/

source /usr/share/modules/init/sh

#for model in CFSV2 CCSM4 CCM4 GEOSV2 GNEMO GFDL; do
for model in CFSV2; do 
 echo $model
  sh job_run_convert_Dyn_FCST_to_postproc.scr ${currentmon} NOAHMP AFRICOM 9 ${currentyear} ${model} ${CSYR} ${CEYR}


done

echo 'Model Processing submitted.'

