#!/bin/sh

ulimit -s unlimited


# 1. User-specified entries:
#  Forecast years and NMME model to check:

   mon=${1}
   iMon=${mon}"01"
   year=${2}
   nmmemod=${3}

    for model in ${nmmemod}; do
          echo ${model}
	  if [ $model == CFSv2 ]; then
            ens_numf=12
	  elif [ $model == GEOSv2 ]; then
            ens_numf=4
	  elif [ $model == CCM4 ]; then
            ens_numf=10
	  elif [ $model == GNEMO ]; then
            ens_numf=10
	  elif [ $model == CCSM4 ]; then
            ens_numf=10
	  elif [ $model == GFDL ]; then
	    ens_numf=15
 	  fi
     done

   DIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM/data/NMME/final/6-Hourly/'


   tcnt=$(( 20*ens_numf ))
   tsiz=$(( 26*ens_numf ))

   cnt=`ls ${DIR}/${model}/${year}/${iMon}/ens*/* | wc -l`
   echo "File count (14xens#("${ens_numf}")):" ${nmmemod}, ${year}, ${mon}, ${cnt}"/"$tcnt

   cd ${DIR}/${model}/${year}
   siz=`du -sm $iMon | cut -f1`
   echo "Folder size (~26MBxens#("${ens_numf}")):" ${nmmemod}, ${year}, ${mon}, ${siz}"MB/~"${tsiz}"MB"
