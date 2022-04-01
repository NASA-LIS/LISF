#!/bin/sh
# Copy files to fill final "6-month" for writing out full average for Jul 1 forecasts
#

  currentmon=${1}
  echo ${currentmon}
  currentyear=${2}
  echo " == Date: "$currentmon", "$currentyear

#

for mo in ${currentmon} ; do

 for yr in ${currentyear} ; do


  year=${yr}
  echo "Year :: "${year}

  if [ ${mo} == jan ]; then
      mo2='09'
      mo3='10'
      let year1=$year
      let year2=$year
  elif [ ${mo} == feb ]; then
      mo2='10'
      mo3='11'
      let year1=$year
      let year2=$year
  elif [ ${mo} == mar ]; then
      mo2='11'
      mo3='12'
      let year1=$year
      let year2=$year
  elif [ ${mo} == apr ]; then
      mo2='12'
      mo3='01'
      let year1=$year
      let year2=year+1
  elif [ ${mo} == may ]; then
      mo2='01'
      mo3='02'
      let year1=year+1
      let year2=year+1
  elif [ ${mo} == jun ]; then
      mo2='02'
      mo3='03'
      let year1=year+1
      let year2=year+1
  elif [ ${mo} == jul ]; then
      mo2='03'
      mo3='04'
      let year1=year+1
      let year2=year+1
  elif [ ${mo} == aug ]; then
      mo2='04'
      mo3='05'
      let year1=year+1
      let year2=year+1
  elif [ ${mo} == sep ]; then
      mo2='05'
      mo3='06'
      let year1=year+1
      let year2=year+1
  elif [ ${mo} == oct ]; then
      mo2='06'
      mo3='07'
      let year1=year+1
      let year2=year+1
  elif [ ${mo} == nov ]; then
      mo2='07'
      mo3='08'
      let year1=year+1
      let year2=year+1
  elif [ ${mo} == dec ]; then
      mo2='08'
      mo3='09'
      let year1=year+1
      let year2=year+1
  fi

 NMME_DATA_DIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM/data/NMME/final/6-Hourly'
 for model in CCSM4 CCM4 GNEMO ; do

   cd ${NMME_DATA_DIR}/${model}/${year}/${mo}01/

   for member in ens1 ens2 ens3 ens4 ens5 ens6 ens7 ens8 ens9 ens10; do

     cd ${member}

     pwd

       echo ${member}" :: PRECTOT."${year1}${mo2}".nc4; PRECTOT."${year2}${mo3}".nc4"


       cp PRECTOT.${year1}${mo2}.nc4 PRECTOT.${year2}${mo3}.nc4


     cd ../
    done
   done



 for model in GEOSv2 ; do


   cd ${NMME_DATA_DIR}/${model}/${year}/${mo}01/

   for member in ens1 ens2 ens3 ens4; do

     cd ${member}

     pwd

        echo ${member}" :: PRECTOT."${year1}${mo2}".nc4; PRECTOT."${year2}${mo3}".nc4"


       cp PRECTOT.${year1}${mo2}.nc4 PRECTOT.${year2}${mo3}.nc4


     cd ../
    done
   done



 for model in CFSv2 ; do


   cd ${NMME_DATA_DIR}/${model}/${year}/${mo}01/

   for member in ens1 ens2 ens3 ens4 ens5 ens6 ens7 ens8 ens9 ens10 ens11 ens12; do

     cd ${member}

     pwd

       echo ${member}" :: PRECTOT."${year1}${mo2}".nc4; PRECTOT."${year2}${mo3}".nc4"

       cp PRECTOT.${year1}${mo2}.nc4 PRECTOT.${year2}${mo3}.nc4


     cd ../
    done
   done


 for model in GFDL ; do


   cd ${NMME_DATA_DIR}/${model}/${year}/${mo}01/

   for member in ens1 ens2 ens3 ens4 ens5 ens6 ens7 ens8 ens9 ens10 ens11 ens12 ens13 ens14 ens15; do

     cd ${member}

     pwd

       echo ${member}" :: PRECTOT."${year1}${mo2}".nc4; PRECTOT."${year2}${mo3}".nc4"

       cp PRECTOT.${year1}${mo2}.nc4 PRECTOT.${year2}${mo3}.nc4


     cd ../
    done
   done



 done

done





