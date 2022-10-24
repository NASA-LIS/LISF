#!/bin/sh
#
#  Title:  CFSv2 Operational Timeseries (Oper_TS) Download Script
#
#  Source:  NOAA NCEI Data Archive Center
#
# -- The initial forecast month date-members (making up first 12-members):
#
#   mon01      MMDD    MMDD    MMDD
#   =====      ====    ====    ====
#  "jan01" : ['1217', '1222', '1227']
#  "feb01" : ['0121', '0126', '0131']
#  "mar01" : ['0215', '0220', '0225']
#  "apr01" : ['0317', '0322', '0327']
#  "may01" : ['0416', '0421', '0426']
#  "jun01" : ['0521', '0526', '0531']
#  "jul01" : ['0620', '0625', '0630']
#  "aug01" : ['0720', '0725', '0730']
#  "sep01" : ['0819', '0824', '0829']
#  "oct01" : ['0918', '0923', '0928']
#  "nov01" : ['1018', '1023', '1028']
#  "dec01" : ['1117', '1122', '1127']
#
#  Author:  Kristi Arsenault, NASA/GSFC;  April 15, 2021
#  Author:  Kristi Arsenault, NASA/GSFC;   Oct   3, 2022
# ________________________________________________________________

# Initial date
# e.g.,
#  year=2022
#  mon=01

  while getopts y:m:c: flag
  do
    case "${flag}" in
       y) year=${OPTARG};;
       m) mon=${OPTARG};;
       c) cfsv2datadir=${OPTARG};;
    esac
  done
  if [[ -z "$year" ]] || [[ -z "$mon" ]] || [[ -z "$cfsv2datadir" ]]; then
    echo "`basename ${0}`: usage: [-y year] [-m month] [-c CFSv2_local_directory]"
    exit 1
  fi
  echo "Year : $year";
  echo "Month: $mon";
  echo 

# Local directory path:
#  cfsv2datadir=/discover/nobackup/projects/usaf_lis/GHI_S2S/CFSv2

# Source data server and directory:
  srcdir=https://www.ncei.noaa.gov/data/climate-forecast-system/access/operational-9-month-forecast/time-series/
 
  # Note: Can set the "cfsv2datadir" in the s2s.config as entries

# _______________________________________________

  echo ""
  echo "CFSv2 local data path  :: "${cfsv2datadir}
  echo "Main CFSv2 source path :: "${srcdir}
  echo ""

  # Make local year directory to download files into:
  yearmo=${year}${mon}
  cd ${cfsv2datadir}
  echo ${year}
  mkdir -p ${year}
  cd ${year}

# Loop over and download each date for given forecast initial date and cycle:

  if [ ${mon} -eq "01" ]; then
#  "jan01" : ['1217', '1222', '1227']
    echo "January ..."
    prevmon=12
    year2=$((year-1))
    day1=17
    day2=22
    day3=27

  elif [ ${mon} -eq "02" ]; then
#  "feb01" : ['0121', '0126', '0131']
    echo "February ..."
    prevmon=01
    year2=${year}
    day1=21
    day2=26
    day3=31

  elif [ ${mon} -eq "03" ]; then
#  "mar01" : ['0215', '0220', '0225']
    echo "March ..."
    prevmon=02
    year2=${year}
    day1=15
    day2=20
    day3=25

  elif [ ${mon} -eq "04" ]; then
#  "apr01" : ['0317', '0322', '0327']
    echo "April ..."
    prevmon=03
    year2=${year}
    day1=17
    day2=22
    day3=27

  elif [ ${mon} -eq "05" ]; then
#  "may01" : ['0416', '0421', '0426']
    echo "May ..."
    prevmon=04
    year2=${year}
    day1=16
    day2=21
    day3=26

  elif [ ${mon} -eq "06" ]; then
#  "jun01" : ['0521', '0526', '0531']
    echo "June ..."
    prevmon=05
    year2=${year}
    day1=21
    day2=26
    day3=31

  elif [ ${mon} -eq "07" ]; then
#  "jul01" : ['0620', '0625', '0630']
    echo "July ..."
    prevmon=06
    year2=${year}
    day1=20
    day2=25
    day3=30

  elif [ ${mon} -eq "08" ]; then
#  "aug01" : ['0720', '0725', '0730']
    echo "August ..."
    prevmon=07
    year2=${year}
    day1=20
    day2=25
    day3=30

  elif [ ${mon} -eq "09" ]; then
#  "sep01" : ['0819', '0824', '0829']
    echo "September ..."
    prevmon=08
    year2=${year}
    day1=19
    day2=24
    day3=29

  elif [ ${mon} -eq "10" ]; then
#  "oct01" : ['0918', '0923', '0928']
    echo "October ..."
    prevmon=09
    year2=${year}
    day1=18
    day2=23
    day3=28

  elif [ ${mon} -eq "11" ]; then
#  "nov01" : ['1018', '1023', '1028']
    echo "November ..."
    prevmon=10
    year2=${year}
    day1=18
    day2=23
    day3=28

  elif [ ${mon} -eq "12" ]; then
#  "dec01" : ['1117', '1122', '1127']
    echo "December ..."
    prevmon=11
    year2=${year}
    day1=17
    day2=22
    day3=27
  fi
  echo "Previous mon,days 1-2-3 :: "${prevmon}", "${day1}"-"${day2}"-"${day3}

  # Initial forecast dates:
  for prevmondays in ${day1} ${day2} ${day3}; do

    icdate=${year2}${prevmon}${prevmondays}
    #echo ${icdate}
    mkdir -p ${icdate}
    cd ${icdate}

    # Loop over variable type:
    for vartype in dlwsfc dswsfc q2m wnd10m prate tmp2m pressfc; do

      # Forecast cycle (00,06,12,18):
      for cycle in 00 06 12 18; do
        #echo "Cycle :: "${cycle}

        # File to be downloaded:
        file=${srcdir}/${year2}/${year2}${prevmon}/${icdate}/${icdate}${cycle}/${vartype}.01.${icdate}${cycle}.daily.grb2
	if [ ! -f "${vartype}.01.${icdate}${cycle}.daily.grb2" ]; then
            wget --no-check-certificate -nc -nv ${file}
	fi
	if [ ! -f "${vartype}.01.${icdate}${cycle}.daily.grb2" ]; then
	    echo "${vartype}.01.${icdate}${cycle}.daily.grb2 is not available!"; exit 1;
	fi

      done
    done
    cd ../
  done
  cd ${cfsv2datadir} 
  echo " -- Done downloading CFSv2 Reforecast files -- "
exit 0
# ____________________________


# Some example website documentation for wget downloads:
#
#  https://oceanobservatories.org/knowledgebase/how-can-i-download-all-files-at-once-from-a-data-request/
#  url=https://www.ncei.noaa.gov/thredds/catalog/model-cfs_refor_6h_9m_flx/2007/200711/20071127/catalog.html
#  wget -r -l2 -nd -nc -np -e robots=off -A.grb2 --no-check-certificate  ${url}
#
