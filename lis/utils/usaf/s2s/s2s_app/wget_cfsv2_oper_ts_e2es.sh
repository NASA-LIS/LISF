#!/bin/sh

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.4
# 
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# -------------------------END NOTICE -- DO NOT EDIT-----------------------

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
#  Author:  Sarith Mahanama, NASA/GSFC;   May 22, 2023
# ________________________________________________________________


# A function to find closest 8 days to the CFSv2 icdate in question
neighb_days(){
    vartype=$1
    icdate=$2
    cycle=$3
    mon=$4
    count_days=1
    while [ $count_days -le 4 ];do
	nwdate=`date -d"$icdate -$count_days Day" +%Y%m%d`
	echo "wget ${srcdir}/cfs.${nwdate}/${cycle}/time_grib_01/${vartype}.01.${nwdate}${cycle}.daily.grb2" >> $CFSV2_LOG
	nwdate=`date -d"$icdate +$count_days Day" +%Y%m%d`
	if [ `date -d"$nwdate" +%m` != ${mon} ];then
	    echo "wget ${srcdir}/cfs.${nwdate}/${cycle}/time_grib_01/${vartype}.01.${nwdate}${cycle}.daily.grb2" >> $CFSV2_LOG
	fi
	((count_days++))
    done
    echo "   " >> $CFSV2_LOG
}

print_message(){
    echo "Note: If all recommended substitutes are also not available, you could try a different forecast hour from any of above dates." >> $CFSV2_LOG
    echo ""  >> $CFSV2_LOG
}

ret_code_pipe=$(mktemp)
main_loop() {
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
		if [ "$download" = 'Y' ] || [ "$download" = 'y' ]; then
		    file=${srcdir}/cfs.${icdate}/${cycle}/time_grib_01/${vartype}.01.${icdate}${cycle}.daily.grb2
		    if [ ! -f "${vartype}.01.${icdate}${cycle}.daily.grb2" ]; then
			wget ${file}
		    fi
		fi
	    
		# File check 1: missing file
		if [ ! -f "${vartype}.01.${icdate}${cycle}.daily.grb2" ]; then
		    have_patch=`grep ${vartype}.01.${icdate}${cycle}.daily.grb2 ${patchfile}`
		    if [[ $have_patch == "" ]]; then
			echo "${vartype}.01.${icdate}${cycle}.daily.grb2:  MISSING " >> $CFSV2_LOG
			echo "Possible substitutes in order of preference are:"      >> $CFSV2_LOG
			neighb_days ${vartype} $icdate ${cycle} ${mon}
			ret_code=1
		    fi
		fi
	    
		# File check 2: corrupted file
		if [ -f "${vartype}.01.${icdate}${cycle}.daily.grb2" ]; then
		    python $LISHDIR/s2s_app/s2s_api.py -i "${vartype}.01.${icdate}${cycle}.daily.grb2" -d $yearmo -c $configfile
		    py_code=$?
		
		    if [ $py_code -gt 0 ]; then
			have_patch=`grep ${vartype}.01.${icdate}${cycle}.daily.grb2 ${patchfile}`
			if [[ $have_patch == "" ]]; then
			    echo "${vartype}.01.${icdate}${cycle}.daily.grb2: CORRUPTED ">> $CFSV2_LOG
			    echo "Possible substitutes in order of preference are:"      >> $CFSV2_LOG
			    neighb_days ${vartype} $icdate ${cycle} ${mon}
			    ret_code=1
			else
			    supfile=`grep ${vartype}.01.${icdate}${cycle}.daily.grb2 ${patchfile}  | cut -d',' -f3 | tr -d ' '`
			    python $LISHDIR/s2s_app/s2s_api.py -i "${patchdir}${supfile}" -d $yearmo -c $configfile
			    py_code=$?
			    if [ $py_code -gt 0 ]; then
				echo "${vartype}.01.${icdate}${cycle}.daily.grb2: Replacement ${supfile} is also CORRUPTED!" >> $CFSV2_LOG
				echo "Try downloading the next file (DON'T forget to update ${patchfile}"                    >> $CFSV2_LOG
				neighb_days ${vartype} $icdate ${cycle} ${mon}
				ret_code=1
			    fi		    
			fi
		    fi
		fi
	    done
	done
	cd ../
    done
    echo $ret_code > $ret_code_pipe
}
# ________________________________________________________________
# Main script
# ________________________________________________________________

# process command line arguments

ret_code=0
while getopts y:m:c:d: flag
do
  case "${flag}" in
     y) year=${OPTARG};;
     m) mon=${OPTARG};;
     c) configfile=${OPTARG};;
     d) download=${OPTARG};;
     *) echo "     "
	echo "USAGE: s2s_app/wget_cfsv2_oper_ts_e2es.sh -y YEAR -m MONTH -c FULL_PATH/CONFIG_FILE -d DOWNLOAD"
	echo "     "
	   echo "where MANDATORY input parameters:"
	   echo "---------------------------------"
	   echo "  YEAR:        forecast start year"
	   echo "  MONTH:       forecast start month [1 to 12]"
	   echo "  CONFIG_FILE: E2ES main config file for forecast with the full path of the E2ES directory"
	   echo "  DOWNLOAD: Download CFSv2 forcings (Y/N). If N only the file check will be performed."
	   exit 1 
	   ;;	 	   
  esac
done
if [[ -z "$year" ]] || [[ -z "$mon" ]] || [[ -z "$configfile" ]] || [[ -z "$download" ]]; then
  echo "`basename ${0}`: usage: [-y year] [-m month ] [-c FULL_PATH/config_file] [-d download (Y/N)]"
  exit 1
fi
echo "Year : $year";
echo "Month: $mon";
echo "Configfile: $configfile";
echo

# Read config file and extract information

export NODE_NAME=`uname -n`
if [[ $NODE_NAME =~ discover* ]] || [[ $NODE_NAME =~ borg* ]]; then
    cfsv2datadir=`grep fcst_download_dir $configfile | cut -d':' -f2 | tr -d "[:space:]"`"/Oper_TS/"
else
    cfsv2datadir=`grep fcst_download_dir $configfile | cut -d':' -f2 | tr -d "[:space:]"`
fi
patchfile=`grep supplementarydir $configfile | cut -d':' -f2 | tr -d "[:space:]"`"/bcsd_fcst/patch_files/patch_files_list.txt"
patchdir=`grep supplementarydir $configfile | cut -d':' -f2 | tr -d "[:space:]"`"/bcsd_fcst/patch_files/"
export LISFDIR=`grep LISFDIR $configfile | cut -d':' -f2 | tr -d "[:space:]"`
export LISHDIR=${LISFDIR}/lis/utils/usaf/s2s/
export LISFMOD=`grep LISFMOD $configfile | cut -d':' -f2 | tr -d "[:space:]"`
export SUPDIR=`grep supplementarydir $configfile | cut -d':' -f2 | tr -d "[:space:]"`
export DATATYPE=`grep DATATYPE  $configfile | cut -d':' -f2 | tr -d "[:space:]"`
export E2ESROOT=`grep E2ESDIR $configfile | cut -d':' -f2 | tr -d "[:space:]"`

if [ $DATATYPE == "hindcast" ]; then
    export E2ESDIR=`grep E2ESDIR $configfile | cut -d':' -f2 | tr -d "[:space:]"`"/hindcast/"
else
    export E2ESDIR=`grep E2ESDIR $configfile | cut -d':' -f2 | tr -d "[:space:]"`
fi

if [[ $NODE_NAME =~ discover* ]] || [[ $NODE_NAME =~ borg* ]]; then
    unset LD_LIBRARY_PATH
    source /etc/profile.d/modules.sh
fi
if [ -e $LISFDIR/env/discover/$LISFMOD ]; then
    # Discover: LISFMOD is in LISF repository
    module use -a $LISFDIR/env/discover/
    module --ignore-cache load $LISFMOD
else
    # non-Discover: LISFMOD is in ${supplementary}/env/
    module use -a $SUPDIR/env/
    module load $LISFMOD
fi

umask 022
ulimit -s unlimited
  
# Source data server and directory (Microsoft planetary_computer):
  srcdir=https://noaacfs.blob.core.windows.net/cfs

# _______________________________________________

  echo ""
  echo "CFSv2 local data path  :: "${cfsv2datadir}
  echo "Main CFSv2 source path :: "${srcdir}
  echo ""

  # Make local year directory to download files into:
  mon=$(expr "$mon" + 0)
  if [ $mon -lt 10 ]; then
      mon="0"$mon
  else
      mon=$mon
  fi
  yearmo=${year}${mon}
  cd ${cfsv2datadir}
  echo ${year}
  mkdir -p ${year}
  cd ${year}

  # open CFSv2 missing/corrupted file info log
  SCRDIR=${E2ESDIR}/scratch/${yearmo}/
  mkdir -p -m 775 ${SCRDIR}/
  CFSV2_LOG=${SCRDIR}/CFSv2_missing_corrupted_files
  /bin/rm -f $CFSV2_LOG

echo " #####################################################################################" >> $CFSV2_LOG
echo "                                  MISSING/INCOMPLETE CFSV2 FILES                      " >> $CFSV2_LOG
echo " #####################################################################################" >> $CFSV2_LOG
echo "                         " >> $CFSV2_LOG
echo "  A replacement file is required for each missing or corrupted file. CFSv2 replacement files are saved in:" >> $CFSV2_LOG
echo "  ${patchdir}, "                   >> $CFSV2_LOG 
echo "  and comma-delimited lines in:  " >> $CFSV2_LOG
echo "  ${patchfile} "                   >> $CFSV2_LOG
echo "  lists the replacement file names for each corrupted file. The table has three columns:  " >> $CFSV2_LOG
echo "  YYYYMMDDHH, bad_file_name, replacement_file_name.     " >> $CFSV2_LOG
echo "                         " >> $CFSV2_LOG
echo " (1) cd ${patchdir}      " >> $CFSV2_LOG
echo " (2) Each problematic file name in the section below is followed by a list of wget commands to download a suitable replacement file in order of preference." >> $CFSV2_LOG
echo "     Download the first suggested replacement file and add a new entry to: " >> $CFSV2_LOG
echo "     ${patchfile} " >> $CFSV2_LOG
echo " (3) Repeat the same procedure to download replacements and update: " >> $CFSV2_LOG
echo "     ${patchfile} " >> $CFSV2_LOG
echo "     for every missing/corrupted file." >> $CFSV2_LOG
echo " (4) Relaunch the forecast: s2s_app/s2s_run.sh -y YEAR -m MONTH -c CONFIGFILE" >> $CFSV2_LOG
echo " (5) If any of the replacement files fail, you will be redirected to this file." >> $CFSV2_LOG
echo "     $CFSV2_LOG" >> $CFSV2_LOG
echo " (6) Repeat steps (2) and (3) using a different replacement file for the original bad file.">> $CFSV2_LOG
echo "                         " >> $CFSV2_LOG

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
echo " "
echo "=================================================================================================="
echo " CFSv2 file checker is running to ensure all forcings files are available and not corrupted......"
echo "=================================================================================================="

# Run the main loop
main_loop &

# Display the rotating hyphen animation
animation="-\|/"
while kill -0 $! >/dev/null 2>&1; do
    for (( i=0; i<${#animation}; i++ )); do
        echo -ne "\rPlease wait... ${animation:$i:1}"
        sleep 0.1
    done
done
ret_code=$(cat $ret_code_pipe)
/bin/rm $ret_code_pipe

if [ $ret_code -gt 0 ]; then
    echo "*** Missing or Incomplete CFSv2 forcing files were found ***."
    echo "Please follow the instructions in:"
    echo $CFSV2_LOG
    print_message
else
    echo "**************************************************************">> $CFSV2_LOG
    echo " SUCCESS ! All CFSv2 forcings files passed the file check."    >> $CFSV2_LOG
    echo "**************************************************************">> $CFSV2_LOG

    echo "**************************************************************"
    echo " SUCCESS ! All CFSv2 forcings files passed the file check."    
    echo "**************************************************************"
    
fi  
echo " -- Done checking (and/or downloading) CFSv2 Forecast files -- "

exit $ret_code
# ____________________________

