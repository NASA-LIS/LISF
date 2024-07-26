#!/bin/bash

MONTH=$1
CFILE=$2
BCSD=$3
NJOBS=`qstat -u $USER | grep $USER | wc -l`
clim_syr=`grep clim_start_year  $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
clim_eyr=`grep clim_end_year  $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
((clim_syr--))

if [[ $MONTH -lt 10 ]]; then
    MM="0"$MONTH
else
    MM=$MONTH
fi

current_year=hindcast/scratch/$MM/current_year_$MM
bcsd_year=hindcast/scratch/$MM/bcsd_year_$MM
NJMAX_BCSD=207
NJMAX_FCST=130

if [[ $BCSD == "BCSD" ]]; then
    if [ -e $bcsd_year ]; then
	echo "last completed year: " `cat $bcsd_year`
    else
	echo "missing " $bcsd_year
	echo $clim_syr > $bcsd_year
    fi
    year=`cat $bcsd_year`
    update_file=$bcsd_year
    NJMAX=$NJMAX_BCSD
else
    if [ -e $current_year ]; then
	echo "last completed year: " `cat $current_year`
    else
	echo "missing " $current_year
	echo $clim_syr > $current_year
    fi
    year=`cat $current_year`
    update_file=$current_year
    NJMAX=$NJMAX_FCST
fi

if [ $NJOBS -gt $NJMAX ]; then
    echo "Currently, " $NJOBS " are in the queue, please try later." 
    exit
fi

# submit jobs

((year++))

if [ $NJOBS -le $NJMAX ] && [ $year -le $clim_eyr ]; then
    # NCCS permits 256 jobs per user at a given time. The s2s_run.sh (hindcast) consists of `100 jobs
    if [[ $BCSD == "BCSD" ]]; then
	#s2s_app/s2s_run.sh -y $year -m $MONTH -c $CFILE -s LDTICS -o Y
	s2s_app/s2s_run.sh -y $year -m $MONTH -c $CFILE -s BCSD -o Y
	echo $year > $update_file
    else
	#s2s_app/s2s_run.sh -y $year -m $MONTH -c $CFILE -s FCST
        s2s_app/s2s_run.sh -y $year -m $MONTH -c $CFILE
	echo $year > $update_file
    fi
fi
chmod 664 $update_file
