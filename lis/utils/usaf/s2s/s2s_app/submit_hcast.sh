#!/bin/bash

MONTH=$1
CFILE=$2
NJOBS=`qstat -u $USER | grep $USER | wc -l`
clim_syr=`grep clim_start_year  $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
clim_eyr=`grep clim_end_year  $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
((clim_syr--))

if [[ $mon -lt 10 ]]; then
    MM="0"$MONTH
else
    MM=$MONTH
fi

if [ -e current_year_$MM ]; then
    echo "last completed year: " `cat current_year_$MM`
else
    echo "missing " current_year_$MM
    echo $clim_syr > current_year_$MM
fi

year=`cat current_year_$MM`
((year++))

if [ $NJOBS -lt 120 ] && [ $year -le $clim_eyr ]; then
    # NCCS permits 256 jobs per user at a given time. The s2s_run.sh (hindcast) consists of `100 jobs  
    s2s_app/s2s_run.sh -y $year -m $MONTH -c $CFILE
    echo $year > current_year_$MM
fi
