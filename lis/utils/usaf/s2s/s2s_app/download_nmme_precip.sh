#!/bin/bash
# 
#  SCRIPT:  download_nmme_hindcasts_sshukla.csh
#
#  PURPOSE:  Download monthly hindcasts of the NMME (version-2) 
#            ensemble of climate forecast models.
#
#  AUTHOR:  Shrad Shukla, UCSB, MAY-2016
#  Edited:  Abheera Hazra, NASA, MAR-2019
#
#  SOURCE WEBSITE:  https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/
#
#  REFERENCE:  
#  Kirtman, B. P., and Coauthors, 2014: The North American Multimodel Ensemble: 
#   Phase-1 seasonal-to-interannual prediction; Phase-2 toward developing 
#   intraseasonal prediction. Bull. Amer. Meteor. Soc., 95, 585–601. 
#   doi: http://dx.doi.org/10.1175/BAMS-D-12-00050.1
#
#__________________________________________________________________


# Specify start and end years to process:
# Default should be 2008 to 2020:

currentyear=$1
currentmon=$2
CFILE=$3

SYR=$currentyear
EYR=$currentyear

# Output target directory (where files are to be written):
TARGETDIR=`grep nmme_download_dir $CFILE | cut -d':' -f2 | tr -d "[:space:]"`

for model in 'NCEP-CFSv2' 'NASA-GEOSS2S' 'CanSIPS-IC3' 'COLA-RSMAS-CCSM4' 'GFDL-SPEAR'
do
    OUTDIR=$TARGETDIR"/"$model
    mkdir -p $OUTDIR

    # Loop over main fields and months:
    for var in 'prec'; do
	for month in $currentmon; do
	    zero=0
	    rm -rf $OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"

	    if [ $model == 'CMC1-CanCM3' ] || [ $model == 'CMC2-CanCM4' ]; then
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi

	    if [ $model == 'GEM-NEMO' ] || [ $model == 'CanCM4i' ]; then
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi

	    if [ $model == 'CanSIPS-IC3' ]; then
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi
	    
	    if [ $model == 'COLA-RSMAS-CCSM3' ] || [ $model == 'COLA-RSMAS-CCSM4' ]; then
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi

	    if [ $model == 'GFDL-SPEAR' ] || [ $model == 'GFDL-CM2p5-FLOR-A06' ] || [ $model == 'GFDL-CM2p5-FLOR-B01' ] || [ $model == 'GFDL-CM2p1-aer04' ]; then 
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi
	    
	    if [ $model == 'NASA-GMAO-062012' ] || [ $model == 'NASA-GEOSS2S' ]; then 
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi

	    if [ $model == 'NASA-GEOSS2S' ]; then
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi
	    
	    if [ $model == 'NCEP-CFSv2' ]; then
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.EARLY_MONTH_SAMPLES/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi
	    echo $str
	    curl -v $str
    done ## month
  done   ## variable
done
