#!/bin/bash
# 
#  SCRIPT:  download_nmme_hindcasts_sshukla.csh
#
#  PURPOSE:  Download monthly hindcasts of the NMME (version-2) 
#            ensemble of climate forecast models.
#
#  AUTHOR:  Shrad Shukla, UCSB, MAY-2016
#  Edited:  Abheera Hazra, NASA, MAR-2019
#  Edited:  Kristi Arsenault, NASA, Sep-2024;
#            Updated data downloads to latest NMME models:
#             CanSIPS-IC4(CanESM5,GEM5.2-NEMO), COLA-RSMAS-CESM1
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


# Specify current year, mon (and in config file) to process:

 echo " == NMME precip (current) == "
 echo " == Enter:  yyyy Mon config_file "

currentyear=$1
currentmon=$2
CFILE=$3

SYR=$currentyear
EYR=$currentyear

# Output target directory (where files are to be written):
TARGETDIR=`grep nmme_download_dir $CFILE | cut -d':' -f2 | tr -d "[:space:]"`

# NMME model ensemble suite, valid upto Aug-2024:
#for model in 'NCEP-CFSv2' 'NASA-GEOSS2S' 'CanSIPS-IC3' 'COLA-RSMAS-CCSM4' 'GFDL-SPEAR'
# NMME model ensemble suite, valid starting Aug-2024:
for model in 'NCEP-CFSv2' 'GFDL-SPEAR' 'NASA-GEOSS2S' 'CanSIPS-IC4' 'COLA-RSMAS-CESM1'
do
    OUTDIR=$TARGETDIR"/"$model
    mkdir -p $OUTDIR

    # Loop over main fields and months:
    for var in 'prec'; do
	for month in $currentmon; do
	    zero=0
	    rm -rf $OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"

            # Older ECCC models: 
	    if [ $model == 'CMC1-CanCM3' ] || [ $model == 'CMC2-CanCM4' ]; then
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi

	    if [ $model == 'GEM-NEMO' ] || [ $model == 'CanCM4i' ]; then
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi

	    if [ $model == 'CanSIPS-IC3' ]; then
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi

            # New ECCC CanSIPS-IC4 models: CanESM5 & GEM5.2-NEMO
	    if [ $model == 'CanSIPS-IC4' ]; then

              # CanESM5:
		strA="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.CanESM5/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var".CanESM5.mon_"$month"."$SYR".nc"
                echo $strA
                curl -v $strA

              # GEM5.2-NEMO:
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.GEM5.2-NEMO/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var".GEM5.2-NEMO.mon_"$month"."$SYR".nc"
	    fi

            # RSMAS-COLA NCAR-based models (latest==CESM1):
	    if [ $model == 'COLA-RSMAS-CCSM3' ] || [ $model == 'COLA-RSMAS-CCSM4' ] || [ $model == 'COLA-RSMAS-CESM1' ]; then
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi

            # GFDL Models (latest==GFDL-SPEAR):
	    if [ $model == 'GFDL-SPEAR' ] || [ $model == 'GFDL-CM2p5-FLOR-A06' ] || [ $model == 'GFDL-CM2p5-FLOR-B01' ] || [ $model == 'GFDL-CM2p1-aer04' ]; then 
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi
	    
            # Older NASA Models:
	    if [ $model == 'NASA-GMAO-062012' ] || [ $model == 'NASA-GEOSS2S' ]; then 
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi

            # Current NASA Model:
	    if [ $model == 'NASA-GEOSS2S' ]; then
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi
	    
            # CFSv2 Model:
	    if [ $model == 'NCEP-CFSv2' ]; then
		str="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.EARLY_MONTH_SAMPLES/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYR"-"$EYR"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"."$SYR".nc"
	    fi

            # Download the model files via curl:
	    echo $str
	    curl -v $str

    done ## month
  done   ## variable
done
