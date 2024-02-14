#!/bin/bash
# 
#  SCRIPT:  download_nmme_hindcasts_sshukla.csh
#
#  PURPOSE:  Download monthly hindcasts of the NMME (version-2) 
#            ensemble of climate forecast models.
#
#  AUTHOR:  Shrad Shukla, UCSB, MAY-2016
#  Edited:  Abheera Hazra, NASA, MAR-2019
#  Edited:  Ryan Zamora, NASA, MAR-2022
#  Edited:  KR Arsenault, NASA, FEB-2024; Added Jan-, Feb-2011 CFSv2 downloads
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
# Default should be 1982 to 2021:


month=$1
CFILE=$2

# Output target directory (where files are to be written):
TARGETDIR=`grep nmme_download_dir $CFILE | cut -d':' -f2 | tr -d "[:space:]"`

# Variable is precipitation
var='prec'

for model in 'NCEP-CFSv2' 'NASA-GEOSS2S' 'CanSIPS-IC3' 'COLA-RSMAS-CCSM4' 'GFDL-SPEAR'
#for model in 'NCEP-CFSv2'
do
    OUTDIR=$TARGETDIR"/"$model
    mkdir -p $OUTDIR

    # rm -rf $OUTDIR"/"$var"."$model".mon_"$month"_"$SYRA"_"$EYRA".nc"

    # CFSv2
    if [ $model == 'NCEP-CFSv2' ]; then
        SYRA=1982
        EYRA=2010

        SYRB=2011
        EYRB=2021

        # Years: Jan-1982 to Dec-2010:
        strA="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.HINDCAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRA"-"$EYRA"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"_"$SYRA"_"$EYRA".nc"

        # Years: Mar-2011 to present:
        strB="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.EARLY_MONTH_SAMPLES/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRB"-"$EYRB"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"_"$SYRB"_"$EYRB".nc"

        # ABOVE DATASETS - MISSING:  Jan- and Feb-2011:
        #  ... So must download separately here:
        if [ $month == 'Jan' -o $month == 'Feb' ]; then
          SYRC=2011
          EYRC=2011

          strC="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.HINDCAST/.PENTAD_SAMPLES/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRC"-"$EYRC"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"_"$SYRC"_"$EYRC".nc"
          echo $strC
          curl -v $strC
        fi

        echo $strA
        curl -v $strA
        echo $strB
        curl -v $strB
    fi

    # GEOSv2
    if [ $model == 'NASA-GEOSS2S' ]; then
        SYRA=1982
        EYRA=2016
        SYRB=2017
        EYRB=2021
        SYRC=1982
        EYRC=2017
        SYRD=2018
        EYRD=2021

        strA="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.HINDCAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRA"-"$EYRA"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"_"$SYRA"_"$EYRA".nc"
        strB="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRB"-"$EYRB"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"_"$SYRB"_"$EYRB".nc"
        strC="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.HINDCAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRC"-"$EYRC"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"_"$SYRC"_"$EYRC".nc"
        strD="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRD"-"$EYRD"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"_"$SYRD"_"$EYRD".nc"
        
        echo $strA
        curl -v $strA
        echo $strB
        curl -v $strB
        echo $strC
        curl -v $strC
        echo $strD
        curl -v $strD
    fi

    # CCM4 & GNEMO5 
    if [ $model == 'CanSIPS-IC3' ]; then
        SYRA=1991
        EYRA=2020
        SYRB=2021
        EYRB=2021

        strA="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.HINDCAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRA"-"$EYRA"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"_"$SYRA"_"$EYRA".nc"
        strB="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRB"-"$EYRB"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"_"$SYRB"_"$EYRB".nc"
        
        echo $strA
        curl -v $strA
        echo $strB
        curl -v $strB
    fi

    # CCSM4 
    if [ $model == 'COLA-RSMAS-CCSM4' ]; then
        SYRA=1982
        EYRA=2021

        strA="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRA"-"$EYRA"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"_"$SYRA"_"$EYRA".nc"
        
        echo $strA
        curl -v $strA
    fi

    # GFDL
    if [ $model == 'GFDL-SPEAR' ]; then 
        SYRA=1991
        EYRA=2020
        SYRB=2021
        EYRB=2021

        strA="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.HINDCAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRA"-"$EYRA"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"_"$SYRA"_"$EYRA".nc"
        strB="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.FORECAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRB"-"$EYRB"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"_"$SYRB"_"$EYRB".nc"
        
        echo $strA
        curl -v $strA
        echo $strB
        curl -v $strB
    fi
done
