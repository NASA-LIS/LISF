#!/bin/bash
# 
#  SCRIPT:  download_nmme_hindcasts.sh
#
#  PURPOSE:  Download monthly hindcasts of the NMME (version-2) 
#            ensemble of climate forecast models.
#
#  AUTHOR:  Shrad Shukla, UCSB, MAY-2016; orig: download_nmme_hindcasts_sshukla.csh
#  Edited:  Abheera Hazra, NASA, MAR-2019
#  Edited:  Ryan Zamora, NASA, MAR-2022
#  Edited:  KR Arsenault, NASA, FEB-2024; Added Jan-, Feb-2011 CFSv2 downloads
#  Edited:  Kristi Arsenault, NASA, Aug-2024;
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
#_______________________________________________________________________________

# Specify start and end years to process:
# Previous default was 1982 to 2021
# Latest default should be 1991 to 2020 
#  (Note: Allowing up to 2021, where models have available data)

month=$1
CFILE=$2

# Output target directory (where files are to be written):
TARGETDIR=`grep nmme_download_dir $CFILE | cut -d':' -f2 | tr -d "[:space:]"`
echo " --- "
echo " --- Where NMME files will be written to:  "${TARGETDIR}
echo " --- "

# Variable is precipitation
var='prec'

# NMME model ensemble suite, valid upto Aug-2024:
#for model in 'NCEP-CFSv2' 'NASA-GEOSS2S' 'CanSIPS-IC3' 'COLA-RSMAS-CCSM4' 'GFDL-SPEAR'
# NMME model ensemble suite, valid starting Aug-2024:
for model in 'NCEP-CFSv2' 'GFDL-SPEAR' 'NASA-GEOSS2S' 'CanSIPS-IC4' 'COLA-RSMAS-CESM1'
do
    OUTDIR=$TARGETDIR"/"$model
    mkdir -p $OUTDIR

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

    # CanSIPS-IC3: CCM4 & GNEMO5 
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

    # CanSIPS-IC4: CanESM5 & GEM5.2-NEMO
    if [ $model == 'CanSIPS-IC4' ]; then
        SYRA=1991
        EYRA=2020

        # CanESM5:
        modelA="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.CanESM5/.HINDCAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRA"-"$EYRA"%29VALUES/data.nc -o "$OUTDIR"/"$var".CanESM5.mon_"$month"_"$SYRA"_"$EYRA".nc"

        # GEM5.2-NEMO:
        modelB="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.GEM5.2-NEMO/.HINDCAST/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRA"-"$EYRA"%29VALUES/data.nc -o "$OUTDIR"/"$var".GEM5.2-NEMO.mon_"$month"_"$SYRA"_"$EYRA".nc"

        # Download CanSIPS-IC4 models:
        echo $modelA
        echo $modelB
        curl -v $modelA
        curl -v $modelB
    fi

    # COLA-RSMAS/NCAR CCSM4 
    if [ $model == 'COLA-RSMAS-CCSM4' ]; then
        SYRA=1982
        EYRA=2021

        strA="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/."$model"/.MONTHLY/."$var"/S/%280000%201%20"$month"%20"$SYRA"-"$EYRA"%29VALUES/data.nc -o "$OUTDIR"/"$var"."$model".mon_"$month"_"$SYRA"_"$EYRA".nc"
        
        echo $strA
        curl -v $strA
    fi

    # COLA-RSMAS/NCAR CESM1
    if [ $model == 'COLA-RSMAS-CESM1' ]; then
        SYRA=1991
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
