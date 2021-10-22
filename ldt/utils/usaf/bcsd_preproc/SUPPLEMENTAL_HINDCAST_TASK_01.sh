#!/bin/sh
# This script takes the monthly global nmme precip that Abheera generated and
# subsets it to the AFRICOM domain.
#

# Load neccesary modules
module load cdo


declare -a month_abbreviations=('jan' 'feb' 'mar' 'apr' 'may' 'jun' 'jul' 'aug' 'sep' 'oct' 'nov' 'dec')

for ((fcst_init_month=1; fcst_init_month<=12; fcst_init_month++)); do

    fcst_init_date=${month_abbreviations[$fcst_init_month - 1]}'01'      # e.g. jan

    for ((fcst_init_year=2008; fcst_init_year<=2020; fcst_init_year++)); do

        for ((ens_num=1; ens_num<=61; ens_num++)); do
            INDIR='/gpfsm/dnb05/projects/p63/FORECASTS/GEOS5/BCSD_Test/EXPERIMENTS/NMME/data/AF/PRECTOT_Monthly/'${fcst_init_year}'/'${fcst_init_date}'/ens'${ens_num}'/nmme'
            OUTDIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM/data/NMME/raw/Monthly/'${fcst_init_date}'/'${fcst_init_year}'/ens'${ens_num}

            mkdir -p ${OUTDIR}

            fcst0=`date -d "${fcst_init_month}/01/${fcst_init_year} + 0 month" +%Y%m`
            fcst1=`date -d "${fcst_init_month}/01/${fcst_init_year} + 1 month" +%Y%m`
            fcst2=`date -d "${fcst_init_month}/01/${fcst_init_year} + 2 month" +%Y%m`
            fcst3=`date -d "${fcst_init_month}/01/${fcst_init_year} + 3 month" +%Y%m`
            fcst4=`date -d "${fcst_init_month}/01/${fcst_init_year} + 4 month" +%Y%m`
            fcst5=`date -d "${fcst_init_month}/01/${fcst_init_year} + 5 month" +%Y%m`
            fcst6=`date -d "${fcst_init_month}/01/${fcst_init_year} + 6 month" +%Y%m`
            fcst7=`date -d "${fcst_init_month}/01/${fcst_init_year} + 7 month" +%Y%m`
            fcst8=`date -d "${fcst_init_month}/01/${fcst_init_year} + 8 month" +%Y%m`

            echo " == FCSTDATE: "$fcst0", "$fcst1", "$fcst2", "$fcst3", "$fcst4", "$fcst5", "$fcst6", "$fcst7", "$fcst8

            cdo sellonlatbox,-19.875,59.875,-39.875,39.875 ${INDIR}/${fcst_init_date}.nmme.monthly.${fcst0}.nc ${OUTDIR}/${fcst_init_date}.nmme.monthly.${fcst0}.nc
            cdo sellonlatbox,-19.875,59.875,-39.875,39.875 ${INDIR}/${fcst_init_date}.nmme.monthly.${fcst1}.nc ${OUTDIR}/${fcst_init_date}.nmme.monthly.${fcst1}.nc
            cdo sellonlatbox,-19.875,59.875,-39.875,39.875 ${INDIR}/${fcst_init_date}.nmme.monthly.${fcst2}.nc ${OUTDIR}/${fcst_init_date}.nmme.monthly.${fcst2}.nc
            cdo sellonlatbox,-19.875,59.875,-39.875,39.875 ${INDIR}/${fcst_init_date}.nmme.monthly.${fcst3}.nc ${OUTDIR}/${fcst_init_date}.nmme.monthly.${fcst3}.nc
            cdo sellonlatbox,-19.875,59.875,-39.875,39.875 ${INDIR}/${fcst_init_date}.nmme.monthly.${fcst4}.nc ${OUTDIR}/${fcst_init_date}.nmme.monthly.${fcst4}.nc
            cdo sellonlatbox,-19.875,59.875,-39.875,39.875 ${INDIR}/${fcst_init_date}.nmme.monthly.${fcst5}.nc ${OUTDIR}/${fcst_init_date}.nmme.monthly.${fcst5}.nc
            cdo sellonlatbox,-19.875,59.875,-39.875,39.875 ${INDIR}/${fcst_init_date}.nmme.monthly.${fcst6}.nc ${OUTDIR}/${fcst_init_date}.nmme.monthly.${fcst6}.nc
            cdo sellonlatbox,-19.875,59.875,-39.875,39.875 ${INDIR}/${fcst_init_date}.nmme.monthly.${fcst7}.nc ${OUTDIR}/${fcst_init_date}.nmme.monthly.${fcst7}.nc
            cdo sellonlatbox,-19.875,59.875,-39.875,39.875 ${INDIR}/${fcst_init_date}.nmme.monthly.${fcst8}.nc ${OUTDIR}/${fcst_init_date}.nmme.monthly.${fcst8}.nc
        done;
    done;
done;
