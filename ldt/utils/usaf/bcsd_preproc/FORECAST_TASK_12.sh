#!/bin/sh
# This is a temporary task that prepares an all zero V10M variable for LIS preparation due to the USAF-LIS observational forcing only including average windspeed
module load nco

fcst_init_date=$1'01'
fcst_init_month=$2
ensf=$3

for ((fcst_init_year=2021; fcst_init_year<=2021; fcst_init_year++)); do
    fcst0=`date -d "${fcst_init_month}/01/${fcst_init_year} + 0 month" +%Y%m`
    fcst1=`date -d "${fcst_init_month}/01/${fcst_init_year} + 1 month" +%Y%m`
    fcst2=`date -d "${fcst_init_month}/01/${fcst_init_year} + 2 month" +%Y%m`
    fcst3=`date -d "${fcst_init_month}/01/${fcst_init_year} + 3 month" +%Y%m`
    fcst4=`date -d "${fcst_init_month}/01/${fcst_init_year} + 4 month" +%Y%m`
    fcst5=`date -d "${fcst_init_month}/01/${fcst_init_year} + 5 month" +%Y%m`
    fcst6=`date -d "${fcst_init_month}/01/${fcst_init_year} + 6 month" +%Y%m`
    fcst7=`date -d "${fcst_init_month}/01/${fcst_init_year} + 7 month" +%Y%m`
    fcst8=`date -d "${fcst_init_month}/01/${fcst_init_year} + 8 month" +%Y%m`
    fcst9=`date -d "${fcst_init_month}/01/${fcst_init_year} + 9 month" +%Y%m`
    echo " == FCSTDATE: "$fcst0", "$fcst1", "$fcst2", "$fcst3", "$fcst4", "$fcst5", "$fcst6", "$fcst7", "$fcst8", "$fcst9

    for ((ens_num=1; ens_num<=$ensf; ens_num++)); do
        INDIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM/data/forecasts/CFSv2_25km/final/6-Hourly/'$fcst_init_year'/'${fcst_init_date}'/ens'$ens_num

        ncap2 -O -s 'V10M=array(0.0,0.0,U10M)' ${INDIR}/CFSv2.${fcst0}.nc4 ${INDIR}/CFSv2.${fcst0}.nc4
        ncap2 -O -s 'V10M=array(0.0,0.0,U10M)' ${INDIR}/CFSv2.${fcst1}.nc4 ${INDIR}/CFSv2.${fcst1}.nc4
        ncap2 -O -s 'V10M=array(0.0,0.0,U10M)' ${INDIR}/CFSv2.${fcst2}.nc4 ${INDIR}/CFSv2.${fcst2}.nc4
        ncap2 -O -s 'V10M=array(0.0,0.0,U10M)' ${INDIR}/CFSv2.${fcst3}.nc4 ${INDIR}/CFSv2.${fcst3}.nc4
        ncap2 -O -s 'V10M=array(0.0,0.0,U10M)' ${INDIR}/CFSv2.${fcst4}.nc4 ${INDIR}/CFSv2.${fcst4}.nc4
        ncap2 -O -s 'V10M=array(0.0,0.0,U10M)' ${INDIR}/CFSv2.${fcst5}.nc4 ${INDIR}/CFSv2.${fcst5}.nc4
        ncap2 -O -s 'V10M=array(0.0,0.0,U10M)' ${INDIR}/CFSv2.${fcst6}.nc4 ${INDIR}/CFSv2.${fcst6}.nc4
        ncap2 -O -s 'V10M=array(0.0,0.0,U10M)' ${INDIR}/CFSv2.${fcst7}.nc4 ${INDIR}/CFSv2.${fcst7}.nc4
        ncap2 -O -s 'V10M=array(0.0,0.0,U10M)' ${INDIR}/CFSv2.${fcst8}.nc4 ${INDIR}/CFSv2.${fcst8}.nc4
    done;
done;
