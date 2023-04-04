#!/bin/sh
#
#  GEFS Operational Product: 
#
#   pgrb2ap5
#   - 0.5 deg
#   - Control run (1-member) + perturbed runs (20-members)
#   - Forecast files > f000 or 03Z and up, include all forcing fields needed.
#   - Analysis and 00Z forecast fields are missing some required forcing fields.
#     * NOTE: This applies to all 4 cycles
#   - 30 members are now available with latest version 12
#
#  K. Arsenault; July 22, 2019;  Initial code
#  K. Arsenault; Jan. 26, 2021;  Automate download of files
# ____________________________________________________________________________

# Local GEFS directory:
  #maindir=/discover/nobackup/projects/lis/MET_FORCING/GEFS/Oper
  maindir=$(pwd)

# HTTPS GEFS data source -- NOMADS server:
  httpsdir=https://nomads.ncep.noaa.gov/pub/data/nccf/com/gens/prod/

# Backup HTTPS GEFS - FTP server:
#  httpsdir=https://ftp.ncep.noaa.gov/data/nccf/com/gens/prod/

# Date
# Check if yesterday's daily file is available for independent check for wget download:
  date=20220102
#  date=`date +%Y%m%d`
  echo ${date}

# Change directory to the GEFS-Operational download directory:
  mkdir -p ${maindir}/${date}
  cd ${maindir}/${date}

# Cycle
  for CC in 00 06 12 18; do
  #for CC in 18; do

    mkdir -p ${maindir}/${date}/${CC}
    cd ${maindir}/${date}/${CC}

    echo " Downloading GEFS files for :: gefs."${date}", for cycle :: "${CC}

    for f in 03 06 27; do
    #for f in 09; do
        
        if [ ${CC} != "00" ] && [ ${f} == "27" ]; then
            continue
        fi

        for ens in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21; do
            if [ ${ens} == "00" ]; then
                gec_url=${httpsdir}/gefs.${date}/${CC}/atmos/pgrb2ap5/gec${ens}.t${CC}z.pgrb2a.0p50.f0${f}
#                echo ${gec_url}

                wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r -c -nH -nd -np -e robots=off --reject idx  -A 'gec*' ${gec_url}
            else
                gep_url=${httpsdir}/gefs.${date}/${CC}/atmos/pgrb2ap5/gep${ens}.t${CC}z.pgrb2a.0p50.f0${f}
                #echo ${gep_url}
 
                wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r -c -nH -nd -np -e robots=off --reject idx  -A 'gep*' ${gep_url}
            fi
            sleep 3s;
        done
    done

    cd ${maindir}/${date}

  done

 exit 0 
# Control member ("gecCC"):
 # Analysis
 wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gens/prod/gefs.20190620/${CC}/pgrb2ap5/gec00.t00z.pgrb2a.0p50.anl

 # Forecast
 wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gens/prod/gefs.20190620/${CC}/pgrb2ap5/gec00.t00z.pgrb2a.0p50.f000

# Perturbed member ("gepNN")
 # Analysis
 wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gens/prod/gefs.20190620/${CC}/pgrb2ap5/gep01.t00z.pgrb2a.0p50.anl
 # Forecast
 wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gens/prod/gefs.20190620/${CC}/pgrb2ap5/gep01.t00z.pgrb2a.0p50.f000

