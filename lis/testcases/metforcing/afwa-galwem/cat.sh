#!/bin/sh

#cat sp.bin t.bin gh.bin r.bin 10u.bin 10v.bin > forcing.bin
for hr in 00 06 12 18
do
   cat 20160411_CY.${hr}_FH.000_DF.GR2_sp.bin  \
       20160411_CY.${hr}_FH.000_DF.GR2_t.bin   \
       20160411_CY.${hr}_FH.000_DF.GR2_gh.bin  \
       20160411_CY.${hr}_FH.000_DF.GR2_r.bin   \
       20160411_CY.${hr}_FH.000_DF.GR2_10u.bin \
       20160411_CY.${hr}_FH.000_DF.GR2_10v.bin > forcing_${hr}.bin
done
