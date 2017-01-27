#!/bin/bash

if [ "$#" -lt 3 ]
then
  echo "ERROR USAGE: $0 <Prefix> <Label 1> <Label 2> ... [Label N]"
  exit 1
fi

echo "*** Combining plots using ImageMagick ***"

PREFIX="$1"
typeset LABELLIST=${@:2}
echo "PREFIX: $PREFIX"
echo "LABELS: $LABELLIST"

# Combine forcing variable plots
typeset FNAMELIST=""
FNAMELIST+="Psurf_f_inst "
FNAMELIST+="Tair_f_inst "
FNAMELIST+="EWind_f_inst "
FNAMELIST+="LWdown_f_inst "
FNAMELIST+="NWind_f_inst "
FNAMELIST+="Rainf_f_tavg "
FNAMELIST+="SWdown_f_inst "
FNAMELIST+="Qair_f_inst"

echo "FORCING VARIABLE NAMES: $FNAMELIST"

for VARNAME in $FNAMELIST; do
typeset FILELIST=""
for LABEL in $LABELLIST; do
  FILELIST+=$LABEL"_forcing_"$VARNAME".gif "
done
convert $FILELIST +append $PREFIX"_forcing_"$VARNAME.gif
done

# Combine output variable plots
typeset ONAMELIST=""
ONAMELIST+="Albedo_inst "
ONAMELIST+="Qg_tavg "
ONAMELIST+="SnowCover_inst "
ONAMELIST+="SnowDepth_inst "
ONAMELIST+="AvgSurfT_tavg "
ONAMELIST+="Qle_tavg "
ONAMELIST+="Qh_tavg "
#ONAMELIST+="effective_mixing_ratio "
ONAMELIST+="SmLiqFrac_inst_Lv1 "
ONAMELIST+="SmLiqFrac_inst_Lv2 "
ONAMELIST+="SmLiqFrac_inst_Lv3 "
ONAMELIST+="SmLiqFrac_inst_Lv4 "
ONAMELIST+="SoilMoist_inst_Lv1 "
ONAMELIST+="SoilMoist_inst_Lv2 "
ONAMELIST+="SoilMoist_inst_Lv3 "
ONAMELIST+="SoilMoist_inst_Lv4 "
ONAMELIST+="SoilTemp_inst_Lv1 "
ONAMELIST+="SoilTemp_inst_Lv2 "
ONAMELIST+="SoilTemp_inst_Lv3 "
ONAMELIST+="SoilTemp_inst_Lv4"

echo "OUTPUT VARIABLE NAMES: $ONAMELIST"

for VARNAME in $ONAMELIST; do
typeset FILELIST=""
for LABEL in $LABELLIST; do
  FILELIST+=$LABEL"_output_"$VARNAME".gif "
done
convert $FILELIST +append $PREFIX"_output_"$VARNAME.gif
done

echo ""

exit
