#!/bin/bash

# usage function
function usage {
  echo "" 2>&1
  echo "------------------------------------------------------------------" 2>&1
  echo "Usage: $myname [options]" 2>&1
  echo "" 2>&1
  echo "Options:" 2>&1
  echo "  -h   : help - print usage and exit" 2>&1
  echo "  -d   : domain list - \"d01,...,dnn\"" 2>&1
  echo "  -t   : datetime list - \"YYYYMMDDHRMN,...,YYYYMMDDHRMN\"" 2>&1
  echo "" 2>&1
  echo "------------------------------------------------------------------" 2>&1
  echo "" 2>&1
}

# convert to single image function
function combine_noah33_plots {

if [ "$#" -lt 3 ]
then
  echo "ERROR USAGE: $0 <Prefix> <Label 1> <Label 2> ... [Label N]"
  exit 1
fi

echo "*** Combining Noah.3.3 plots using ImageMagick ***"

PREFIX="$1"
typeset LABELLIST=${@:2}
echo "PREFIX: $PREFIX"
echo "LABELS: $LABELLIST"

# Combine variable plots
typeset VARLIST=""
VARLIST+="forcing_Psurf_f_inst "
VARLIST+="forcing_Tair_f_inst "
VARLIST+="forcing_EWind_f_inst "
VARLIST+="forcing_LWdown_f_inst "
VARLIST+="forcing_NWind_f_inst "
VARLIST+="forcing_Rainf_f_tavg "
VARLIST+="forcing_SWdown_f_inst "
VARLIST+="forcing_Qair_f_inst "
VARLIST+="output_Albedo_inst "
VARLIST+="output_Qg_tavg "
VARLIST+="output_SnowCover_inst "
VARLIST+="output_SnowDepth_inst "
VARLIST+="output_AvgSurfT_tavg "
VARLIST+="output_Qle_tavg "
VARLIST+="output_Qh_tavg "
#VARLIST+="output_effective_mixing_ratio "
VARLIST+="output_SmLiqFrac_inst_Lv1 "
VARLIST+="output_SmLiqFrac_inst_Lv2 "
VARLIST+="output_SmLiqFrac_inst_Lv3 "
VARLIST+="output_SmLiqFrac_inst_Lv4 "
VARLIST+="output_SoilMoist_inst_Lv1 "
VARLIST+="output_SoilMoist_inst_Lv2 "
VARLIST+="output_SoilMoist_inst_Lv3 "
VARLIST+="output_SoilMoist_inst_Lv4 "
VARLIST+="output_SoilTemp_inst_Lv1 "
VARLIST+="output_SoilTemp_inst_Lv2 "
VARLIST+="output_SoilTemp_inst_Lv3 "
VARLIST+="output_SoilTemp_inst_Lv4"

echo "FORCING and OUTPUT variable file list:"
echo "$VARLIST"

for VARNAME in $VARLIST; do
  typeset FILELIST=""
  for LABEL in $LABELLIST; do
    FILELIST+=$LABEL"_"$VARNAME".gif "
  done
  convert $FILELIST +append $PREFIX"_"$VARNAME.gif
done

echo ""

}

myname=$(basename $0) #name of this script
fscript=ferret_LIS_HIST.jnl

ferret -version >/dev/null 2>&1 || \
  { echo "ERROR: Please load module ferret"; exit 1; }

if [ ! -e $fscript ]; then
  echo "ERROR: Ferret script file not found. [$fscript]"
  exit 1
fi

# Generate default domain and datetime lists
flist=`find . -name "LIS_HIST_*\.d*\.nc"`
if [ "$flist" == "" ]; then
  echo "ERROR: No LIS_HIST files found."
  exit 1
fi
flist=${flist//.\/LIS_HIST_/}
flist=${flist//.nc/}
dlist=""
maxtime="000000000000"
mintime="999912312359"
for item in $flist; do
 dlist+=${item/*./}" "
 if [ ${item/.*/} -gt $maxtime ]; then 
   maxtime=${item/.*/}; fi
 if [ ${item/.*/} -lt $mintime ]; then 
   mintime=${item/.*/}; fi
done
dlist=$(echo $dlist | tr ' ' '\n' | sort -u)
tlist=$mintime" "$maxtime

# Parse options
while getopts "hd:t:" opt
do
  case $opt in
    h ) usage ; exit 0 ;;
    d ) dlist=$OPTARG ;;
    t ) tlist=$OPTARG ;;
    \?) usage ; exit 1 ;;
  esac
done
shift $(($OPTIND - 1))

# Replace commas with spaces
dlist=${dlist//,/ }
tlist=${tlist//,/ }

# Generate plots then combine plots from same domains
for domain in $dlist; do
flist=""
plist=""
for time in $tlist; do
 flist+=" LIS_HIST_$time.$domain.nc"
 plist+=" plot_$time.$domain"
done
ferret -gif -script $fscript $flist
combine_noah33_plots seq_plot.$domain $plist
done

# Create archive
tar -czf seq_plots.tar.gz seq_plot*.gif

exit

