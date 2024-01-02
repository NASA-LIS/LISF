#!/bin/sh 

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

#-----------------------------------------------------------------------

ls_offline=`find ./src/offline/ -name *.txt -print -o -name *.F90 -print`
if [ -n "$ls_offline" ]; then
perl protex -b $ls_offline > doc/offline.tex
else
touch doc/offline.tex
fi

#-----------------------------------------------------------------------

ls_core1=`find ./src/core/ -name *.txt -print -o -name *.F90 -print`
echo $ls_core1
if [ -n "$ls_core1" ]; then
perl protex -b $ls_core1 > doc/core1.tex
else
touch doc/core1.tex
fi
ls_core2=`find ./src/core/ -name *.c -print`
echo $ls_core2
if [ -n "$ls_core2" ]; then
perl protex -C -b $ls_core2 > doc/core2.tex
else
touch doc/core2.tex
fi

#-----------------------------------------------------------------------

ls_plugins=`find ./src/plugins/ -name *.txt -print -o -name *.F90 -print`
echo $ls_plugins
if [ -n "$ls_plugins" ]; then
perl protex -b $ls_plugins > doc/plugins.tex
else
touch doc/plugins.tex
fi

#-----------------------------------------------------------------------

ls_interp=`find ./src/interp/ -name *.txt -print -o -name *.F90 -print`
echo $ls_interp
if [ -n "$ls_interp" ]; then
perl protex -b $ls_interp > doc/interp.tex
else
touch doc/interp.tex
fi

#-----------------------------------------------------------------------

ls_runmodes=`find ./src/runmodes/ -name *.txt -print -o -name *.F90 -print`
echo $ls_runmodes
if [ -n "$ls_runmodes" ]; then
perl protex -b $ls_runmodes > doc/runmodes.tex
else
touch doc/runmodes.tex
fi

#-----------------------------------------------------------------------

ls_domains=`find ./src/domains/ -name *.txt -print -o -name *.F90 -print`
echo $ls_domains
if [ -n "$ls_domains" ]; then
perl protex -b $ls_domains > doc/domains.tex
else
touch doc/domains.tex
fi

#-----------------------------------------------------------------------

ls_params=`find ./src/params/ -name *.txt -print -o -name *.F90 -print`
echo $ls_params
if [ -n "$ls_params" ]; then
perl protex -b $ls_params > doc/params.tex
else
touch doc/params.tex
fi

#-----------------------------------------------------------------------

ls_metforcing=`find ./src/metforcing/ -name rdhm356 -prune -o -name *.txt -print -o -name *.F90 -print`
echo $ls_metforcing
if [ -n "$ls_metforcing" ]; then
perl protex -b $ls_metforcing > doc/metforcing.tex
else
touch doc/metforcing.tex
fi

ls_metforcing=`find ./src/metforcing/rdhm356 -name *.txt -print -o -name *.F90 -print`
echo $ls_metforcing
if [ -n "$ls_metforcing" ]; then
perl protex -b $ls_metforcing >> doc/metforcing.tex
fi

ls_metforcing=`find ./src/metforcing/rdhm356 -name *.c -print`
echo $ls_metforcing
if [ -n "$ls_metforcing" ]; then
perl protex -C -b $ls_metforcing >> doc/metforcing.tex
fi

#-----------------------------------------------------------------------

ls_lsms=`find ./src/surfacemodels/land/cable -name *.txt -print -o -name *.F90 -print;      \
         find ./src/surfacemodels/land/clm2 -name *.txt -print -o -name *.F90 -print;       \
         find ./src/surfacemodels/land/clsm.f2.5 -name *.txt -print -o -name *.F90 -print;  \
         find ./src/surfacemodels/land/geowrsi.2 -name *.txt -print -o -name *.F90 -print;  \
         find ./src/surfacemodels/land/hyssib -name *.txt -print -o -name *.F90 -print;     \
         find ./src/surfacemodels/land/mosaic -name *.txt -print -o -name *.F90 -print;     \
         find ./src/surfacemodels/land/noah.2.7.1 -name *.txt -print -o -name *.F90 -print; \
         find ./src/surfacemodels/land/noah.3.2 -name *.txt -print -o -name *.F90 -print;   \
         find ./src/surfacemodels/land/noah.3.3 -name *.txt -print -o -name *.F90 -print;   \
         find ./src/surfacemodels/land/noah.3.6 -name *.txt -print -o -name *.F90 -print;   \
         find ./src/surfacemodels/land/noah.3.9 -name *.txt -print -o -name *.F90 -print;   \
         find ./src/surfacemodels/land/noahmp.3.6 -name *.txt -print -o -name *.F90 -print;   \
         find ./src/surfacemodels/land/noahmp.4.0.1 -name *.txt -print -o -name *.F90 -print;   \
         find ./src/surfacemodels/land/rdhm.3.5.6 -name *.txt -print -o -name *.F90 -print; \
         find ./src/surfacemodels/land/sib2 -name *.txt -print -o -name *.F90 -print;       \
         find ./src/surfacemodels/land/template -name *.txt -print -o -name *.F90 -print;`

echo $ls_lsms
if [ -n "$ls_lsms" ]; then
perl protex -b $ls_lsms > doc/lsms.tex
else
touch doc/lsms.tex
fi

ls_lsms=`find ./src/surfacemodels/land/vic.4.1.1 -name *.txt -print -o -name *.F90 -print;`
echo $ls_lsms
if [ -n "$ls_lsms" ]; then
perl protex -b $ls_lsms >> doc/lsms.tex
fi

ls_lsms=`find ./src/surfacemodels/land/vic.4.1.1 -name *.c -print`
echo $ls_lsms
if [ -n "$ls_lsms" ]; then
perl protex -C -b $ls_lsms >> doc/lsms.tex
fi

ls_lsms=`find ./src/surfacemodels/land/vic.4.1.2.l -name *.txt -print -o -name *.F90 -print;`
echo $ls_lsms
if [ -n "$ls_lsms" ]; then
perl protex -b $ls_lsms >> doc/lsms.tex
fi

ls_lsms=`find ./src/surfacemodels/land/vic.4.1.2.l -name *.c -print`
echo $ls_lsms
if [ -n "$ls_lsms" ]; then
perl protex -C -b $ls_lsms >> doc/lsms.tex
fi

#-----------------------------------------------------------------------

ls_daalgs=`find ./src/dataassim/algorithm/ -name *.txt -print -o -name *.F90 -print`
echo $ls_daalgs
if [ -n "$ls_daalgs" ]; then
perl protex -b $ls_daalgs > doc/daalgs.tex
else
touch doc/daalgs.tex
fi

#-----------------------------------------------------------------------

ls_daperts=`find ./src/dataassim/perturb/ -name *.txt -print -o -name *.F90 -print`
echo $ls_daperts
if [ -n "$ls_daperts" ]; then
perl protex -b $ls_daperts > doc/daperts.tex
else
touch doc/daperts.tex
fi

#-----------------------------------------------------------------------

ls_dabiasest=`find ./src/dataassim/biasEstimation/ -name *.txt -print -o -name *.F90 -print`
echo $ls_dabiasest
if [ -n "$ls_dabiasest" ]; then
perl protex -b $ls_dabiasest > doc/dabiasest.tex
else
touch doc/dabiasest.tex
fi

ls_daobs=`find ./src/dataassim/obs/ -name *.txt -print -o -name *.F90 -print`
echo $ls_daobs
if [ -n "$ls_daobs" ]; then
perl protex -b $ls_daobs > doc/daobs.tex
else
touch doc/daobs.tex
fi

#-----------------------------------------------------------------------

ls_pealgs=`find ./src/optUE/algorithm/ -name *.txt -print -o -name *.F90 -print`
echo $ls_pealgs
if [ -n "$ls_pealgs" ]; then
perl protex -b $ls_pealgs > doc/pealgs.tex
else
touch doc/pealgs.tex
fi

#-----------------------------------------------------------------------

ls_petypes=`find ./src/optUE/type/ -name *.txt -print -o -name *.F90 -print`
echo $ls_petypes
if [ -n "$ls_petypes" ]; then
perl protex -b $ls_petypes > doc/petypes.tex
else
touch doc/petypes.tex
fi

#-----------------------------------------------------------------------

ls_rtms=`find ./src/rtms/ -name *.txt -print -o -name *.F90 -print`
echo $ls_rtms
if [ -n "$ls_rtms" ]; then
perl protex -b $ls_rtms > doc/rtms.tex
else
touch doc/rtms.tex
fi

#-----------------------------------------------------------------------
