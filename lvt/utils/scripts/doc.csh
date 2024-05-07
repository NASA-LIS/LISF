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

ls_main=`find ./src/main/ -name *.F90 -print`
perl protex -b $ls_main > doc/main.tex

ls_core1=`find ./src/core/ -name *.F90 -print`
echo $ls_core1
perl protex -b $ls_core1 > doc/core1.tex
ls_core2=`find ./src/core/ -name *.c -print`
echo $ls_core2
perl protex -C -b $ls_core2 > doc/core2.tex

ls_plugins=`find ./src/plugins/ -name *.F90 -print`
echo $ls_plugins
perl protex -b $ls_plugins > doc/plugins.tex

ls_interp=`find ./src/interp/ -name *.F90 -print`
echo $ls_interp
perl protex -b $ls_interp > doc/interp.tex

ls_runmodes=`find ./src/runmodes/ -name *.F90 -print`
echo $ls_runmodes
perl protex -b $ls_runmodes > doc/runmodes.tex

ls_obs=`find ./src/obs/ -name *.F90 -print`
echo $ls_obs
perl protex -b $ls_obs > doc/obs.tex

ls_domains=`find ./src/domains/ -name *.F90 -print`
echo $ls_domains
perl protex -b $ls_domains > doc/domains.tex

ls_params=`find ./src/params/ -name *.F90 -print`
echo $ls_params
perl protex -b $ls_params > doc/params.tex



