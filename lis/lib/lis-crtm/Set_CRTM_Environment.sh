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

# Set current directory as new CRTM root.
export CRTM_ROOT=${PWD}

export PATH=$PATH:$CRTM_ROOT/scripts/shell/Utility

# Construct the new CRTM environment variables
export CRTM_SOURCE_ROOT="${CRTM_ROOT}/src"
export CRTM_FIXFILE_ROOT="${CRTM_ROOT}/fix"
export CRTM_TEST_ROOT="${CRTM_ROOT}/test"
export CRTM_EXTERNALS_ROOT="${CRTM_ROOT}/externals"
export CRTM_SCRIPTS_ROOT="${CRTM_ROOT}/scripts"
export CRTM_DOC_ROOT="${CRTM_ROOT}/doc"
export CRTM_VALIDATION_ROOT="${CRTM_ROOT}/validation"

alias crtmsrc="cd $CRTM_SOURCE_ROOT"
alias crtmfix="cd $CRTM_FIXFILE_ROOT"
alias crtmtest="cd $CRTM_TEST_ROOT"
alias crtmscripts="cd $CRTM_SCRIPTS_ROOT"
alias crtmext="cd $CRTM_EXTERNALS_ROOT"
alias crtmdoc="cd $CRTM_DOC_ROOT"
alias crtmval="cd $CRTM_VALIDATION_ROOT"

echo "All CRTM environment variables now rooted at ${CRTM_ROOT}"
