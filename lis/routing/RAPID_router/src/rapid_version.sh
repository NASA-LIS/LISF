#!/bin/bash
#*******************************************************************************
#rapid_version.sh
#*******************************************************************************

#Purpose:
#This script allows determining the version of RAPID that is being used, but  
#only if git is installed and if the RAPID git repository is present.  Otherwise
#'unknown' is used.  
#Author:
#Cedric H. David, 2016-2020.


#*******************************************************************************
#Check if a program exists and perform tasks
#*******************************************************************************
if type 'git' > /dev/null; then
     #git is installed
     if git rev-parse --git-dir > /dev/null 2>&1; then
          #this is a git repository
          git describe 
     else
          #this is not a git repository
          echo "unknown, NOT a git repository"
     fi
else
     #git is not installed
     echo "unknown, git NOT installed"
fi


#*******************************************************************************
#end
#*******************************************************************************
