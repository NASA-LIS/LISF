//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <stdio.h>
#include <string.h>
#include "ftn.h"
//BOP
//
// !ROUTINE: lis_log_msgC.c
//
// !DESCRIPTION:
//   This file provides methods for diagnostic and error logging in LIS 
//   for codes written in the C language. The diagnostic messages are 
//   appended with the appropriate time-date and processor stamps. 
//EOP
void FTN(lis_log_msg)(char *, int);
void FTN(lis_log_blocked_msg)(char *, int);

//BOP
// !ROUTINE: lis_log_msgC
// \label{lis_log_msgC}
//
// !DESCRIPTION:
//  This routine is a C wrapper to the lis\_log\_msg routine.
//
// !REVISION HISTORY:
//  12 Mar 2004: James Geiger; Initial version
//
// !INTERFACE:
void lis_log_msgC(char * string)
//EOP
{
  printf("[log]: %s\n",string);
}

