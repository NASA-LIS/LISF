//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "ftn.h"
double vic411_svp(double);
//BOP
//
// !ROUTINE: fsvp
// \label{fsvp411}
// 
// !REVISION HISTORY: 
//  12 Dec 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
double FTN(vic411_fsvp)(double * temp)
//
// !DESCRIPTION: 
//  This routine returns the saturated vapor pressure (vic411\_svp) corresponding to
//  the given temperature by calling VIC's vic411\_svp routine.
//
//  The arguments are:
//  \begin{description}
//  \item[temp] temperature, in Celcius
//  \end{description}
//
//EOP
{
   return vic411_svp(*temp);
}
