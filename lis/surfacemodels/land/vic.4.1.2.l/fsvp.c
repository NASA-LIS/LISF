//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.3
//
// Copyright (c) 2020 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "ftn.h"
double svp(double);
//BOP
//
// !ROUTINE: fsvp
// \label{fsvp}
// 
// !REVISION HISTORY: 
//  12 Dec 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
double FTN(fsvp)(double * temp)
//
// !DESCRIPTION: 
//  This routine returns the saturated vapor pressure (svp) corresponding to
//  the given temperature by calling VIC's svp routine.
//
//  The arguments are:
//  \begin{description}
//  \item[temp] temperature, in Celcius
//  \end{description}
//
//EOP
{
   return svp(*temp);
}
