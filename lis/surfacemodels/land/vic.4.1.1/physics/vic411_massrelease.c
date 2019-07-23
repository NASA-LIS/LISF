/*
 * SUMMARY:      vic411_MassRelease.c - Calculates mass release of snow from canopy
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Brian Connelly and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:     6-Oct-1996 at 15:42:13
 * LAST-MOD: Mon Sep 28 16:21:39 1998 by VIC Administrator <vicadmin@u.washington.edu>
 * DESCRIPTION:  Calculates mass release of snow from canopy
 * DESCRIP-END.
 * FUNCTIONS:    vic411_MassRelease()
 * COMMENTS:     
 */

#include <stdio.h>
#include <vic411_vicNl.h>
  
static char vcid[] = "$Id: massrelease.c,v 3.1 1999/02/16 18:02:07 vicadmin Exp $";

/*****************************************************************************
  Function name: vic411_MassRelease()

  Purpose      : Calculates mass release of snow from canopy

  Required     :
    float *InterceptedSnow
    float *TempInterceptionStorage
    float *ReleasedMass
    float *Drip 

  Returns      : none

  Modifies     : see under required (i.e. all variables are modified)

  Comments     :
*****************************************************************************/
void vic411_MassRelease(double *InterceptedSnow, double *TempInterceptionStorage, 
                 double *ReleasedMass, double *Drip ) 
{
  double TempDrip;
  double TempReleasedMass;
  double Threshold;
  double MaxRelease;
  
  /* If the amount of snow in the canopy is greater than some minimum
     value, MIN_INTERCEPTION_STORAGE, then calculte mass release and Drip */
  
  if (*InterceptedSnow > MIN_INTERCEPTION_STORAGE) {
    Threshold  = 0.10 * *InterceptedSnow;
    MaxRelease = 0.17 * *InterceptedSnow;
    
    /* If the amount of vic411_snow_melt after interception, vic411_snow_melt, is >= the
       theshhold then there is mass release.  If vic411_snow_melt is < the treshhold
       then there is no mass release but that water remains in
       *TempInterceptionStorage which will be augmented during the next
       compute period */
    
    if ((*TempInterceptionStorage) >= Threshold) {
      
      
      *Drip += Threshold;
      *InterceptedSnow -= Threshold;
      *TempInterceptionStorage -= Threshold;
      if (*InterceptedSnow < MIN_INTERCEPTION_STORAGE)
        TempReleasedMass = 0.0;
      else
        TempReleasedMass = min((*InterceptedSnow - MIN_INTERCEPTION_STORAGE),
                               MaxRelease); 
      *ReleasedMass += TempReleasedMass;
      *InterceptedSnow -= TempReleasedMass;
      vic411_MassRelease(InterceptedSnow, TempInterceptionStorage, ReleasedMass,
                  Drip); 
    }

    else {
      TempDrip = min(*TempInterceptionStorage, *InterceptedSnow);
      *Drip += TempDrip;
      *InterceptedSnow -= TempDrip;
    }
  }

  /* (*InterceptedSnow < MIN_INTERCEPTION_STORAGE) If the amount of snow in
     the canopy is less than some minimum value, MIN_INTERCEPTION_STORAGE,
     then only melt can occur and there is no mass release. */

  else {
    TempDrip = min(*TempInterceptionStorage, *InterceptedSnow);
    *Drip += TempDrip;
    *InterceptedSnow -= TempDrip;
    *TempInterceptionStorage = 0.0;
  }
}
