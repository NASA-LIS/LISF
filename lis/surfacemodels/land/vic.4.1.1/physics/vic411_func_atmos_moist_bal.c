#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_func_atmos_moist_bal.c,v 5.1 2001/08/15 22:42:20 cherkaue Exp $";

double vic411_func_atmos_moist_bal(double VPcanopy, va_list ap) {
/**********************************************************************
  vic411_func_atmos_moist_bal.c      Keith Cherkauer        March 2, 2001

  This routine solves the atmospheric exchange moisture balance.

**********************************************************************/

  double  InLatentHeat;
  double  Lv;
  double  Ra;
  double  atmos_density;
  double  gamma;
  double  vp; // atmospheric vapor pressure

  double *LatentHeat;
 
  // internal routine variables
  double  vic411_Error;

  // extract variables from va_arg
  InLatentHeat  = (double)  va_arg(ap, double);
  Lv            = (double)  va_arg(ap, double);
  Ra            = (double)  va_arg(ap, double);
  atmos_density = (double)  va_arg(ap, double);
  gamma         = (double)  va_arg(ap, double);
  vp            = (double)  va_arg(ap, double);

  LatentHeat    = (double *)va_arg(ap, double *);

  // compute sensible heat flux between canopy and atmosphere
  (*LatentHeat) = Lv * atmos_density * Cp * (vp - VPcanopy) / ( gamma * Ra );

  // compute energy balance error
  vic411_Error = InLatentHeat - (*LatentHeat);

  return ( vic411_Error );

}
