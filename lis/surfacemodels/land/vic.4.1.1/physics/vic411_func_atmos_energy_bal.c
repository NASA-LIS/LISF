#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_func_atmos_energy_bal.c,v 5.1 2001/08/15 22:42:20 cherkaue Exp $";

double vic411_func_atmos_energy_bal(double Tcanopy, va_list ap) {
/**********************************************************************
  vic411_func_atmos_energy_bal.c      Keith Cherkauer        February 6, 2001

  This routine solves the atmospheric exchange energy balance.

**********************************************************************/

  double  LatentHeat;
  double  NetRadiation;
  double  Ra;
  double  Tair;
  double  atmos_density;
  double  InSensible;

  double *SensibleHeat;
 
  // internal routine variables
  double  vic411_Error;

  // extract variables from va_arg
  LatentHeat    = (double)  va_arg(ap, double);
  NetRadiation  = (double)  va_arg(ap, double);
  Ra            = (double)  va_arg(ap, double);
  Tair          = (double)  va_arg(ap, double);
  atmos_density = (double)  va_arg(ap, double);
  InSensible    = (double)  va_arg(ap, double);

  SensibleHeat  = (double *)va_arg(ap, double *);

  // compute sensible heat flux between canopy and atmosphere
  (*SensibleHeat) = atmos_density * Cp * (Tair - Tcanopy) / Ra;

  // compute energy balance error
  //vic411_Error = NetRadiation + LatentHeat + (*SensibleHeat);
  vic411_Error = InSensible - (*SensibleHeat);

  return ( vic411_Error );

}
