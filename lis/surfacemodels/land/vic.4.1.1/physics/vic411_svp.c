#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_svp.c,v 5.1 2001/08/15 23:44:57 cherkaue Exp $";

double vic411_svp(double temp)
/**********************************************************************
  This routine computes the saturated vapor pressure using Handbook
  of Hydrology eqn 4.2.2

  Pressure in Pa

**********************************************************************/
{
  double SVP;
  
  SVP = A_SVP * exp((B_SVP * temp)/(C_SVP+temp));

  if(temp<0) SVP *= 1.0 + .00972 * temp + .000042 * temp * temp;

  return (SVP*1000.);
}

double vic411_svp_slope(double temp)
/**********************************************************************
  This routine computes the gradient of d(vic411_svp)/dT using Handbook
  of Hydrology eqn 4.2.3

  returned value in Pa
**********************************************************************/
{
  return (B_SVP * C_SVP) / ((C_SVP + temp) * (C_SVP + temp)) * vic411_svp(temp);
}


 
