#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_advected_sensible_heat.c,v 5.2 2004/05/06 22:29:38 tbohn Exp $";

double vic411_advected_sensible_heat(double SnowCoverFract, 
			      double AirDens, 
			      double Tair, 
			      double TGrnd, 
			      double Ra) {
/**************************************************************************
  vic411_advected_sensible_heat       Keith Cherkauer            October 6, 2000

  This routine computes the sensible heat flux advected from bare soil
  patches to snow covered patches.  Based on the work of Marsh, P.,
  J. W. Pomeroy, and N. Neumann, Sensible heat flux and local advection 
  over a heterogeneous landscape at an Arctic tundra site during snowmelt,
  Annals of Glaciology, 25, 1997. pp 132-136.

**************************************************************************/

  double Qbare;
  double Qadv;
  double Fs;
  double val;

  /* Compute sensible heat flux from bare patches */
  Qbare = AirDens * Cp * (Tair - TGrnd) / Ra;

  /* Compute fraction of sensible heat that contributes to snowmelt */
  if ( SnowCoverFract > 0.6 ) Fs = 1.0;
  else if ( SnowCoverFract > 0.5 ) {
    val = 3. * SnowCoverFract - 1.8;
    Fs = pow(10, val);
  }
  else if ( SnowCoverFract > 0.2 ) {
    val = 5.6667 * SnowCoverFract - 3.1333;
    Fs = pow(10, val);
  }
  else Fs = 0.01;

  /* Compute advected sensible heat flux */
  Qadv = ( Qbare * (1.0 - SnowCoverFract) ) / SnowCoverFract * Fs;

  return (Qadv);

}

