#include <stdio.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_calc_cloud_cover_fraction.c,v 4.2 2004/05/06 22:29:38 tbohn Exp $";

void vic411_calc_cloud_cover_fraction(vic411_atmos_data_struct *atmos,
			       vic411_dmy_struct        *dmy,
			       int                nrecs,
			       int                Ndays,
			       int                stepspday,
			       double            *tskc) {
/********************************************************************
  vic411_calc_cloud_cover_fraction     Keith Cherkauer    January 12, 2000

  This routine is designed to estimate cloud cover when observations
  of shortwave radiation are available.

*********************************************************************/

  vic411_nrerror("The function to estimate cloud cover from observed solar radiation does not yeat work.");

}
