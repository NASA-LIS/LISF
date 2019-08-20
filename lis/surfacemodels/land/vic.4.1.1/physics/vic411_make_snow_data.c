#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>
 
static char vcid[] = "$Id: vic411_make_snow_data.c,v 3.1 1999/02/16 18:02:07 vicadmin Exp $";

vic411_snow_data_struct **vic411_make_snow_data(int nveg)
/**********************************************************************
	vic411_make_snow_data	Keith Cherkauer		January 22, 1997

  This routine makes an array of snow cover data structures, one 
  for each vegetation type plus bare soil.

  modifications:
  07-09-98 modified to make te make a two dimensional array which 
           also accounts for a variable number of snow elevation
           bands                                               KAC

**********************************************************************/
{
  extern vic411_option_struct vic411_options;

  int                i;
  vic411_snow_data_struct **temp;

  temp = (vic411_snow_data_struct **) calloc(nveg, 
				      sizeof(vic411_snow_data_struct *));

  for(i=0;i<nveg;i++) {
    temp[i] = (vic411_snow_data_struct *) calloc(vic411_options.SNOW_BAND, 
					  sizeof(vic411_snow_data_struct));
  }
    
  return temp;
}
