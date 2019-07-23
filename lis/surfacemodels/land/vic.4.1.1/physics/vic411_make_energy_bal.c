#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>
 
static char vcid[] = "$Id: vic411_make_energy_bal.c,v 4.2 2004/11/02 03:08:44 vicadmin Exp $";

vic411_energy_bal_struct **vic411_make_energy_bal(int nveg)
/**********************************************************************
	vic411_make_energy_bal	Keith Cherkauer		May 26, 1996

  This routine makes an array of frozen soil data structures, one 
  for each vegetation type and bare soil.

  Modifications:
  01-Nov-04 Removed modification of Nnodes, as this was preventing
	    correct reading/writing of state files for QUICK_FLUX
	    =TRUE.						TJB

**********************************************************************/
{
  extern vic411_option_struct vic411_options;

  int i, j;
  vic411_energy_bal_struct **temp;

  temp = (vic411_energy_bal_struct**) calloc(nveg, 
				      sizeof(vic411_energy_bal_struct*));

  /** Initialize all records to unfrozen conditions */
  for(i = 0; i < nveg; i++) {
    temp[i] = (vic411_energy_bal_struct*) calloc(vic411_options.SNOW_BAND, 
					  sizeof(vic411_energy_bal_struct));
    for(j = 0; j < vic411_options.SNOW_BAND; j++) {
      temp[i][j].frozen = FALSE;
    }
  }

  return temp;
}
