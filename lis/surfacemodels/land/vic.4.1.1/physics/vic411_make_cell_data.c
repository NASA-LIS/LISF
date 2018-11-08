#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>
 
static char vcid[] = "$Id: vic411_make_cell_data.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

vic411_cell_data_struct **vic411_make_cell_data(int veg_type_num, int Nlayer)
/**********************************************************************
	vic411_make_cell_data	Keith Cherkauer		July 9, 1997

  This subroutine makes an array of type cell, which contains soil
  column variables for a single grid cell.

**********************************************************************/
{
  extern vic411_option_struct vic411_options;

  int i;
  vic411_cell_data_struct **temp;

  temp = (vic411_cell_data_struct**) calloc(veg_type_num, 
                                  sizeof(vic411_cell_data_struct*));
  for(i=0;i<veg_type_num;i++) {
    temp[i] = (vic411_cell_data_struct*) calloc(vic411_options.SNOW_BAND, 
					 sizeof(vic411_cell_data_struct));
/*     for(j=0;j<vic411_options.SNOW_BAND;j++) { */
/*       temp[i][j].layer  */
/* 	= (vic411_layer_data_struct*)calloc(Nlayer, */
/* 				     sizeof(vic411_layer_data_struct)); */
      
/*     } */
  }
  return temp;
}
