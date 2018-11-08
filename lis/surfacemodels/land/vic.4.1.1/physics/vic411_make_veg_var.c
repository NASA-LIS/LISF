#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>
 
static char vcid[] = "$Id: vic411_make_veg_var.c,v 3.1 1999/02/16 18:02:07 vicadmin Exp $";

vic411_veg_var_struct **vic411_make_veg_var(int veg_type_num)
/**********************************************************************
	vic411_make_veg_var	Dag Lohman		January 1996

  This routine makes an array of vegitation variables for each vegitation
  type.

  Modifications:
  07-13-98 modified to add structure definitions for all defined 
           elevation bands                                       KAC

**********************************************************************/
{
  extern vic411_option_struct vic411_options;
  
  int              i;
  vic411_veg_var_struct **temp;

  temp = (vic411_veg_var_struct **) calloc(veg_type_num, 
				    sizeof(vic411_veg_var_struct *));
  for(i=0;i<veg_type_num;i++)
    temp[i] = (vic411_veg_var_struct *) calloc(vic411_options.SNOW_BAND, 
					sizeof(vic411_veg_var_struct));
  return temp;
}
