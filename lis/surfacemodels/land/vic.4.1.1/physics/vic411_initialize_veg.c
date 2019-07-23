#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_initialize_veg.c,v 4.2 2003/04/22 18:46:33 cherkaue Exp $";

void vic411_initialize_veg(vic411_veg_var_struct      **veg_var,
		    vic411_veg_con_struct       *veg_con,
		    vic411_global_param_struct   *gp,
		    int                    Nveg)
/**********************************************************************
  vic411_initialize_veg		Dag Lohmann	 January 1996

  This routine initailizes the vegetation variable array.

  Modifications:
  07-13-98 modified to initialize vegetation structure for all 
           defined elevation bands                                 KAC
  11-18-02 modified to get the maximum number of vegetation types
           passed to it.  This allows the maximum number of vegetation
           types to include the wetland vegetation fraction when the 
           lake model is active.                                  LCB

**********************************************************************/
{
  extern vic411_option_struct   vic411_options;

  int i, j;

  for ( i = 0 ; i < Nveg ; i++) {
    for ( j = 0 ; j < vic411_options.SNOW_BAND ; j++ ) {
      veg_var[i][j].Wdew = 0.0;
      veg_var[i][j].throughfall = 0.0;
    }
  }
}
