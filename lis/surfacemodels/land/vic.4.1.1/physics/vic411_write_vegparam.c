#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_write_vegparam.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

void vic411_write_vegparam(vic411_veg_con_struct *veg_con)
/**********************************************************************
	vic411_write_vegparam		Dag Lohmann	January 1996

  This routine writes vegetation parameters to stdout, used primarily 
  for debugging, and making sure the model is reading the proper 
  parameters..

  Modifications:
  5/21/96	Routine was modified to allow for variable
		number of layers				KAC
  4-12-98  Updated for new standard vegetation parameters       KAC

**********************************************************************/
{
  extern vic411_veg_lib_struct *vic411_veg_lib;
  extern vic411_option_struct vic411_options;

  int i, j, l;
  int vegclass;

  printf("Vegetation Parameters:\n");
  printf("\tvegetat_type_num = %d\n",  veg_con[0].vegetat_type_num);

  for (i = 0; i < veg_con[0].vegetat_type_num; i++) {
    vegclass = veg_con[i].veg_class;
    printf("\n\tveg_class            = %d\n",  vic411_veg_lib[vegclass].veg_class);
    printf("\tCv                   = %f\n", veg_con[i].Cv);
    if(vic411_veg_lib[vegclass].overstory)
      printf("\tOverstory            = TRUE\n");
    else 
      printf("\tOverstory            = FALSE\n");
    printf("\trarc                 = %f s/m\n", vic411_veg_lib[vegclass].rarc);
    printf("\trmin                 = %f s/m\n", vic411_veg_lib[vegclass].rmin);
    for(l=0;l<vic411_options.ROOT_ZONES;l++)
      printf("\tzone_depth _fract%d   = %f %f\n",l+1,
	     veg_con[i].zone_depth[l],veg_con[i].zone_fract[l]);
    for(l=0;l<vic411_options.Nlayer;l++)
      printf("\troot_percent%d        = %f\n",l+1,veg_con[i].root[l]);
    for (j = 0; j < 12; j++) 
      printf("\tLAI[%02d]             = %f\n",j,vic411_veg_lib[vegclass].LAI[j]);
    for (j = 0; j < 12; j++) 
      printf("\talbedo[%02d]          = %f\n",j,vic411_veg_lib[vegclass].albedo[j]);
    for (j = 0; j < 12; j++) 
      printf("\tdisplacement[%02d]    = %f m\n",j,
	     vic411_veg_lib[vegclass].displacement[j]);
    for (j = 0; j < 12; j++) 
      printf("\troughness[%02d]       = %f m\n",j,vic411_veg_lib[vegclass].roughness[j]);
    printf("\twind_h                  = %f s/m\n", vic411_veg_lib[vegclass].wind_h);
  }
}

