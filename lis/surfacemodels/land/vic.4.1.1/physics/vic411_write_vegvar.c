#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_write_vegvar.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

void vic411_write_vegvar(vic411_veg_var_struct *veg, 
		  int             n)
/**********************************************************************
  vic411_write_vegvar		Keith Cherkauer		May 29, 1996

  This routine writes vegetation variables to stdout.  Used primarily
  for debugging purposes.

  Modifications:
  5/21/96	Routine was modified to allow for variable
		number of layers				KAC

**********************************************************************/
{
  printf("Vegetation Variables: vegtype %i\n",n);
  printf("\tcanopyevap  = %f\n", veg->canopyevap);
  printf("\tWdew        = %f\n", veg->Wdew);
  printf("\tthroughfall = %f\n", veg->throughfall);
}

