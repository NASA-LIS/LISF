#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>
 
static char vcid[] = "$Id: vic411_write_atmosdata.c,v 4.2 2004/05/11 20:35:10 tbohn Exp $";

void vic411_write_atmosdata(vic411_atmos_data_struct *atmos, int nrecs)
/**********************************************************************
	vic411_write_atmosdata		Dag Lohmann	Januray 1996

  This routine writes atmospheric data to the screen.

  Modifications:
    28-Aug-99 Changed to reflect the new vic411_atmos_data_struct.     Bart Nijssen
    07-May-04 No longer close the debug file, since the next cell
	      must write to it.					TJB
**********************************************************************/
{
#if LINK_DEBUG
  extern vic411_debug_struct debug;

  int i;
  int j;

  /*  first write all the SNOW_STEP data  - only write if the modelstep !=
      SNOWSTEP */
  if (vic411_NR > 0) {
    for (i = 0; i < nrecs; i++) {
      for (j = 0; j < vic411_NF; j++) {
	fprintf(debug.fg_snowstep_atmos,"%d\t%d",  i, j);
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].prec[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].air_temp[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].wind[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].vpd[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].vp[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].pressure[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].density[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].shortwave[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].longwave[j]); 
	fprintf(debug.fg_snowstep_atmos,"\n");
      }
    }
  /* don't close the debug output file, as we need to write to it for the next cell as well */
/*  fclose(debug.fg_snowstep_atmos);*/
  }
  
  /* then write all the dt data */
  for (i = 0; i < nrecs; i++) {
    fprintf(debug.fg_modelstep_atmos,"%d",  i);
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].prec[vic411_NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].air_temp[vic411_NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].wind[vic411_NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].vpd[vic411_NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].vp[vic411_NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].pressure[vic411_NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].density[vic411_NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].shortwave[vic411_NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].longwave[vic411_NR]); 
    fprintf(debug.fg_modelstep_atmos,"\n");
  }
  fclose(debug.fg_modelstep_atmos);
#endif

}


