#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>
#include <string.h>
 
static char vcid[] = "$Id: vic411_read_forcing_data.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

double **vic411_read_forcing_data(FILE                **infile,
			   vic411_global_param_struct   vic411_global_param)
/**********************************************************************
  vic411_read_forcing_data    Keith Cherkauer      January 10, 2000

  This subroutine controls the order and number of forcing variables
  read from the forcing data files.  Two forcing files are allowed, 
  variables, time step and file format must be defined in the global
  control file.

**********************************************************************/
{
  extern vic411_option_struct    vic411_options;
  extern vic411_param_set_struct vic411_param_set;
  extern int              vic411_NR, vic411_NF;

  char                 errorstr[MAXSTRING];
  int                  i;
  double             **forcing_data;

  /** Allocate data arrays for input forcing data **/
  forcing_data = (double **)calloc(N_FORCING_TYPES,sizeof(double*));
  for(i=0;i<N_FORCING_TYPES;i++) 
    if (vic411_param_set.TYPE[i].SUPPLIED) 
      forcing_data[i] = (double *)calloc((vic411_global_param.nrecs * vic411_NF),
			   sizeof(double));

  /** Read First Forcing Data File **/
  if(vic411_param_set.FORCE_DT[0] > 0) {
    vic411_read_atmos_data(infile[0], vic411_global_param, 0, vic411_global_param.forceskip[0],
		    forcing_data);
  }
  else {
    sprintf(errorstr,"ERROR: File time step must be defined for at least the first forcing file (FILE_DT).\n");
    vic411_vicerror(errorstr);
  }

  /** Read Second Forcing Data File **/
  if(vic411_param_set.FORCE_DT[1] > 0) {
    vic411_read_atmos_data(infile[1], vic411_global_param, 1, vic411_global_param.forceskip[1], 
		    forcing_data);
  }

  return(forcing_data);

}
