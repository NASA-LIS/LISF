#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_open_state_file.c,v 4.10 2006/10/18 20:58:59 vicadmin Exp $";


FILE *vic411_open_state_file(vic411_global_param_struct *global,
		      vic411_filenames_struct     filenames,
		      int                  Nlayer,
		      int                  Nnodes) 
/*********************************************************************
  vic411_open_state_file      Keith Cherkauer           April 15, 2000

  This subroutine opens the model state file for output.

  Modifications:
  04-10-03 Modified to open and write to a binary state file.    KAC
  06-03-03 modified to handle both ASCII and BINARY state files.  KAC
  2005-11-29 SAVE_STATE is set in global param file, not in vic411_user_def.h GCT
  2005-12-06 Moved setting of statename to vic411_get_global_param     GCT
  2006-08-23 Changed order of fread/fwrite statements from ...1, sizeof...
             to ...sizeof, 1,... GCT
  2006-Oct-16 Merged infiles and outfiles structs into vic411_filep_struct;
	      This included moving global->statename to filenames->statefile. TJB

*********************************************************************/
{
  extern vic411_option_struct vic411_options;

  FILE   *statefile;
  char    filename[MAXSTRING];
  double  Nsum;

  /* open state file */
  sprintf(filename,"%s", filenames.statefile);
  if ( vic411_options.BINARY_STATE_FILE )
    statefile = vic411_open_file(filename,"wb");
  else
    statefile = vic411_open_file(filename,"w");

  /* Write save state date information */
  if ( vic411_options.BINARY_STATE_FILE ) {
    fwrite( &global->stateyear, sizeof(int), 1, statefile );
    fwrite( &global->statemonth, sizeof(int), 1, statefile );
    fwrite( &global->stateday, sizeof(int), 1, statefile );
  }
  else {
    fprintf(statefile,"%i %i %i\n", global->stateyear, 
	    global->statemonth, global->stateday);
  }

  /* Write simulation flags */
  if ( vic411_options.BINARY_STATE_FILE ) {
    fwrite( &Nlayer, sizeof(int), 1, statefile );
    fwrite( &Nnodes, sizeof(int), 1, statefile );
  }
  else {
    fprintf(statefile,"%i %i\n", Nlayer, Nnodes);
  }

  return(statefile);

}

