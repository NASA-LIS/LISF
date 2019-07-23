#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_vicerror.c,v 4.3 2006/10/18 20:59:01 vicadmin Exp $";

void vic411_vicerror(char error_text[])
/**********************************************************************
	vic411_vicerror.c	Keith Cherkauer		April 23, 1997

  This subroutine was written to handle numerical errors within the
  VIC model.  This will flush all file buffers so that all records 
  that have been run will be written to disk before the model is exited.

  Modifications:
  2006-Sep-23 Implemented flexible output configuration; uses the new
              out_data and out_data_files structures. TJB
  2006-Oct-16 Merged infiles and outfiles structs into vic411_filep_struct. TJB

**********************************************************************/
{
        extern vic411_option_struct vic411_options;
	extern vic411_Error_struct vic411_Error;
#if LINK_DEBUG
        extern vic411_debug_struct debug;
#endif

        vic411_filenames_struct fnames;
	void _exit();

        vic411_options.COMPRESS=FALSE;	/* turn off compression of last set of files */

	fprintf(stderr,"VIC model run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now writing output files...\n");
        vic411_close_files(&(vic411_Error.filep), vic411_Error.out_data_files, &fnames);
	fprintf(stderr,"...now exiting to system...\n");
        fflush(stdout);
        fflush(stderr);
	_exit(1);
}
