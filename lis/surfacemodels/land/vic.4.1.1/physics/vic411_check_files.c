#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_check_files.c,v 5.4.2.1 2009/01/24 01:15:42 vicadmin Exp $";

void vic411_check_files(vic411_filep_struct     *filep, 
		 vic411_filenames_struct *fnames)
/**********************************************************************
	vic411_check_files		Dag Lohmann		January 1996

  This routine opens files for soil, vegetation, and global parameters.

  Modifcations:
  02-27-01 Added controls for lake model parameter file    KAC
  2005-Apr-13 Added logic for OUTPUT_FORCE option.			TJB
  2006-Oct-16 Merged infiles and outfiles structs into vic411_filep_struct.	TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB

**********************************************************************/
{
  extern vic411_option_struct  vic411_options;
  extern FILE          *vic411_open_file(char string[], char type[]);

  filep->soilparam   = vic411_open_file(fnames->soil, "r");
#if !OUTPUT_FORCE
  filep->veglib      = vic411_open_file(fnames->veglib, "r");
  filep->vegparam    = vic411_open_file(fnames->veg, "r");
  if(vic411_options.SNOW_BAND>1)
    filep->snowband    = vic411_open_file(fnames->snowband, "r");
  if ( vic411_options.LAKES )
    filep->lakeparam = vic411_open_file(fnames->lakeparam,"r");
#endif /* !OUTPUT_FORCE */

}


