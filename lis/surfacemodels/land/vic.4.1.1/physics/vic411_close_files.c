#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vic411_vicNl.h>
 
static char vcid[] = "$Id: vic411_close_files.c,v 5.8 2006/10/18 20:58:57 vicadmin Exp $";

void vic411_close_files(vic411_filep_struct         *filep,
                 vic411_out_data_file_struct *out_data_files,
                 vic411_filenames_struct     *fnames)
/**********************************************************************
	vic411_close_files	Dag Lohmann		January 1996

  This routine closes all forcing data files, and output files.

  Modifications:
  7-19-96  Files are now gzipped when they are closed.  This
	   was added to save space when using large volumes
	   of data.						KAC
  02-27-01 Now closes files opened for lake model applications  KAC
  11-18-02 Now closes lake debugging file.                      LCB
  29-Oct-03 Distinguishing between input lakeparam file and output
	    lake file.						TJB
  2005-Mar-24 Added support for ALMA output files.		TJB
  2005-Apr-10 Added logic for OUTPUT_FORCE option.		TJB
  2006-Sep-23 Implemented flexible output configuration; uses new
	      out_data_files structure. TJB
  2006-Oct-16 Merged infiles and outfiles structs into vic411_filep_struct. TJB

**********************************************************************/
{
  extern vic411_option_struct vic411_options;
#if LINK_DEBUG
  extern vic411_debug_struct debug;
#endif
  int filenum;

  /**********************
    Close All Input Files
    **********************/

  fclose(filep->forcing[0]);
  if(vic411_options.COMPRESS) vic411_compress_files(fnames->forcing[0]);
  if(filep->forcing[1]!=NULL) {
    fclose(filep->forcing[1]);
    if(vic411_options.COMPRESS) vic411_compress_files(fnames->forcing[1]);
  }

  /*******************
    Close Output Files
    *******************/
/*
  for (filenum=0; filenum<vic411_options.Noutfiles; filenum++) {
    fclose(out_data_files[filenum].fh);
    if(vic411_options.COMPRESS) vic411_compress_files(out_data_files[filenum].filename);
  }
*/
#if !OUTPUT_FORCE

  /*******************************
    Close All Used Debugging Files
    *******************************/ 

#if LINK_DEBUG
  if(debug.DEBUG || debug.PRT_TEMP) {
    fclose(debug.fg_temp);
  }
  if(debug.DEBUG || debug.PRT_MOIST) {
    fclose(debug.fg_moist);
  }
  if(debug.DEBUG || debug.PRT_KAPPA) {
    fclose(debug.fg_kappa);
  }
  if(debug.DEBUG || debug.PRT_LAKE) {
    fclose(debug.fg_lake);
  }
  if(debug.DEBUG || debug.PRT_BALANCE) {
    fclose(debug.fg_balance);
  }
  if(debug.DEBUG || debug.PRT_FLUX) {
    fclose(debug.fg_energy);
  }
  if(debug.DEBUG || debug.PRT_SNOW) {
    fclose(debug.fg_snow);
  }
  if(debug.DEBUG || debug.PRT_GRID) {
    fclose(debug.fg_grid);
  }
#endif /* LINK_DEBUG */
#endif /* !OUTPUT_FORCE */

}
