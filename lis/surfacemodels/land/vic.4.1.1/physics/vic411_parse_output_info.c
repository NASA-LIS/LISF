#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>
#include <string.h>
 
static char vcid[] = "$Id: vic411_parse_output_info.c,v 5.2.2.5 2009/05/17 21:22:08 vicadmin Exp $";

void vic411_parse_output_info(vic411_filenames_struct      *names,
                       FILE                  *gp,
                       vic411_out_data_file_struct  **out_data_files,
                       vic411_out_data_struct       *out_data)
/**********************************************************************
  vic411_parse_output_info	Ted Bohn	            September 10 2006

  This routine reads the VIC model global control file, getting
  information for output variables list (if any).

  Modifications:
  2006-Nov-07 Changed default precision from %.1f to %.4f.	TJB
  2007-Jan-15 Modified to expect "OUT_TYPE_" at beginning of
	      output data type strings.				TJB
  2007-Apr-21 Added initialization for format, outfilenum, and
	      outvarnum.					TJB
  2008-Feb-15 Added check on number of output files defined vs.
	      N_OUTFILES.					TJB
  2009-Feb-09 Sets PRT_SNOW_BAND to FALSE if N_OUTFILES has been
	      specified.					TJB
  2009-Mar-15 Added default values for format, typestr, and
	      multstr, so that they can be omitted from global
	      param file.					TJB
**********************************************************************/
{
  extern vic411_option_struct    vic411_options;

  char cmdstr[MAXSTRING];
  char optstr[MAXSTRING];
  int  i;
  int  outfilenum;
  int  fn;
  char varname[20];
  int  outvarnum;
  char format[10];
  char typestr[20];
  int  type;
  char multstr[20];
  float mult;
  int  tmp_noutfiles;
  char ErrStr[MAXSTRING];

  strcpy(format,"*");

  /** Read through global control file to find output info **/

  fgets(cmdstr,MAXSTRING,gp);

  outfilenum = -1;
  outvarnum = 0;
  while(!feof(gp)) {
    if(cmdstr[0]!='#' && cmdstr[0]!='\n' && cmdstr[0]!='\0') {

      sscanf(cmdstr,"%s",optstr);

      if(strcasecmp("N_OUTFILES",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&tmp_noutfiles);
        vic411_free_out_data_files(out_data_files);
        vic411_options.Noutfiles = tmp_noutfiles;
        *out_data_files = (vic411_out_data_file_struct *)calloc(vic411_options.Noutfiles, sizeof(vic411_out_data_file_struct));
        outfilenum = -1;
        vic411_init_output_list(out_data, FALSE, "%.4f", OUT_TYPE_FLOAT, 1);
        // PRT_SNOW_BAND is ignored if N_OUTFILES has been specified
        vic411_options.PRT_SNOW_BAND = FALSE;
      }
      else if(strcasecmp("OUTFILE",optstr)==0) {
        outfilenum++;
        if (!vic411_options.Noutfiles) {
          vic411_nrerror("vic411_Error in global param file: \"N_OUTFILES\" must be specified before you can specify \"OUTFILE\".");
        }
        if (outfilenum >= vic411_options.Noutfiles) {
          sprintf(ErrStr, "vic411_Error in global param file: number of output files specified in N_OUTFILES (%d) is less than actual number of output files defined in the global param file.",vic411_options.Noutfiles);
          vic411_nrerror(ErrStr);
        }
        sscanf(cmdstr,"%*s %s %d",(*out_data_files)[outfilenum].prefix,&((*out_data_files)[outfilenum].nvars));
        (*out_data_files)[outfilenum].varid = (int *)calloc((*out_data_files)[outfilenum].nvars, sizeof(int));
        outvarnum = 0;
      }
      else if(strcasecmp("OUTVAR",optstr)==0) {
        if (outfilenum < 0) {
          vic411_nrerror("vic411_Error in global param file: \"OUTFILE\" must be specified before you can specify \"OUTVAR\".");
        }
        strcpy(format,"");
        strcpy(typestr,"");
        strcpy(multstr,"");
        sscanf(cmdstr,"%*s %s %s %s %s",varname, format, typestr, multstr);
        if (strcasecmp("",format) == 0) {
          strcpy(format,"*");
          type = OUT_TYPE_DEFAULT;
          mult = 0; // 0 means default multiplier
        }
        else {
          if (strcasecmp("OUT_TYPE_USINT", typestr)==0)
            type = OUT_TYPE_USINT;
          else if (strcasecmp("OUT_TYPE_SINT", typestr)==0)
            type = OUT_TYPE_SINT;
          else if (strcasecmp("OUT_TYPE_FLOAT", typestr)==0)
            type = OUT_TYPE_FLOAT;
          else if (strcasecmp("OUT_TYPE_DOUBLE", typestr)==0)
            type = OUT_TYPE_DOUBLE;
          else
            type = OUT_TYPE_DEFAULT;
          if (strcmp("*", multstr)==0)
            mult = 0; // 0 means use default multiplier
          else
            mult = (float)atof(multstr);
        }
        if (vic411_set_output_var((*out_data_files), TRUE, outfilenum, out_data, varname, outvarnum, format, type, mult) != 0) {
          vic411_nrerror("vic411_Error in global param file: Invalid output variable specification.");
        }
        strcpy(format,"");
        outvarnum++;
      }

    }
    fgets(cmdstr,MAXSTRING,gp);
  }
  fclose(gp);

}
