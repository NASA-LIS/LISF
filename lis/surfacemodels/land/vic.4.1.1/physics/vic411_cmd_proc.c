#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_cmd_proc.c,v 4.3.2.1 2006/11/17 01:06:20 vicadmin Exp $";

vic411_filenames_struct vic411_cmd_proc(int argc, char *argv[]) 
/**********************************************************************
  vic411_cmd_proc                  Keith Cherkauer                1997

  This routine checks the command line for valid program vic411_options.  If
  no vic411_options are found, or an invalid combination of them appear, the
  routine calls vic411_usage() to print the model vic411_usage to the screen, before
  exiting execution.

  Modifications:
  11-18-98  Added comment block to vic411_cmd_proc() and fixed routine so
            that it will exit if global command file is not defined
            using the "-g" vic411_flag.                                KAC
  30-Oct-03 Added -v option to display vic411_version information.	TJB

**********************************************************************/
{
  extern vic411_option_struct vic411_options;
#if LINK_DEBUG
  extern vic411_debug_struct debug;
#endif
  extern int getopt();
  extern char *optarg;
  extern char *vic411_optstring;

  vic411_filenames_struct names;
  int              optchar;
  char             GLOBAL_SET;
  
  if(argc==1) {
    vic411_usage(argv[0]);
    exit(1);
  }
  
  GLOBAL_SET = FALSE;

  while((optchar = getopt(argc, argv, vic411_optstring)) != EOF) {
    switch((char)optchar) {
    case 'v':
      /** Version information **/
      vic411_display_current_settings(DISP_VERSION,(vic411_filenames_struct*)NULL,(vic411_global_param_struct*)NULL);
      exit(0);
      break;
    case 'o':
      /** Compile-time vic411_options information **/
      vic411_display_current_settings(DISP_COMPILE_TIME,(vic411_filenames_struct*)NULL,(vic411_global_param_struct*)NULL);
      exit(0);
      break;
    case 'g':
      /** Global Parameters File **/
      strcpy(names.global, optarg);
      GLOBAL_SET = TRUE;
      break;
    default:
      /** Print Usage if Invalid Command Line Arguments **/
      vic411_usage(argv[0]);
      exit(1);
      break;
    }
  }

  if(!GLOBAL_SET) {
    fprintf(stderr,"ERROR: Must set global control file using the '-g' vic411_flag\n");
    vic411_usage(argv[0]);
    exit(1);
  }

  return names;
}


void vic411_usage(char *temp)
/**********************************************************************
	vic411_usage		Keith Cherkauer		May 27, 1996

  This routine prints out vic411_usage details.

**********************************************************************/
{
  fprintf(stderr,"Usage: %s [-v | -o | -g<global_parameter_file>]\n",temp);
  fprintf(stderr,"  v: display vic411_version information\n");
  fprintf(stderr,"  o: display compile-time vic411_options settings (set in vic411_user_def.h)\n");
  fprintf(stderr,"  g: read model parameters from <global_parameter_file>.\n");
  fprintf(stderr,"       <global_parameter_file> is a file that contains all needed model\n");
  fprintf(stderr,"       parameters as well as model option flags, and the names and\n");
  fprintf(stderr,"       locations of all other files.\n");
}
