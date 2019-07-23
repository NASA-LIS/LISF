#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>
#include <string.h>
 
static char vcid[] = "$Id: vic411_get_global_param.c,v 5.22.2.27 2009/09/20 02:32:07 vicadmin Exp $";

/********************************************************************/
/*			GLOBAL VARIABLES                            */
/********************************************************************/
int vic411_NR;		      /* array index for atmos struct that indicates
			 the model step avarage or sum */
int vic411_NF;		      /* array index loop counter limit for atmos
			 struct that indicates the SNOW_STEP values */
 
vic411_global_param_struct vic411_get_global_param(vic411_filenames_struct *names,
                                     FILE             *gp)
/**********************************************************************
  vic411_get_global_param	Keith Cherkauer	            March 1998

  This routine reads the VIC model global control file, getting
  values for global parameters, model vic411_options, and debugging controls.

  NOTE: any additions or removals of parameters in this file must also
  be made in vic411_display_current_settings.c.

  Modifications:
  7-19-96 Modified to read time step		        KAC
  4-5-98  Modified to read model vic411_options and debugging
          controls from a single file                   KAC
  01-20-00 modified to work with new radiation estimation routines,
           new simplified frozen soil moisture, and new new open
           format forcing file rad routines.              KAC
  02-27-01 added reads for lake model parameters          KAC
  04-21-03 added parameters for blowing snow algorithm, printing
           lake variables during debugging and reading Bart's 
           new Arno parameters.                           KAC
  11-May-04 Modified to display compile-time and run-time vic411_options
	    if VERBOSE is set to TRUE.						TJB
  13-Oct-04 Added validation for GRND_FLUX option.              		TJB
  01-Nov-04 Added validation for Nnodes with QUICK_FLUX option, as
	    part of fix for QUICK_FLUX state file compatibility.		TJB
  2005-03-08 Added EQUAL_AREA option.						TJB
  2005-03-24 Added ALMA_OUTPUT option.						TJB
  2005-04-07 Fixed state file warning check.					TJB
  2005-Apr-13 Added logic for OUTPUT_FORCE option.				TJB
  2005-Apr-23 Changed ARNO_PARAMS to NIJSSEN2001_BASEFLOW			TJB
  2005-11-29 SAVE_STATE is set in vic411_global_param (not at compile time)		GCT
  2005-12-06 Moved setting of statename from vic411_open_state_file to here.		GCT
  2005-12-07 Added checks for range of STATEMONTH and STATEDAY			GCT
  2005-12-07 Allow user to use NO_FLUX in addition to NOFLUX for NOFLUX in
             global.param.file							GCT
  2006-09-13 Replaced NIJSSEN2001_BASEFLOW with BASEFLOW option.		TJB/GCT
  2006-Sep-23 Implemented flexible output configuration; removed the
              OPTIMIZE and LDAS_OUTPUT vic411_options; implemented aggregation of
	      output variables.							TJB
  2006-Oct-16 Merged infiles and outfiles structs into vic411_filep_struct;
	      This included moving global->statename to filenames->statefile;
	      also added f_path_pfx to store forcing file path and prefix.	TJB
  2006-Nov-07 Removed LAKE_MODEL option.					TJB
  2007-Jan-03 Added ALMA_INPUT option.						TJB
  2007-Jan-15 Added PRT_HEADER option.						TJB
  2007-Apr-03 Added CONTINUEONERROR option.					GCT
  2007-Apr-24 Added EXP_TRANS option.						JCA
  2007-Apr-24 Added IMPLICIT option.						JCA
  2007-Apr-23 Added initialization of global parameters.			TJB
  2007-Apr-23 Added check for FULL_ENERGY if lake model is run.			TJB
  2007-Aug-08 Added EXCESS_ICE option.						JCA
  2007-Sep-14 Added initialization of names->soil_dir.				TJB
  2007-Oct-10 Added validation of dt, start date, end date, and nrecs.		TJB
  2007-Oct-31 Added validation of input/output files.				TJB
  2008-Jan-25 Removed setting of SNOW_STEP = global.dt for
	      OUTPUT_FORCE == TRUE.						TJB
  2008-Jan-28 Added check that end date falls AFTER start date.			TJB
  2008-Mar-12 Relocated code validating IMPLICIT and EXCESS_ICE vic411_options.	TJB
  2008-Apr-21 Added SNOW_ALBEDO option.						KAC via TJB
  2008-Apr-21 Added SNOW_DENSITY option.					TJB
  2009-Jan-12 Added COMPUTE_TREELINE and JULY_TAVG_SUPPLIED vic411_options.		TJB
  2009-Jan-16 Modified aero_resist_used and Ra_used to become arrays of
	      two elements (surface and overstory); added
	      vic411_options.AERO_RESIST_CANSNOW.					TJB
  2009-May-17 Added AR_406_LS to vic411_options.AERO_RESIST_CANSNOW.			TJB
  2009-May-17 Added vic411_options.MIN_LIQ.						TJB
  2009-May-18 Added vic411_options.PLAPSE.						TJB
  2009-May-20 Added vic411_options.GRND_FLUX_TYPE.					TJB
  2009-May-22 Added TFALLBACK value to vic411_options.CONTINUEONERROR.			TJB
  2009-Jun-15 Changed order of vic411_options to match global parameter file.		TJB
  2009-Aug-29 Now handles commented lines that are indented.			TJB
  2009-Sep-19 Moved TFALLBACK to its own separate option.			TJB
**********************************************************************/
{
  extern vic411_option_struct    vic411_options;
  extern vic411_param_set_struct vic411_param_set;
#if LINK_DEBUG
  extern vic411_debug_struct     debug;
#endif
  extern int              vic411_NF, vic411_NR;

  char cmdstr[MAXSTRING];
  char optstr[MAXSTRING];
  char flgstr[MAXSTRING];
  char ErrStr[MAXSTRING];
  int  file_num;
  int  field;
  int  i;
  int  tmpstartdate;
  int  tmpenddate;
  int  lastvalidday;
  int  lastday[] = {
            31, /* JANUARY */
            28, /* FEBRUARY */
            31, /* MARCH */
            30, /* APRIL */
            31, /* MAY */
            30, /* JUNE */
            31, /* JULY */
            31, /* AUGUST */
            30, /* SEPTEMBER */
            31, /* OCTOBER */
            30, /* NOVEMBER */
            31, /* DECEMBER */
        } ;
  vic411_global_param_struct global;

  /** Initialize global parameters (that aren't part of the vic411_options struct) **/
  global.dt            = MISSING;
  global.nrecs         = MISSING;
  global.startyear     = MISSING;
  global.startmonth    = MISSING;
  global.startday      = MISSING;
  global.starthour     = MISSING;
  global.endyear       = MISSING;
  global.endmonth      = MISSING;
  global.endday        = MISSING;
  global.resolution    = MISSING;
  global.MAX_SNOW_TEMP = 0;
  global.MIN_RAIN_TEMP = 0;
  global.measure_h     = 2.0;
  global.wind_h        = 10.0;
  for(i = 0; i < 2; i++) {
    global.forceyear[i]  = MISSING;
    global.forcemonth[i] = 1;
    global.forceday[i]   = 1;
    global.forcehour[i]  = 0;
    global.forceskip[i]  = 0;
    strcpy(names->f_path_pfx[i],"MISSING");
  }
  file_num             = 0;
  global.skipyear      = 0;
  strcpy(names->init_state,   "MISSING");
  global.stateyear     = MISSING;
  global.statemonth    = MISSING;
  global.stateday      = MISSING;
  strcpy(names->statefile,    "MISSING");
  strcpy(names->soil,         "MISSING");
  strcpy(names->soil_dir,     "MISSING");
  strcpy(names->veg,          "MISSING");
  strcpy(names->veglib,       "MISSING");
  strcpy(names->snowband,     "MISSING");
  strcpy(names->lakeparam,    "MISSING");
  strcpy(names->result_dir,   "MISSING");
  global.out_dt        = MISSING;


  /** Read through global control file to find parameters **/

  fgets(cmdstr,MAXSTRING,gp);

  while(!feof(gp)) {
    if(cmdstr[0]!='#' && cmdstr[0]!='\n' && cmdstr[0]!='\0') {

      sscanf(cmdstr,"%s",optstr);

      /* Handle case of comment line in which '#' is indented */
      if (optstr[0] == '#') {
        fgets(cmdstr,MAXSTRING,gp);
        continue;
      }

      /*************************************
       Get Model Global Parameters
      *************************************/
      if(strcasecmp("NLAYER",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&vic411_options.Nlayer);
      }
      else if(strcasecmp("NODES",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&vic411_options.Nnode);
      }
      else if(strcasecmp("TIME_STEP",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.dt);
      }
      else if(strcasecmp("SNOW_STEP",optstr)==0) {
	sscanf(cmdstr,"%*s %d",&vic411_options.SNOW_STEP);
      }
      else if(strcasecmp("STARTYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.startyear);
      }
      else if(strcasecmp("STARTMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.startmonth);
      }
      else if(strcasecmp("STARTDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.startday);
      }
      else if(strcasecmp("STARTHOUR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.starthour);
      }
      else if(strcasecmp("NRECS",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.nrecs);
      }
      else if(strcasecmp("ENDYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.endyear);
      }
      else if(strcasecmp("ENDMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.endmonth);
      }
      else if(strcasecmp("ENDDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.endday);
      }
      else if(strcasecmp("FULL_ENERGY",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) {
	  vic411_options.FULL_ENERGY=TRUE;
 	  vic411_options.QUICK_FLUX=TRUE;
	  vic411_options.GRND_FLUX=TRUE;
	}
	else vic411_options.FULL_ENERGY = FALSE;
      }
      else if(strcasecmp("FROZEN_SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) {
	  vic411_options.FROZEN_SOIL=TRUE;
	  vic411_options.QUICK_FLUX=FALSE;
	  vic411_options.GRND_FLUX=TRUE;
	}
        else vic411_options.FROZEN_SOIL = FALSE;
      }
      else if(strcasecmp("GRND_FLUX",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.GRND_FLUX=TRUE;
        else vic411_options.GRND_FLUX = FALSE;
      }
      else if(strcasecmp("QUICK_FLUX",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.QUICK_FLUX=TRUE;
        else vic411_options.QUICK_FLUX = FALSE;
      }
      else if(strcasecmp("QUICK_SOLVE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.QUICK_SOLVE=TRUE;
        else vic411_options.QUICK_SOLVE = FALSE;
      }
      else if( (strcasecmp("NOFLUX",optstr)==0) || (strcasecmp("NO_FLUX",optstr)==0) ) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.NOFLUX=TRUE;
        else vic411_options.NOFLUX = FALSE;
      }
      else if(strcasecmp("IMPLICIT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.IMPLICIT=TRUE;
        else vic411_options.IMPLICIT = FALSE;
      }
      else if(strcasecmp("EXP_TRANS",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.EXP_TRANS=TRUE;
        else vic411_options.EXP_TRANS = FALSE;
      }
      else if (strcasecmp("SNOW_ALBEDO", optstr)==0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if(strcasecmp("SUN1999",flgstr)==0) vic411_options.SNOW_ALBEDO=SUN1999;
        else vic411_options.SNOW_ALBEDO = USACE;
      }
      else if (strcasecmp("SNOW_DENSITY", optstr)==0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if(strcasecmp("DENS_SNTHRM",flgstr)==0) vic411_options.SNOW_DENSITY=DENS_SNTHRM;
        else vic411_options.SNOW_DENSITY = DENS_BRAS;
      }
      else if(strcasecmp("BLOWING",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.BLOWING=TRUE;
        else vic411_options.BLOWING = FALSE;
      }
      else if(strcasecmp("DIST_PRCP",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.DIST_PRCP=TRUE;
        else vic411_options.DIST_PRCP = FALSE;
      }
      else if(strcasecmp("PREC_EXPT",optstr)==0) {
	sscanf(cmdstr,"%*s %f",&vic411_options.PREC_EXPT);
      }
      else if(strcasecmp("CORRPREC",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.CORRPREC=TRUE;
        else vic411_options.CORRPREC = FALSE;
      }
      else if(strcasecmp("MIN_WIND_SPEED",optstr)==0) {
	sscanf(cmdstr,"%*s %f",&vic411_options.MIN_WIND_SPEED);
      }
      else if(strcasecmp("MIN_RAIN_TEMP",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&global.MIN_RAIN_TEMP);
      }
      else if(strcasecmp("MAX_SNOW_TEMP",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&global.MAX_SNOW_TEMP);
      }
      else if(strcasecmp("CONTINUEONERROR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.CONTINUEONERROR=TRUE;
        else vic411_options.CONTINUEONERROR = FALSE;
      }
      else if(strcasecmp("COMPUTE_TREELINE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("FALSE",flgstr)==0) vic411_options.COMPUTE_TREELINE=FALSE;
        else {
          vic411_options.COMPUTE_TREELINE = TRUE;
          vic411_options.AboveTreelineVeg = atoi( flgstr );
        }
      }
      else if(strcasecmp("EQUAL_AREA",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.EQUAL_AREA=TRUE;
        else vic411_options.EQUAL_AREA = FALSE;
      }
      else if(strcasecmp("RESOLUTION",optstr)==0) {
        sscanf(cmdstr,"%*s %f",&global.resolution);
      }
      else if (strcasecmp("AERO_RESIST_CANSNOW", optstr)==0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if(strcasecmp("AR_406",flgstr)==0) vic411_options.AERO_RESIST_CANSNOW=AR_406;
        else if(strcasecmp("AR_406_LS",flgstr)==0) vic411_options.AERO_RESIST_CANSNOW=AR_406_LS;
        else if(strcasecmp("AR_406_FULL",flgstr)==0) vic411_options.AERO_RESIST_CANSNOW=AR_406_FULL;
        else if(strcasecmp("AR_410",flgstr)==0) vic411_options.AERO_RESIST_CANSNOW=AR_410;
        else if(strcasecmp("AR_COMBO",flgstr)==0) vic411_options.AERO_RESIST_CANSNOW=AR_COMBO;
      }
      else if (strcasecmp("GRND_FLUX_TYPE", optstr)==0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if(strcasecmp("GF_406",flgstr)==0) vic411_options.GRND_FLUX_TYPE=GF_406;
        else if(strcasecmp("GF_410",flgstr)==0) vic411_options.GRND_FLUX_TYPE=GF_410;
        else if(strcasecmp("GF_FULL",flgstr)==0) vic411_options.GRND_FLUX_TYPE=GF_FULL;
      }
      else if(strcasecmp("MIN_LIQ",optstr)==0) {
        sscanf(cmdstr,"%*s %s", flgstr);
        if(strcasecmp("FALSE", flgstr) == 0) vic411_options.MIN_LIQ = FALSE;
        else {
	  vic411_options.MIN_LIQ = TRUE;
	}
      }
      else if(strcasecmp("PLAPSE",optstr)==0) {
        sscanf(cmdstr,"%*s %s", flgstr);
        if(strcasecmp("FALSE", flgstr) == 0) vic411_options.PLAPSE = FALSE;
        else {
	  vic411_options.PLAPSE = TRUE;
	}
      }
      else if(strcasecmp("TFALLBACK",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.TFALLBACK=TRUE;
        else vic411_options.TFALLBACK = FALSE;
      }

      /*************************************
       Define state files
      *************************************/
      else if(strcasecmp("INIT_STATE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("FALSE",flgstr)==0) vic411_options.INIT_STATE=FALSE;
        else {
	  vic411_options.INIT_STATE = TRUE;
	  strcpy(names->init_state,flgstr);
	}
      }
      else if(strcasecmp("STATENAME",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->statefile);
        vic411_options.SAVE_STATE = TRUE;
      }
      else if(strcasecmp("STATEYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.stateyear);
      }
      else if(strcasecmp("STATEMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.statemonth);
      }
      else if(strcasecmp("STATEDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.stateday);
      }
      else if(strcasecmp("BINARY_STATE_FILE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("FALSE",flgstr)==0) vic411_options.BINARY_STATE_FILE=FALSE;
	else vic411_options.BINARY_STATE_FILE=TRUE;
      }

      /*************************************
       Define forcing files
      *************************************/
      else if(strcasecmp("FORCING1",optstr)==0) {
	if ( strcmp( names->f_path_pfx[0], "MISSING" ) != 0 ) 
	  vic411_nrerror("Tried to define FORCING1 twice, if you want to use two forcing files, the second must be defined as FORCING2");
        sscanf(cmdstr,"%*s %s", names->f_path_pfx[0]);
	file_num = 0;
	field=0;
      }
      else if(strcasecmp("FORCING2",optstr)==0) {
        sscanf(cmdstr,"%*s %s", names->f_path_pfx[1]);
        if (strcasecmp("FALSE",names->f_path_pfx[1])==0)
          strcpy(names->f_path_pfx[1],"MISSING");
	file_num = 1;
	field=0;
      }
      else if (strcasecmp("FORCE_FORMAT",optstr)==0) {
	sscanf(cmdstr, "%*s %s", flgstr);
	if (strcasecmp(flgstr, "BINARY") == 0)
	  vic411_param_set.FORCE_FORMAT[file_num] = BINARY;
	else if (strcasecmp(flgstr, "ASCII") == 0)
	  vic411_param_set.FORCE_FORMAT[file_num] = ASCII;
	else
	  vic411_nrerror("FORCE_FORMAT must be either ASCII or BINARY.");
      }
      else if (strcasecmp("FORCE_ENDIAN",optstr)==0) {
	sscanf(cmdstr, "%*s %s", flgstr);
	if (strcasecmp(flgstr, "LITTLE") == 0)
	  vic411_param_set.FORCE_ENDIAN[file_num] = LITTLE;
	else if (strcasecmp(flgstr, "BIG") == 0)
	  vic411_param_set.FORCE_ENDIAN[file_num] = BIG;
	else
	  vic411_nrerror("FORCE_ENDIAN must be either BIG or LITTLE.");
      }
      else if(strcasecmp("N_TYPES",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&vic411_param_set.N_TYPES[file_num]);
      }
      else if(strcasecmp("FORCE_TYPE",optstr)==0) {
	vic411_get_force_type(cmdstr,file_num,&field);
      }
      else if(strcasecmp("FORCE_DT",optstr)==0) {
	sscanf(cmdstr,"%*s %d ", &vic411_param_set.FORCE_DT[file_num]);
      }
      else if(strcasecmp("FORCEYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.forceyear[file_num]);
      }
      else if(strcasecmp("FORCEMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.forcemonth[file_num]);
      }
      else if(strcasecmp("FORCEDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.forceday[file_num]);
      }
      else if(strcasecmp("FORCEHOUR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.forcehour[file_num]);
      }
      else if(strcasecmp("GRID_DECIMAL",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&vic411_options.GRID_DECIMAL);
      }
      else if(strcasecmp("WIND_H",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&global.wind_h);
      }
      else if(strcasecmp("MEASURE_H",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&global.measure_h);
      }
      else if(strcasecmp("ALMA_INPUT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.ALMA_INPUT=TRUE;
        else vic411_options.ALMA_INPUT = FALSE;
      }

      /*************************************
       Define parameter files
      *************************************/

      else if(strcasecmp("SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->soil);
      }
      else if(strcasecmp("ARC_SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.ARC_SOIL=TRUE;
        else vic411_options.ARC_SOIL = FALSE;
      }
      else if(strcasecmp("SOIL_DIR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->soil_dir);
      }
      else if (strcasecmp("ARNO_PARAMS", optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) {
          vic411_nrerror("Please change \"ARNO_PARAMS  TRUE\" to \"BASEFLOW  NIJSSEN2001\" in your global parameter file.");
        }
        else {
          vic411_nrerror("Please change \"ARNO_PARAMS  FALSE\" to \"BASEFLOW  ARNO\" in your global parameter file.");
        }
      }
      else if (strcasecmp("NIJSSEN2001_BASEFLOW", optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) {
          vic411_nrerror("Please change \"NIJSSEN2001_BASEFLOW  TRUE\" to \"BASEFLOW  NIJSSEN2001\" in your global parameter file.");
        }
        else {
          vic411_nrerror("Please change \"NIJSSEN2001_BASEFLOW  FALSE\" to \"BASEFLOW  ARNO\" in your global parameter file.");
        }
      }
      else if (strcasecmp("BASEFLOW", optstr)==0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if(strcasecmp("NIJSSEN2001",flgstr)==0) vic411_options.BASEFLOW=NIJSSEN2001;
        else vic411_options.BASEFLOW = ARNO;
      }
      else if(strcasecmp("JULY_TAVG_SUPPLIED",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("FALSE",flgstr)==0) vic411_options.JULY_TAVG_SUPPLIED=FALSE;
	else vic411_options.JULY_TAVG_SUPPLIED=TRUE;
      }
      else if(strcasecmp("VEGPARAM",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->veg);
      }
      else if(strcasecmp("ROOT_ZONES",optstr)==0) {
	sscanf(cmdstr,"%*s %d",&vic411_options.ROOT_ZONES);
      }
      else if(strcasecmp("GLOBAL_LAI",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.GLOBAL_LAI=TRUE;
        else vic411_options.GLOBAL_LAI = FALSE;
      }
      else if(strcasecmp("VEGLIB",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->veglib);
      }
      else if(strcasecmp("SNOW_BAND",optstr)==0) {
	sscanf(cmdstr,"%*s %d %s",&vic411_options.SNOW_BAND,names->snowband);
      }
      else if(strcasecmp("LAKES",optstr)==0) {
        sscanf(cmdstr,"%*s %s", flgstr);
        if(strcasecmp("FALSE", flgstr) == 0) vic411_options.LAKES = FALSE;
        else {
	  vic411_options.LAKES = TRUE;
	  strcpy(names->lakeparam, flgstr);
	}
      }
      else if(strcasecmp("LAKE_PROFILE",optstr)==0) {
        sscanf(cmdstr,"%*s %s", flgstr);
        if(strcasecmp("FALSE", flgstr) == 0) vic411_options.LAKE_PROFILE = FALSE;
        else {
	  vic411_options.LAKE_PROFILE = TRUE;
	}
      }

      /*************************************
       Define output files
      *************************************/
      else if(strcasecmp("RESULT_DIR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->result_dir);
      }
      else if(strcasecmp("OUT_STEP",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.out_dt);
      }
      else if(strcasecmp("SKIPYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global.skipyear);
      }
      else if(strcasecmp("COMPRESS",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.COMPRESS=TRUE;
        else vic411_options.COMPRESS = FALSE;
      }
      else if(strcasecmp("BINARY_OUTPUT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.BINARY_OUTPUT=TRUE;
        else vic411_options.BINARY_OUTPUT = FALSE;
      }
      else if(strcasecmp("ALMA_OUTPUT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.ALMA_OUTPUT=TRUE;
        else vic411_options.ALMA_OUTPUT = FALSE;
      }
      else if(strcasecmp("MOISTFRACT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.MOISTFRACT=TRUE;
        else vic411_options.MOISTFRACT = FALSE;
      }
      else if(strcasecmp("PRT_HEADER",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.PRT_HEADER=TRUE;
        else vic411_options.PRT_HEADER = FALSE;
      }
      else if(strcasecmp("PRT_SNOW_BAND",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) vic411_options.PRT_SNOW_BAND=TRUE;
        else vic411_options.PRT_SNOW_BAND = FALSE;
      }

      /*************************************
       Define output file contents
      *************************************/
      else if(strcasecmp("N_OUTFILES",optstr)==0) {
        ; // do nothing
      }
      else if(strcasecmp("OUTFILE",optstr)==0) {
        ; // do nothing
      }
      else if(strcasecmp("OUTVAR",optstr)==0) {
        ; // do nothing
      }

      /******************************
        Get Model Debugging Options
	****************************/

#if LINK_DEBUG
      else if(strcasecmp("PRT_FLUX",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_FLUX=TRUE;
        else debug.PRT_FLUX = FALSE;
      }
      else if(strcasecmp("PRT_BALANCE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_BALANCE=TRUE;
        else debug.PRT_BALANCE = FALSE;
      }
      else if(strcasecmp("PRT_SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_SOIL=TRUE;
        else debug.PRT_SOIL = FALSE;
      }
      else if(strcasecmp("PRT_VEGE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_VEGE=TRUE;
        else debug.PRT_VEGE = FALSE;
      }
      else if(strcasecmp("PRT_GLOBAL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_GLOBAL=TRUE;
        else debug.PRT_GLOBAL = FALSE;
      }
      else if(strcasecmp("PRT_ATMOS",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_ATMOS=TRUE;
        else debug.PRT_ATMOS = FALSE;
      }
      else if(strcasecmp("PRT_SNOW",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_SNOW=TRUE;
        else debug.PRT_SNOW = FALSE;
      }
      else if(strcasecmp("PRT_MOIST",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_MOIST=TRUE;
        else debug.PRT_MOIST = FALSE;
      }
      else if(strcasecmp("PRT_LAKE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_LAKE=TRUE;
        else debug.PRT_LAKE = FALSE;
      }
      else if(strcasecmp("PRT_TEMP",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_TEMP=TRUE;
        else debug.PRT_TEMP = FALSE;
      }
      else if(strcasecmp("DEBUG_DIR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",debug.debug_dir);
      }
#endif

      /***********************************
        Unrecognized Global Parameter Flag
        ***********************************/
      else {
	fprintf(stderr,"WARNING: Unrecognized option in the global parameter file:\n\t%s is unknown - check your spelling\n", optstr);
      }
    }
    fgets(cmdstr,MAXSTRING,gp);
  }

  /******************************************
    Check for undefined required parameters
  ******************************************/

  // Validate model time step
  if (global.dt == MISSING)
    vic411_nrerror("Model time step has not been defined.  Make sure that the global file defines TIME_STEP.");
  else if (global.dt < 1) {
    sprintf(ErrStr,"The specified model time step (%d) < 1 hour.  Make sure that the global file defines a positive number of hours for TIME_STEP.",global.dt);
    vic411_nrerror(ErrStr);
  }

  // Validate the output step
  if (global.out_dt == 0 || global.out_dt == MISSING) {
    global.out_dt = global.dt;
  }
  else if (global.out_dt < global.dt || global.out_dt > 24 || (float)global.out_dt/(float)global.dt != (float)(global.out_dt/global.dt)){
    vic411_nrerror("Invalid output step specified.  Output step must be an integer multiple of the model time step; >= model time step and <= 24");
  }

  // Validate SNOW_STEP and set vic411_NR and vic411_NF
  if (global.dt < 24 && global.dt != vic411_options.SNOW_STEP)
    vic411_nrerror("If the model step is smaller than daily, the snow model should run\nat the same time step as the rest of the model.");
  if (global.dt % vic411_options.SNOW_STEP != 0 || vic411_options.SNOW_STEP > global.dt)
    vic411_nrerror("SNOW_STEP should be <= TIME_STEP and divide TIME_STEP evenly ");
  vic411_NF = global.dt/vic411_options.SNOW_STEP;
  if (vic411_NF == 1)
    vic411_NR = 0;
  else
    vic411_NR = vic411_NF;

  // Validate simulation start date
  if (global.startyear == MISSING)
    vic411_nrerror("Simulation start year has not been defined.  Make sure that the global file defines STARTYEAR.");
  else if (global.startyear < 0) {
    sprintf(ErrStr,"The specified simulation start year (%d) < 0.  Make sure that the global file defines a positive integer for STARTYEAR.",global.startyear);
    vic411_nrerror(ErrStr);
  }
  if (global.startmonth == MISSING)
    vic411_nrerror("Simulation start month has not been defined.  Make sure that the global file defines STARTMONTH.");
  else if (global.startmonth < 0) {
    sprintf(ErrStr,"The specified simulation start month (%d) < 0.  Make sure that the global file defines a positive integer for STARTMONTH.",global.startmonth);
    vic411_nrerror(ErrStr);
  }
  if (global.startday == MISSING)
    vic411_nrerror("Simulation start day has not been defined.  Make sure that the global file defines STARTDAY.");
  else if (global.startday < 0) {
    sprintf(ErrStr,"The specified simulation start day (%d) < 0.  Make sure that the global file defines a positive integer for STARTDAY.",global.startday);
    vic411_nrerror(ErrStr);
  }
  if (global.starthour == MISSING) {
    if (global.dt == 24)
      global.starthour = 0;
    else
      vic411_nrerror("Simulation start hour has not been defined, yet model time step is less than 24 hours.  Make sure that the global file defines STARTHOUR.");
  }
  else if (global.starthour < 0) {
    sprintf(ErrStr,"The specified simulation start hour (%d) < 0.  Make sure that the global file defines a positive integer for STARTHOUR.",global.starthour);
    vic411_nrerror(ErrStr);
  }

  // Validate simulation end date and/or number of timesteps
  if (global.nrecs == MISSING && global.endyear == MISSING && global.endmonth == MISSING && global.endday == MISSING)
    vic411_nrerror("The model global file MUST define EITHER the number of records to simulate (NRECS), or the year (ENDYEAR), month (ENDMONTH), and day (ENDDAY) of the last full simulation day");
  else if (global.nrecs == MISSING) {
    if (global.endyear == MISSING)
      vic411_nrerror("Simulation end year has not been defined.  Make sure that the global file defines ENDYEAR.");
    else if (global.endyear < 0) {
      sprintf(ErrStr,"The specified simulation end year (%d) < 0.  Make sure that the global file defines a positive integer for ENDYEAR.",global.endyear);
      vic411_nrerror(ErrStr);
    }
    if (global.endmonth == MISSING)
      vic411_nrerror("Simulation end month has not been defined.  Make sure that the global file defines ENDMONTH.");
    else if (global.endmonth < 0) {
      sprintf(ErrStr,"The specified simulation end month (%d) < 0.  Make sure that the global file defines a positive integer for ENDMONTH.",global.endmonth);
      vic411_nrerror(ErrStr);
    }
    if (global.endday == MISSING)
      vic411_nrerror("Simulation end day has not been defined.  Make sure that the global file defines ENDDAY.");
    else if (global.endday < 0) {
      sprintf(ErrStr,"The specified simulation end day (%d) < 0.  Make sure that the global file defines a positive integer for ENDDAY.",global.endday);
      vic411_nrerror(ErrStr);
    }
    tmpstartdate = global.startyear*10000 + global.startmonth*100 + global.startday;
    tmpenddate = global.endyear*10000 + global.endmonth*100 + global.endday;
    if (tmpenddate < tmpstartdate) {
      sprintf(ErrStr,"The specified simulation end date (%04d-%02d-%02d) is EARLIER than the specified start date (%04d-%02d-%02d).",global.endyear,global.endmonth,global.endday,global.startyear,global.startmonth,global.startday);
      vic411_nrerror(ErrStr);
    }
  }
  else if (global.nrecs < 1) {
    sprintf(ErrStr,"The specified duration of simulation (%d) < 1 time step.  Make sure that the global file defines a positive integer for NRECS.",global.nrecs);
    vic411_nrerror(ErrStr);
  }

  // Validate forcing files and variables
  if ( strcmp ( names->f_path_pfx[0], "MISSING" ) == 0 )
    vic411_nrerror("No forcing file has been defined.  Make sure that the global file defines FORCING1.");
  for(i=0;i<2;i++) {
    if ( i == 0 || (i == 1 && vic411_param_set.N_TYPES[i] != MISSING) ) {
      if (vic411_param_set.N_TYPES[i] == MISSING) {
        sprintf(ErrStr,"Need to specify the number forcing variables types in forcing file %d.", i);
        vic411_nrerror(ErrStr);
      }
      if (vic411_param_set.FORCE_FORMAT[i] == MISSING) {
        sprintf(ErrStr,"Need to specify the INPUT_FORMAT (ASCII or BINARY) for forcing file %d.",i);
        vic411_nrerror(ErrStr);
      }
      if (vic411_param_set.FORCE_INDEX[i][vic411_param_set.N_TYPES[i]-1] == MISSING) {
        sprintf(ErrStr,"Did not define enough forcing variables in forcing file %d.",i);
        vic411_nrerror(ErrStr);
      }
      if(vic411_param_set.FORCE_DT[i] == MISSING ) {
        sprintf(ErrStr,"Must define time steps (FORCE_DT <dt>) in control file for focing file %d.",file_num);
        vic411_nrerror(ErrStr);
      }
    }
  }
  if(vic411_param_set.N_TYPES[1] != MISSING && global.forceyear[1] == MISSING) {
    global.forceyear[1] = global.forceyear[0];
    global.forcemonth[1] = global.forcemonth[0];
    global.forceday[1] = global.forceday[0];
    global.forcehour[1] = global.forcehour[0];
    global.forceskip[1] = 0;
  }

  // Validate result directory
  if ( strcmp ( names->result_dir, "MISSING" ) == 0 )
    vic411_nrerror("No results directory has been defined.  Make sure that the global file defines the result directory on the line that begins with \"RESULT_DIR\".");

  // Validate soil parameter file information
  if ( strcmp ( names->soil, "MISSING" ) == 0 )
    vic411_nrerror("No soil parameter file has been defined.  Make sure that the global file defines the soil parameter file on the line that begins with \"SOIL\".");
  if (vic411_options.ARC_SOIL && strcmp ( names->soil_dir, "MISSING" ) == 0)
    vic411_nrerror("\"ARC_SOIL\" was specified as TRUE, but no soil parameter directory (\"SOIL_DIR\") has been defined.  Make sure that the global file defines the soil parameter directory on the line that begins with \"SOIL_DIR\".");

  /*******************************************************************************
    Validate parameters required for normal simulations but NOT for OUTPUT_FORCE
  *******************************************************************************/

#if !OUTPUT_FORCE

  // Validate veg parameter information
  if ( strcmp ( names->veg, "MISSING" ) == 0 )
    vic411_nrerror("No vegetation parameter file has been defined.  Make sure that the global file defines the vegetation parameter file on the line that begins with \"VEGPARAM\".");
  if ( strcmp ( names->veglib, "MISSING" ) == 0 )
    vic411_nrerror("No vegetation library file has been defined.  Make sure that the global file defines the vegetation library file on the line that begins with \"VEGLIB\".");
  if(vic411_options.ROOT_ZONES<0)
    vic411_nrerror("ROOT_ZONES must be defined to a positive integer greater than 0, in the global control file.");

  // Validate the elevation band file information
  if(vic411_options.SNOW_BAND > 1) {
    if ( strcmp ( names->snowband, "MISSING" ) == 0 ) {
      sprintf(ErrStr, "\"SNOW_BAND\" was specified with %d elevation bands, but no elevation band file has been defined.  Make sure that the global file defines the elevation band file on the line that begins with \"SNOW_BAND\" (after the number of bands).", vic411_options.SNOW_BAND);
      vic411_nrerror(ErrStr);
    }
    if(vic411_options.SNOW_BAND > MAX_BANDS) {
      sprintf(ErrStr,"Global file wants more snow bands (%d) than are defined by MAX_BANDS (%d).  Edit vic411_user_def.h and recompile.",vic411_options.SNOW_BAND,MAX_BANDS);
      vic411_nrerror(ErrStr);
    }
  }
  else if (vic411_options.SNOW_BAND <= 0) {
    sprintf(ErrStr,"Invalid number of elevation bands specified in global file (%d).  Number of bands must be >= 1.",vic411_options.SNOW_BAND);
    vic411_nrerror(ErrStr);
  }

  // Validate the input state file information
  if( vic411_options.INIT_STATE ) {
    if ( strcmp ( names->init_state, "MISSING" ) == 0 )
      vic411_nrerror("\"INIT_STATE\" was specified, but no input state file has been defined.  Make sure that the global file defines the inputstate file on the line that begins with \"INIT_STATE\".");
  }

  // Validate the output state file information
  if( vic411_options.SAVE_STATE ) {
    if ( strcmp ( names->statefile, "MISSING" ) == 0)
      vic411_nrerror("\"SAVE_STATE\" was specified, but no output state file has been defined.  Make sure that the global file defines the output state file on the line that begins with \"SAVE_STATE\".");
    if ( global.stateyear == MISSING || global.statemonth == MISSING || global.stateday == MISSING )  {
      sprintf(ErrStr,"Incomplete specification of the date to save state for state file (%s).\nSpecified date (yyyy-mm-dd): %04d-%02d-%02d\nMake sure STATEYEAR, STATEMONTH, and STATEDAY are set correctly in your global parameter file.\n", names->statefile, global.stateyear, global.statemonth, global.stateday);
      vic411_nrerror(ErrStr);
    }
    // Check for month, day in range
    lastvalidday = lastday[global.statemonth - 1];
    if ( global.statemonth == 2 ) {
      if ( (global.stateyear % 4) == 0 && ( (global.stateyear % 100) != 0 || (global.stateyear % 400) == 0 ) ){
        lastvalidday = 29;
      }
    }
    if ( global.stateday > lastvalidday || global.statemonth > 12 || global.statemonth < 1 || global.stateday > 31 || global.stateday < 1 ){
      sprintf(ErrStr,"Unusual specification of the date to save state for state file (%s).\nSpecified date (yyyy-mm-dd): %04d-%02d-%02d\nMake sure STATEYEAR, STATEMONTH, and STATEDAY are set correctly in your global parameter file.\n", names->statefile, global.stateyear, global.statemonth, global.stateday);
      vic411_nrerror(ErrStr);
    }
  }
  // Set the statename here to be able to compare with INIT_STATE name
  if( vic411_options.SAVE_STATE ) {
    sprintf(names->statefile,"%s_%04i%02i%02i", names->statefile,
          global.stateyear, global.statemonth, global.stateday);
  }
  if( vic411_options.INIT_STATE && vic411_options.SAVE_STATE && (strcmp( names->init_state, names->statefile ) == 0))  {
      sprintf(ErrStr,"The save state file (%s) has the same name as the initialize state file (%s).  The initialize state file will be destroyed when the save state file is opened.", names->statefile, names->init_state);
      vic411_nrerror(ErrStr);
  }

  // Validate soil parameter/simulation mode combinations
  if(vic411_options.Nlayer > MAX_LAYERS) {
    sprintf(ErrStr,"Global file wants more soil moisture layers (%d) than are defined by MAX_LAYERS (%d).  Edit vic411_user_def.h and recompile.",vic411_options.Nlayer,MAX_LAYERS);
    vic411_nrerror(ErrStr);
  }
  if(vic411_options.Nnode > MAX_NODES) {
    sprintf(ErrStr,"Global file wants more soil thermal nodes (%d) than are defined by MAX_NODES (%d).  Edit vic411_user_def.h and recompile.",vic411_options.Nnode,MAX_NODES);
    vic411_nrerror(ErrStr);
  }
  if(vic411_options.QUICK_FLUX) {
    if((vic411_options.FULL_ENERGY || vic411_options.FROZEN_SOIL) && vic411_options.Nnode != 3) {
      sprintf(ErrStr,"To run the model in FULL_ENERGY or FROZEN_SOIL modes with QUICK_FLUX=TRUE, you must define exactly 3 soil thermal nodes.  Currently Nnodes is set to  %d.",vic411_options.Nnode);
      vic411_nrerror(ErrStr);
    }
    else if (!vic411_options.FULL_ENERGY && !vic411_options.FROZEN_SOIL && vic411_options.Nnode != 1) {
      sprintf(ErrStr,"To run the model with FULL_ENERGY=FALSE, FROZEN_SOIL=FALSE, and QUICK_FLUX=TRUE, you must define exactly 1 soil thermal node.  Currently Nnodes is set to  %d.",vic411_options.Nnode);
      vic411_nrerror(ErrStr);
    }
  }
  else {
    if(vic411_options.Nnode < 4) {
      sprintf(ErrStr,"To run the model with QUICK_FLUX=FALSE, you must define at least 4 soil thermal nodes.  Currently Nnodes is set to %d.",vic411_options.Nnode);
      vic411_nrerror(ErrStr);
    }
  }
  if((vic411_options.FULL_ENERGY || vic411_options.FROZEN_SOIL) && vic411_options.Nlayer<3) {
    sprintf(ErrStr,"You must define at least 3 soil moisture layers to run the model in FULL_ENERGY or FROZEN_SOIL modes.  Currently Nlayers is set to  %d.",vic411_options.Nlayer);
    vic411_nrerror(ErrStr);
  }
  if((!vic411_options.FULL_ENERGY && !vic411_options.FROZEN_SOIL) && vic411_options.Nlayer<1) {
    sprintf(ErrStr,"You must define at least 1 soil moisture layer to run the model.  Currently Nlayers is set to  %d.",vic411_options.Nlayer);
    vic411_nrerror(ErrStr);
  }
  if(!vic411_options.FULL_ENERGY && !vic411_options.FROZEN_SOIL && vic411_options.GRND_FLUX) {
    sprintf(ErrStr,"Both FULL_ENERGY and FROZEN_SOIL are FALSE, but GRND_FLUX is TRUE.\nThis combination of vic411_options is not recommended.  Unless you intend\nto use this combination of vic411_options, we recommend commenting out\nthe GRND_FLUX entry in your global file.  To do this,place a \"#\"\nat the beginning of the line containing \"GRND_FLUX\".\n");
    vic411_nrerror(ErrStr);
  }
  if(vic411_options.FULL_ENERGY && !vic411_options.GRND_FLUX) {
    sprintf(ErrStr,"FULL_ENERGY is TRUE, but GRND_FLUX is explicitly set to FALSE.\nThis combination of vic411_options is not recommended.  Unless you intend\nto use this combination of vic411_options, we recommend commenting out\nthe GRND_FLUX entry in your global file.  To do this,place a \"#\"\nat the beginning of the line containing \"GRND_FLUX\".\n");
    vic411_nrerror(ErrStr);
  }
  if(vic411_options.IMPLICIT)  {
    if ( QUICK_FS ) 
      fprintf(stderr,"WARNING: IMPLICIT and QUICK_FS are both TRUE.\n\tThe QUICK_FS option is ignored when IMPLICIT=TRUE\n");
  }
  if( EXCESS_ICE ) {
    if ( !vic411_options.FULL_ENERGY )
      vic411_nrerror("set FULL_ENERGY = TRUE to run EXCESS_ICE option.");
    if ( !vic411_options.FROZEN_SOIL )
      vic411_nrerror("set FROZEN_SOIL = TRUE to run EXCESS_ICE option.");
    if ( !vic411_options.GRND_FLUX )
      vic411_nrerror("set GRND_FLUX = TRUE to run EXCESS_ICE option.");
    if ( vic411_options.QUICK_SOLVE ) {
      fprintf(stderr,"WARNING: QUICK_SOLVE and EXCESS_ICE are both TRUE.\n\tThis is an incompatible combination.  Setting QUICK_SOLVE to FALSE.\n");
      vic411_options.QUICK_SOLVE=FALSE;  
    }    
    if ( QUICK_FS ) 
      vic411_nrerror("QUICK_FS = TRUE and EXCESS_ICE = TRUE are incompatible vic411_options.");
  }

  // Validate lake parameter information
  if (vic411_options.LAKES) {
    if (!vic411_options.FULL_ENERGY) {
      sprintf(ErrStr, "FULL_ENERGY must be TRUE if the lake model is to be run.");
      vic411_nrerror(ErrStr);
    }
    if ( strcmp ( names->lakeparam, "MISSING" ) == 0 )
      vic411_nrerror("\"LAKES\" was specified, but no lake parameter file has been defined.  Make sure that the global file defines the lake parameter file on the line that begins with \"LAKES\".");
    if (global.resolution == 0) {
      sprintf(ErrStr, "The model grid cell resolution (RESOLUTION) must be defined in the global control file when the lake model is active.");
      vic411_nrerror(ErrStr);
    }
    if (global.resolution > 360 && !vic411_options.EQUAL_AREA) {
      sprintf(ErrStr, "For EQUAL_AREA=FALSE, the model grid cell resolution (RESOLUTION) must be set to the number of lat or lon degrees per grid cell.  This cannot exceed 360.");
      vic411_nrerror(ErrStr);
    }
    if (vic411_options.COMPUTE_TREELINE) {
      sprintf(ErrStr, "LAKES = TRUE and COMPUTE_TREELINE = TRUE are incompatible vic411_options.");
      vic411_nrerror(ErrStr);
    }
  }

  /*********************************
    Output major vic411_options to stderr
  *********************************/
#if VERBOSE
  vic411_display_current_settings(DISP_ALL,names,&global);
#else
  vic411_display_current_settings(DISP_VERSION,names,&global);
#endif

#if VERBOSE
  fprintf(stderr,"Time Step = %d hour(s)\n",global.dt);
  fprintf(stderr,"Simulation start date = %02i/%02i/%04i\n",
	  global.startday, global.startmonth, global.startyear);
  if ( global.nrecs > 0 )
    fprintf(stderr,"Number of Records = %d\n\n",global.nrecs);
  else 
    fprintf(stderr,"Simulation end date = %02i/%02i/%04i\n\n",
	    global.endday, global.endmonth, global.endyear);
  fprintf(stderr,"Full Energy...................(%d)\n",vic411_options.FULL_ENERGY);
  fprintf(stderr,"Use Distributed Precipitation.(%d)\n",vic411_options.DIST_PRCP);
  if(vic411_options.DIST_PRCP)
    fprintf(stderr,"..Using Precipitation Exponent of %f\n",vic411_options.PREC_EXPT);
  if ( vic411_options.GRND_FLUX ) {
    fprintf(stderr,"Ground heat flux will be estimated ");
    if ( vic411_options.QUICK_FLUX ) 
      fprintf(stderr,"using Liang, Wood and Lettenmaier (1999).\n");
    else 
      fprintf(stderr,"using Cherkauer and Lettenmaier (1999).\n");
  }
  else
    fprintf(stderr,"Ground heat flux not computed (no energy balance).\n");
  fprintf(stderr,"Use Frozen Soil Model.........(%d)\n",vic411_options.FROZEN_SOIL);
  if( vic411_options.IMPLICIT ) 
    fprintf(stderr,".... Using the implicit solution for the soil heat equation.\n");
  else
    fprintf(stderr,".... Using the explicit solution for the soil heat equation.\n");
  if( vic411_options.EXP_TRANS )
    fprintf(stderr,".... Thermal nodes are exponentially distributed with depth.\n");
  else
    fprintf(stderr,".... Thermal nodes are linearly distributed with depth (except top two nodes).\n");
  if( EXCESS_ICE )
    fprintf(stderr,".... Excess ground ice is being considered.\n\t\tTherefore, ground ice (as a volumetric fraction) must be initialized for each\n\t\t   soil layer in the soil file.\n\t\tCAUTION: When excess ice melts, subsidence occurs.\n\t\t  Therefore, soil layer depths, damping depth, thermal node depths,\n\t\t     bulk densities, porosities, and other properties are now dynamic!\n\t\t  EXERCISE EXTREME CAUTION IN INTERPRETING MODEL OUTPUT.\n\t\t  It is recommended to add OUT_SOIL_DEPTH to your list of output variables.\n");
  if ( QUICK_FS ){
    fprintf(stderr,".... Using linearized UFWC curve with %d temperatures.\n", QUICK_FS_TEMPS);
  }
  fprintf(stderr,"Run Snow Model Using a Time Step of %d hours\n", 
	  vic411_options.SNOW_STEP);
  fprintf(stderr,"Compress Output Files.........(%d)\n",vic411_options.COMPRESS);
  fprintf(stderr,"Correct Precipitation.........(%d)\n",vic411_options.CORRPREC);
  fprintf(stderr,"\n");
  fprintf(stderr,"Using %d Snow Bands\n",vic411_options.SNOW_BAND);
  fprintf(stderr,"Using %d Root Zones\n",vic411_options.ROOT_ZONES);
  if ( vic411_options.SAVE_STATE )
    fprintf(stderr,"Model state will be saved on = %02i/%02i/%04i\n\n",
	    global.stateday, global.statemonth, global.stateyear);
  if ( vic411_options.BINARY_OUTPUT ) 
    fprintf(stderr,"Model output is in standard BINARY format.\n");
  else 
    fprintf(stderr,"Model output is in standard ASCII format.\n");
  if ( LINK_DEBUG ) 
    fprintf(stderr,"Debugging code has been included in the executable.\n");
  else 
    fprintf(stderr,"Debugging code has not been compiled.\n");
#endif

#endif // !OUTPUT_FORCE

  return global;

}
