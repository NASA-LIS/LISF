#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>
#include <string.h>

static char vcid[] = "$Id: vic411_read_vegparam.c,v 4.5.2.11 2009/10/01 21:12:16 vicadmin Exp $";

vic411_veg_con_struct *vic411_read_vegparam(FILE *vegparam,
                              int   gridcel,
                              int   Nveg_type)
/**********************************************************************
  vic411_read_vegparam.c    Keith Cherkauer and Dag Lohmann       1997

  This routine reads in vegetation parameters for the current grid cell.
  It also relates each vegetation class in the cell to the appropriate
  parameters in the vegetation library.

  Modifications:
  09-24-98  Modified to read root zone distribution information so
           that soil layer root fractions can be computed for new 
	   soil layer depths - see vic411_calc_root_fractions.c           KAC
  07-15-99 Modified to read LAI values from a new line in the vegetation
           parameter file.  Added specifically to work with the new
	   global LAI files.
  11-18-02 Added code to read in blowing snow parameters.          LCB
  03-27-03 Modified code to update Wdmax based on LAI values read in
           for the current grid cell.  If LAI is not obtained from this
           function, then the values cacluated in vic411_read_veglib.c are
           left unchanged.						DP & KAC
  2006-Nov-07 Allocates MaxVeg+1 veg tiles.				TJB
  2007-May-11 Changed some 'fscanf' statements to 'fgets' and 'sscanf' 
	      to count rootzone and BLOWING fields. Also tests for
	      fetch < 1.						GCT
  2007-Oct-31 Added missing brackets in if(vic411_options.GLOBAL_LAI) block.	TJB
  2008-Oct-23 Added blocks to free vegarr[].				LCB via TJB
  2009-Jan-16 Added logic for COMPUTE_TREELINE option.			TJB
  2009-Jun-09 Modified to use extension of vic411_veg_lib structure to contain
	      bare soil information.					TJB
  2009-Jun-17 Modified to understand both tabs and spaces as delimiters.TJB
  2009-Jun-17 Fixed incorrect placement of free vegarr[] for case of
	      GLOBAL_LAI==FALSE.					TJB
  2009-Jul-26 Allocate extra veg tile for COMPUTE_TREELINE case.	TJB
  2009-Jul-31 Removed extra veg tile for lake/wetland case.		TJB
  2009-Sep-14 Made error messages clearer.				TJB
  2009-Oct-01 Added error message for case of LAI==0 and overstory==1.	TJB
**********************************************************************/
{

  void vic411_ttrim( char *string );
  extern vic411_veg_lib_struct *vic411_veg_lib;
  extern vic411_option_struct   vic411_options;
#if LINK_DEBUG
  extern vic411_debug_struct    debug;
#endif

  vic411_veg_con_struct *temp;
  int             vegcel, i, j, k, vegetat_type_num, skip, veg_class;
  int             MaxVeg;
  int             Nfields, NfieldsMax;
  int             NoOverstory;
  float           depth_sum;
  float           sum;
  char            str[500];
  char            ErrStr[MAXSTRING];
  char            line[MAXSTRING];
  char            tmpline[MAXSTRING];
  const char      delimiters[] = " \t";
  char            *token;
  char            *vegarr[500];
  size_t	  length;

  if(vic411_options.GLOBAL_LAI) skip=2;
  else skip=1;

  NoOverstory = 0;

#if !NO_REWIND
  rewind(vegparam);
#endif  

  while ( ( fscanf(vegparam, "%d %d", &vegcel, &vegetat_type_num) == 2 ) && vegcel != gridcel ){
    if (vegetat_type_num < 0) {
      sprintf(ErrStr,"ERROR number of vegetation tiles (%i) given for cell %i is < 0.\n",vegetat_type_num,vegcel);
      vic411_nrerror(ErrStr);
    }
    for (i = 0; i <= vegetat_type_num * skip; i++){
      if ( fgets(str, 500, vegparam) == NULL ){
        sprintf(ErrStr,"ERROR unexpected EOF for cell %i while reading root zones and LAI\n",vegcel);
        vic411_nrerror(ErrStr);
      }
    }
  }
  fgets(str, 500, vegparam); // read newline at end of veg class line to advance to next line
  if (vegcel != gridcel) {
    fprintf(stderr, "vic411_Error in vegetation file.  Grid cell %d not found\n",
            gridcel);
    exit(99);
  }
  if(vegetat_type_num >= MAX_VEG) {
    sprintf(ErrStr,"Vegetation parameter file wants more vegetation tiles in grid cell %i (%i) than are allowed by MAX_VEG (%i) [NOTE: bare soil class is assumed].  Edit vic411_vicNl_def.h and recompile.",gridcel,vegetat_type_num+1,MAX_VEG);
    vic411_nrerror(ErrStr);
  }

  // Make sure to allocate extra memory for bare soil tile
  // and optionally an above-treeline veg tile
  MaxVeg = vegetat_type_num+1;
  if ( vic411_options.AboveTreelineVeg >= 0 )
    MaxVeg++;

  /** Allocate memory for vegetation grid cell parameters **/
  temp = (vic411_veg_con_struct*) calloc( MaxVeg, sizeof(vic411_veg_con_struct));
  temp[0].Cv_sum = 0.0;

  for (i = 0; i < vegetat_type_num; i++) {
    temp[i].zone_depth = calloc(vic411_options.ROOT_ZONES,sizeof(float));
    temp[i].zone_fract = calloc(vic411_options.ROOT_ZONES,sizeof(float));
    temp[i].vegetat_type_num = vegetat_type_num;

    // Read the root zones line
    if ( fgets( line, MAXSTRING, vegparam ) == NULL ){
      sprintf(ErrStr,"ERROR unexpected EOF for cell %i while reading vegetat_type_num %d\n",vegcel,vegetat_type_num);
      vic411_nrerror(ErrStr);
    }
    strcpy(tmpline, line);
    vic411_ttrim( tmpline );
    token = strtok (tmpline, delimiters);    /*  token => veg_class, move 'line' pointer to next field */  
    Nfields = 0;
    vegarr[Nfields] = calloc( 500, sizeof(char));
    strcpy(vegarr[Nfields],token);
    Nfields++;

    token = strtok (NULL, delimiters);
    while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
    while ( token != NULL ) {
      vegarr[Nfields] = calloc( 500, sizeof(char));      
      strcpy(vegarr[Nfields],token);
      Nfields++;
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
    }

    NfieldsMax = 2 + 2 * vic411_options.ROOT_ZONES;  /* Number of expected fields this line */
    if( vic411_options.BLOWING ){
      NfieldsMax += 3;
    }
    if ( Nfields != NfieldsMax ) {
      sprintf(ErrStr,"ERROR - cell %d - expecting %d fields but found %d in veg line %s\n",gridcel,NfieldsMax, Nfields, line);
      vic411_nrerror(ErrStr);
    }

    temp[i].LAKE = 0;
    temp[i].veg_class = atoi( vegarr[0] );
    temp[i].Cv = atof( vegarr[1] );
    depth_sum = 0;
    sum = 0.;
    for(j=0;j<vic411_options.ROOT_ZONES;j++) {
      temp[i].zone_depth[j] = atof( vegarr[2 + j*2] );
      temp[i].zone_fract[j] = atof( vegarr[3 + j*2] );
      depth_sum += temp[i].zone_depth[j];
      sum += temp[i].zone_fract[j];
    }
    if(depth_sum <= 0) {
      sprintf(str,"Root zone depths must sum to a value greater than 0.");
      vic411_nrerror(str);
    }
    if(sum != 1.) {
      fprintf(stderr,"WARNING: Root zone fractions sum to more than 1 ( = %f), normalizing fractions.  If the sum is large, check that your vegetation parameter file is in the form - <zone 1 depth> <zone 1 fract> <zone 2 depth> <zone 2 fract> ...\n", sum);
      for(j=0;j<vic411_options.ROOT_ZONES;j++) {
	temp[i].zone_fract[j] /= sum;
      }
    }

    if(vic411_options.BLOWING) {
      j = 2 * vic411_options.ROOT_ZONES;
      temp[i].sigma_slope = atof( vegarr[2 + j] );
      temp[i].lag_one = atof( vegarr[3 + j] );
      temp[i].fetch = atof( vegarr[4 + j]) ;
      if( temp[i].sigma_slope <= 0. || temp[i].lag_one <= 0.) {
        sprintf(str,"Deviation of terrain slope must be greater than 0.");
        vic411_nrerror(str);
      }
      if( temp[i].fetch < 1.0  ) {
	sprintf(str,"ERROR - BLOWING parameter fetch should be >> 1 but cell %i has fetch = %.2f\n", gridcel, temp[i].fetch );
        vic411_nrerror(str);
      }
    }

    veg_class = MISSING;
    for(j=0;j<Nveg_type;j++)
      if(temp[i].veg_class == vic411_veg_lib[j].veg_class)
	veg_class = j;
    if(veg_class == MISSING) {
      sprintf(ErrStr,"The vegetation class id %i in vegetation tile %i from cell %i is not defined in the vegetation library file.", temp[i].veg_class, i, gridcel);
      vic411_nrerror(ErrStr);
    }
    else
      temp[i].veg_class = veg_class;

    temp[0].Cv_sum += temp[i].Cv;

    for(k=0; k<Nfields; k++)
      free(vegarr[k]);

    if ( vic411_options.GLOBAL_LAI ) {
      // Read the LAI line
      if ( fgets( line, MAXSTRING, vegparam ) == NULL ){
        sprintf(ErrStr,"ERROR unexpected EOF for cell %i while reading LAI for vegetat_type_num %d\n",vegcel,vegetat_type_num);
        vic411_nrerror(ErrStr);
      }      
      Nfields = 0;
      vegarr[Nfields] = calloc( 500, sizeof(char));      
      strcpy(tmpline, line);
      vic411_ttrim( tmpline );
      token = strtok (tmpline, delimiters); 
      strcpy(vegarr[Nfields],token);
      Nfields++;
 
      while( ( token = strtok (NULL, delimiters)) != NULL ){
        vegarr[Nfields] = calloc( 500, sizeof(char));      
        strcpy(vegarr[Nfields],token);
        Nfields++;
      }
      NfieldsMax = 12; /* For LAI */
      if ( Nfields != NfieldsMax ) {
        sprintf(ErrStr,"ERROR - cell %d - expecting %d LAI values but found %d in line %s\n",gridcel, NfieldsMax, Nfields, line);
        vic411_nrerror(ErrStr);
      }

      for ( j = 0; j < 12; j++ ) {
        vic411_veg_lib[temp[i].veg_class].LAI[j] = atof( vegarr[j] );
        if (vic411_veg_lib[temp[i].veg_class].overstory && vic411_veg_lib[temp[i].veg_class].LAI[j] == 0) {
          sprintf(ErrStr,"ERROR: cell %d, veg tile %d: the specified veg class (%d) is listed as an overstory class in the veg LIBRARY, but the LAI given in the veg PARAM FILE for this tile for month %d is 0.\n",gridcel, i+1, temp[i].veg_class+1, j+1);
          vic411_nrerror(ErrStr);
        }
        vic411_veg_lib[temp[i].veg_class].Wdmax[j] = 
	  LAI_WATER_FACTOR * vic411_veg_lib[temp[i].veg_class].LAI[j];
      }
      for(k=0; k<Nfields; k++)
        free(vegarr[k]);
    }

    // Determine if cell contains non-overstory vegetation
    if ( vic411_options.COMPUTE_TREELINE && !vic411_veg_lib[temp[i].veg_class].overstory )
      NoOverstory++;

  }

  // Determine if we have bare soil
  if(temp[0].Cv_sum>1.0){
    fprintf(stderr,"WARNING: Cv exceeds 1.0 at grid cell %d, fractions being adjusted to equal 1\n", gridcel);
    for(j=0;j<vegetat_type_num;j++)
      temp[j].Cv = temp[j].Cv / temp[0].Cv_sum;
    temp[0].Cv_sum = 1.;
  }
  else if(temp[0].Cv_sum>0.99 && temp[0].Cv_sum<1.0){
    fprintf(stderr,"WARNING: Cv > 0.99 and Cv < 1.0 at grid cell %d, model assuming that bare soil is not to be run - fractions being adjusted to equal 1\n", gridcel);
    for(j=0;j<vegetat_type_num;j++)
      temp[j].Cv = temp[j].Cv / temp[0].Cv_sum;
    temp[0].Cv_sum = 1.;
  }

  // Handle veg above the treeline
  if ( vic411_options.SNOW_BAND > 1 && vic411_options.COMPUTE_TREELINE
       && ( !NoOverstory && temp[0].Cv_sum == 1. ) ) {

    // All vegetation in the current cell is defined with overstory.
    // Add default non-overstory vegetation so that snow bands above treeline
    // can be sucessfully simulated.

    if ( vic411_options.AboveTreelineVeg < 0 ) {

      // Above treeline snowband should be treated as bare soil
      for ( j = 0; j < vegetat_type_num; j++ )
        temp[j].Cv -= ( 0.001 / (float)vegetat_type_num );
      temp[0].Cv_sum -= 0.001;

    }
    else {

      // Above treeline snowband should use the defined vegetation
      // add vegetation to typenum
      // check that veg type exists in library and does not have overstory
      if(vegetat_type_num > 0) {

        for ( j = 0; j < vegetat_type_num; j++ ) {
          temp[j].Cv -= ( 0.001 / (float)vegetat_type_num );
          temp[j].vegetat_type_num++;
        }

        temp[vegetat_type_num].Cv         = 0.001;
        temp[vegetat_type_num].veg_class  = vic411_options.AboveTreelineVeg;
        temp[vegetat_type_num].Cv_sum     = temp[vegetat_type_num-1].Cv_sum;
        temp[vegetat_type_num].zone_depth = calloc( vic411_options.ROOT_ZONES,
                                                  sizeof(float));
        temp[vegetat_type_num].zone_fract = calloc( vic411_options.ROOT_ZONES,
                                                  sizeof(float));
        temp[vegetat_type_num].vegetat_type_num = vegetat_type_num+1;

        // Since root zones are not defined they are copied from the last
        // vegetation type.
        for ( j = 0; j < vic411_options.ROOT_ZONES; j++ ) {
          temp[vegetat_type_num].zone_depth[j]
            = temp[vegetat_type_num-1].zone_depth[j];
          temp[vegetat_type_num].zone_fract[j]
            = temp[vegetat_type_num-1].zone_fract[j];
        }

      }

      // Identify current vegetation class
      veg_class = MISSING;
      for ( j = 0; j < Nveg_type; j++ ) {
        if(temp[vegetat_type_num].veg_class == vic411_veg_lib[j].veg_class) {
          veg_class = j;
          break;
        }
      }
      if ( veg_class == MISSING ) {
        sprintf(ErrStr,"The vegetation class id %i defined for above-treeline from cell %i is not defined in the vegetation library file.", temp[vegetat_type_num].veg_class, gridcel);
        vic411_nrerror(ErrStr);
      }
      else {
        temp[vegetat_type_num].veg_class = veg_class;
      }

      if ( vic411_veg_lib[veg_class].overstory ) {
        sprintf(ErrStr,"Vegetation class %i is defined to have overstory, so it cannot be used as the default vegetation type for above canopy snow bands.", vic411_veg_lib[veg_class].veg_class );
        vic411_nrerror(ErrStr);
      }

    }

  }

  // Bare soil tile
  if (temp[0].Cv_sum < 1.) {
    j = vegetat_type_num;
    temp[j].veg_class = Nveg_type; // Create a veg_class ID for bare soil, which is not mentioned in the veg library
    temp[j].Cv = 1.0 - temp[0].Cv_sum;
    // Don't allocate any root-zone-related arrays
    if(vic411_options.BLOWING) {
      if (vegetat_type_num > 0) {
        temp[j].sigma_slope = temp[0].sigma_slope;
        temp[j].lag_one = temp[0].lag_one;
        temp[j].fetch = temp[0].fetch;
      }
      else {
        temp[j].sigma_slope = 0.005;
        temp[j].lag_one = 0.95;
        temp[j].fetch = 2000;
      }
    }
  }

  return temp;
} 





#if 0
/* trim trailing newlines */

#define END '\0'
#define NEW '\n'

void vic411_ttrim( char *c ) 
{
  while( (*c++ != END) );
    --c;
  for( ; *--c == NEW; *c = END );

}
#endif
