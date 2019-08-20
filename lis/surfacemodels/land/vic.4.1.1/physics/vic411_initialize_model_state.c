#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>
#include <string.h>

static char vcid[] = "$Id: vic411_initialize_model_state.c,v 5.8.2.18 2009/09/28 22:12:48 vicadmin Exp $";

int vic411_initialize_model_state(vic411_dist_prcp_struct    *prcp,
			   vic411_dmy_struct           dmy,
			   vic411_global_param_struct *vic411_global_param,
			   vic411_filep_struct         filep,
			   int                  cellnum,
			   int                  Nveg,
			   int                  Nnodes,
			   int                  Ndist,
			   double               surf_temp, 
			   vic411_soil_con_struct     *soil_con,
			   vic411_veg_con_struct      *veg_con,
			   vic411_lake_con_struct      lake_con,
			   char               **init_STILL_STORM,
			   int                **init_DRY_TIME,
			   vic411_save_data_struct    *save_data)
/**********************************************************************
  vic411_initialize_model_state      Keith Cherkauer	    April 17, 2000

  This routine initializes the model state (energy balance, water balance,
  and snow components).  If a state file is provided to the model than its
  contents are checked to see if it agrees with the current simulation
  set-up, if so it is used to initialize the model state.  If no state
  file is provided the model initializes all variables with defaults and
  the user should expect to throw out the beginning of the simulation 
  period as model start-up.

  UNITS: (m, s, kg, C, moisture in mm) unless otherwise specified

  Modifications:
  4-17-00 Modified from initialize_energy_bal.c and vic411_initialize_snow.c
          to provide a single controlling routine for initializing the
          model state.
  9-00    Fixed bug where initialization of soil node temperatures 
          and moitures was within two vegetation loops, thus only
          the first vegetation type was properly initialized.     KAC
  2-19-03 Modified to initialize soil and vegetation parameters for
          the dry grid cell fraction, if distributed precipitation
          is activated.                                           KAC
  11-18-02 Modified to initialize lake and wetland algorithms 
          variables.                                              LCB
  2-10-03 Fixed looping problem with initialization of soil moisture. KAC
  3-12-03 Modified so that soil layer ice content is only calculated 
          when frozen soil is implemented and active in the current 
          grid cell.                                                KAC
  04-10-03 Modified to read storm parameters from model state file.  KAC
  04-25-03 Modified to work with vegetation type specific storm 
           parameters.                                              KAC
  07-May-04 Initialize soil_con->dz_node[Nnodes] to 0.0, since it is
	    accessed in vic411_set_node_parameters().				TJB
  01-Nov-04 Added support for state files containing SPATIAL_FROST
	    and LAKE_MODEL state variables.				TJB
  2006-Apr-21 Replaced Cv (uninitialized) with lake_con.Cl[0] in
	      surfstor calculation.					TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
              save_data structure to track changes in moisture storage
              over each time step; this needs initialization here.	TJB
  2006-Oct-10 Added snow[veg][band].snow_canopy to save_data.swe.	TJB
  2006-Oct-16 Merged infiles and outfiles structs into vic411_filep_struct;
	      This included removing the unused init_snow file.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2007-Apr-24 Added EXP_TRANS option.					JCA
  2007-Apr-24 Zsum_node loaded into soil_con structure for later use
	      without having to recalculate.				JCA
  2007-Aug-09 Added features for EXCESS_ICE option.			JCA
  2007-Aug-21 Return value of ErrorFlag if error in
	      vic411_distribute_node_moisture_properties.			JCA
  2007-Sep-18 Check for soil moist exceeding max moist moved from
	      vic411_read_initial_model_state to here.				JCA
  2007-Oct-24 Modified vic411_initialize_lake() to return ErrorFlag.		TJB
  2008-Mar-01 Reinserted missing logic for QUICK_FS in calls to
	      vic411_distribute_node_moisture_properties() and
	      vic411_estimate_layer_ice_content().				TJB
  2009-Feb-09 Removed dz_node from call to
	      vic411_distribute_node_moisture_properties.			KAC via TJB
  2009-Feb-09 Removed dz_node from call to find_0_degree_front.		KAC via TJB
  2009-Mar-15 Modified to not call vic411_estimate_layer_ice_content() if
	      not modeling frozen soil.					KAC via TJB
  2009-Mar-16 Added resid_moist to argument list of
	      vic411_estimate_layer_ice_content().  This allows computation
	      of min_liq, the minimum allowable liquid water content
	      in each layer as a function of temperature.		TJB
  2009-Jun-09 Modified to use extension of vic411_veg_lib structure to contain
	      bare soil information.					TJB
  2009-Jul-26 Added initial estimate of incoming longwave at surface
	      (LongUnderOut) for use in canopy snow T iteration.	TJB
  2009-Jul-31 Removed extra lake/wetland veg tile.			TJB
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.		TJB
  2009-Sep-19 Made initialization of Tfoliage more accurate for snow bands.	TJB
  2009-Sep-28 Added initialization of energy structure.			TJB
**********************************************************************/
{
  extern vic411_option_struct vic411_options;
  extern vic411_veg_lib_struct *vic411_veg_lib;
#if LINK_DEBUG
  extern vic411_debug_struct debug;
#endif
#if QUICK_FS
  extern double temps[];
#endif

  char     ErrStr[MAXSTRING];
  char     FIRST_VEG;
  int      i, j, ii, veg, index, dist;
  int      lidx;
  int      dry;
  int      band;
#if SPATIAL_FROST
  int      frost_area;
#endif
  int      ErrorFlag;
  double   Cv;
  double   Zsum, dp;
  double   tmpdp, tmpadj, Bexp;
  double   Tair;
  double   tmp;
  double  *M;
  double   moist[MAX_VEG][MAX_BANDS][MAX_LAYERS];
#if SPATIAL_FROST
  double   ice[MAX_VEG][MAX_BANDS][MAX_LAYERS][FROST_SUBAREAS];
#else
  double   ice[MAX_VEG][MAX_BANDS][MAX_LAYERS];
#endif // SPATIAL_FROST
#if QUICK_FS
  double   Aufwc, Bufwc;
#endif
  double   Clake;
  double   mu;
  double   TreeAdjustFactor[MAX_BANDS];
#if EXCESS_ICE
  double   sum_mindepth, sum_depth_pre, sum_depth_post, tmp_mindepth;
#endif

  vic411_cell_data_struct     ***cell;
  vic411_energy_bal_struct     **energy;
  vic411_lake_var_struct        *lake_var;
  vic411_snow_data_struct      **snow;
  vic411_veg_var_struct       ***veg_var;

  cell    = prcp->cell;
  energy  = prcp->energy;
  lake_var = &prcp->lake_var;
  snow    = prcp->snow;
  veg_var = prcp->veg_var;

  // Initialize soil depths
  dp = soil_con->dp;

  FIRST_VEG = TRUE;

  // increase initial soil surface temperature if air is very cold
  Tair = surf_temp;
  if ( surf_temp < -1. ) surf_temp = -1.;
  
  // initialize storm parameters to start a new simulation
  (*init_STILL_STORM) = (char *)malloc((Nveg+1)*sizeof(char));
  (*init_DRY_TIME)    = (int *)malloc((Nveg+1)*sizeof(int));
  for ( veg = 0 ; veg <= Nveg ; veg++ )
    (*init_DRY_TIME)[veg] = -999;
  
  /********************************************
    Initialize all snow pack variables 
    - some may be reset if state file present
  ********************************************/

  vic411_initialize_snow(snow, Nveg, cellnum);

  /********************************************
    Initialize all soil layer variables 
    - some may be reset if state file present
  ********************************************/

  vic411_initialize_soil(cell[WET], soil_con, veg_con, Nveg);
  if ( vic411_options.DIST_PRCP )
    vic411_initialize_soil(cell[DRY], soil_con, veg_con, Nveg);

  /********************************************
    Initialize all vegetation variables 
    - some may be reset if state file present
  ********************************************/

  vic411_initialize_veg(veg_var[WET], veg_con, vic411_global_param, Nveg);
  if ( vic411_options.DIST_PRCP )
    vic411_initialize_veg(veg_var[DRY], veg_con, vic411_global_param, Nveg);

  /********************************************
    Initialize all lake variables 
  ********************************************/

  if ( vic411_options.LAKES && lake_con.Cl[0] > 0) {
    ErrorFlag = vic411_initialize_lake(lake_var, lake_con, soil_con, surf_temp);
    if (ErrorFlag == ERROR) return(ErrorFlag);
  }

  /********************************************
    Initialize all spatial frost variables 
  ********************************************/

#if SPATIAL_FROST
  for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
    if ( FROST_SUBAREAS == 1 ) soil_con->frost_fract[frost_area] = 1.;
    else if (FROST_SUBAREAS == 2 ) soil_con->frost_fract[frost_area] = 0.5;
    else {
      soil_con->frost_fract[frost_area] = 1. / (FROST_SUBAREAS - 1);
      if ( frost_area == 0 || frost_area == FROST_SUBAREAS-1 ) 
	soil_con->frost_fract[frost_area] /= 2.;
    }
  }
#endif // SPATIAL_FROST

  /********************************************************
    Compute grid cell fractions for all subareas used in 
    spatial distribution of soil frost routines.
  ********************************************************/

#if QUICK_FS
  if(vic411_options.FROZEN_SOIL) {

    /***********************************************************
      Prepare table of maximum unfrozen water content values
      - This linearizes the equation for maximum unfrozen water
        content, reducing computation time for the frozen soil
        model.
    ***********************************************************/

    for(lidx=0;lidx<vic411_options.Nlayer;lidx++) { 
      for(ii=0;ii<QUICK_FS_TEMPS;ii++) {
	Aufwc = vic411_maximum_unfrozen_water(temps[ii], 1.0, 
				       soil_con->bubble[lidx], 
				       soil_con->expt[lidx]);
	Bufwc = vic411_maximum_unfrozen_water(temps[ii+1], 1.0, 
				       soil_con->bubble[lidx], 
				       soil_con->expt[lidx]);
	soil_con->ufwc_table_layer[lidx][ii][0] 
	  = vic411_linear_interp(0., temps[ii], temps[ii+1], Aufwc, Bufwc);
	soil_con->ufwc_table_layer[lidx][ii][1] 
	  = (Bufwc - Aufwc) / (temps[ii+1] - temps[ii]);
      }
    }
  }  
#endif // QUICK_FS

  /************************************************************************
    CASE 1: Not using quick ground heat flux, and initial conditions files 
    provided
  ************************************************************************/

  if(vic411_options.INIT_STATE) {

#if EXCESS_ICE
    sum_mindepth = 0;
    sum_depth_pre = 0;
    for( lidx = 0; lidx < vic411_options.Nlayer; lidx++ ){
      tmp_mindepth = (float)(int)(soil_con->min_depth[lidx] * 1000 + 0.5) / 1000;	
      sum_mindepth += tmp_mindepth;
      sum_depth_pre += soil_con->depth[lidx];
    }
#endif

    vic411_read_initial_model_state(filep.init_state, prcp, vic411_global_param,  
			     Nveg, vic411_options.SNOW_BAND, cellnum, soil_con,
			     Ndist, *init_STILL_STORM, *init_DRY_TIME, lake_con);



#if EXCESS_ICE
    // calculate dynamic soil and veg properties if excess_ice is present
    sum_depth_post = 0;
    for( lidx = 0; lidx < vic411_options.Nlayer; lidx++ )
      sum_depth_post += soil_con->depth[lidx];
    if( sum_depth_post != sum_depth_pre) {
      /*update soil_con properties*/
      for( lidx = 0; lidx < vic411_options.Nlayer; lidx++ ) {
	soil_con->bulk_density[lidx] = (1.0-soil_con->effective_porosity[lidx])*soil_con->soil_density[lidx]; 
	soil_con->max_moist[lidx] = soil_con->depth[lidx] * soil_con->effective_porosity[lidx] * 1000.;	
      } //loop for each soil layer      
      
      /********update remaining soil_con properties**********/
      /* update Maximum Infiltration for Upper Layers */
      if(vic411_options.Nlayer==2)
	soil_con->max_infil = (1.0+soil_con->b_infilt)*soil_con->max_moist[0];
      else
	soil_con->max_infil = (1.0+soil_con->b_infilt)*(soil_con->max_moist[0]+soil_con->max_moist[1]);
      
      /* Soil Layer Critical and Wilting Point Moisture Contents */
      for(lidx=0;lidx<vic411_options.Nlayer;lidx++) {//soil layer
	soil_con->Wcr[lidx]  = soil_con->Wcr_FRACT[lidx] * soil_con->max_moist[lidx];
	soil_con->Wpwp[lidx] = soil_con->Wpwp_FRACT[lidx] * soil_con->max_moist[lidx];
	if(soil_con->Wpwp[lidx] > soil_con->Wcr[lidx]) {
	  sprintf(ErrStr,"Updated wilting point moisture (%f mm) is greater than updated critical point moisture (%f mm) for layer %d.\n\tIn the soil parameter file, Wpwp_FRACT MUST be <= Wcr_FRACT.\n",
		  soil_con->Wpwp[lidx], soil_con->Wcr[lidx], lidx);
	  vic411_nrerror(ErrStr);
	}
	if(soil_con->Wpwp[lidx] < soil_con->resid_moist[lidx] * soil_con->depth[lidx] * 1000.) {
	  sprintf(ErrStr,"Updated wilting point moisture (%f mm) is less than updated residual moisture (%f mm) for layer %d.\n\tIn the soil parameter file, Wpwp_FRACT MUST be >= resid_moist / (1.0 - bulk_density/soil_density).\n",
		  soil_con->Wpwp[lidx], soil_con->resid_moist[lidx] * soil_con->depth[lidx] * 1000., lidx);
	  vic411_nrerror(ErrStr);
	}
      }      
      
      /* If BASEFLOW = NIJSSEN2001 then convert ARNO baseflow
	 parameters d1, d2, d3, and d4 to Ds, Dsmax, Ws, and c */
      if(vic411_options.BASEFLOW == NIJSSEN2001) {
	lidx = vic411_options.Nlayer-1;
	soil_con->Dsmax = soil_con->Dsmax_orig * 
	  pow((double)(1./(soil_con->max_moist[lidx]-soil_con->Ws_orig)), -soil_con->c) +
	  soil_con->Ds_orig * soil_con->max_moist[lidx];
	soil_con->Ds = soil_con->Ds_orig * soil_con->Ws_orig / soil_con->Dsmax_orig;
	soil_con->Ws = soil_con->Ws_orig/soil_con->max_moist[lidx];
      }
      
      /*********** update root fractions ***************/
      vic411_calc_root_fractions(veg_con, soil_con);
      
#if VERBOSE
      /* write changes to screen */
      fprintf(stderr,"Soil properties initialized from state file:\n");
      for(lidx=0;lidx<vic411_options.Nlayer;lidx++) {//soil layer
	fprintf(stderr,"\tFor layer %d:\n",lidx+1);
	fprintf(stderr,"\t\tDepth of soil layer = %.2f m.\n",soil_con->depth[lidx]);
	fprintf(stderr,"\t\tEffective porosity = %.2f.\n",soil_con->effective_porosity[lidx]);
	fprintf(stderr,"\t\tBulk density = %.2f kg/m^3.\n",soil_con->bulk_density[lidx]);
      }
      fprintf(stderr,"\tDamping depth = %.2f m.\n",soil_con->dp);
      if(sum_depth_post == sum_mindepth)
	fprintf(stderr,"\tExcess ice is no longer present in the soil column.\n");
#endif //VERBOSE

    }//updated initial conditions due to state file
#endif //EXCESS_ICE

    /******Check that soil moisture does not exceed maximum allowed************/
    for ( veg = 0 ; veg <= Nveg ; veg++ ) {
      for( band = 0; band < vic411_options.SNOW_BAND; band++ ) {
	for( lidx = 0; lidx < vic411_options.Nlayer; lidx++ ) {	  
	  for ( dist = 0; dist < Ndist; dist ++ ) {
	    if ( cell[dist][veg][band].layer[lidx].moist > soil_con->max_moist[lidx] ) {
	      fprintf( stderr, "WARNING: Maximum soil moisture exceeded in layer %d for veg type %d and snow band %d.  Value of %f reset to maximum (%f mm).\n", lidx, veg, band, cell[dist][veg][band].layer[lidx].moist, soil_con->max_moist[lidx] );
	      cell[dist][veg][band].layer[lidx].moist = soil_con->max_moist[lidx];
	    }	    
	  }
	}
      }      
    }

    /****** initialize moist and ice ************/
    for ( veg = 0 ; veg <= Nveg ; veg++ ) {
      // Initialize soil for existing vegetation types
      Cv = veg_con[veg].Cv;
      
      if ( Cv > 0 ) {
	for( band = 0; band < vic411_options.SNOW_BAND; band++ ) {
	  for( lidx = 0; lidx < vic411_options.Nlayer; lidx++ ) {
	    moist[veg][band][lidx] = cell[0][veg][band].layer[lidx].moist;

#if SPATIAL_FROST
	    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
	      ice[veg][band][lidx][frost_area] 
		= cell[0][veg][band].layer[lidx].ice[frost_area];
#else
	    ice[veg][band][lidx] = cell[0][veg][band].layer[lidx].ice;
#endif
	  }
	}
      }
    }

  }
  
  /************************************************************************
    CASE 2: Initialize soil if using quick heat flux, and no initial
    soil properties file given
  ************************************************************************/
    
  else if(vic411_options.QUICK_FLUX) {

    /* Initialize soil node thicknesses */
    soil_con->dz_node[0]   = soil_con->depth[0];
    soil_con->dz_node[1]   = soil_con->depth[0];
    soil_con->dz_node[2]   = 2. * (dp - 1.5 * soil_con->depth[0]);    
    soil_con->Zsum_node[0] = 0;
    soil_con->Zsum_node[1] = soil_con->depth[0];
    soil_con->Zsum_node[2] = dp;

    for ( veg = 0 ; veg <= Nveg ; veg++ ) {
      // Initialize soil for existing vegetation types
      Cv = veg_con[veg].Cv;
      
      if ( Cv > 0 ) {
	for( band = 0; band < vic411_options.SNOW_BAND; band++ ) {

	  /* Initialize soil node temperatures */
	  energy[veg][band].T[0] = surf_temp;
	  energy[veg][band].T[1] = surf_temp;
	  energy[veg][band].T[2] = soil_con->avg_temp;

	  /* Initialize soil layer thicknesses */
	  for ( lidx = 0; lidx < vic411_options.Nlayer; lidx++ ) {
	    moist[veg][band][lidx] = cell[0][veg][band].layer[lidx].moist;
#if SPATIAL_FROST
	    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
	      ice[veg][band][lidx][frost_area] = 0.;
#else
	    ice[veg][band][lidx] = 0.;
#endif
	  }
	}
      }
    }
  }

  /*****************************************************************
    CASE 3: Initialize Energy Balance Variables if not using quick
    ground heat flux, and no Initial Condition File Given 
  *****************************************************************/
  else if(!vic411_options.QUICK_FLUX) {
    for ( veg = 0 ; veg <= Nveg ; veg++ ) {
      // Initialize soil for existing vegetation types
      Cv = veg_con[veg].Cv;
      
      if ( Cv > 0 ) {
	for( band = 0; band < vic411_options.SNOW_BAND; band++ ) {
	  
	  if(!vic411_options.EXP_TRANS){  
	    /* Initialize soil node temperatures and thicknesses 
	       Nodes set at surface, the depth of the first layer,
	       twice the depth of the first layer, and at the
	       damping depth.  Extra nodes are placed equal distance
	       between the damping depth and twice the depth of the
	       first layer. */
	    
	    energy[veg][band].T[0] = surf_temp;
	    soil_con->dz_node[0] = soil_con->depth[0];
	    soil_con->dz_node[1] = soil_con->depth[0];
	    soil_con->dz_node[2] = soil_con->depth[0];
	    energy[veg][band].T[Nnodes-1] = soil_con->avg_temp;
	    energy[veg][band].T[1] = vic411_exp_interp(soil_con->depth[0], 0., dp, 
						surf_temp, soil_con->avg_temp);
	    energy[veg][band].T[2] = vic411_exp_interp(2. * soil_con->depth[0], 0., dp, 
						surf_temp, soil_con->avg_temp);
	    
	    soil_con->Zsum_node[0] = 0;
	    soil_con->Zsum_node[1] = soil_con[0].depth[0];
	    Zsum   = 2. * soil_con[0].depth[0];
	    soil_con->Zsum_node[2] = Zsum;
	    tmpdp  = dp - soil_con[0].depth[0] * 2.5;
	    tmpadj = 3.5;
	    for ( index = 3; index < Nnodes-1; index++ ) {
	      if ( FIRST_VEG ) {
		soil_con->dz_node[index] = tmpdp/(((double)Nnodes-tmpadj));
	      }
	      Zsum += (soil_con->dz_node[index]
		       +soil_con->dz_node[index-1])/2.;
	      soil_con->Zsum_node[index] = Zsum;
	      energy[veg][band].T[index] = vic411_exp_interp(Zsum,0.,soil_con[0].dp,
						      surf_temp,
						      soil_con[0].avg_temp);
	    }
	    if ( FIRST_VEG ) {
	      FIRST_VEG = FALSE;
	      soil_con->dz_node[Nnodes-1] = (dp - Zsum 
					     - soil_con->dz_node[Nnodes-2] 
					     / 2. ) * 2.;
	      Zsum += (soil_con->dz_node[Nnodes-2]
		       +soil_con->dz_node[Nnodes-1])/2.;
	      soil_con->Zsum_node[Nnodes-1] = Zsum;
	      if((int)(Zsum*1000+0.5) != (int)(dp*1000+0.5)) {
		sprintf(ErrStr,"Sum of thermal node thicknesses (%f) in vic411_initialize_model_state do not equal dp (%f), check initialization procedure",Zsum,dp);
		vic411_nrerror(ErrStr);
	      }
	    }
	  }
	  else{ /* exponential grid transformation, EXP_TRANS = TRUE*/
	    
	    /*calculate exponential function parameter */
	    if ( FIRST_VEG ) {
	      Bexp = logf(dp+1.)/(double)(Nnodes-1); //to force Zsum=dp at bottom node
	      for ( index = 0; index <= Nnodes-1; index++ )
		soil_con->Zsum_node[index] = expf(Bexp*index)-1.;
	    }	    
	    
	    //top node	  
	    index=0;
	    if ( FIRST_VEG )
	      soil_con->dz_node[index] = soil_con->Zsum_node[index+1]-soil_con->Zsum_node[index];
	    energy[veg][band].T[index] = surf_temp;
	    //middle nodes
	    for ( index = 1; index < Nnodes-1; index++ ) {
	      if ( FIRST_VEG ) {
		soil_con->dz_node[index] = (soil_con->Zsum_node[index+1]-soil_con->Zsum_node[index])/2.+(soil_con->Zsum_node[index]-soil_con->Zsum_node[index-1])/2.;
	      }
	      energy[veg][band].T[index] = vic411_exp_interp(soil_con->Zsum_node[index],0.,soil_con[0].dp,
						      surf_temp,soil_con[0].avg_temp);
	    }
	    //bottom node
	    index=Nnodes-1;
	    if ( FIRST_VEG )
	      soil_con->dz_node[index] = soil_con->Zsum_node[index]-soil_con->Zsum_node[index-1];
	    energy[veg][band].T[index] = soil_con[0].avg_temp;
	  }
	  
	  //initialize moisture and ice for each soil layer
	  for ( lidx = 0; lidx < vic411_options.Nlayer; lidx++ ) {
	    moist[veg][band][lidx] = cell[0][veg][band].layer[lidx].moist;
#if SPATIAL_FROST
	    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
	      ice[veg][band][lidx][frost_area] = 0.;
#else
	    ice[veg][band][lidx] = 0.;
#endif
	  }
	}
      }
    }
  }

  /*********************************
    CASE 4: Unknown option
  *********************************/
  else {
    for ( veg = 0 ; veg <= Nveg ; veg++ ) {
      // Initialize soil for existing vegetation types
      Cv = veg_con[veg].Cv;

      if ( Cv > 0 ) {
	for( band = 0; band < vic411_options.SNOW_BAND; band++ ) {
	  // Initialize soil for existing snow elevation bands
	  if ( soil_con->AreaFract[band] > 0. ) {	  
	    for ( index = 0; index < vic411_options.Nlayer; index++ ) {
	      soil_con->dz_node[index] = 1.;
	    }
	  }
	}
      }
    }
  }

  /********************************************
    Initialize subsidence 
  ********************************************/

#if EXCESS_ICE
  for ( lidx = 0; lidx < vic411_options.Nlayer; lidx++ ) 
    soil_con->subsidence[lidx] = 0.0;
    
#endif // EXCESS_ICE

  /******************************************
    Initialize soil thermal node properties 
  ******************************************/

  if ( vic411_options.GRND_FLUX ) {
    
    FIRST_VEG = TRUE;
    for ( veg = 0 ; veg <= Nveg ; veg++) {
      // Initialize soil for existing vegetation types
      Cv = veg_con[veg].Cv;

      if ( Cv > 0 ) {
	for( band = 0; band < vic411_options.SNOW_BAND; band++ ) {
	  // Initialize soil for existing snow elevation bands
	  if ( soil_con->AreaFract[band] > 0. ) {
	    
	    /** Set soil properties for all soil nodes **/
	    if(FIRST_VEG) {
	      FIRST_VEG = FALSE;
	      vic411_set_node_parameters(soil_con->dz_node, soil_con->Zsum_node, soil_con->max_moist_node,
				  soil_con->expt_node, soil_con->bubble_node,
				  soil_con->alpha, soil_con->beta,
				  soil_con->gamma, soil_con->depth,
				  soil_con->max_moist, soil_con->expt, 
				  soil_con->bubble, soil_con->quartz, 
#if QUICK_FS
				  soil_con->ufwc_table_node,
#endif // QUICK_FS
#if EXCESS_ICE
				  soil_con->porosity, soil_con->effective_porosity,
				  soil_con->porosity_node, soil_con->effective_porosity_node,
#endif // EXCESS_ICE
				  Nnodes, vic411_options.Nlayer, soil_con->FS_ACTIVE);	  
	    }
	
	    /* set soil moisture properties for all soil thermal nodes */
	    ErrorFlag = vic411_distribute_node_moisture_properties(energy[veg][band].moist,
						energy[veg][band].ice,
						energy[veg][band].kappa_node,
						energy[veg][band].Cs_node,
						soil_con->Zsum_node,
						energy[veg][band].T,
						soil_con->max_moist_node,
#if QUICK_FS
						soil_con->ufwc_table_node,
#else
						soil_con->expt_node,
						soil_con->bubble_node,
#endif // QUICK_FS
#if EXCESS_ICE
						soil_con->porosity_node,
						soil_con->effective_porosity_node,
#endif // EXCESS_ICE
						moist[veg][band], 
						soil_con->depth,
						soil_con->soil_density,
						soil_con->bulk_density,
						soil_con->quartz,
						Nnodes, vic411_options.Nlayer,
						soil_con->FS_ACTIVE);
	    if ( ErrorFlag == ERROR )
	      return ( ErrorFlag );
	    
	    /* initialize layer moistures and ice contents */
	    for ( dry = 0; dry < Ndist; dry++ ) {
	      for ( lidx = 0; lidx < vic411_options.Nlayer; lidx++ ) {
		cell[dry][veg][band].layer[lidx].moist 
		  = moist[veg][band][lidx];
#if SPATIAL_FROST
		for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )

		  cell[dry][veg][band].layer[lidx].ice[frost_area] = ice[veg][band][lidx][frost_area];
#else
		cell[dry][veg][band].layer[lidx].ice = ice[veg][band][lidx];
#endif
	      }
	      if ( !( vic411_options.LAKES && veg == lake_con.lake_idx ) && vic411_options.FROZEN_SOIL && soil_con->FS_ACTIVE ) {
		ErrorFlag = vic411_estimate_layer_ice_content(cell[dry][veg][band].layer,
						       soil_con->Zsum_node,
						       energy[veg][band].T,
						       soil_con->max_moist_node,
#if QUICK_FS
						       soil_con->ufwc_table_node,
#else
						       soil_con->expt_node,
						       soil_con->bubble_node,
#endif // QUICK_FS
						       soil_con->depth,
						       soil_con->max_moist,
#if QUICK_FS
						       soil_con->ufwc_table_layer,
#else
						       soil_con->expt,
						       soil_con->bubble,
#endif // QUICK_FS
#if SPATIAL_FROST
						       soil_con->frost_fract, 
						       soil_con->frost_slope, 
#endif // SPATIAL_FROST
#if EXCESS_ICE
						       soil_con->porosity,
						       soil_con->effective_porosity,
#endif // EXCESS_ICE
						       soil_con->bulk_density,
						       soil_con->soil_density,
						       soil_con->quartz,
						       soil_con->resid_moist,
						       Nnodes, vic411_options.Nlayer, 
						       soil_con->FS_ACTIVE);
		
	      }
	    }
	    
	    /* Find freezing and thawing front depths */
	    if(!vic411_options.QUICK_FLUX && soil_con->FS_ACTIVE) 
	      vic411_find_0_degree_fronts(&energy[veg][band], soil_con->Zsum_node,
				   energy[veg][band].T, Nnodes);
	  }
	}
      }
    }	
  }

  // initialize miscellaneous energy balance terms
  for ( veg = 0 ; veg <= Nveg ; veg++) {
    for ( band = 0; band < vic411_options.SNOW_BAND; band++ ) {
      /* Set fluxes to 0 */
      energy[veg][band].advected_sensible = 0.0;
      energy[veg][band].advection         = 0.0;
      energy[veg][band].AtmosError        = 0.0;
      energy[veg][band].AtmosLatent       = 0.0;
      energy[veg][band].AtmosLatentSub    = 0.0;
      energy[veg][band].AtmosSensible     = 0.0;
      energy[veg][band].canopy_advection  = 0.0;
      energy[veg][band].canopy_latent     = 0.0;
      energy[veg][band].canopy_latent_sub = 0.0;
      energy[veg][band].canopy_refreeze   = 0.0;
      energy[veg][band].canopy_sensible   = 0.0;
      energy[veg][band].deltaCC           = 0.0;
      energy[veg][band].deltaH            = 0.0;
      energy[veg][band].error             = 0.0;
      energy[veg][band].fusion            = 0.0;
      energy[veg][band].grnd_flux         = 0.0;
      energy[veg][band].latent            = 0.0;
      energy[veg][band].latent_sub        = 0.0;
      energy[veg][band].longwave          = 0.0;
      energy[veg][band].LongOverIn        = 0.0;
      energy[veg][band].LongUnderIn       = 0.0;
      energy[veg][band].LongUnderOut      = 0.0;
      energy[veg][band].melt_energy       = 0.0;
      energy[veg][band].NetLongAtmos      = 0.0;
      energy[veg][band].NetLongOver       = 0.0;
      energy[veg][band].NetLongUnder      = 0.0;
      energy[veg][band].NetShortAtmos     = 0.0;
      energy[veg][band].NetShortGrnd      = 0.0;
      energy[veg][band].NetShortOver      = 0.0;
      energy[veg][band].NetShortUnder     = 0.0;
      energy[veg][band].out_long_canopy   = 0.0;
      energy[veg][band].out_long_surface  = 0.0;
      energy[veg][band].refreeze_energy   = 0.0;
      energy[veg][band].sensible          = 0.0;
      energy[veg][band].shortwave         = 0.0;
      energy[veg][band].ShortOverIn       = 0.0;
      energy[veg][band].ShortUnderIn      = 0.0;
      energy[veg][band].snow_flux         = 0.0;
      /* Initial estimate of LongUnderOut for use by vic411_snow_intercept() */
      tmp = energy[veg][band].T[0] + KELVIN;
      energy[veg][band].LongUnderOut = STEFAN_B * tmp * tmp * tmp * tmp;
      energy[veg][band].Tfoliage     = Tair + soil_con->Tfactor[band];
    }
  }

  // initialize Tfallback counters
  for ( veg = 0 ; veg <= Nveg ; veg++) {
    for ( band = 0; band < vic411_options.SNOW_BAND; band++ ) {
      energy[veg][band].Tfoliage_fbcount = 0;
      energy[veg][band].Tcanopy_fbcount = 0;
      energy[veg][band].Tsurf_fbcount = 0;
      for ( index = 0; index < Nnodes-1; index++ ) {
	energy[veg][band].T_fbcount[index] = 0;
      }
    }
  }

  // Compute treeline adjustment factors
  for ( band = 0; band < vic411_options.SNOW_BAND; band++ ) {
    if ( soil_con->AboveTreeLine[band] ) {
      Cv = 0;
      for ( veg = 0 ; veg < veg_con[0].vegetat_type_num ; veg++ ) {
        if ( vic411_veg_lib[veg_con[veg].veg_class].overstory )
          Cv += veg_con[veg].Cv;
      }
      TreeAdjustFactor[band] = 1. / ( 1. - Cv );
    }
    else TreeAdjustFactor[band] = 1.;
  }

  // Save initial moisture storages for differencing at end of time step
  save_data->total_soil_moist = 0;
  save_data->swe = 0;
  save_data->wdew = 0;
  save_data->surfstor = 0;
  for( veg = 0; veg <= Nveg; veg++ ) {

    Cv = veg_con[veg].Cv;

    if ( Cv > 0 ) {

      for ( dist = 0; dist < Ndist; dist++ ) {

        if(dist==0)
          mu = prcp[0].mu[veg];
        else
          mu = 1. - prcp[0].mu[veg];

        for( band = 0; band < vic411_options.SNOW_BAND; band++ ) {

          if (soil_con->AreaFract[band] > 0.
              && ( veg == veg_con[0].vegetat_type_num
                  || ( !soil_con->AboveTreeLine[band]
                      || (soil_con->AboveTreeLine[band] && !vic411_veg_lib[veg_con[veg].veg_class].overstory)))) {

            for(lidx=0;lidx<vic411_options.Nlayer;lidx++) {
              save_data->total_soil_moist +=
                (cell[dist][veg][band].layer[lidx].moist) // layer[].moist appears to contain both liquid and ice
                * Cv * mu * soil_con->AreaFract[band] * TreeAdjustFactor[band];
            }
            save_data->wdew += veg_var[dist][veg][band].Wdew * Cv * mu * soil_con->AreaFract[band] * TreeAdjustFactor[band];

          }

        }

      }

      for( band = 0; band < vic411_options.SNOW_BAND; band++ ) {
        save_data->swe += (snow[veg][band].swq + snow[veg][band].snow_canopy)
                         * Cv * soil_con->AreaFract[band] * TreeAdjustFactor[band];
      }

    }

  }
//<devel -- bug.patch>
// lake_con.basin is uninitialized.
//  Clake = lake_var->sarea/lake_con.basin[0];
//  save_data->surfstor = (lake_var->volume/lake_var->sarea)*1000. * Clake * lake_con.Cl[0] * soil_con->AreaFract[0] * TreeAdjustFactor[0];
  save_data->surfstor = 0.;
//</devel -- bug.patch>

  return(0);
}


int vic411_update_thermal_nodes(vic411_dist_prcp_struct    *prcp,
			  int                  Nveg,
			  int                  Nnodes,
			  int                  Ndist,
			  vic411_soil_con_struct     *soil_con,
			  vic411_veg_con_struct      *veg_con)
/**********************************************************************
  vic411_update_thermal_nodes           Jennifer Adam        August 16, 2007

  This routine is run after subsidence occurs (used only for EXCESS_ICE option).
  This routine updates the node depths and interpolates the current
  node temperatures to the new depths, then recalculates the nodal
  thermal properties.  Much of this routine is taken directly from
  vic411_initialize_model_state.

  Modifications:
  2009-Feb-09 Removed dz_node from call to
	      vic411_distribute_node_moisture_properties.			KAC via TJB
  2009-Feb-09 Removed dz_node from call to find_0_degree_front.		KAC via TJB
**********************************************************************/
{
  extern vic411_option_struct vic411_options;
  extern vic411_veg_lib_struct *vic411_veg_lib;
#if LINK_DEBUG
  extern vic411_debug_struct debug;
#endif

  char     ErrStr[MAXSTRING];
  char     FIRST_VEG;
  int      veg, index, dist;
  int      lidx;
  int      dry;
  int      band;
  int      ErrorFlag;
  double   Cv;
  double   Zsum, dp;
  double   tmpdp, tmpadj, Bexp;
  double   moist[MAX_VEG][MAX_BANDS][MAX_LAYERS];

  vic411_cell_data_struct     ***cell;
  vic411_energy_bal_struct     **energy;

  double Tnode_prior[MAX_NODES];
  double Zsum_prior[MAX_NODES];

  cell    = prcp->cell;
  energy  = prcp->energy;
  
  dp = soil_con->dp;

  FIRST_VEG = TRUE;

  /*****************************************************************
    Update soil thermal node depths, thicknesses, and temperatures.
    CASE 3: Initialize Energy Balance Variables if not using quick
    ground heat flux, and no Initial Condition File Given 
    (Currently this is the only case that works with EXCESS_ICE.)
  *****************************************************************/

  /*****************************************************************
    Update soil thermal node depths and thicknesses.
  *****************************************************************/
  //set previous Zsum
  for ( index = 0; index < Nnodes; index++ ) 
    Zsum_prior[index] = soil_con->Zsum_node[index];

  if(!vic411_options.EXP_TRANS){  
    /* Nodes set at surface, the depth of the first layer,
       twice the depth of the first layer, and at the
       damping depth.  Extra nodes are placed equal distance
       between the damping depth and twice the depth of the
       first layer. */
    
    soil_con->dz_node[0] = soil_con->depth[0];
    soil_con->dz_node[1] = soil_con->depth[0];
    soil_con->dz_node[2] = soil_con->depth[0];	  
    soil_con->Zsum_node[0] = 0;
    soil_con->Zsum_node[1] = soil_con[0].depth[0];
    Zsum   = 2. * soil_con[0].depth[0];
    soil_con->Zsum_node[2] = Zsum;
    tmpdp  = dp - soil_con[0].depth[0] * 2.5;
    tmpadj = 3.5;
    for ( index = 3; index < Nnodes-1; index++ ) {
      soil_con->dz_node[index] = tmpdp/(((double)Nnodes-tmpadj));
      Zsum += (soil_con->dz_node[index]
	       +soil_con->dz_node[index-1])/2.;
      soil_con->Zsum_node[index] = Zsum;
    }
    soil_con->dz_node[Nnodes-1] = (dp - Zsum 
				   - soil_con->dz_node[Nnodes-2] 
				   / 2. ) * 2.;
    Zsum += (soil_con->dz_node[Nnodes-2]
	     +soil_con->dz_node[Nnodes-1])/2.;
    soil_con->Zsum_node[Nnodes-1] = Zsum;
    if((int)(Zsum*1000+0.5) != (int)(dp*1000+0.5)) {
      sprintf(ErrStr,"Sum of thermal node thicknesses (%f) in vic411_initialize_model_state do not equal dp (%f), check initialization procedure",Zsum,dp);
      vic411_nrerror(ErrStr);
    }
  }
  else{ /* exponential grid transformation, EXP_TRANS = TRUE*/
    
    /*calculate exponential function parameter */
    Bexp = logf(dp+1.)/(double)(Nnodes-1); //to force Zsum=dp at bottom node
    for ( index = 0; index <= Nnodes-1; index++ )
      soil_con->Zsum_node[index] = expf(Bexp*index)-1.;
    
    //top node	  
    index=0;
    soil_con->dz_node[index] = soil_con->Zsum_node[index+1]-soil_con->Zsum_node[index];
    //middle nodes
    for ( index = 1; index < Nnodes-1; index++ ) {
      soil_con->dz_node[index] = (soil_con->Zsum_node[index+1]-soil_con->Zsum_node[index])/2.+(soil_con->Zsum_node[index]-soil_con->Zsum_node[index-1])/2.;
    }
    //bottom node
    index=Nnodes-1;
    soil_con->dz_node[index] = soil_con->Zsum_node[index]-soil_con->Zsum_node[index-1];
  }
#if VERBOSE
  fprintf(stderr,"More updated parameters in soil_con: dz_node and Zsum_node.\n");
#endif

  /******************************************
    Update soil thermal node temperatures via linear interpolation.
  ******************************************/
  for ( veg = 0 ; veg <= Nveg ; veg++ ) {
    // Initialize soil for existing vegetation types
    Cv = veg_con[veg].Cv;
    
    if ( Cv > 0 ) {
      for( band = 0; band < vic411_options.SNOW_BAND; band++ ) {
	if ( soil_con->AreaFract[band] > 0. ) {
	  //set previous temperatures
	  for ( index = 0; index < Nnodes; index++ ) 
	    Tnode_prior[index] = energy[veg][band].T[index];
	  //top node: no need to update surface temperature
	  //remaining nodes
	  for ( index = 1; index < Nnodes; index++ ) {
	    energy[veg][band].T[index] = vic411_linear_interp(soil_con->Zsum_node[index],Zsum_prior[index-1],Zsum_prior[index],Tnode_prior[index-1],Tnode_prior[index]);	
	  }//node
	}	
      }//band
    }
  }//veg

  /******************************************
    Update soil thermal node properties 
  ******************************************/  
  if ( vic411_options.GRND_FLUX ) {
    
    FIRST_VEG = TRUE;
    for ( veg = 0 ; veg <= Nveg ; veg++) {
      // Initialize soil for existing vegetation types
      Cv = veg_con[veg].Cv;

      if ( Cv > 0 ) {
	for( band = 0; band < vic411_options.SNOW_BAND; band++ ) {
	  // Initialize soil for existing snow elevation bands
	  if ( soil_con->AreaFract[band] > 0. ) {
	    /** Set soil properties for all soil nodes **/
	    if(FIRST_VEG) {
	      FIRST_VEG = FALSE;
	      vic411_set_node_parameters(soil_con->dz_node, soil_con->Zsum_node, soil_con->max_moist_node,
				  soil_con->expt_node, soil_con->bubble_node,
				  soil_con->alpha, soil_con->beta,
				  soil_con->gamma, soil_con->depth,
				  soil_con->max_moist, soil_con->expt, 
				  soil_con->bubble, soil_con->quartz, 
#if QUICK_FS
				  soil_con->ufwc_table_node,
#endif // QUICK_FS
#if EXCESS_ICE
				  soil_con->porosity, soil_con->effective_porosity,
				  soil_con->porosity_node, soil_con->effective_porosity_node,
#endif // EXCESS_ICE
				  Nnodes, vic411_options.Nlayer, soil_con->FS_ACTIVE);	  
	    }

	    for ( lidx = 0; lidx < vic411_options.Nlayer; lidx++ ) 
	      moist[veg][band][lidx] = cell[0][veg][band].layer[lidx].moist;

	    /* set soil moisture properties for all soil thermal nodes */
	    if ( !( vic411_options.LAKES && veg_con->LAKE != 0 ) ) {
	      ErrorFlag = vic411_distribute_node_moisture_properties(energy[veg][band].moist,
						  energy[veg][band].ice,
						  energy[veg][band].kappa_node,
						  energy[veg][band].Cs_node,
						  soil_con->Zsum_node,
						  energy[veg][band].T,
						  soil_con->max_moist_node,
#if QUICK_FS
						  soil_con->ufwc_table_node,
#else
						  soil_con->expt_node,
						  soil_con->bubble_node,
#endif // QUICK_FS
#if EXCESS_ICE
						  soil_con->porosity_node,
						  soil_con->effective_porosity_node,
#endif // EXCESS_ICE
						  moist[veg][band],
						  soil_con->depth,
						  soil_con->soil_density,
						  soil_con->bulk_density,
						  soil_con->quartz,
						  Nnodes, vic411_options.Nlayer,
						  soil_con->FS_ACTIVE);
	      if ( ErrorFlag == ERROR )
		return ( ErrorFlag );
	    }

	    /* initialize layer moistures and ice contents */
	    for ( dry = 0; dry < Ndist; dry++ ) {	      
	      if ( !( vic411_options.LAKES && veg_con->LAKE != 0 ) ) 
		ErrorFlag = vic411_estimate_layer_ice_content(cell[dry][veg][band].layer,
						       soil_con->Zsum_node,
						       energy[veg][band].T,
						       soil_con->max_moist_node,
#if QUICK_FS
						       soil_con->ufwc_table_node,
#else
						       soil_con->expt_node,
						       soil_con->bubble_node,
#endif // QUICK_FS
						       soil_con->depth,
						       soil_con->max_moist,
#if QUICK_FS
						       soil_con->ufwc_table_layer,
#else
						       soil_con->expt,
						       soil_con->bubble,
#endif // QUICK_FS
#if SPATIAL_FROST
						       soil_con->frost_fract, 
						       soil_con->frost_slope, 
#endif // SPATIAL_FROST
#if EXCESS_ICE
						       soil_con->porosity,
						       soil_con->effective_porosity,
#endif // EXCESS_ICE
						       soil_con->bulk_density,
						       soil_con->soil_density,
						       soil_con->quartz,
						       soil_con->resid_moist,
						       Nnodes, vic411_options.Nlayer, 
						       soil_con->FS_ACTIVE);	      
	    }
	    
	    /* Find freezing and thawing front depths */
	    if(!vic411_options.QUICK_FLUX && soil_con->FS_ACTIVE) 
	      if ( !( vic411_options.LAKES && veg_con->LAKE != 0 ) ) 
		vic411_find_0_degree_fronts(&energy[veg][band], soil_con->Zsum_node,
				     energy[veg][band].T, Nnodes);
	  }
	}//band
      }
    }//veg	
  }

  return(0);  
}
