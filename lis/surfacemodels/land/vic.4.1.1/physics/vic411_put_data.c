#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_put_data.c,v 5.14.2.26 2009/10/13 19:47:21 vicadmin Exp $";

int  vic411_put_data(vic411_dist_prcp_struct  *prcp,
	      vic411_atmos_data_struct *atmos,
              vic411_soil_con_struct   *soil_con,
	      vic411_veg_con_struct    *veg_con,
              vic411_lake_con_struct   *lake_con,
              vic411_out_data_file_struct   *out_data_files,
              vic411_out_data_struct   *out_data,
              vic411_save_data_struct  *save_data,
	      vic411_dmy_struct        *dmy,
              int                rec)
/**********************************************************************
	vic411_put_data.c	Dag Lohmann		January 1996

  This routine converts data units, and stores finalized values
  in an array for later output to the output files.

  modifications:
  06-24-98  modified for new distributed presipitation data structures KAC
  01-20-00 modified to deal with simplified frozen soil moisture layers
           and frost depth / thaw depth accounting                 KAC
  03-08-00 modified to eliminate extra lines for storing bare
           soil variables.                                         KAC
  6-8-2000 modified to handle spatially distribute frozen soil     KAC
  10-6-2000 modified to handle partial snow cover                  KAC
  02-27-01 modified to output lake model variables                 KAC
  11-18-02 updated output of lake variables to reflect algorithm
           changes.  Also added output variables for blowing snow
           algorithm.                                              LCB
  03-12-03 modified to add additional energy balance variable storage
           when output of snow bands is selected.                  KAC
  03-12-03 Modifed to add AboveTreeLine to vic411_soil_con_struct so that
           the model can make use of the computed treeline.     KAC
  30-Oct-03 Snow_flux was incorrectly set to Tcanopy.  Fixed.   TJB
  25-Aug-04 Sub_snow was incorrectly set to blowing_flux.  Now it is
            set to vapor_flux.                                  TJB
  28-Sep-04 Now out_data->aero_resist stores the aerodynamic resistance
            used in flux calculations.                          TJB
  2005-Mar-24 Modified to compute ALMA output variables.        TJB
  2005-Apr-23 Now aero_cond is aggregated instead of aero_resist.       TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
              out_data and out_data_files structures; removed the
              OPTIMIZE and LDAS_OUTPUT vic411_options; uses the
	      new save_data structure; implemented aggregation.		TJB
  2006-Oct-10 Shortened the names of output variables whose names were
	      too long; fixed typos in others; created new OUT_IN_LONG
	      variable.							TJB
  2006-Nov-07 Added OUT_SOIL_TNODE.					TJB
  2006-Nov-07 Assigned value to overstory.				TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2006-Nov-30 Added OUT_DELSURFSTOR.					TJB
  2006-Nov-30 Convert pressure and vapor pressure to kPa for output.	TJB
  2006-Dec-20 Changed OUT_SURF_TEMP from average of T[0] and T[1] to
	      direct assignment of T[0].				TJB
  2007-Apr-21 Moved initialization of tmp_fract to immediately before the
	        #if SPATIAL_FROST
	      block, so that it would be initialized in all cases.	TJB
  2007-Aug-17 Added EXCESS_ICE output variables.			JCA
  2007-Aug-22 Added OUT_WATER_ERROR as output variable.			JCA
  2007-Nov-06 Lake area is now the larger of lake.areai and lake.sarea.
	      Added wetland canopyevap and canopy_vapor_flux to grid
	      cell flux aggregation.					LCB via TJB
  2008-Apr-21 Made computation of out_data[OUT_SURFSTOR] more robust.	TJB
  2008-Sep-09 Calculate sarea in order to output lake surface area at
	      the end of the time step.  The stored variable
	      lake->sarea represents the sarea from the beginning of
	      the time step, not the updated value from the end of the
	      time step.						LCB via TJB
  2008-Sep-09 Added SOIL_TNODE_WL as an output variable, the soil
	      temperature in the wetland fraction of the grid cell.	LCB via TJB
  2008-Sep-09 Allow output of wetland frost/thaw depths even if Clake
	      is 1.0 since wetland energy balance is always computed.	LCB via TJB
  2008-Sep-09 Lake depth assignment moved up to precede sarea
	      assignment.						LCB via TJB
  2008-Sep-09 Check to make sure that area > 0.0 when checking to see
	      if ice area > sarea.					LCB via TJB
  2008-Oct-23 Changed data type of vic411_put_data() to be int so that it
	      can return ErrorFlag.					TJB
  2009-Jan-12 Added a final return of (0) since the data type of vic411_put_data()
	      is int rather than void.					TJB
  2009-Jan-16 Modified aero_resist_used and Ra_used to become arrays of
	      two elements (surface and overstory); added
	      vic411_options.AERO_RESIST_CANSNOW.				TJB
  2009-Jan-16 Added AERO_COND1&2 and AERO_RESIST1&2 to track
	      surface and overstory values; changed AERO_COND
	      and AERO_RESIST to track "scene" values.			TJB
  2009-Feb-09 Removed checks on PRT_SNOW_BAND option.			TJB
  2009-Feb-22 Added OUT_VPD.						TJB
  2009-May-17 Added OUT_ASAT.						TJB
  2009-Jun-09 Modified to use extension of vic411_veg_lib structure to contain
	      bare soil information.					TJB
  2009-Jun-09 Added OUT_PET_*, potential evap computed for various
	      reference land cover types.				TJB
  2009-Jun-09 Cell_data structure now only stores final aero_resist
	      values (called "aero_resist").  Preliminary uncorrected
	      aerodynamic resistances for current vegetation and various
	      reference land cover types for use in potential evap
	      calculations is stored in temporary array aero_resist.	TJB
  2009-Jun-19 Added T vic411_flag to indicate whether TFALLBACK occurred.	TJB
  2009-Jul-31 Modified so that wetland veg is now included in vic411_main loop
	      over veg tiles and aggregated the same way as all other
	      veg tiles.						TJB
  2009-Aug-28 OUT_LAKE_ICE_TEMP and OUT_LAKE_SURF_TEMP are [C].		TJB
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.		TJB
  2009-Sep-28 Created vic411_collect_wb_terms and vic411_collect_eb_terms to handle
	      summing of storages and fluxes from upland veg tiles,
	      wetland veg tile, and lake.  Added logic to handle an
	      initial (pre-simulation) call for purpose of initializing
	      water and energy balance checks.				TJB
  2009-Sep-30 Miscellaneous fixes for lake model.			TJB
  2009-Oct-05 Modifications for taking changes in lake area into account.	TJB
  2009-Oct-08 Extended T fallback scheme to snow and ice T.		TJB
**********************************************************************/
{
  extern vic411_global_param_struct vic411_global_param;
  extern vic411_veg_lib_struct  *vic411_veg_lib;
  extern vic411_option_struct    vic411_options;
#if LINK_DEBUG
  extern vic411_debug_struct     debug;
#endif

  int                     veg;
  int                     index;
  int                     Ndist;
  int                     dist;
  int                     band;
  int                     Nbands;
  int                     overstory;
  int                     HasVeg;
  int                     IsWet;
  char              *AboveTreeLine;
  double            *AreaFract;
  double            *depth;
  double            *dz;
#if SPATIAL_FROST
  double            *frost_fract;
  double             frost_slope;
#endif // SPATIAL_FROST
  double             dp;
  int                skipyear;
  double                  Cv;
  double                  Clake;
  double                  Cv_save;
  double                  mu;
  double                  cv_baresoil;
  double                  cv_veg;
  double                  cv_overstory;
  double                  cv_snow;
  double                  inflow;
  double                  outflow;
  double                  storage;
  double                  TreeAdjustFactor[MAX_BANDS];
  double                  ThisAreaFract;
  double                  ThisTreeAdjust;
  int                     n;
  int                     v;
  int                     i;
  int                     dt_sec;
  int                     out_dt_sec;
  int                     out_step_ratio;
  static int              step_count;
  int                     ErrorFlag;
  static int              Tfoliage_fbcount_total;
  static int              Tcanopy_fbcount_total;
  static int              Tsnowsurf_fbcount_total;
  static int              Tsurf_fbcount_total;
  static int              Tsoil_fbcount_total;

  vic411_cell_data_struct     ***cell;
  vic411_energy_bal_struct     **energy;
  vic411_lake_var_struct         lake_var;
  vic411_snow_data_struct      **snow;
  vic411_veg_var_struct       ***veg_var;

  AboveTreeLine = soil_con->AboveTreeLine;
  AreaFract = soil_con->AreaFract;
  depth = soil_con->depth;
  dz = soil_con->dz_node;
#if SPATIAL_FROST
  frost_fract = soil_con->frost_fract;
  frost_slope = soil_con->frost_slope;
#endif // SPATIAL_FROST
  dp = soil_con->dp;
  skipyear = vic411_global_param.skipyear;
  dt_sec = vic411_global_param.dt*SECPHOUR;
  out_dt_sec = vic411_global_param.out_dt*SECPHOUR;
  out_step_ratio = (int)(out_dt_sec/dt_sec);
  if (rec >= 0) step_count++;
  if (rec == 0) {
    Tsoil_fbcount_total = 0;
    Tsurf_fbcount_total = 0;
    Tsnowsurf_fbcount_total = 0;
    Tcanopy_fbcount_total = 0;
    Tfoliage_fbcount_total = 0;
  }

  if(vic411_options.DIST_PRCP) 
    Ndist = 2;
  else 
    Ndist = 1;

  // Compute treeline adjustment factors
  for ( band = 0; band < vic411_options.SNOW_BAND; band++ ) {
    if ( AboveTreeLine[band] ) {
      Cv = 0;
      for ( veg = 0 ; veg < veg_con[0].vegetat_type_num ; veg++ ) {
	if ( vic411_veg_lib[veg_con[veg].veg_class].overstory ) {
          if (vic411_options.LAKES && veg_con[veg].LAKE) {
            if (band == 0) {
              // Fraction of tile that is flooded
              if (lake_var.new_ice_area > lake_var.surface[0])
                Clake = lake_var.new_ice_area/lake_con->basin[0];
              else
                Clake = lake_var.surface[0]/lake_con->basin[0];
	      Cv += veg_con[veg].Cv*(1-Clake);
            }
          }
          else {
	    Cv += veg_con[veg].Cv;
          }
        }
      }
      TreeAdjustFactor[band] = 1. / ( 1. - Cv );
    }
    else TreeAdjustFactor[band] = 1.;
    if ( TreeAdjustFactor[band] != 1 && rec == 0 )
      fprintf( stderr, "WARNING: Tree adjust factor for band %i is equal to %f.\n", band, TreeAdjustFactor[band] );
  }

  cv_baresoil = 0;
  cv_veg = 0;
  cv_overstory = 0;
  cv_snow = 0;

  // Initialize output data to zero
  vic411_zero_output_list(out_data);

  // Set output versions of input forcings
  out_data[OUT_AIR_TEMP].data[0]  = atmos->air_temp[vic411_NR];
  out_data[OUT_DENSITY].data[0]   = atmos->density[vic411_NR];
  out_data[OUT_LONGWAVE].data[0]  = atmos->longwave[vic411_NR];
  out_data[OUT_PREC].data[0]      = atmos->out_prec;
  out_data[OUT_PRESSURE].data[0]  = atmos->pressure[vic411_NR]/kPa2Pa;
  out_data[OUT_QAIR].data[0]      = EPS * atmos->vp[vic411_NR]/atmos->pressure[vic411_NR];
  out_data[OUT_RAINF].data[0]     = atmos->out_rain;
  out_data[OUT_REL_HUMID].data[0] = 100.*atmos->vp[vic411_NR]/(atmos->vp[vic411_NR]+atmos->vpd[vic411_NR]);
  out_data[OUT_SHORTWAVE].data[0] = atmos->shortwave[vic411_NR];
  out_data[OUT_SNOWF].data[0]     = atmos->out_snow;
  out_data[OUT_VP].data[0]        = atmos->vp[vic411_NR]/kPa2Pa;
  out_data[OUT_VPD].data[0]       = atmos->vpd[vic411_NR]/kPa2Pa;
  out_data[OUT_WIND].data[0]      = atmos->wind[vic411_NR];
 
  cell    = prcp->cell;
  energy  = prcp->energy;
  lake_var = prcp->lake_var;
  snow    = prcp->snow;
  veg_var = prcp->veg_var;

  /****************************************
    Store Output for all Vegetation Types (except lakes)
  ****************************************/
  for ( veg = 0 ; veg <= veg_con[0].vegetat_type_num ; veg++) {

    Cv = veg_con[veg].Cv;
    Clake = 0;
    Nbands = vic411_options.SNOW_BAND;
    IsWet = 0;
    // Check if this is lake/wetland tile
    if (vic411_options.LAKES && veg_con[veg].LAKE) {
      // Only consider non-flooded portion of wetland tile
      if (lake_var.new_ice_area > lake_var.surface[0])
        Clake = lake_var.new_ice_area/lake_con->basin[0];
      else
        Clake = lake_var.surface[0]/lake_con->basin[0];
      Cv_save = Cv;
      Cv *= (1-Clake);
      Nbands = 1;
      IsWet = 1;
    }

    if (veg < veg_con[0].vegetat_type_num)
      HasVeg = 1;
    else
      HasVeg = 0;

    if ( Cv > 0  || Clake > 0) {

      overstory = vic411_veg_lib[veg_con[veg].veg_class].overstory;

      /*********************************
        Store Output for all Bands 
      *********************************/
      for(band=0;band<Nbands;band++) {

        ThisAreaFract = AreaFract[band];
        ThisTreeAdjust = TreeAdjustFactor[band];
        if (IsWet) {
          ThisAreaFract = 1;
          ThisTreeAdjust = 1;
        }

        if(ThisAreaFract > 0. && ( veg == veg_con[0].vegetat_type_num
           || ( !AboveTreeLine[band] || (AboveTreeLine[band] && !overstory)))) {

          /*******************************************************
            Store Output for Wet and Dry Fractions
          *******************************************************/
          for ( dist = 0; dist < Ndist; dist++ ) {
            if(dist==0) 
              mu = prcp[0].mu[veg];
            else 
              mu = 1. - prcp[0].mu[veg];

            /** compute running totals of various landcovers **/
            if (HasVeg)
              cv_veg += Cv * mu * ThisAreaFract * ThisTreeAdjust;
            else
              cv_baresoil += Cv * mu * ThisAreaFract * ThisTreeAdjust;
            if (overstory)
              cv_overstory += Cv * mu * ThisAreaFract * ThisTreeAdjust;
            if (snow[veg][band].swq> 0.0)
              cv_snow += Cv * mu * ThisAreaFract * ThisTreeAdjust;

	    /*********************************
              Record Water Balance Terms 
	    *********************************/
            vic411_collect_wb_terms(cell[dist][veg][band],
                             veg_var[dist][veg][band],
                             snow[veg][band],
                             lake_var,
                             mu,
                             Cv,
                             ThisAreaFract,
                             ThisTreeAdjust,
                             HasVeg,
                             0,
                             overstory,
                             depth,
                             out_data);

	  } // End wet/dry loop

	  /**********************************
	    Record Energy Balance Terms
	  **********************************/
          vic411_collect_eb_terms(energy[veg][band],
                           snow[veg][band],
                           cell[WET][veg][band],
                           &Tsoil_fbcount_total,
                           &Tsurf_fbcount_total,
                           &Tsnowsurf_fbcount_total,
                           &Tcanopy_fbcount_total,
                           &Tfoliage_fbcount_total,
                           Cv,
                           ThisAreaFract,
                           ThisTreeAdjust,
                           HasVeg,
                           0,
                           overstory,
                           band,
                           depth,
                           dz,
#if SPATIAL_FROST
                           frost_fract,
                           frost_slope,
#endif // SPATIAL_FROST
                           out_data);

          // Store Wetland-Specific Variables

          if (IsWet) {
            // Wetland soil temperatures
            for(i=0;i<vic411_options.Nnode;i++) {
              out_data[OUT_SOIL_TNODE_WL].data[i] = energy[veg][band].T[i];
            }
          }

          /**********************************
            Record Lake Variables
          **********************************/
          if (IsWet) {

            // Recover the total lake/wetland tile fraction
            Cv = Cv_save;

            // Override some variables of soil under lake with those of wetland
            // This is for those variables whose lake values shouldn't be included
            // in grid cell average
            // Note: doing this for eb terms will lead to reporting of eb errors 
            // this should be fixed when we implement full thermal solution beneath lake
            for (i=0; i<MAX_FRONTS; i++) {
              lake_var.energy.fdepth[i]      = energy[veg][band].fdepth[i];
              lake_var.energy.tdepth[i]      = energy[veg][band].fdepth[i];
            }
            for (i=0; i<vic411_options.Nnode; i++) {
              lake_var.energy.ice[i]         = energy[veg][band].ice[i];
              lake_var.energy.T[i]           = energy[veg][band].T[i];
            }
            for (i=0; i<N_PET_TYPES; i++) {
              lake_var.soil.pot_evap[i]      = cell[WET][veg][band].pot_evap[i];
            }
            lake_var.soil.rootmoist          = cell[WET][veg][band].rootmoist;
            lake_var.energy.deltaH           = energy[veg][band].deltaH;
            lake_var.energy.fusion           = energy[veg][band].fusion;
            lake_var.energy.grnd_flux        = energy[veg][band].grnd_flux;


  	    /*********************************
              Record Water Balance Terms 
	    *********************************/
            vic411_collect_wb_terms(lake_var.soil,
                             veg_var[WET][0][0],
                             lake_var.snow,
                             lake_var,
                             1.0,
                             Cv*Clake,
                             ThisAreaFract,
                             ThisTreeAdjust,
                             0,
                             1,
                             overstory,
                             depth,
                             out_data);

	    /**********************************
	      Record Energy Balance Terms
	    **********************************/
            vic411_collect_eb_terms(lake_var.energy,
                             lake_var.snow,
                             lake_var.soil,
                             &Tsoil_fbcount_total,
                             &Tsurf_fbcount_total,
                             &Tsnowsurf_fbcount_total,
                             &Tcanopy_fbcount_total,
                             &Tfoliage_fbcount_total,
                             Cv*Clake,
                             ThisAreaFract,
                             ThisTreeAdjust,
                             0,
                             1,
                             overstory,
                             band,
                             depth,
                             dz,
#if SPATIAL_FROST
                             frost_fract,
                             frost_slope,
#endif // SPATIAL_FROST
                             out_data);

            // Store Lake-Specific Variables

            // Lake ice
            if (lake_var.new_ice_area > 0.0) {
              out_data[OUT_LAKE_ICE].data[0]   = (lake_var.ice_water_eq/lake_var.new_ice_area) * ice_density / RHO_W;
              out_data[OUT_LAKE_ICE_TEMP].data[0]   = lake_var.tempi;
              out_data[OUT_LAKE_ICE_HEIGHT].data[0] = lake_var.hice;
            }
            else {
              out_data[OUT_LAKE_ICE].data[0]   = 0.0;
              out_data[OUT_LAKE_ICE_TEMP].data[0]   = 0.0;
              out_data[OUT_LAKE_ICE_HEIGHT].data[0]   = 0.0;
            }

            // Lake dimensions
            out_data[OUT_LAKE_DEPTH].data[0] = lake_var.ldepth;
            if(lake_var.new_ice_area >= lake_var.surface[0]) {
              out_data[OUT_LAKE_ICE_FRACT].data[0]  = 1.0;
              out_data[OUT_LAKE_SURF_AREA].data[0]  = lake_var.new_ice_area;
            }
            else {
              if (out_data[OUT_LAKE_SURF_AREA].data[0] > 0)
                out_data[OUT_LAKE_ICE_FRACT].data[0]  = lake_var.new_ice_area/lake_var.surface[0];
              else
                out_data[OUT_LAKE_ICE_FRACT].data[0]  = 0;
              out_data[OUT_LAKE_SURF_AREA].data[0]  = lake_var.surface[0];
            }
            out_data[OUT_LAKE_VOLUME].data[0]     = lake_var.volume;

            // Other lake characteristics
            out_data[OUT_LAKE_SURF_TEMP].data[0]  = lake_var.temp[0];
            if (out_data[OUT_LAKE_SURF_AREA].data[0] > 0) {
              out_data[OUT_LAKE_MOIST].data[0]      = (lake_var.volume / soil_con->cell_area) * 1000.; // mm over gridcell
              out_data[OUT_SURFSTOR].data[0]        = (lake_var.volume / soil_con->cell_area) * 1000.; // same as OUT_LAKE_MOIST
            }
            else {
              out_data[OUT_LAKE_MOIST].data[0] = 0;
              out_data[OUT_SURFSTOR].data[0] = 0;
            }

          } // End if vic411_options.LAKES etc.

	} // End if ThisAreaFract etc.

      } // End loop over bands

    } // End if Cv > 0

  } // End loop over veg
 

  /*****************************************
    Aggregation of Dynamic Soil Properties      
   *****************************************/
#if EXCESS_ICE
  for(index=0;index<vic411_options.Nlayer;index++) {
    out_data[OUT_SOIL_DEPTH].data[index]  = soil_con->depth[index];
    out_data[OUT_SUBSIDENCE].data[index]  = soil_con->subsidence[index];
    out_data[OUT_POROSITY].data[index]    = soil_con->effective_porosity[index];
  }  
  for(index=0;index<vic411_options.Nnode;index++) 
    out_data[OUT_ZSUM_NODE].data[index]   = soil_con->Zsum_node[index];
#endif // EXCESS_ICE

  /*****************************************
    Finish aggregation of special-case variables
   *****************************************/
  // Normalize quantities that aren't present over entire grid cell
  if (cv_baresoil > 0) {
    out_data[OUT_BARESOILT].data[0] /= cv_baresoil;
  }
  if (cv_veg > 0) {
    out_data[OUT_VEGT].data[0] /= cv_veg;
  }
  if (cv_overstory > 0) {
    out_data[OUT_AERO_COND2].data[0] /= cv_overstory;
  }
  if (cv_snow > 0) {
    out_data[OUT_SALBEDO].data[0] /= cv_snow;
    out_data[OUT_SNOW_SURF_TEMP].data[0] /= cv_snow;
    out_data[OUT_SNOW_PACK_TEMP].data[0] /= cv_snow;
  }

  // Radiative temperature
  out_data[OUT_RAD_TEMP].data[0] = pow(out_data[OUT_RAD_TEMP].data[0],0.25);

  // Aerodynamic conductance and resistance
  if (out_data[OUT_AERO_COND1].data[0] > SMALL) {
    out_data[OUT_AERO_RESIST1].data[0] = 1 / out_data[OUT_AERO_COND1].data[0];
  }
  else {
    out_data[OUT_AERO_RESIST1].data[0] = HUGE_RESIST;
  }
  if (out_data[OUT_AERO_COND2].data[0] > SMALL) {
    out_data[OUT_AERO_RESIST2].data[0] = 1 / out_data[OUT_AERO_COND2].data[0];
  }
  else {
    out_data[OUT_AERO_RESIST2].data[0] = HUGE_RESIST;
  }
  if (out_data[OUT_AERO_COND].data[0] > SMALL) {
    out_data[OUT_AERO_RESIST].data[0] = 1 / out_data[OUT_AERO_COND].data[0];
  }
  else {
    out_data[OUT_AERO_RESIST].data[0] = HUGE_RESIST;
  }

  /*****************************************
    Compute derived variables
   *****************************************/
  // Water balance terms
  out_data[OUT_DELSOILMOIST].data[0] = 0;
  for (index=0; index<vic411_options.Nlayer; index++) {
    out_data[OUT_SOIL_MOIST].data[index] = out_data[OUT_SOIL_LIQ].data[index]+out_data[OUT_SOIL_ICE].data[index];
    out_data[OUT_DELSOILMOIST].data[0] += out_data[OUT_SOIL_MOIST].data[index];
    out_data[OUT_SMLIQFRAC].data[index] = out_data[OUT_SOIL_LIQ].data[index]/out_data[OUT_SOIL_MOIST].data[index];
    out_data[OUT_SMFROZFRAC].data[index] = 1 - out_data[OUT_SMLIQFRAC].data[index];
  }
  out_data[OUT_DELSOILMOIST].data[0] -= save_data->total_soil_moist;
  out_data[OUT_DELSWE].data[0] = out_data[OUT_SWE].data[0] + out_data[OUT_SNOW_CANOPY].data[0] - save_data->swe;
  out_data[OUT_DELINTERCEPT].data[0] = out_data[OUT_WDEW].data[0] - save_data->wdew;
  out_data[OUT_DELSURFSTOR].data[0] = out_data[OUT_SURFSTOR].data[0] - save_data->surfstor;

  // Energy terms
  out_data[OUT_REFREEZE].data[0] = (out_data[OUT_RFRZ_ENERGY].data[0]/Lf)*dt_sec;
  out_data[OUT_R_NET].data[0] = out_data[OUT_NET_SHORT].data[0] + out_data[OUT_NET_LONG].data[0];

  // Save current moisture state for use in next time step
  save_data->total_soil_moist = 0;
  for (index=0; index<vic411_options.Nlayer; index++) {
    save_data->total_soil_moist += out_data[OUT_SOIL_MOIST].data[index];
  }
  save_data->surfstor = out_data[OUT_SURFSTOR].data[0];
  save_data->swe = out_data[OUT_SWE].data[0] + out_data[OUT_SNOW_CANOPY].data[0];
  save_data->wdew = out_data[OUT_WDEW].data[0];

  /********************
    Check Water Balance 
    ********************/
  inflow  = out_data[OUT_PREC].data[0];
  outflow = out_data[OUT_EVAP].data[0] + out_data[OUT_RUNOFF].data[0] + out_data[OUT_BASEFLOW].data[0];
  storage = 0.;
  for(index=0;index<vic411_options.Nlayer;index++)
    if(vic411_options.MOISTFRACT)
      storage += (out_data[OUT_SOIL_LIQ].data[index] + out_data[OUT_SOIL_ICE].data[index]) 
	* depth[index] * 1000;
    else
      storage += out_data[OUT_SOIL_LIQ].data[index] + out_data[OUT_SOIL_ICE].data[index];
  storage += out_data[OUT_SWE].data[0] + out_data[OUT_SNOW_CANOPY].data[0] + out_data[OUT_WDEW].data[0] + out_data[OUT_SURFSTOR].data[0];
  out_data[OUT_WATER_ERROR].data[0] = vic411_calc_water_balance_error(rec,inflow,outflow,storage);
  
  /********************
    Check Energy Balance 
  ********************/
  if(vic411_options.FULL_ENERGY)
    vic411_calc_energy_balance_error(rec, out_data[OUT_NET_SHORT].data[0] + out_data[OUT_NET_LONG].data[0],
			      out_data[OUT_LATENT].data[0]+out_data[OUT_LATENT_SUB].data[0],
			      out_data[OUT_SENSIBLE].data[0]+out_data[OUT_ADV_SENS].data[0],
			      out_data[OUT_GRND_FLUX].data[0]+out_data[OUT_DELTAH].data[0]+out_data[OUT_FUSION].data[0],
			      out_data[OUT_ADVECTION].data[0] - out_data[OUT_DELTACC].data[0]
			      - out_data[OUT_SNOW_FLUX].data[0] + out_data[OUT_RFRZ_ENERGY].data[0]);



  /******************************************************************************************
    Return to parent function if this was just an initialization of wb and eb storage terms
  ******************************************************************************************/
  if (rec < 0) return(0);



  /********************
    Report T Fallback Occurrences
  ********************/
  if (rec == vic411_global_param.nrecs-1) {
    fprintf(stdout,"Total number of fallbacks in Tfoliage: %d\n", Tfoliage_fbcount_total);
    fprintf(stdout,"Total number of fallbacks in Tcanopy: %d\n", Tcanopy_fbcount_total);
    fprintf(stdout,"Total number of fallbacks in Tsnowsurf: %d\n", Tsnowsurf_fbcount_total);
    fprintf(stdout,"Total number of fallbacks in Tsurf: %d\n", Tsurf_fbcount_total);
    fprintf(stdout,"Total number of fallbacks in soil T profile: %d\n", Tsoil_fbcount_total);
  }

  /********************
    Temporal Aggregation 
    ********************/
  for (v=0; v<N_OUTVAR_TYPES; v++) {
    if (out_data[v].aggtype == AGG_TYPE_END) {
      for (i=0; i<out_data[v].nelem; i++) {
        out_data[v].aggdata[i] = out_data[v].data[i];
      }
    }
    else if (out_data[v].aggtype == AGG_TYPE_SUM) {
      for (i=0; i<out_data[v].nelem; i++) {
        out_data[v].aggdata[i] += out_data[v].data[i];
      }
    }
    else if (out_data[v].aggtype == AGG_TYPE_AVG) {
      for (i=0; i<out_data[v].nelem; i++) {
        out_data[v].aggdata[i] += out_data[v].data[i]/out_step_ratio;
      }
    }
  }
  out_data[OUT_AERO_RESIST].aggdata[0] = 1/out_data[OUT_AERO_COND].aggdata[0];
  out_data[OUT_AERO_RESIST1].aggdata[0] = 1/out_data[OUT_AERO_COND1].aggdata[0];
  //<devel -- bug.patch>
  //handle case when out_data[OUT_AERO_COND2].aggdata[0] is 0
  if ( out_data[OUT_AERO_COND2].aggdata[0] > SMALL )
  {
  out_data[OUT_AERO_RESIST2].aggdata[0] = 1/out_data[OUT_AERO_COND2].aggdata[0];
  }
  else
  {
  out_data[OUT_AERO_RESIST2].aggdata[0] = HUGE_RESIST;
  }
  //</devel -- bug.patch>

  /********************
    Output procedure
    (only execute when we've completed an output interval)
    ********************/
  if (step_count == out_step_ratio) {

    /***********************************************
      Change of units for ALMA-compliant output
    ***********************************************/
    if (vic411_options.ALMA_OUTPUT) {
      out_data[OUT_BASEFLOW].aggdata[0] /= out_dt_sec;
      out_data[OUT_EVAP].aggdata[0] /= out_dt_sec;
      out_data[OUT_EVAP_BARE].aggdata[0] /= out_dt_sec;
      out_data[OUT_EVAP_CANOP].aggdata[0] /= out_dt_sec;
      out_data[OUT_EVAP_LAKE].aggdata[0] /= out_dt_sec;
      out_data[OUT_INFLOW].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_ICE_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_LAKE_SURF_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_PREC].aggdata[0] /= out_dt_sec;
      out_data[OUT_RAINF].aggdata[0] /= out_dt_sec;
      out_data[OUT_REFREEZE].aggdata[0] /= out_dt_sec;
      out_data[OUT_RUNOFF].aggdata[0] /= out_dt_sec;
      out_data[OUT_SNOW_MELT].aggdata[0] /= out_dt_sec;
      out_data[OUT_SNOWF].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_BLOWING].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_CANOP].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_SNOW].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_SNOW].aggdata[0] += out_data[OUT_SUB_CANOP].aggdata[0];
      out_data[OUT_SUB_SURFACE].aggdata[0] /= out_dt_sec;
      out_data[OUT_TRANSP_VEG].aggdata[0] /= out_dt_sec;
      out_data[OUT_BARESOILT].aggdata[0] += KELVIN;
      out_data[OUT_SNOW_PACK_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_SNOW_SURF_TEMP].aggdata[0] += KELVIN;
      for (index=0; index<vic411_options.Nlayer; index++) {
        out_data[OUT_SOIL_TEMP].aggdata[index] += KELVIN;
      }
      for (index=0; index<vic411_options.Nnode; index++) {
        out_data[OUT_SOIL_TNODE].aggdata[index] += KELVIN;
        out_data[OUT_SOIL_TNODE_WL].aggdata[index] += KELVIN;
      }
      out_data[OUT_SURF_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_VEGT].aggdata[0] += KELVIN;
      out_data[OUT_FDEPTH].aggdata[0] /= 100;
      out_data[OUT_TDEPTH].aggdata[0] /= 100;
      out_data[OUT_DELTACC].aggdata[0] *= out_dt_sec;
      out_data[OUT_DELTAH].aggdata[0] *= out_dt_sec;
      out_data[OUT_AIR_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_PRESSURE].aggdata[0] *= 1000;
      out_data[OUT_VP].aggdata[0] *= 1000;
      out_data[OUT_VPD].aggdata[0] *= 1000;
    }

    /*************
      Write Data
    *************/
    if(rec >= skipyear) {
      if (vic411_options.BINARY_OUTPUT) {
        for (v=0; v<N_OUTVAR_TYPES; v++) {
          for (i=0; i<out_data[v].nelem; i++) {
            out_data[v].aggdata[i] *= out_data[v].mult;
          }
        }
      }
      /* commented by Shugong Wang. No need to write output by VIC itself */
      /* vic411_write_data(out_data_files, out_data, dmy, vic411_global_param.out_dt); */
    }

    // Reset the step count
    step_count = 0;

    // Reset the agg data
    for (v=0; v<N_OUTVAR_TYPES; v++) {
      for (i=0; i<out_data[v].nelem; i++) {
        out_data[v].aggdata[i] = 0;
      }
    }

  } // End of output procedure

  return (0);

}

void vic411_collect_wb_terms(vic411_cell_data_struct  cell,
                      vic411_veg_var_struct    veg_var,
                      vic411_snow_data_struct  snow,
                      vic411_lake_var_struct   lake_var,
                      double            mu,
                      double            Cv,
                      double            AreaFract,
                      double            TreeAdjustFactor,
                      int               HasVeg,
                      int               IsWet,
                      int               overstory,
                      double           *depth,
                      vic411_out_data_struct  *out_data)
{

  extern vic411_option_struct    vic411_options;
  double AreaFactor;
  double tmp_evap;
  double tmp_cond1;
  double tmp_cond2;
  double tmp_moist;
  double tmp_ice;
  int index;
#if SPATIAL_FROST
  int                     frost_area;
#endif

  AreaFactor = Cv * mu * AreaFract * TreeAdjustFactor;

  /** record evaporation components **/
  tmp_evap = 0.0;
  for(index=0;index<vic411_options.Nlayer;index++)
    tmp_evap += cell.layer[index].evap;
  if (HasVeg)
    out_data[OUT_TRANSP_VEG].data[0] += tmp_evap * AreaFactor;
  else 
    out_data[OUT_EVAP_BARE].data[0] += tmp_evap * AreaFactor;
  tmp_evap += snow.vapor_flux * 1000.;
  out_data[OUT_SUB_SNOW].data[0] += snow.vapor_flux * 1000. * AreaFactor; 
  out_data[OUT_SUB_SURFACE].data[0] += snow.surface_flux * 1000. * AreaFactor; 
  out_data[OUT_SUB_BLOWING].data[0] += snow.blowing_flux * 1000. * AreaFactor; 
  if (HasVeg) {
    tmp_evap += snow.canopy_vapor_flux * 1000.;
    out_data[OUT_SUB_CANOP].data[0] += snow.canopy_vapor_flux * 1000. * AreaFactor; 
  }
  if (HasVeg) {
    tmp_evap += veg_var.canopyevap;
    out_data[OUT_EVAP_CANOP].data[0] += veg_var.canopyevap * AreaFactor; 
  }
  if (IsWet)  {
    tmp_evap += lake_var.evapw;
    out_data[OUT_EVAP_LAKE].data[0] += lake_var.evapw * AreaFactor; // mm over gridcell
  }
  out_data[OUT_EVAP].data[0] += tmp_evap * AreaFactor; 

  /** record potential evap **/
  out_data[OUT_PET_SATSOIL].data[0] += cell.pot_evap[0] * AreaFactor;
  out_data[OUT_PET_H2OSURF].data[0] += cell.pot_evap[1] * AreaFactor;
  out_data[OUT_PET_SHORT].data[0] += cell.pot_evap[2] * AreaFactor;
  out_data[OUT_PET_TALL].data[0] += cell.pot_evap[3] * AreaFactor;
  out_data[OUT_PET_NATVEG].data[0] += cell.pot_evap[4] * AreaFactor;
  out_data[OUT_PET_VEGNOCR].data[0] += cell.pot_evap[5] * AreaFactor;

  /** record saturated area fraction **/
  out_data[OUT_ASAT].data[0] += cell.asat * AreaFactor; 

  /** record vic411_runoff **/
  out_data[OUT_RUNOFF].data[0]   += cell.vic411_runoff * AreaFactor;

  /** record baseflow **/
  out_data[OUT_BASEFLOW].data[0] += cell.baseflow * AreaFactor; 

  /** record inflow **/
  if (HasVeg) 
    out_data[OUT_INFLOW].data[0] += (cell.inflow + veg_var.canopyevap) * AreaFactor;
  else 
    out_data[OUT_INFLOW].data[0] += (cell.inflow) * AreaFactor;
 
  /** record canopy interception **/
  if (HasVeg) 
    out_data[OUT_WDEW].data[0] += veg_var.Wdew * AreaFactor;

  /** record aerodynamic conductance and resistance **/
  if (cell.aero_resist[0] > SMALL) {
    tmp_cond1 = (1/cell.aero_resist[0]) * AreaFactor;
  }
  else {
    tmp_cond1 = 0;
  }
  out_data[OUT_AERO_COND1].data[0] += tmp_cond1;
  if (overstory) {
    if (cell.aero_resist[1] > SMALL) {
      tmp_cond2 = (1/cell.aero_resist[1]) * AreaFactor;
    }
    else {
      tmp_cond2 = 0;
    }
    out_data[OUT_AERO_COND2].data[0] += tmp_cond2;
  }
  if (overstory) {
    out_data[OUT_AERO_COND].data[0] += tmp_cond2;
  }
  else {
    out_data[OUT_AERO_COND].data[0] += tmp_cond1;
  }

  /** record layer moistures **/
  for(index=0;index<vic411_options.Nlayer;index++) {
    tmp_moist = cell.layer[index].moist;
#if SPATIAL_FROST
    tmp_ice = 0;
    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
      tmp_ice  += (cell.layer[index].ice[frost_area] * frost_fract[frost_area]);
#else
    tmp_ice   = cell.layer[index].ice;
#endif
    tmp_moist -= tmp_ice;
    if(vic411_options.MOISTFRACT) {
      tmp_moist /= depth[index] * 1000.;
      tmp_ice /= depth[index] * 1000.;
    }
    out_data[OUT_SOIL_LIQ].data[index] += tmp_moist * AreaFactor;
    out_data[OUT_SOIL_ICE].data[index] += tmp_ice * AreaFactor;
  }
  out_data[OUT_SOIL_WET].data[0] += cell.wetness * AreaFactor;
  out_data[OUT_ROOTMOIST].data[0] += cell.rootmoist * AreaFactor;

  /** record layer temperatures **/
  for(index=0;index<vic411_options.Nlayer;index++) {
    out_data[OUT_SOIL_TEMP].data[index] += cell.layer[index].T * AreaFactor;
  }

  /*****************************
    Record Snow Pack Variables 
  *****************************/
  
  /** record snow water equivalence **/
  out_data[OUT_SWE].data[0] += snow.swq * AreaFactor * 1000.;
  
  /** record snowpack depth **/
  out_data[OUT_SNOW_DEPTH].data[0] += snow.depth * AreaFactor * 100.;
  
  /** record snowpack albedo, temperature **/
  if (snow.swq> 0.0) {
    out_data[OUT_SALBEDO].data[0] += snow.albedo * AreaFactor;
    out_data[OUT_SNOW_SURF_TEMP].data[0] += snow.surf_temp * AreaFactor;
    out_data[OUT_SNOW_PACK_TEMP].data[0] += snow.pack_temp * AreaFactor;
  }

  /** record canopy intercepted snow **/
  if (HasVeg)
    out_data[OUT_SNOW_CANOPY].data[0] += (snow.snow_canopy) * AreaFactor * 1000.;

  /** record snowpack melt **/
  out_data[OUT_SNOW_MELT].data[0] += snow.melt * AreaFactor;

  /** record snow cover fraction **/
  out_data[OUT_SNOW_COVER].data[0] += snow.coverage * AreaFactor;

}

void vic411_collect_eb_terms(vic411_energy_bal_struct energy,
                      vic411_snow_data_struct  snow,
                      vic411_cell_data_struct  cell_wet,
                      int              *Tsoil_fbcount_total,
                      int              *Tsurf_fbcount_total,
                      int              *Tsnowsurf_fbcount_total,
                      int              *Tcanopy_fbcount_total,
                      int              *Tfoliage_fbcount_total,
                      double            Cv,
                      double            AreaFract,
                      double            TreeAdjustFactor,
                      int               HasVeg,
                      int               IsWet,
                      int               overstory,
                      int               band,
                      double           *depth,
                      double           *dz,
#if SPATIAL_FROST
                      double           *frost_fract,
                      double            frost_slope,
#endif // SPATIAL_FROST
                      vic411_out_data_struct  *out_data)
{

  extern vic411_option_struct    vic411_options;
  double AreaFactor;
  double tmp_fract;
  double rad_temp;
  double surf_temp;
  int index;
#if SPATIAL_FROST
  int    frost_area;
#endif // SPATIAL_FROST

  AreaFactor = Cv * AreaFract * TreeAdjustFactor;

  /**********************************
    Record Frozen Soil Variables
  **********************************/

  /** record freezing and thawing front depths **/
  if(vic411_options.FROZEN_SOIL) {
    for(index = 0; index < MAX_FRONTS; index++) {
      if(energy.fdepth[index] != MISSING)
        out_data[OUT_FDEPTH].data[index] += energy.fdepth[index] * AreaFactor * 100.;
      if(energy.tdepth[index] != MISSING)
        out_data[OUT_TDEPTH].data[index] += energy.tdepth[index] * AreaFactor * 100.;
    }
  }

  tmp_fract = 0;
#if SPATIAL_FROST
  for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
    if ( cell_wet.layer[0].ice[frost_area] )
      tmp_fract  += frost_fract[frost_area];
#else
  if ( cell_wet.layer[0].ice > 0 )
    tmp_fract   = 1.;
#endif
  out_data[OUT_SURF_FROST_FRAC].data[0] += tmp_fract * AreaFactor;

  tmp_fract = 0;
#if SPATIAL_FROST
  if ( (energy.T[0] + frost_slope / 2.) > 0 ) {
    if ( (energy.T[0] - frost_slope / 2.) <= 0 )
      tmp_fract += vic411_linear_interp( 0, (energy.T[0] + frost_slope / 2.), (energy.T[0] - frost_slope / 2.), 1, 0) * AreaFactor;
  }
  else
    tmp_fract += 1 * AreaFactor;
#else
  if ( energy.T[0] <= 0 )
    tmp_fract = 1 * AreaFactor;
#endif

  /**********************************
    Record Energy Balance Variables
  **********************************/

  /** record surface radiative temperature **/
  if ( overstory && snow.snow && !(vic411_options.LAKES && IsWet)) {
    rad_temp = energy.Tcanopy + KELVIN;
  }
  else
    rad_temp = energy.Tsurf + KELVIN;

  /** record surface skin temperature **/
  surf_temp = energy.Tsurf;

  /** record landcover temperature **/
  if(HasVeg) {
    // landcover is bare soil
    out_data[OUT_BARESOILT].data[0] += (rad_temp-KELVIN) * AreaFactor;
  }
  else {
    // landcover is vegetation
    if ( overstory && !snow.snow )
      // here, rad_temp will be wrong since it will pick the understory temperature
      out_data[OUT_VEGT].data[0] += energy.Tfoliage * AreaFactor;
    else
      out_data[OUT_VEGT].data[0] += (rad_temp-KELVIN) * AreaFactor;
  }

  /** record mean surface temperature [C]  **/
  out_data[OUT_SURF_TEMP].data[0] += surf_temp * AreaFactor;
  
  /** record thermal node temperatures **/
  for(index=0;index<vic411_options.Nnode;index++) {
    out_data[OUT_SOIL_TNODE].data[index] += energy.T[index] * AreaFactor;
  }
  if (IsWet) {
    for(index=0;index<vic411_options.Nnode;index++) {
      out_data[OUT_SOIL_TNODE_WL].data[index] = energy.T[index];
    }
  }

  /** record temperature flags  **/
  out_data[OUT_SURFT_FBFLAG].data[0] += energy.Tsurf_fbflag * AreaFactor;
  *Tsurf_fbcount_total += energy.Tsurf_fbcount;
  for (index=0; index<vic411_options.Nnode; index++) {
    out_data[OUT_SOILT_FBFLAG].data[index] += energy.T_fbflag[index] * AreaFactor;
    *Tsoil_fbcount_total += energy.T_fbcount[index];
  }
  out_data[OUT_SNOWT_FBFLAG].data[0] += snow.surf_temp_fbflag * AreaFactor;
  *Tsnowsurf_fbcount_total += snow.surf_temp_fbcount;
  out_data[OUT_TFOL_FBFLAG].data[0] += energy.Tfoliage_fbflag * AreaFactor;
  *Tfoliage_fbcount_total += energy.Tfoliage_fbcount;
  out_data[OUT_TCAN_FBFLAG].data[0] += energy.Tcanopy_fbflag * AreaFactor;
  *Tcanopy_fbcount_total += energy.Tcanopy_fbcount;

  /** record net shortwave radiation **/
  out_data[OUT_NET_SHORT].data[0] += energy.NetShortAtmos * AreaFactor;

  /** record net longwave radiation **/
  out_data[OUT_NET_LONG].data[0]  += energy.NetLongAtmos * AreaFactor;

  /** record incoming longwave radiation at ground surface (under veg) **/
  if ( snow.snow && overstory )
    out_data[OUT_IN_LONG].data[0] += energy.LongOverIn * AreaFactor;
  else
    out_data[OUT_IN_LONG].data[0] += energy.LongUnderIn * AreaFactor;

  /** record albedo **/
  if ( snow.snow && overstory )
    out_data[OUT_ALBEDO].data[0]    += energy.AlbedoOver * AreaFactor;
  else
    out_data[OUT_ALBEDO].data[0]    += energy.AlbedoUnder * AreaFactor;

  /** record latent heat flux **/
  out_data[OUT_LATENT].data[0]    -= energy.AtmosLatent * AreaFactor;

  /** record latent heat flux from sublimation **/
  out_data[OUT_LATENT_SUB].data[0] -= energy.AtmosLatentSub * AreaFactor;

  /** record sensible heat flux **/
  out_data[OUT_SENSIBLE].data[0]  -= energy.AtmosSensible * AreaFactor;

  /** record ground heat flux (+ heat storage) **/
  out_data[OUT_GRND_FLUX].data[0] -= energy.grnd_flux * AreaFactor;

  /** record heat storage **/
  out_data[OUT_DELTAH].data[0]    -= energy.deltaH * AreaFactor;

  /** record heat of fusion **/
  out_data[OUT_FUSION].data[0]    -= energy.fusion * AreaFactor;

  /** record energy balance error **/
  out_data[OUT_ENERGY_ERROR].data[0] += energy.error * AreaFactor;

  /** record radiative effective temperature [K], 
      emissivities set = 1.0  **/
  out_data[OUT_RAD_TEMP].data[0] += ((rad_temp) * (rad_temp) * (rad_temp) * (rad_temp)) * AreaFactor;
  
  /** record snowpack cold content **/
  out_data[OUT_DELTACC].data[0] += energy.deltaCC * AreaFactor;
  
  /** record snowpack advection **/
  if (snow.snow && overstory)
    out_data[OUT_ADVECTION].data[0] += energy.canopy_advection * AreaFactor;
  out_data[OUT_ADVECTION].data[0] += energy.advection * AreaFactor;
  
  /** record snow energy flux **/
  out_data[OUT_SNOW_FLUX].data[0] += energy.snow_flux * AreaFactor;
  
  /** record refreeze energy **/
  if (snow.snow && overstory)
    out_data[OUT_RFRZ_ENERGY].data[0] += energy.canopy_refreeze * AreaFactor;
  out_data[OUT_RFRZ_ENERGY].data[0] += energy.refreeze_energy * AreaFactor;

  /** record melt energy **/
  out_data[OUT_MELT_ENERGY].data[0] += energy.melt_energy * AreaFactor;

  /** record advected sensible heat energy **/
  if ( !overstory )
    out_data[OUT_ADV_SENS].data[0] -= energy.advected_sensible * AreaFactor;
 
  /**********************************
    Record Band-Specific Variables
  **********************************/

  /** record band snow water equivalent **/
  out_data[OUT_SWE_BAND].data[band] += snow.swq * Cv * 1000.;

  /** record band snowpack depth **/
  out_data[OUT_SNOW_DEPTH_BAND].data[band] += snow.depth * Cv * 100.;

  /** record band canopy intercepted snow **/
  if (HasVeg)
    out_data[OUT_SNOW_CANOPY_BAND].data[band] += (snow.snow_canopy) * Cv * 1000.;

  /** record band snow melt **/
  out_data[OUT_SNOW_MELT_BAND].data[band] += snow.melt * Cv;

  /** record band snow coverage **/
  out_data[OUT_SNOW_COVER_BAND].data[band] += snow.coverage * Cv;

  /** record band cold content **/
  out_data[OUT_DELTACC_BAND].data[band] += energy.deltaCC * Cv;
    
  /** record band advection **/
  out_data[OUT_ADVECTION_BAND].data[band] += energy.advection * Cv;
    
  /** record band snow flux **/
  out_data[OUT_SNOW_FLUX_BAND].data[band] += energy.snow_flux * Cv;
    
  /** record band refreeze energy **/
  out_data[OUT_RFRZ_ENERGY_BAND].data[band] += energy.refreeze_energy * Cv;
    
  /** record band melt energy **/
  out_data[OUT_MELT_ENERGY_BAND].data[band] += energy.melt_energy * Cv;

  /** record band advected sensble heat **/
  out_data[OUT_ADV_SENS_BAND].data[band] -= energy.advected_sensible * Cv;

  /** record surface layer temperature **/
  out_data[OUT_SNOW_SURFT_BAND].data[band] += snow.surf_temp * Cv;

  /** record pack layer temperature **/
  out_data[OUT_SNOW_PACKT_BAND].data[band] += snow.pack_temp * Cv;

  /** record latent heat of sublimation **/
  out_data[OUT_LATENT_SUB_BAND].data[band] += energy.latent_sub * Cv;

  /** record band net downwards shortwave radiation **/
  out_data[OUT_NET_SHORT_BAND].data[band] += energy.NetShortAtmos * Cv;

  /** record band net downwards longwave radiation **/
  out_data[OUT_NET_LONG_BAND].data[band] += energy.NetLongAtmos * Cv;

  /** record band albedo **/
  if (snow.snow && overstory)
    out_data[OUT_ALBEDO_BAND].data[band] += energy.AlbedoOver * Cv;
  else
    out_data[OUT_ALBEDO_BAND].data[band] += energy.AlbedoUnder * Cv;

  /** record band net latent heat flux **/
  out_data[OUT_LATENT_BAND].data[band] -= energy.latent * Cv;

  /** record band net sensible heat flux **/
  out_data[OUT_SENSIBLE_BAND].data[band] -= energy.sensible * Cv;

  /** record band net ground heat flux **/
  out_data[OUT_GRND_FLUX_BAND].data[band] -= energy.grnd_flux * Cv;

}
