#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: lakes.eb.c,v 5.11.2.27 2009/10/09 16:59:41 vicadmin Exp $";

int vic411_solve_lake(double             snowfall,
	       double             rainfall,
	       double             tair, 
	       double             wind, 
	       double             vp, 
	       double             shortin, 
	       double             longin, 
	       double             vpd, 
	       double             pressure, 
	       double             air_density, 
	       vic411_lake_var_struct   *lake, 
	       vic411_lake_con_struct    lake_con, 
	       vic411_soil_con_struct    soil_con, 
	       int                dt, 
	       int                rec, 
	       double             wind_h,
	       vic411_dmy_struct         dmy,
	       double             fracprv)
{

/**********************************************************************
* This subroutine solves the energy budget for open water bodies.
*
* Parameter description :
*
  lat  	Latitiude of the lake (degrees).
  lng 	Longitude of the lake (degrees).
  tair	Air temperature (C).
  vp	Air vapor pressure (kPa).
  pressure      	Air pressure (kPa).
  wind		Wind speed (m s-1).
  longwave		Downwelling long wave radiation (W m-2).
  shortwave		Incoming short wave radiation (W m-2).
  prec	         Precipitation (mm).
  lake_energy.latent	 Latent heat flux (W m-2).
  lake_energy.sensible	 Sensible heat flux (W m-2).
  lake_con.eta_a	 Decline of solar radiation input with depth (m-1).  
  lake_con.surface[numnod]	Area of the lake (m2).
  lake.temp[numnod]	Temperature of the lake water (C).
  lake.tempi                Temperature of the lake ice (C).
  lake.hice      	         Height of the lake ice (m).
  lake.fraci	         Fractional coverage of ice (-).
  lake_snow.depth	         Height of the snow on top of the lake (m).
  lake.numnod        	Number of nodes in the lake (-).
  dt		         Time step size (hrs).

  Modifications:
  2006-Oct-16 Now set mixdepth=0 for case of complete ice cover; this
	      guarantees that it is initialized for all cases.		TJB
  2006-Nov-15 Convert swq and surf_water from mm over lake to mm over ice
	      fraction at beginning of function; this was needed to avoid
	      a water budget error since swq and surf_water were being
	      converted to mm over lake at end of the function.		TJB
  2007-Apr-21 Added initialization of energy_ice_melt_bot and qf for case
	      in which fracprv >= FRACLIM but hice == 0.0.		TJB
  2007-Apr-23 Added initialization of lake_energy->Tsurf.		TJB
  2007-Nov-06 Lake snow physics now consistent with land snow pack.
	      Lake ice now limited by available lake water.		LCB via TJB
  2008-Apr-21 Corrected omissions in previous mods.  Also added code to
	      handle previously-unhandled case in deltaH computation.	TJB
  2008-Oct-23 Fixed problems with vic411_usage of uninitialized variables.	LCB via TJB
  2008-Oct-23 Deleted call to vic411_ice_depth() in section 11, since it
	      wasn't working correctly and wasn't necessary.  Now,
	      the reported ice depth is the average over the lake.	LCB via TJB
  2009-Jul-31 Removed lakemain().					TJB
  2009-Sep-28 Removed initialize_prcp() and update_prcp().  Changed
	      final units of runoff_out, baseflow_out, and evapw to
	      mm over lake (to be consistent with changes in vic411_put_data).	TJB
**********************************************************************/

  double LWnetw,LWneti;
  double sw_water, sw_ice;
  double T[MAX_LAKE_NODES];   /* temp of the water column, open fraction. */
  double Tnew[MAX_LAKE_NODES];  
  double Ti[MAX_LAKE_NODES];   /* temp of the water column, ice fraction. */
  double water_density[MAX_LAKE_NODES], water_cp[MAX_LAKE_NODES];
  float albi, albw;
  double Ts, tempalbs;
  double Tcutoff, Tcutk;  /* Lake freezing temperature (K). */
  double Qhw, Qhi;
  double Qgw, Qgi;
  double Qew, Qei;
  double eflux, eadd;
  double Tki, Tkw;     /* Surface temp. of ice and water in Kelvin. */
  double qbot, qw;
  int freezeflag;
  int mixdepth;
  double new_ice_area; /* Ice area formed by freezing in the open water portion. */
  int k, i, ErrorFlag;
  double tw1, tw2, r1, r2, rtu, rnet;
  double Ls, Le;
  double sumjoula, sumjoulb, sumjouli;
  double temp;
  double water_energy;
  double temphw, temphi;
  double SW_under_ice;
  int last_snow;
  double cc;
  double ice_energy;
  double error_over_water;
  double sw_underice_visible;
  double sw_underice_nir;
  double T_lower, T_upper;
  double Tsurf;
  double new_ice_height;
  double windw, windi;
  double energy_ice_formation;
  double energy_ice_melt_bot;
  double melt, deltaCC_ice, Qnet_ice;
  double joule_intermediate;
  double energy_out_bottom;
  double energy_out_bottom_ice;    
  double sw_out;
  int index;
  double tempdepth, in;
  double qf;
  double inputs, outputs, internal, phasechange;
  double new_ice_water_eq;
  double temp_refreeze_energy;
  vic411_energy_bal_struct *lake_energy;
  vic411_snow_data_struct  *lake_snow;
  
  /**********************************************************************
   * 1. Initialize variables.
   **********************************************************************/
 
  lake_energy = &(lake->energy);
  lake_snow = &(lake->snow);

  lake_energy->advection=0.0;
  lake_energy->deltaCC = 0.0;
  lake_energy->grnd_flux = 0.0;
  lake_energy->snow_flux = 0.0;
  lake->snowmlt = 0.0;
  qbot=qw=0.0;
  new_ice_height = new_ice_area = new_ice_water_eq = 0.0;
  lake->evapw=0.0;
  energy_ice_formation = 0.0;
  energy_out_bottom = energy_out_bottom_ice = 0.0;
  lake_snow->vapor_flux=0.0;
  lake_energy->Tsurf = lake->temp[0];
  temp_refreeze_energy = 0.0;

  if(lake->activenod > 0 || lake->areai > 0.0) {

    /* --------------------------------------------------------------------
     * Calculate the water freezing point.
     * ------------------------------------------------------------------- */

    vic411_rhoinit(&Tcutoff, pressure);          /* degrees C */
 
   /* --------------------------------------------------------------------
     * Initialize liquid water and ice temperature profiles to be used
     * later in the module.
     * -------------------------------------------------------------------- */

    for ( k = 0; k < lake->activenod; k++ ) {
      T[k] = lake->temp[k];
      Ti[k] = lake->temp[k];
      water_density[k] = vic411_calc_density(T[k]);
      water_cp[k] = vic411_specheat(T[k]);
    }

    vic411_energycalc( lake->temp, &sumjoulb, lake->activenod, lake->dz, 
		lake->surfdz,lake->surface, water_cp, water_density);

    /**********************************************************************
     * 2. Calculate added precipitation and total snow height.
     **********************************************************************/

    // Convert swq from mm/(lake area) to mm/(ice area)
    if (fracprv > 0.0 && lake_snow->swq > 0.0)
      lake_snow->swq /= fracprv;

    if ( fracprv >= 1.0) {  /* areai is relevant */

      // If there is no snow, add the rain over ice directly to the lake.
      if(lake_snow->swq <=0.0 && rainfall > 0.0) {
        lake->volume += (rainfall/1000.)*lake->areai;
        rainfall = 0.0;
      }
    }
    else if ( fracprv > FRACLIM && fracprv < 1.0 ) { /* sarea is relevant */

      /* Precip over open water directly increases lake volume. */
      lake->volume += ( (snowfall/1000. + rainfall/1000.)*(1-fracprv) * lake->sarea );

      // If there is no snow, add the rain over ice directly to the lake.
      if(lake_snow->swq <= 0.0 && rainfall > 0.0) {
        lake->volume += (rainfall/1000.)*fracprv*lake->sarea;
        rainfall = 0.0; /* Because do not want it added to snow->surf_water */
      }
    }
    else {
      lake->volume += ((rainfall+snowfall)/1000.)*lake->sarea;
      rainfall = 0.0;
      snowfall = 0.0;
    }

    /**********************************************************************
     * 3. Calculate incoming solar radiation over water and ice.
     **********************************************************************/

    /* --------------------------------------------------------------------
     * Calculate the albedo of the lake for ice, snow and liquid water.
     * To be consistent with VIC snow model, tempalbs = NEW_SNOW_ALB if snow
     * is falling, but aging routing does not get reset.
     * -------------------------------------------------------------------- */

    vic411_alblake(Tcutoff, tair, &lake->SAlbedo, &tempalbs, &albi, &albw, snowfall,
	    lake_snow->coldcontent, dt, &lake_snow->last_snow, 
	    lake_snow->swq, lake_snow->depth, &lake_snow->MELTING, dmy.day_in_year );

    /* --------------------------------------------------------------------
     * Calculate the incoming solar radiaton for both the ice fraction
     * and the liquid water fraction of the lake.
     * -------------------------------------------------------------------- */

    if (lake_snow->swq > SNOWCRIT *RHOSNOW/RHO_W) {
      sw_ice = shortin * (1.-tempalbs);
      lake_energy->AlbedoLake = (fracprv) * tempalbs + (1. - fracprv) * albw;
    }
    else if (lake_snow->swq > 0. && lake_snow->swq <= SNOWCRIT*RHOSNOW/RHO_W) {
      sw_ice = shortin * (1.- (albi+tempalbs)/2.);
      lake_energy->AlbedoLake = (fracprv) * (albi+tempalbs)/2. + (1. - fracprv) * albw;
    }
    else if (fracprv > 0. && lake_snow->swq <= 0.) {
      sw_ice = shortin * (1.-albi);
      lake_energy->AlbedoLake = (fracprv) * albi + (1. - fracprv) * albw;
    }
    else /* lake->fraci = 0 */ {
      sw_ice=0.0;
      lake_energy->AlbedoLake = albw;
    }

    sw_water = shortin * (1.-albw);

    /**********************************************************************
     * 4. Calculate initial energy balance over water.
     **********************************************************************/
   
    if( (1.-fracprv) > SMALL && lake->activenod > 0) {
      freezeflag=1;         /* Calculation for water, not ice. */
      windw = wind*log((2. + ZWATER)/ZWATER)/log(wind_h/ZWATER);

      ErrorFlag = vic411_water_energy_balance( lake->activenod, lake->surface, &lake->evapw, 
					dt, freezeflag, lake->dz, lake->surfdz,
					(double)soil_con.lat, Tcutoff, tair, windw, 
					pressure, vp, air_density, longin, sw_water, 
					sumjoulb, wind_h, &Qhw, &Qew, &LWnetw, T, 
					water_density, &lake_energy->deltaH, 
					&energy_ice_formation, fracprv, 
					&new_ice_area, water_cp, &new_ice_height,
					&energy_out_bottom, &new_ice_water_eq,
					lake->volume-lake->ice_water_eq);
      if ( ErrorFlag == ERROR ) return (ERROR);

      /* --------------------------------------------------------------------
       * Do the convective mixing of the lake water.
       * -------------------------------------------------------------------- */

      mixdepth = 0;        /* Set to zero for this time step. */
      vic411_tracer_mixer( T, &mixdepth, freezeflag, lake->surface, 
		    lake->activenod, lake->dz, lake->surfdz, water_cp );

      lake_energy->AtmosLatent      = ( 1. - fracprv ) * Qew;
      lake_energy->AtmosSensible    = ( 1. - fracprv ) * Qhw;
      lake_energy->NetLongAtmos     = ( 1. - fracprv ) * LWnetw;
      lake_energy->NetShortAtmos    = ( 1. - fracprv ) * sw_water;
      lake_energy->refreeze_energy  = energy_ice_formation*(1.-fracprv);
      lake_energy->deltaH          *= ( 1. - fracprv );
      lake_energy->grnd_flux        = -1.*(energy_out_bottom*(1.-fracprv));
      lake_energy->Tsurf            = (1. - fracprv)*T[0];

    }          /* End of water fraction calculations. */
    else {
      // ice covers 100% of lake, reset open water fluxes
      mixdepth = 0;
      LWnetw = 0;
      Qew    = 0;
      Qhw    = 0;

      lake_energy->AtmosLatent      = 0.0;
      lake_energy->AtmosSensible    = 0.0;
      lake_energy->NetLongAtmos     = 0.0;
      lake_energy->NetShortAtmos    = 0.0;
      lake_energy->refreeze_energy  = 0.0;
      lake_energy->deltaH           = 0.0;
      lake_energy->grnd_flux        = 0.0;
      lake_energy->Tsurf            = 0.0;

    }
 
    /**********************************************************************
     *  6. Calculate initial energy balance over ice.
     **********************************************************************/

    if ( fracprv >= FRACLIM ) {

      freezeflag = 0;         /* Calculation for ice. */
      Le = (677. - 0.07 * tair) * JOULESPCAL * GRAMSPKG; /* ice*/
      windi = ( wind * log((2. + soil_con.snow_rough) / soil_con.snow_rough) 
		/ log(wind_h/soil_con.snow_rough) );
      if ( windi < 1.0 ) windi = 1.0;
      lake->aero_resist = (log((2. + soil_con.snow_rough) / soil_con.snow_rough)
			   * log(wind_h/soil_con.snow_rough) / (von_K*von_K)) / windi;

      /* Calculate snow/ice temperature and change in ice thickness from 
         surface melting. */
      ErrorFlag = vic411_ice_melt( wind_h+soil_con.snow_rough, lake->aero_resist, &(lake->aero_resist),
			    Le, lake_snow, lake, dt,  0.0, soil_con.snow_rough, 1.0, 
			    rainfall, snowfall,  windi, Tcutoff, tair, sw_ice, 
			    longin, air_density, pressure,  vpd,  vp, &lake->snowmlt, 
			    &lake_energy->advection, &lake_energy->deltaCC, 
			    &lake_energy->snow_flux, &Qei, &Qhi, &Qnet_ice, 
			    &temp_refreeze_energy, &LWneti, fracprv);
      if ( ErrorFlag == ERROR ) return (ERROR);

      lake_energy->refreeze_energy += temp_refreeze_energy*fracprv;
      lake->tempi = lake_snow->surf_temp;  

      /**********************************************************************
       *  7. Adjust temperatures of water column in ice fraction.
       **********************************************************************/

      /* --------------------------------------------------------------------
       * Calculate inputs to vic411_temp_area..
       * -------------------------------------------------------------------- */

      if (lake->activenod > 0) {
        ErrorFlag = vic411_water_under_ice( freezeflag, sw_ice, wind, Ti, water_density, 
				     (double)soil_con.lat, lake->activenod, lake->dz, lake->surfdz,
				     Tcutoff, &qw, lake->surface, &temphi, water_cp, 
				     mixdepth, lake->hice, lake_snow->swq*RHO_W/RHOSNOW,
				     (double)dt, &energy_out_bottom_ice);     
        if ( ErrorFlag == ERROR ) return (ERROR);
      }
      else
        temphi = -sumjoulb;
	
      /**********************************************************************
       *   8.  Calculate change in ice thickness and fraction
       *    within fraction that already has ice.
       **********************************************************************/
      /* Check to see if ice has already melted (from the top) in this time step. */
      if(lake->ice_water_eq > 0.0) {
        ErrorFlag = vic411_lakeice(&lake->tempi, Tcutoff, sw_ice, Ti[0], 
			    fracprv, dt, lake_energy->snow_flux, qw, 
			    &energy_ice_melt_bot, lake_snow->swq*RHO_W/RHOSNOW, 
			    lake_energy->deltaCC, rec, dmy, &qf,
			    &lake->ice_water_eq, lake->volume-new_ice_water_eq, lake->surface[0]);
        if ( ErrorFlag == ERROR ) return (ERROR);
      }

      lake_energy->AtmosLatent += fracprv * Qei;
      lake_energy->advection       *= fracprv;
      lake_energy->AtmosSensible   += fracprv * Qhi;
      lake_energy->NetLongAtmos    += fracprv * LWneti;
      lake_energy->NetShortAtmos   += fracprv * sw_ice;
      lake_energy->deltaH          += fracprv*temphi;
      lake_energy->grnd_flux  += -1.*(energy_out_bottom_ice*fracprv);
//      lake_energy->refreeze_energy += energy_ice_melt_bot;
      lake_energy->refreeze_energy += energy_ice_melt_bot*fracprv;
      lake_energy->Tsurf += fracprv*lake_snow->surf_temp;

    }
    else {
      /* No Lake Ice Fraction */
      LWneti = 0;
      Qei    = 0.;
      Qhi    = 0;
      qf=0.0;
      temphi =0.0;
      lake_energy->refreeze_energy=0.0;
      if(fracprv > 0.0) {
        energy_ice_melt_bot = (lake->hice*RHOICE + (snowfall/1000.)*RHO_W)*Lf/(dt*SECPHOUR);
        lake->areai = 0.0;
        lake->hice = 0.0;
        lake->ice_water_eq = 0.0;

      }
      else {
        energy_ice_melt_bot = 0.0;
        lake->areai = 0.0;
        lake->hice = 0.0;
        lake->ice_water_eq = 0.0;
      }
    }
  
    /**********************************************************************
     * 9. Average water temperature.
     **********************************************************************/

    if(lake->activenod > 0) {
      // Average ice-covered and non-ice water columns.
      vic411_colavg( lake->temp, T, Ti, fracprv,lake->density, lake->activenod,
              lake->dz, lake->surfdz);

      // Calculate depth average temperature of the lake
      lake->tempavg = 0.0;
      for(i=0; i< lake->activenod; i++) {
        lake->tempavg += lake->temp[i]/lake->activenod;
      }
    }
    else
      lake->tempavg = -99;

    /**********************************************************************
     * 10. Calculate the final water heat content and energy balance.
     **********************************************************************/

    /* Incoming energy. */ 
    inputs = (sw_ice + LWneti + lake_energy->advection + Qhi + Qei);
    outputs = energy_out_bottom_ice;
    internal = temphi;
    phasechange = -1*(lake_energy->refreeze_energy)-1.*energy_ice_melt_bot;
    
    lake_energy->error = inputs - outputs - internal - phasechange;
    
    lake_energy->snow_flux       = 0.0;

    // Sign convention
    lake_energy->deltaH *= -1;

    lake_energy->error            = ( lake_energy->NetShortAtmos 
				      + lake_energy->NetLongAtmos 
				      + lake_energy->AtmosSensible 
				      + lake_energy->AtmosLatent 		   
				      + lake_energy->deltaH 
				      + lake_energy->grnd_flux
				      + lake_energy->refreeze_energy 
				      + lake_energy->advection );

    temphw=temphi=0.0;

    /**********************************************************************
     * 11. Final accounting for passing variables back to VIC.
     **********************************************************************/
 
    /* Adjust water and ice evaporation to be representative of 
       entire lake. */
    lake->evapw           *= ( (1. - fracprv ) * dt * SECPHOUR ); // in mm
    lake_snow->vapor_flux *= fracprv; // in meters 
    lake_snow->blowing_flux *= fracprv; // in meters 
    lake_snow->surface_flux *= fracprv; // in meters 	

    // Adjust snow variables to represent average depth over lakefraction, based on
    // ice fraction at the beginning of the time step. Converted to average over
    // wetland land tile in update_prcp.
    
    lake_snow->swq        *= fracprv;
    lake_snow->surf_water *= fracprv;
    lake_snow->pack_water *= fracprv;
    lake_snow->depth = lake_snow->swq * RHO_W / RHOSNOW;
    lake_snow->coldcontent *= fracprv;

   /* Update ice area to include new ice growth in water fraction. */

    lake->new_ice_area  = lake->areai;

    if(new_ice_area > 0.0 ) {
      lake->new_ice_area += new_ice_area;
      lake->ice_water_eq += new_ice_water_eq;
    }
    
    if(lake->ice_water_eq > 0.0) {
      lake->hice = (lake->ice_water_eq/lake->new_ice_area)*RHO_W/RHOICE;
    }
    else
      lake->hice = 0.0;
    
    /* Change area of ice-covered fraction if ice has thinned. */
    if(lake->hice <= 0.0) {
      lake->new_ice_area = 0.0;
      lake->hice = 0.0;
    }
    else if(lake->hice < FRACMIN) {
      lake->new_ice_area = (lake->new_ice_area * lake->hice) / FRACMIN;
      lake->hice = FRACMIN;
    }
  } /* End of if activenods > 0 */

 

  return (0);

}   /* End of vic411_solve_lake function. */

void vic411_latsens (double Tsurf, double Tcutk, double hice, double tair, double wind,
	      double pressure, double vp, double air_density, 
	      double *evap, double *qsen, double wind_h)
{
/**********************************************************************
 * Calculate the partitioning of the energy balance into latent and
 * sensible heat fluxes.
 *
 * Parameters :
 *
 * Tsurf	Lake surface temperature (K).
 * Tcutk        Melting temperature (K)
 * hice		Ice height (m).
 * tair		Air temperature (C).
 * wind		Wind speed (m s-1).
 * pressure	Air pressure (kPa).
 * density	Air density (kg/m3).
 * vp           Vapor pressure (kPa). 
 * evap		Evaporation rate (m s-1).
 * qsen		Sensible heat flux (W m-2).
 * eog		Vapor pressure at lake level (kPa).
 **********************************************************************/

    float dragcoeff;
    double eog,delT;
    double qair, qlake;  /* Specific humidity of the atmosphere and lake, respectively. */
    double delq;		/* Difference in absolute humidity between the lake */
			/* surface and higher up (-). */

/**********************************************************************
 * Calculate the drag coefficient.
 **********************************************************************/

    if(hice > 0.)
      dragcoeff = vic411_lkdrag (Tsurf, tair + KELVIN, wind, ZSNOW, wind_h);
    else
      dragcoeff = vic411_lkdrag (Tsurf, tair + KELVIN, wind, ZWATER, wind_h);
       
/**********************************************************************
 * Determine the coefficients to be used in the calculation of
 * the vapor pressure at lake level depending on whether the lake is 
 * covered with ice or not.
 **********************************************************************/

	if ( (hice <= 0.) && (Tsurf > Tcutk) ) {

/* --------------------------------------------------------------------
 * Lake is not covered with ice (Handbook of Hydrology eq. 4.2.2).
 * eog in kPa.
 * -------------------------------------------------------------------- */
	  eog=.611*exp(17.269*(Tsurf-KELVIN)/(Tsurf+237.3-KELVIN));
	}
       
	else {

/* --------------------------------------------------------------------
 * Lake is covered with ice.
 * -------------------------------------------------------------------- */
       
	  eog=.611*exp(21.874*(Tsurf-KELVIN)/(Tsurf-7.66));
	}


/**********************************************************************
 * Calculate the specific humidity at lake level and the difference
 * between this humidity and the air humidity at measurement height.
 **********************************************************************/

      qlake=0.622*(eog/(pressure-0.378*eog));
      qair = 0.622*(vp/(pressure - 0.378*vp));
      delq = qair - qlake;
 
      /**********************************************************************
       * Calculate the evaporation rate.  Eq. 4 in Hostetler (1991), 
       * after Brutsaert (1982).  evap in mm s-1
       **********************************************************************/

      *evap=-1*dragcoeff*wind*air_density*delq;

/**********************************************************************
 * Calculate the difference in lake surface temperature and the
 * air temperature at measurement height and consequently the
 * sensible heat flux.  Hostetler (1991) eq 5, in W/m^2.
 **********************************************************************/

      delT=tair+KELVIN-Tsurf;
      *qsen=dragcoeff*wind*air_density*Cp;
      *qsen=*qsen*delT;
}


void vic411_alblake (double  Tcutoff, 
	      double  Tair, 
	      double *snowalbedo, 
	      double *albs, 
	      float  *albi, 
	      float  *albw, 
	      double  newsnow, 
	      double  coldcontent, 
	      int     dt, 
	      int    *last_snow, 
	      double  swq,
	      double  depth,
	      char   *MELTING,
	      int     day_in_year)
{
/**********************************************************************
 * Calculate the albedo of snow, ice and water of the lake.
 *
 * Parameters :
 *
 * tair		Air temperature (C).
 * albs		Snow albedo (-).
 * albi		Ice albedo (-).
 * albw		Water albedo (-).
 *
 * Modifications:
 * 2007-Nov-06 Lake snow physics now consistent with land snow pack.	LCB via TJB
 * 2008-Apr-21 Corrected some omissions from the previous mods.		LCB via TJB
 * 2008-Apr-21 Updated to be compatible with new vic411_snow_albedo().		KAC via TJB
 **********************************************************************/
  double albgl, albgs;

  if( (Tair-Tcutoff) > 0.0) {
    if( (Tair-Tcutoff) < 20.) {
      albgl=0.4 - 0.011*(Tair-Tcutoff);
      albgs=0.6 - 0.0245*(Tair-Tcutoff);
    }
    else {	
      albgl=0.4 - 0.011*20.;
      albgs = 0.6 - 0.0245*20.;
    }
  }
  else {
    albgl=0.4;
    albgs=0.6;
  }

  *albi=0.5*albgs + 0.5*albgl; 

  // update number of days since last significant snowfall
  if ( newsnow > TraceSnow )
    *last_snow = 1;
  else if ( swq == 0 )
    *last_snow = 0;
  else
    *last_snow +=1;

  /** Record if snowpack is melting this time step **/
  if (swq > 0.0) {
    if ( coldcontent >= 0 && day_in_year > 60 // ~ March 1
       && day_in_year < 273 // ~ October 1
       ) {
      *MELTING = TRUE;
    }
    else {
      *MELTING = FALSE;
    }
  }
  else {
    *MELTING = FALSE;
  }

  if ( *MELTING && newsnow > TraceSnow )
    *MELTING = FALSE;

  // compute snow surface albedo
  if(swq > 0.0)
    *snowalbedo = vic411_snow_albedo(newsnow, swq, depth, *snowalbedo, coldcontent, dt, *last_snow, *MELTING);
  else if(swq == 0.0 && newsnow > 0.0)
    *snowalbedo = NEW_SNOW_ALB;
  else
    *snowalbedo = 0.0;

  if(newsnow > 0.0)
    *albs = NEW_SNOW_ALB;
  else
    *albs = *snowalbedo;

  *albw = 0.15;

}


void vic411_colavg (double *finaltemp, double *T,double *Ti,float lakeprv,double *density, int numnod, double dz, double surfdz)
  {

/**********************************************************************
 * Calculate the temperature of the lake water at the different node
 * levels and calculate the water density.
 *
 * Parameters :
 * 
 * finaltemp    Averaged water column temperature (C).
 * T		Water temperature for water fraction(C).
 * Ti		Water temperature for ice fraction(C).
 * lakeprv	Fraction of lake covered with ice (-).
 * density  	Water density at each node (kg/m3).
 * numnod	Number of nodes in the lake (-).
 * dz           Thickness of each water layer (m). 
 * surfdz       Thickness of the surface water layer.
 **********************************************************************/
    int j;
    double water_densityw, water_densityi;
    double temp;
    float z;
  
    //   fprintf(stdout, "%d\t",numnod);
      for(j=0; j<numnod; j++)
	{

/* --------------------------------------------------------------------
 * Calculate the densities of the ice and water fractions.
 * -------------------------------------------------------------------- */

         water_densityw = vic411_calc_density(T[j]);
         water_densityi = vic411_calc_density(Ti[j]);
         water_densityw = water_densityw+1000.;
	 water_densityi = water_densityi+1000.; 

/* --------------------------------------------------------------------
 * Calculate the depth differences.
 * -------------------------------------------------------------------- */

         z=dz;
         if (j == 0)
	   z=surfdz;

/* --------------------------------------------------------------------
 * Calculate the lake (water) temperature as a weight of ice and water
 * temperatures.
 * -------------------------------------------------------------------- */

	 /*       temp=((1.-lakeprv)*T[j]*z*water_densityw*sh_water+
              lakeprv*Ti[j]*z*water_densityi*sh_ice)/
              ((z*water_densityw*sh_water+z*water_densityi*sh_ice)*0.5);
	 */
	 temp=((1.-lakeprv)*T[j]*z*water_densityw+
	       lakeprv*Ti[j]*z*water_densityi)/
	   ((1.-lakeprv)*z*water_densityw+lakeprv*z*water_densityi);
	 
   finaltemp[j]=temp;
 /* --------------------------------------------------------------------
 * Recalculate the water density at the different nodes.
 * -------------------------------------------------------------------- */

	 density[j] = vic411_calc_density(finaltemp[j]);
	}
  }

double vic411_calc_density(double ts)
  {

/**********************************************************************
 * Calculate the water density.
 *
 * Parameters :
 *
 * ts		Temperature at the node (C).
 *
 * Returns:
 *
 * rhostps	Water density at the node (kg/m3).
 **********************************************************************/

      double rhow,t,rhostps;

      t=ts;

/* --------------------------------------------------------------------
 * Calculate the water density as a function of temperature.
 * -------------------------------------------------------------------- */

      rhow=999.842594 + .06793952*t - .009095290*t*t +
	.0001001685*t*t*t - .000001120083*t*t*t*t + .000000006536332*t*t*t*t*t;

     

/* --------------------------------------------------------------------
 * Return the difference between 1000 kg/m3 and the calculated
 * density, not the actual value.
 * -------------------------------------------------------------------- */

      rhostps=rhow-1000.;
      return (rhostps);
  }

void vic411_eddy (int freezeflag, double wind, double *T, double *water_density, 
	   double *de, double lat, int numnod, double dz, double surfdz)
   {

/**********************************************************************
 * Calculate the vic411_eddy diffusivity.
 *
 * Parameters :
 *
 * freezeflag	= 1 if liquid water, 0 if ice.
 * wind		Wind speed at 2 meter (m s-1).
 * T		Water temperature at the different nodes (C).
 * density		Water density at the different nodes (C).
 * de		Eddy diffusivity at each node (m2/d).
 * lat	         Latitude of the pixel (degrees).
 * numnod	         Number of nodes in the lake (-).
 **********************************************************************/

      double ks, N2, ws, radmax, Po;
      double dpdz, rad;
      double zhalf[MAX_LAKE_NODES];  
      int k;
      double z;
      double Ri; /* Richardson's number. */      

/**********************************************************************
 * Calculate the density at all nodes and the distance between all nodes.
 **********************************************************************/

      for(k=0; k< numnod; k++)
	{
         zhalf[k] = dz;
	}

/**********************************************************************
 * Calculate the distance between the first node and the surface node.
 **********************************************************************/

      zhalf[0] = (surfdz + dz) * 0.5;

/**********************************************************************
 * If there is ice only molecular diffusivity is taken in account,
 * no vic411_eddy diffusivity.
 **********************************************************************/

      if (freezeflag != 1) {
	for(k=0; k< numnod; k++)
	  {
            de[k]=DM;
	  }
      }

      /**********************************************************************
       * Avoid too low wind speeds for computational stability.
       **********************************************************************/
     
      else {
	
	if(wind < 1.0)
	  wind = 1.0;

      /**********************************************************************
       *Determine the latitudinaly dependent parameter of the Ekman
       * profile, the surface value of the friction velocity and the neutral
       * value of the Prandtl number (Hostetler and Bartlein eq. 6 and 7).
       **********************************************************************/
      
	ks=6.6*pow(sin((double)fabs(lat)*PI/180.),0.5)*pow(wind,-1.84);
	ws=0.0012*wind;
	Po=1.0;

      /**********************************************************************
       * Determine the maximum value of the second term in the calculation
       * of the Richardson number for computational stability.
       **********************************************************************/

	radmax=4.e4;
      
	for(k= 0; k<numnod-1; k++) {

	/**********************************************************************
         * Calculate the vic411_eddy diffusivity for each node.
	 **********************************************************************/

	/* --------------------------------------------------------------------
	 * Calculate the derivative of density with depth, the Brunt-Vaisala
	 * frequency and the depth of the node (Hostetler and Bartlein eq. 9).
	 * -------------------------------------------------------------------- */
	
	  dpdz=(water_density[k+1]-water_density[k])/zhalf[k];
	  N2=(dpdz/(1.e3+water_density[k]))*9.8;   
	  z=surfdz+((float)k)*dz;	 
	
	/* --------------------------------------------------------------------
	 * Calculate the second term in the calculation of the Richardson
	 * number, make sure this number does not get too large for
	 * computational stability.  Also make sure this term does not go
	 * below 1 to avoid negative Richardson numbers (Hostetler and Bartlein eq. 8).
	 * -------------------------------------------------------------------- */
 
	  if ((z*exp(ks*z)/ws) > 1.e8) 
	    rad = radmax;
	
	  else 
	    {
	      rad=1.+40.*N2*(von_K*z)*(von_K*z)/(ws*ws*exp(-2.*ks*z));
	      
	      if (rad > radmax) 
		rad=radmax;
	    
	      if (rad < 1.0) 
		rad=1.0;
	    }

	/* --------------------------------------------------------------------
	 * Calculate the Richardson number and the vic411_eddy diffusivity.
	 * (Hostetler and Bartlein eq. 8 and 5.
	 * -------------------------------------------------------------------- */

	  Ri=(-1.0+sqrt(rad))/20.0;
	  de[k]=DM+(von_K*ws*z/Po)*exp(-ks*z)
	    /(1.0+37.0*Ri*Ri);
	}

      /* --------------------------------------------------------------------
       * The vic411_eddy diffusivity of the last node is assumed to equal the
       * vic411_eddy diffusivity of the second last node.
       * -------------------------------------------------------------------- */
      
	de[numnod-1]=de[numnod-2];
      }
      
   }

void vic411_iceform (double *qfusion, 
	      double *T, 
	      double  Tcutoff, 
	      double  fracprv,
	      double *areaadd, 
	      int     numnod, 
	      int     dt,  
	      double  dz, 
	      double  surfdz,	
	      double *cp, 
	      double *surface, 
	      double *new_ice_height, 
	      double *water_density,
	      double *new_ice_water_eq,
	      double  lvolume)
   {

/**********************************************************************
 * Calculate the form of new ice in the lake as long as the fractional
 * coverage of the ice is not 1 yet.
 *
 * Parameters :
 *
 * T		Water temperatures for all the nodes (C).
 * Tcutoff	Temperature at which water freezes (C).
 * fracprv	Fractional coverage of ice before this calculation (-).
 * fracadd	Added fractional ice coverage (-).
 * lake.fraci	New fractional ice coverage (-).
 * qfusion      Heat flux absorbed into the ice (W m-2).
 * hice		Ice height (m).
 * dt		Time step (s).
 * numnod	Number of nodes in the lake (-).
 * cp           Specific heat (J/Kg K) 
 *
 * Modifications:
 * 2007-Nov-06 Ice formation is now limited by available liquid water.
 *             ice_water_eq is now the state variable for lake ice water
 *             storage (used to be hice).				LCB via TJB
 **********************************************************************/
     double sum, extra;
     int j;
     double xfrac, di;

/**********************************************************************
 * Calculate the newly added ice to the ice layer.
 **********************************************************************/
     *qfusion = 0.0;
     sum = 0.0;
     *new_ice_water_eq = 0.0;

      for(j=0; j<numnod;j++) {
	  if (T[j]  < Tcutoff) {

/* --------------------------------------------------------------------
 * Calculate the ice growth (extra in J).
 * -------------------------------------------------------------------- */
	    
            if (j == 0) 
	      extra = ( (Tcutoff-T[j]) * surfdz * RHO_W
			* cp[j]  * (1.0-fracprv) * (surface[j]+surface[j+1])/2.);
	    else if (j == numnod-1)
	      extra = ( (Tcutoff-T[j]) * dz * RHO_W
			* cp[j]  * (1.0-fracprv) * surface[j]);
	    else
	      extra = ( (Tcutoff-T[j]) * dz * RHO_W
		        * cp[j]  * (1.0-fracprv) * (surface[j]+surface[j+1])/2.);
/* --------------------------------------------------------------------
 * If ice, the temperature of the water is the freezing temperature.
 * -------------------------------------------------------------------- */

            T[j]=Tcutoff;
	 
            sum+=extra;
	    }
	}

/**********************************************************************
 * Calculate the heat flux absorbed into the ice (in W m-2) and the 
 * thickness of the new ice.
 **********************************************************************/

      *new_ice_water_eq = sum/(RHO_W*Lf);

      if(lvolume > *new_ice_water_eq) {

        *qfusion=(sum/(dt*SECPHOUR*surface[0]*(1.0-fracprv))); /*W m-2*/

        di=sum/(Lf*RHOICE);    /* m^3 of ice formed */
      }
      else if (lvolume > 0.0 ){
        *new_ice_water_eq = lvolume;

        di = *new_ice_water_eq*RHO_W/RHOICE;

        // NEED TO CHANGE ICE TEMPERATURE TO ACCOUNT FOR EXTRA qfusion
        *qfusion=(*new_ice_water_eq*RHO_W/RHOICE)/(dt*SECPHOUR*surface[0]*(1.0-fracprv)); /*W m-2*/
        //      *deltaCC=(sum - *new_ice_water_eq*RHO_W/RHOICE)/(dt*SECPHOUR*surface[0]*(1.0-fracprv)); /*W m-2*/
      }
      else {
        *new_ice_water_eq = 0.0;

        di = 0.0;

        *qfusion=0.0; /*W m-2*/
      }

      /**********************************************************************
       * Calculate the added fractional coverage of ice, make sure the
       * total fractional coverage does not exceed 1.
       **********************************************************************/

      *new_ice_height = FRACMIN;
      *areaadd=di/(FRACMIN);

      if ( *areaadd > (1.0 - fracprv)*surface[0]) {
        *new_ice_height = di/((1.-fracprv)*surface[0]);
        *areaadd=(1.-fracprv)*surface[0];
      }

   }

void vic411_icerad (double  sw,
	     double  hi,
	     double  hs,
	     double *avgcond,
	     double *SWnet,
	     double *SW_under_ice)
{

/**********************************************************************
 * Calculate the radiation balance over ice.
 *
 * Paramterers :
 *
 * sw		Net solar radiation at top of snowpack (W m-2).
 * hi		Ice depth (m).
 * hs		Snow depth (m).
 * avgcond	Thermal conductivity of the ice and snow pack
 *		combined (W/mK). 
 * SWnet	Net short wave radiation at the top of the lake ice.
 * SW_under_ice	Incoming short wave radiation at the bottom of
 *		the snow-ice layer.
 **********************************************************************/
 
  double a, b, c, d;

/**********************************************************************
 * Calculate the thermal conductivity of the combined snow and ice
 * layer.
 **********************************************************************/

  *avgcond=(hs*CONDI+hi*CONDS)/(CONDI*CONDS);

/**********************************************************************
 * Calculate the incoming radiation at different levels in the
 * ice-snow layer.
 **********************************************************************/

/* --------------------------------------------------------------------
 * Calculation of constants. Patterson and Hamblin eq, 7
 * -------------------------------------------------------------------- */

  a=-1.*(1.-exp(-lamssw*hs))/(CONDS*lamssw);
  b=-1.*exp(-lamssw*hs)*(1-exp(-lamisw*hi))/(CONDI*lamisw);
  c=-1.*(1.-exp(-lamslw*hs))/(CONDS*lamslw);
  d=-1.*exp(-lamslw*hs)*(1-exp(-lamilw*hi))/(CONDI*lamilw);

/* --------------------------------------------------------------------
 * Solar radiation at bottom of snow pack. RHS of Patterson and Hamblin, eq. 7
 * -------------------------------------------------------------------- */

  *SWnet=sw*a1*(a+b)+sw*a2*(c+d);

/* --------------------------------------------------------------------
 * Solar radiation at bottom of snow/ice layer.
 * Patterson and Hamblin eq. 8 (qf-qo)
 * -------------------------------------------------------------------- */

  *SW_under_ice = ( a1*sw*(1-exp(-(lamssw*hs+lamisw*hi)))
		    +a2*sw*(1-exp(-(lamslw*hs+lamilw*hi))) );
}
  
int vic411_lakeice (double *tempi, double Tcutoff, double sw_ice,
	     double twater, double fracice, int dt, double snowflux,
	     double qw, double *energy_ice_melt_bot, double sdepth,
	     double SWabsorbed, int rec, vic411_dmy_struct dmy, double *qf,
	     double *ice_water_eq, double volume, double sarea)
/**********************************************************************
 * Calculate the growth and decrease in the lake ice cover.
 * Changed from original model since vic411_ice_melt() (based on VIC/DHSVM vic411_snow_melt())
 * now handles melt of snow and ice from the surface, for consistency with
 * the rest of VIC.  This routine now handles melt/freeze at the
 * bottom of the ice pack.
 *
 * Parameters :
 *
 * tempi      Temperature of the ice (K).
 * Tcutoff    Temperature at which water freezes (K).
 * sw_ice     Net short wave radiation over ice (W m-2).
 * hice               Ice height (m).
 * twater       Temperature of the top water layer (K).
 * *fracice   Revised fraction of the lake covered by ice due to growth/decrease of exisitng cover (-).
 * dt         Time step (s).
 * *qbot      Incoming short wave radiation at bottom of the ice (W m-2).
 * *qw                Heat storage in the lake ice (J/m3).
 * ds           Amount of snowmelt (m)
 * qnetice;   Heat flux absorbed into the ice (W m-2).
 *
 * Modifications:
 * 2007-Nov-06 Changed from original model since vic411_ice_melt() (based on
 *             VIC/DHSVM vic411_snow_melt()) now handles melt of snow and ice
 *             from the surface, for consistency with the rest of VIC.
 *             This routine now handles melt/freeze at the bottom of
 *             the ice pack.						LCB via TJB
 **********************************************************************/
{
  
  double condqw, evapl; 
  double qmelts;           /* Energy available to melt ice. ? */
  double dibot;            /* change in ice surface at the bottom of the pack. */
  double delta_ice_frac;   /* Change in ice covered fraction. */
  double delta_ice_depth;
  double xfrac;
  double tprev;
  double RefrozenWater;
  double new_water_eq;
     
  /**********************************************************************
   * Calculate fluxes at the base of the ice.
   **********************************************************************/
  
  /* --------------------------------------------------------------------
   * Flux of heat in the ice at the ice/water interface (P & H, eq. 8)
   * -------------------------------------------------------------------- */

  *qf = snowflux + sw_ice - SWabsorbed;  

  /* --------------------------------------------------------------------
   * Amount of heat used to melt the ice (positive means freezing).
   * -------------------------------------------------------------------- */
  
  *energy_ice_melt_bot = *qf - qw; // from Hostetler 1991

  /* --------------------------------------------------------------------
   * Calculate the growth of the ice pack at the bottom (in meters).
   * -------------------------------------------------------------------- */

  dibot=(*energy_ice_melt_bot/(RHOICE*Lf))*dt*SECPHOUR;

  /* --------------------------------------------------------------------
   * Calculate the water equivalent of the ice volume (in cubic meters).
   * Freezing occurs over surface area of water, not ice, if ice area exceeds
   * water area.
   * -------------------------------------------------------------------- */

  new_water_eq = dibot * sarea * (fracice) * RHOICE/RHO_W;
 
  /**********************************************************************
   * Calculate the new height of the ice pack and the fractional ice
   * cover.
   **********************************************************************/
  if(dibot > 0.0) { /*Freezing; check if enough unfrozen water is available. */
    if(volume - *ice_water_eq >= new_water_eq ) {
      *ice_water_eq += new_water_eq;
    }
    else { /* Freezing is restricted by available water. */
      if( (volume - *ice_water_eq) > 0.0) {
        dibot = (volume - *ice_water_eq) / (sarea*(fracice)*RHOICE/RHO_W);
        *ice_water_eq = volume;
      }
      else {
        dibot = 0.0;
      }
    }
  }
  else { /* Melt */

    *ice_water_eq += new_water_eq;

    // check that ice not completely melted in the current time step
    if ( *ice_water_eq <= 0.0 ) {
      *ice_water_eq = 0.0;
    }

  }
  // *energy_ice_melt_bottom is not currently adjusted if there is not enough water to freeze or ice to melt.
  // This energy should go to ? and warming the water?, respectively.

  return (0);

}

float vic411_lkdrag (float Tsurf, double Tair, double wind, double roughness, double Z1)
{

  /**********************************************************************
   * Calculate the lake drag coefficient.
   *
   * Parameter :
   *
   * Tsurf   	Lake surface temperature (K).
   * Tair		Air temperature (K).
   * wind		Wind speed (m s-1).
   * dragcoeff	Drag coefficient (-).
   **********************************************************************/

     double cdrn, ribn, ribd, rib;
     double cdr, cdrmin;
      
/**********************************************************************
 * Calculate the Richardson number.
 **********************************************************************/

     cdrn=(von_K/log(Z1/roughness))*(von_K/log(Z1/roughness)); /* dimensionless */

     ribn=Z1*G*(1.-Tsurf/Tair);           /* m2/s2 */

      if ((Tsurf/Tair) <= 1.0) {
	ribd=wind*wind+0.1*0.1;
      }
      else {
	ribd=wind*wind+1.0*1.0;
	  }

      rib=ribn/ribd;                /* dimensionless */

/**********************************************************************
 * Calculate the drag coefficient using the Richardson number.
 **********************************************************************/

      if (rib < 0.) {

         cdr=cdrn*(1.0+24.5*sqrt(-cdrn*rib));
      }

      else {
         cdr=cdrn/(1.0+11.5*rib);
	   }

      if((.25*cdrn) > 6.e-4)
	cdrmin=.25*cdrn;
      else
	cdrmin=6.e-4;

      if (cdr < cdrmin) 
	cdr=cdrmin;

      return(cdr);

   }


void vic411_rhoinit(double *tfsp, 
	     double pressure)
{

/**********************************************************************
 * Calculate the temperature at which water 
 * freezes depending on salinity and air pressure.
 *
 * Paramters :
 *
 * tfsp		Lake freezing point (C).
 * pressure     Air pressure (Pa)
 **********************************************************************/

  double salinity;

/**********************************************************************
 * Salinity is assumed to be zero.
 **********************************************************************/

  salinity = 0.;

/**********************************************************************
 * Calculate the lake freezing temperature (C). Pressure in bars.
 **********************************************************************/

  *tfsp = ( -0.0575 * salinity + 1.710523e-3 * pow(salinity,1.5) 
	    - 2.154996e-4 * salinity * salinity - 7.53e-3 * (pressure) / 100. );
      

}
 
double vic411_specheat (double t)
    {

/**********************************************************************
 * Calculate the specific heat of the water depending on water
 * temperature.  Salinity is assumed to be zero.
 *
 * Paramterers :
 *
 * t	Water temperature (C).
 **********************************************************************/

      double cpt;               /* Specific heat (J/Kg K) */

      cpt=4217.4 - 3.720283*t + 0.1412855*t*t - 2.654387e-3*t*t*t +
	2.093236e-5*t*t*t*t;

      return cpt;
    }

void vic411_temp_area(double sw_visible, double sw_nir, double surface_force, 
	       double *T, double *Tnew, double *water_density, double *de,
	       int dt, double *surface, int numnod, double dz, double surfdz, double *temph, 
	       double *cp, double *energy_out_bottom)
    {
/********************************************************************** 				       
  Calculate the water temperature for different levels in the lake.
 
  Parameters :
 
  sw_visible   Shortwave rad in visible band entering top of water column
  sw_nir       Shortwave rad in near infrared band entering top of water column
  surface_force The remaining rerms i nthe top layer energy balance 
  T		Lake water temperature at different levels (K).
  water_density		Water density at different levels (kg/m3).
  de		Diffusivity of water (or ice) (m2/d).
  dt		Time step size (s).
  surface	Area of the lake at different levels (m2).
  numnod	Number of nodes in the lake (-).
  dz        Thickness of the lake layers. 

  Modifications:
  2007-Apr-23 Added initialization of temph.				TJB
  2007-Oct-24 Modified by moving closing bracket for if ( numnod==1 ) up
	      so that the code actually calls vic411_energycalc() even if the
	      lake is represented by only one node.			KAC via TJB

 **********************************************************************/

      double z[MAX_LAKE_NODES], zhalf[MAX_LAKE_NODES];
      double a[MAX_LAKE_NODES], b[MAX_LAKE_NODES], c[MAX_LAKE_NODES];
      double d[MAX_LAKE_NODES];
      double told[MAX_LAKE_NODES];
     
      double dist12;
      int k;
      double surface_1, surface_2, surface_avg, T1;
      double cnextra;
      double swtop; /* The solar radiation at the top of the water column. */
      double top, bot; /* The depth of the top and the bottom of the current water layer. */
      double water_density_new;
      double jouleold, joulenew;
      double energyinput;
      double energymixed;
      double term1, term2;
      double totalenergy;

/**********************************************************************
 * Calculate the distance between the centers of the surface and first
 * lake layers.
 **********************************************************************/

/**********************************************************************
 * Initialize the water density at all the nodes and the depth of all
 * and distance between all nodes.
 **********************************************************************/

      for(k=0; k<numnod; k++) {
	 if(k==0)
	   z[k] = surfdz;
	 else
	   z[k]=dz;
         zhalf[k]=dz;
      }
      zhalf[0]=0.5*(z[0]+z[1]);
      energyinput=0.0;
      energymixed = 0.0;

/**********************************************************************
 * Calculate the right hand side vector in the tridiagonal matrix system
 * of equations.
 **********************************************************************/

	
      surface_1 = surface[0];
      surface_2 = surface[1];
      surface_avg = (surface_1 + surface_2)/2.;

       T1 = (sw_visible*(1*surface_1-surface_2*exp(-lamwsw*surfdz)) + 
      	    sw_nir*(1*surface_1-surface_2*exp(-lamwlw*surfdz)))/surface_avg 
      	+ (surface_force*surface_1)/surface_avg;          /* W m-2 */

       totalenergy = sw_visible  + sw_nir + surface_force;
       energyinput +=T1*surface_avg;
       cnextra = 0.5*(surface_2/surface_avg)*(de[0]/zhalf[0])*((T[1]-T[0])/z[0]);
      
       energymixed += cnextra;

       *temph=0.0;
       
       if(numnod==1)
	 Tnew[0] = T[0]+(T1*dt*SECPHOUR)/((1.e3+water_density[0])*cp[0]*z[0]);
       else {	
	 
	 /* --------------------------------------------------------------------
	  * First calculate d for the surface layer of the lake.
	  * -------------------------------------------------------------------- */
	 
	 d[0]= T[0]+(T1*dt*SECPHOUR)/((1.e3+water_density[0])*cp[0]*z[0])+cnextra*dt*SECPHOUR;
	 
	 *energy_out_bottom = (surface_1 - surface_2)*(sw_visible*exp(-lamwsw*surfdz) + 
						       sw_nir*exp(-lamwlw*surfdz));
	 
	 
	 
/* --------------------------------------------------------------------
 * Calculate d for the remainder of the column.
 * --------------------------------------------------------------------*/

/* ....................................................................
 * All nodes but the deepest node.
 * ....................................................................*/

      for(k=1; k<numnod-1; k++) {
	top = (surfdz+(k-1)*dz);
        bot = (surfdz+(k)*dz);

         surface_1 = surface[k]; 
	 surface_2 = surface[k+1];
         surface_avg =( surface[k]  + surface[k+1]) / 2.;


	  T1 = (sw_visible*(surface_1*exp(-lamwsw*top)-surface_2*exp(-lamwsw*bot)) + 
	       sw_nir*(surface_1*exp(-lamwlw*top)-surface_2*exp(-lamwlw*bot)))/surface_avg; 

	  energyinput +=T1*surface_avg;
	  term1 = 0.5 *(1./surface_avg)*((de[k]/zhalf[k])*((T[k+1]-T[k])/z[k]))*surface_2;
	  term2 = 0.5 *(-1./surface_avg)*((de[k-1]/zhalf[k-1])*((T[k]-T[k-1])/z[k]))*surface_1;
	 
	  cnextra = term1 + term2;
	
	 energymixed += (term1+term2);
         d[k]= T[k]+(T1*dt*SECPHOUR)/((1.e3+water_density[k])*cp[k]*z[k])+cnextra*dt*SECPHOUR;

	 	 *energy_out_bottom += (surface_1 - surface_2)*(sw_visible*exp(-lamwsw*bot) + 
	 						sw_nir*exp(-lamwlw*bot));

      }
/* ....................................................................
 * Calculation for the deepest node.
 * ....................................................................*/
       k=numnod-1;
      surface_1 = surface[k];
      surface_2 = surface[k];
      surface_avg = surface[k];
     
      top = (surfdz+(k-1)*dz);
      bot = (surfdz+(k)*dz);

      //       T1 = (sw_visible*(exp(-lamwsw*top)) + 
      //	    sw_nir*(exp(-lamwlw*top)))*surface_1/surface_avg; 

      T1 = (sw_visible*(surface_1*exp(-lamwsw*top)-surface_2*exp(-lamwsw*bot)) + 
	    sw_nir*(surface_1*exp(-lamwlw*top)-surface_2*exp(-lamwlw*bot)))/surface_avg; 

      energyinput +=T1*surface_avg;
      energyinput /= surface[0];

      cnextra = 0.5 * (-1.*surface_1/surface_avg)*((de[k-1]/zhalf[k-1])*((T[k]-T[k-1])/z[k]));

            *energy_out_bottom = 0.;
      *energy_out_bottom += surface_2*(sw_visible*exp(-lamwsw*bot) + sw_nir*exp(-lamwlw*bot));
      *energy_out_bottom /= surface[0];

      energymixed += cnextra;
    
      d[k] = T[k]+(T1*dt*SECPHOUR)/((1.e3+water_density[k])*cp[k]*z[k])+cnextra*dt*SECPHOUR;

/**********************************************************************
 * Calculate arrays for tridiagonal matrix.
 **********************************************************************/

/* --------------------------------------------------------------------
 * Top node of the column.
 * --------------------------------------------------------------------*/

     
      surface_2 = surface[1];
      surface_avg = (surface[0] + surface[1] ) / 2.;

      b[0] = -0.5 * ( de[0] / zhalf[0] ) *
	(  dt*SECPHOUR / z[0] ) * surface_2/surface_avg;
      a[0] = 1. - b[0];

/* --------------------------------------------------------------------
 * Second to second last node of the column.
 * --------------------------------------------------------------------*/

      for(k=1;k<numnod-1;k++) {
         surface_1 = surface[k];
         surface_2 = surface[k+1];
	 surface_avg = ( surface[k]  + surface[k+1]) / 2.;

         b[k] = -0.5 * ( de[k] / zhalf[k] ) *
	   (  dt*SECPHOUR / z[k] )*surface_2/surface_avg;
         c[k] = -0.5 * ( de[k-1] / zhalf[k-1] ) *
	   (  dt*SECPHOUR / z[k] )*surface_1/surface_avg;
         a[k] = 1. - b[k] - c[k];

	   }
/* --------------------------------------------------------------------
 * Deepest node of the column.
 * --------------------------------------------------------------------*/

      surface_1 = surface[numnod-1];
      surface_avg = surface[numnod-1];
      c[numnod-1] = -0.5 * ( de[numnod-1] / zhalf[numnod-1] ) *
	(  dt*SECPHOUR / z[numnod-1] ) * surface_1/surface_avg;
      a[numnod-1] = 1. - c[numnod-1];

/**********************************************************************
 * Solve the tridiagonal matrix.
 **********************************************************************/
     
      vic411_tridia(numnod,c,a,b,d,Tnew);

    }

/**********************************************************************
 * Adjust energy fluxes for change in density -> should be fixed by
 * moving to lagrangian scheme
 **********************************************************************/
    
    vic411_energycalc(Tnew, &joulenew, numnod,dz, surfdz, surface, cp, water_density);
     
    *temph=0.0;

    *temph = joulenew;

}

void vic411_tracer_mixer (double *T, 
                   int    *mixdepth,
		   int     freezeflag,
                   double *surface,
                   int     numnod, 
                   double  dz, 
		   double surfdz, 
                   double *cp) {
/**********************************************************************
 * Simulate the convective mixing in the lake.
 *
 * Paramters :
 *
 * T		Water temperatures (K).
 * water_density		Water densities (kg/m3).
 * mixdepth	         Top depth of local instability (node number, -).
 * freezeflag	         0 for ice, 1 for liquid water.
 * surface	Area of the lake per node number (m2).
 * numnod	         Number of nodes in the lake (-).
 **********************************************************************/

  int    k,j,m;             /* Counter variables. */
  int    mixprev;
  double avet, avev;
  double vol; /* Volume of the surface layer (formerly vol_tr). */
  double heatcon; /*( Heat content of the surface layer (formerly vol). */ 
  double Tav, densnew;
  double rho_max;
  double water_density[MAX_LAKE_NODES];
  
  for ( k = 0; k < numnod; k++ )
    water_density[k] = vic411_calc_density(T[k]);

/**********************************************************************
 * Initialize the top depth of local instability.
 **********************************************************************/

  mixprev = 0;
      
  for ( k = 0; k < numnod-1; k++ ) {

/**********************************************************************
 * Check for instability at each slice in water column.
 **********************************************************************/
    avet=0.0;
    avev=0.0;
	
    if ( water_density[k] > water_density[k+1] ) {

/* --------------------------------------------------------------------
 * If there is instability apply the mixing scheme.
 * --------------------------------------------------------------------*/

      if (mixprev == 0 && (k+1) > *mixdepth) {

/* ....................................................................
 * Correct the top depth of local instability.
 * ....................................................................*/
	*mixdepth = k+1;

      }

/*----------------------------------------------------------------------
 * Mix from mixprev to k+1  
 *----------------------------------------------------------------------*/
	   
      for ( m = mixprev; m <= k+1; m++ ) {

/* --------------------------------------------------------------------
 * Apply the mixing scheme from the previous depth to the instability
 * up to this node.
 * --------------------------------------------------------------------*/

	if (m == 0) {
	  /* Calculate the heat content and volume of the surface layer. */
	  heatcon = surfdz*(1.e3+water_density[m])*cp[m]*surface[m];
	  vol = surfdz * surface[m];
	}
	else {
	  /* Calculate the heat content and volume of all layers but 
	     the surface layer. */
	  heatcon = dz*(1.e3+water_density[m])*cp[m]*surface[m];
	  vol = dz * surface[m];
	}
	
	/* Calculate the volumetric weighed average lake temperature. */
	avet = avet+T[m]*heatcon;
	avev = avev+heatcon;
	
      }
      
      Tav = avet / avev;

/* --------------------------------------------------------------------
 * Calculate the density of the surface layer.
 * --------------------------------------------------------------------*/

      densnew = vic411_calc_density(Tav);

/* --------------------------------------------------------------------
 * Calculate the maximum density up to the local level.
 * --------------------------------------------------------------------*/

      rho_max = 0.0;
	
      for ( j = 0; j < mixprev; j++ ) {
	if ( (1000.+water_density[j]) > rho_max ) {
	  rho_max=1000.+water_density[j];
	}
      }

/* --------------------------------------------------------------------
 * Adjust temperatures and density in the mixed part of column.
 * --------------------------------------------------------------------*/

      for ( j = mixprev; j <= k+1; j++ ) {
	T[j]=Tav;
	water_density[j] = densnew;
      }

/* --------------------------------------------------------------------
 * Check to make sure that the mixing has not generated new instabilities
 * above the previous depth to the instabilities.
 * --------------------------------------------------------------------*/

      if (rho_max > (1000.+densnew)) {
	
	/* If there are still instabilities iterate again..*/
	mixprev = 0;
	k=-1;
      }
    }

    else {
/**********************************************************************
 * If there are no instabilities up to now then the depth to the
 * instability has to be increased by 1 node.
 **********************************************************************/

      mixprev=k+1;
    }
  }

/**********************************************************************
 * Recalculate the water density.
 **********************************************************************/

  for ( k = 0; k < numnod; k++ ) {
    water_density[k] = vic411_calc_density(T[k]);
  }
}      

void vic411_tridia (int ne, 
	     double *a, 
	     double *b, 
	     double *c, 
	     double *y, 
	     double *x) {
/**********************************************************************
 * Solve a tridiagonal system of equations.
 * Parameters :
 * ns		The number of systems to be solved.
 * nd		First dimension of arrays (larger than or equal to ns).
 * ne		The number of unknowns in each system. This must 
 *		be larger than or equal to 2.
 * a		The sub diagonals of the matrices are stored in locations
 *		a(j,2) through a(j,ne).
 * b		The vic411_main diagonals of the matrices are stored in
 *		locations b(j,1) through b(j,ne).
 * c		The super diagonals of the matrices are stored in
 *		locations c(j,1) through c(j,ne-1).
 * y		The right hand side of the equations is stored in
 *		y(j,1) through y(j,ne).
 * x		The solutions of the systems are returned in
 *		locations x(j,1) through x(j,ne).
 *
 * History : Based on a streamlined vic411_version of the old NCAR ULIB subroutine
 * TRDI used in the PHOENIX climate model of Schneider and Thompson
 * (J.G.R., 1981). Revised by Starley Thompson to solve multiple systems
 * and vectorize well on the CRAY-1. Later revised to include a PARAMETER
 * statement to define loop limits and thus enable Cray short vector
 * loops.
 *
 * Algorithm:  LU decomposition followed by solution.  NOTE: This
 * subroutine executes satisfactorily if the input matrix is diagonally
 * dominant and non-singular.  The diagonal elements are used to pivot, and
 * no tests are made to determine singularity. If a singular or numerically
 * singular matrix is used as input a divide by zero or floating point
 * overflow will result.
 **********************************************************************/

     double alpha[MAX_LAKE_NODES], gamma[MAX_LAKE_NODES]; /* Work arrays dimensioned (nd,ne).*/

     int nm1, i, ib;
  
     nm1 = ne-1;

/**********************************************************************
 * Obtain the LU decompositions.
 **********************************************************************/

	  alpha[0] = 1./b[0];
	  gamma[0] = c[0]*alpha[0];


      for(i=1; i<nm1;i++) {
	  alpha[i] = 1./(b[i]-a[i]*gamma[i-1]);
	  gamma[i] = c[i]*alpha[i];
      }

/**********************************************************************
 * Solve the system.
 **********************************************************************/

      x[0] = y[0]*alpha[0];
	
      for(i=1; i<nm1; i++) {
	  x[i] = (y[i]-a[i]*x[i-1])*alpha[i];
      }

    
      x[nm1] = (y[nm1]-a[nm1]*x[nm1-1])/
	   (b[nm1]-a[nm1]*gamma[nm1-1]);

      for(i=nm1-1; i>=0; i--) {
            x[i] = x[i]-gamma[i]*x[i+1];
      }

}  


void vic411_energycalc (double *finaltemp, double *sumjoule, int numnod, double dz, double surfdz, double *surface, double *cp, double *density)
{
  
  /**********************************************************************
   * Calculate the thermal energy in the column.
   * finaltemp[numnodes]   Average temp for each node of the water column.
   * sumjoule              Thermal energy of water column, in joules 
   * numnod   
   **********************************************************************/
  
  double energy;
  int k;
  
  *sumjoule=0.0;
  
  for(k=0; k<numnod;k++)
    {
      // density = vic411_calc_density(finaltemp[k]);
      
      if(k==0)
	energy=(finaltemp[k]+KELVIN)*surfdz*(1000.+density[k])*cp[k]*(surface[k]+surface[k+1])/2.;
      else if (k==numnod-1)
	energy=(finaltemp[k]+KELVIN)*dz*(1000.+density[k])*cp[k]*surface[k]/2.;
      else
	energy=(finaltemp[k]+KELVIN)*dz*(1000.+density[k])*cp[k]*(surface[k]+surface[k+1])/2.;
      
      *sumjoule+=energy;
    }
}

int vic411_water_balance (vic411_lake_var_struct *lake, vic411_lake_con_struct lake_con, int dt, vic411_dist_prcp_struct *prcp,
		    int rec, int iveg,int band, double lakefrac, vic411_soil_con_struct soil_con, vic411_veg_con_struct veg_con,
#if EXCESS_ICE
		    int SubsidenceUpdate, double total_meltwater,
#endif
		    double deltasnow, double vapor_flux)
/**********************************************************************
 * This routine calculates the water balance of the lake
 
  Modifications:
  2007-Aug-09 Added features for EXCESS_ICE option.				JCA
  2007-Oct-24 Added rec to call for vic411_water_balance so that error and warning
	      messages can report the model record number.			KAC via TJB
  2007-Oct-24 Modified loop for vic411_get_sarea to include the active node.		KAC via TJB
  2007-Nov-06 Ice water content now tracked via ice_water_eq.  Lake area is
	      now the larger of lake->areai and lake->sarea.  Drainage is now
	      modeled as flow over a broad-crested wier.			LCB via TJB
  2008-Apr-21 Corrected some omissions from the previous mods.			LCB via TJB
  2009-Mar-16 Inserted missing logic for SPATIAL_FROST and replaced resid_moist
	      with min_liq.							TJB
  2009-Sep-30 Miscellaneous fixes for lake model.				TJB
  2009-Oct-05 Modified to update/rescale lake and wetland storages and fluxes
	      to account for changes in lake area.				TJB
**********************************************************************/
{
  extern vic411_option_struct   vic411_options;
  int isave_n;
  double d_area, d_level, d_volume;
  double runoff_volume;
  double oldvolume;
  double tempvolume,surfacearea,ldepth;
  double remain;
  double i;
  double m;
  float index;
  double newdepth;
  int j,k, frost_area;
  double in;
  double Tnew[MAX_LAKE_NODES];
  double tempdepth;
  double wetland_runoff;
  double bpercent;
  double inflec, midrate;
  double circum;
  double snowmltvolume, evapvolume;
  double runoff_out_volume, baseflow_out_volume;
  double newfraction, Recharge;
  int ErrorFlag;
  vic411_cell_data_struct    ***cell;
  vic411_veg_var_struct      ***veg_var;
  vic411_snow_data_struct     **snow;
  vic411_energy_bal_struct    **energy;
  int lindex;
  double frac;
  double Dsmax, min_liq, liq, rel_moist;
  double *frost_frac;
  double volume_save;
  int dist, Ndist;
  double *delta_moist;
  double *moist;
  double max_newfraction;
  double depth_in_save;

  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  snow    = prcp->snow;
  energy  = prcp->energy;

#if SPATIAL_FROST
  frost_fract = soil_con.frost_fract;
#endif

  delta_moist = (double*)calloc(vic411_options.Nlayer,sizeof(double));
  moist = (double*)calloc(vic411_options.Nlayer,sizeof(double));

  if (vic411_options.DIST_PRCP)
    Ndist = 2;
  else
    Ndist = 1;

  /**********************************************************************
   * 1. convert vic411_runoff input/output volume to rate
   **********************************************************************/
  
  isave_n = lake->activenod;   /* save initial no. of nodes for later */
  
  /* Assume baseflow enters lake in same proportion as vic411_runoff. */
  bpercent = lake_con.rpercent;

  /* Convert from mm over lake area to m3. */
#if EXCESS_ICE
  if(SubsidenceUpdate > 0 ) 
    runoff_volume = ( (lake->runoff_in + lake->baseflow_in + total_meltwater) / 1000. ) * soil_con.cell_area*lake_con.Cl[0]*lakefrac;
  else
#endif
    runoff_volume = ( (lake->runoff_in + lake->baseflow_in) / 1000. ) * soil_con.cell_area*lake_con.Cl[0]*lakefrac;

  /**********************************************************************
   * 2. calculate change in lake level for lake outflow calculation
   *     grid cell vic411_runoff (m3/TS for non-lake area)
   *     snow meltwater (mm)
   *     open water evaporation (mm)
   *     precip (added in solvelake)  
   **********************************************************************/
  
  oldvolume = lake->volume;
  evapvolume = (lake->evapw/1000.) * lake->sarea;
  snowmltvolume = (lake->snowmlt/1000.) * lake->areai;
 
  // Add vic411_runoff from rest of grid cell and wetland to lake, remove evaporation
  // Evaporation is not allowed to exceed the liquid volume of the lake (after incoming vic411_runoff & baseflow are added)
  if (fabs(evapvolume) > SMALL && evapvolume > ((lake->volume-lake->ice_water_eq) + runoff_volume + snowmltvolume)) {
    evapvolume = (lake->volume-lake->ice_water_eq) + runoff_volume + snowmltvolume;
    lake->volume = lake->ice_water_eq;
  }
  else {
    lake->volume += ( runoff_volume + snowmltvolume-evapvolume);
  }

  // Estimate new lake depth and surface area of the unfrozen portion for recharge and outflow calculations
  volume_save = lake->volume;
  if(lake->ice_water_eq > (lake->volume - lake->ice_water_eq)) /* Ice is not buoyant. */
    ErrorFlag = vic411_get_depth(lake_con, lake->volume-lake->ice_water_eq, &ldepth);  
  else
    ErrorFlag = vic411_get_depth(lake_con, lake->volume, &ldepth);
  if ( ErrorFlag == ERROR ) {
    fprintf(stderr, "Something went wrong in vic411_get_depth; record = %d, volume = %f, depth = %e\n",rec,lake->volume,ldepth);
    return ( ErrorFlag );
  }
  ErrorFlag = vic411_get_sarea(lake_con, ldepth, &surfacearea);
  if ( ErrorFlag == ERROR ) {
    fprintf(stderr, "Something went wrong in vic411_get_sarea; record = %d, depth = %f, sarea = %e\n",rec,ldepth,surfacearea);
    return ( ErrorFlag );
  }

  // Estimate the new lake fraction (before recharge)
  if(lake->new_ice_area > surfacearea)
    newfraction = lake->new_ice_area/lake_con.basin[0];
  else
    newfraction = surfacearea/lake_con.basin[0];

  // Save this estimate of the new lake fraction for use later
  max_newfraction = newfraction;

  /**********************************************************************
   * 3. calculate recharge to wetland
   **********************************************************************/
  
  // Based on the above initial estimate of lake area, the lake
  // will either inundate some of the wetland or will recede and
  // expose new wetland.  In the case of inundation, we will
  // take all above-ground moisture storage (snowpack, etc) and
  // give it to the lake, but at the same time we will take water
  // from the lake and fill the inundated soil to saturation.
  // In the case of a receding lake, newly-exposed wetland is
  // assumed to have saturated soil and 0 above-ground storage.

  // The redistribution of moisture within the wetland will happen
  // during the final step of this function.

  // Note:
  // This does not account for phase changes and temperature
  // changes, so will generate some energy balance errors.
  // In addition, the moisture exchange with the wetland
  // will make our initial estimate of the new lake area
  // somewhat inaccurate (the initial estimate will be an
  // upper bound on the new lake area, resulting in some of
  // the wetland being "splashed" by the lake, causing the
  // snow to run off and the soil to moisten, but not being
  // inundated).  However, lake dimensions will be recalculated
  // after vic411_runoff and baseflow are subtracted from the lake.

  if(newfraction > lakefrac) {

    Recharge = 0.0;
    for(j=0; j<vic411_options.Nlayer; j++) {
      Recharge += ((soil_con.max_moist[j]-cell[WET][iveg][band].layer[j].moist ) / 1000. ) * (newfraction-lakefrac) * lake_con.basin[0];
    }
    Recharge -= (veg_var[WET][iveg][band].Wdew/1000. + snow[iveg][band].snow_canopy + snow[iveg][band].swq) * (newfraction-lakefrac) * lake_con.basin[0];

    // Fill the soil to saturation if possible in inundated area
    // Subtract the Recharge (which may be negative) from the lake
    if(lake->volume-lake->ice_water_eq > Recharge) {
      lake->volume -= Recharge;
      lake->recharge = Recharge*1000/(lake_con.basin[0]*lakefrac); // mm over initial lake area

      /* Update wetland soil moisture storage terms. */
      for(j=0; j<vic411_options.Nlayer; j++) {
        for ( dist = 0; dist < Ndist; dist++ ) {
          delta_moist[j] = (soil_con.max_moist[j]-cell[dist][iveg][band].layer[j].moist);
        }
      }
    }
    else {

      Recharge = (lake->volume-lake->ice_water_eq)* 1000./((1.-lakefrac)*lake_con.basin[0]) - (veg_var[WET][iveg][band].Wdew + snow[iveg][band].snow_canopy*1000. + snow[iveg][band].swq*1000.)*(newfraction-lakefrac)/(1.-lakefrac);
      lake->volume = lake->ice_water_eq;
      lake->recharge = Recharge*(1-lakefrac)/lakefrac; // mm over initial lake area

      for(j=0; j<vic411_options.Nlayer; j++) {

        if(Recharge > (soil_con.max_moist[j]-cell[WET][iveg][band].layer[j].moist)) {
          Recharge -= (soil_con.max_moist[j]-cell[WET][iveg][band].layer[j].moist);
          for ( dist = 0; dist < Ndist; dist++ ) {
            delta_moist[j] = (soil_con.max_moist[j]-cell[WET][iveg][band].layer[j].moist);
          }
        }
        else {
          for ( dist = 0; dist < Ndist; dist++ ) {
            delta_moist[j] = Recharge;
          }
          Recharge = 0.0;
        }
      }

    }

  }
  else {
    lake->recharge = 0.0;
    for(j=0; j<vic411_options.Nlayer; j++) {
      delta_moist[j] = 0.0;
    }
  }

  /**********************************************************************
   * 4. Calculate outflow from lake.  Runoff estimate is based on the
   *    equation for flow over a broad crested weir.  Baseflow estimate
   *    is the ARNO formulation, using the ice content of the adjacent
   *    wetland.  Outgoing vic411_runoff and baseflow are in meters over lake surface.
   **********************************************************************/

  Dsmax = soil_con.Dsmax / 24.;
  lindex = vic411_options.Nlayer-1;
#if SPATIAL_FROST
  min_liq = 0;
  liq = 0;
  for (frost_area=0; frost_area<FROST_SUBAREAS; frost_area++) {
    min_liq += cell[WET][iveg][band].layer[lindex].min_liq[frost_area] * soil_con.depth[lindex] * 1000. * frost_fract[frost_area];
    liq += (soil_con.max_moist[lindex] - cell[WET][iveg][band].layer[lindex].ice[frost_area])*frost_fract[frost_area];
  }
#else
  min_liq = cell[WET][iveg][band].layer[lindex].min_liq * soil_con.depth[lindex] * 1000.;
  liq = soil_con.max_moist[lindex] - cell[WET][iveg][band].layer[lindex].ice;
#endif

  /** Compute relative moisture **/
  rel_moist = (liq-min_liq) / (soil_con.max_moist[lindex]-min_liq);

  /** Compute baseflow as function of relative moisture **/
  frac = Dsmax * soil_con.Ds / soil_con.Ws;
  lake->baseflow_out = frac * rel_moist;
  if (rel_moist > soil_con.Ws) {
    frac = (rel_moist - soil_con.Ws) / (1 - soil_con.Ws);
    lake->baseflow_out += Dsmax * (1 - soil_con.Ds / soil_con.Ws) * pow(frac,soil_con.c);
  }	    
  if(lake->baseflow_out < 0) 
    lake->baseflow_out = 0;

  // extract baseflow volume from the lake m^3
  baseflow_out_volume = lake->baseflow_out*(lake_con.Cl[0]*lakefrac* soil_con.cell_area)/1000.;
  if(lake->volume -lake->ice_water_eq >= baseflow_out_volume)
    lake->volume -= baseflow_out_volume;
  else {
    baseflow_out_volume = lake->volume - lake->ice_water_eq;
    lake->volume -= baseflow_out_volume;
  }

  // Find new lake depth for vic411_runoff calculations
  if(lake->ice_water_eq > (lake->volume - lake->ice_water_eq)) /* Ice is not buoyant. */
    ErrorFlag = vic411_get_depth(lake_con, lake->volume-lake->ice_water_eq, &ldepth);  
  else
    ErrorFlag = vic411_get_depth(lake_con, lake->volume, &ldepth);
  if ( ErrorFlag == ERROR ) {
    fprintf(stderr, "Something went wrong in vic411_get_depth; record = %d, volume = %f, depth = %e\n",rec,lake->volume,ldepth);
    return ( ErrorFlag );
  }

  // Compute vic411_runoff volume in m^3 and extract vic411_runoff volume from lake
  if(ldepth <= lake_con.mindepth )
    runoff_out_volume = 0.0;
  else {
    circum=2*PI*pow(surfacearea/PI,0.5);
    runoff_out_volume = lake_con.wfrac*circum*SECPHOUR*((double)dt)*1.6*pow(ldepth-lake_con.mindepth, 1.5);
    if((lake->volume - lake->ice_water_eq) >= runoff_out_volume) { 
      /*liquid water is available */
      if( (lake->volume - runoff_out_volume) < lake_con.minvolume ) 
	runoff_out_volume = lake->volume - lake_con.minvolume;
      lake->volume -= runoff_out_volume;
    }
    else {
      runoff_out_volume = lake->volume - lake->ice_water_eq;
      if( (lake->volume - runoff_out_volume) < lake_con.minvolume ) 
	runoff_out_volume = lake->volume - lake_con.minvolume;
      lake->volume -= runoff_out_volume;
    }
  }

  // Check that lake volume does not exceed our earlier estimate.
  // This will prevent runaway lake growth for the case in which
  // the lake recharge to wetland is negative and large enough to
  // more than compensate for vic411_runoff and baseflow out of the lake.
  if (lake->volume > volume_save) {
    runoff_out_volume += lake->volume - volume_save;
    lake->volume = volume_save;
  }

  // check that lake volume does not exceed its maximum
  if (lake->volume - lake_con.maxvolume > SMALL) {
    if(lake->ice_water_eq > lake_con.maxvolume) {
      runoff_out_volume += (lake->volume - lake->ice_water_eq);
      lake->volume = lake->ice_water_eq;
    }
    else {
      runoff_out_volume += (lake->volume - lake_con.maxvolume);
      lake->volume = lake_con.maxvolume;
    }
  }
  else if (lake->volume < SMALL)
    lake->volume = 0.0;

  /**********************************************************************/
  /* End of vic411_runoff calculation */
  /**********************************************************************/

  // adjust from m^3 to mm over lake
  lake->runoff_out = 1000. * runoff_out_volume / (soil_con.cell_area * lake_con.Cl[0] * lakefrac);
  lake->baseflow_out = 1000. * baseflow_out_volume / (soil_con.cell_area * lake_con.Cl[0] * lakefrac);
  lake->evapw = 1000.* evapvolume / (soil_con.cell_area * lake_con.Cl[0] * lakefrac);

  // Recalculate lake depth
  if(lake->ice_water_eq > (lake->volume - lake->ice_water_eq)) /* Ice is not buoyant. */
    ErrorFlag = vic411_get_depth(lake_con, lake->volume-lake->ice_water_eq, &(lake->ldepth));
  else
    ErrorFlag = vic411_get_depth(lake_con, lake->volume, &(lake->ldepth));
 
  if ( ErrorFlag == ERROR ) {
    fprintf(stderr, "Something went wrong in vic411_get_depth; record = %d, volume = %f, depth = %e\n",rec,lake->volume,lake->ldepth);
    return ( ErrorFlag );
  }

  /**********************************************************************
   *  4. Adjust the activenodes and lake area array. 
   **********************************************************************/
 
  if(lake->ldepth > MAX_SURFACE_LAKE && lake->ldepth < 2*MAX_SURFACE_LAKE) {
    /* Not quite enough for two full layers. */
    lake->surfdz = lake->ldepth/2.;
    lake->dz = lake->ldepth/2.;
    lake->activenod = 2;
  }
  else if(lake->ldepth >= 2* MAX_SURFACE_LAKE) {
    /* More than two layers. */	
    lake->surfdz = MAX_SURFACE_LAKE;
    lake->activenod = (int) (lake->ldepth/MAX_SURFACE_LAKE);
    if(lake->activenod > MAX_LAKE_NODES)
      lake->activenod = MAX_LAKE_NODES;
    lake->dz = (lake->ldepth-lake->surfdz)/((float)(lake->activenod-1));
  }	
  else if(lake->ldepth > SMALL) {
    lake->surfdz = lake->ldepth;
    lake->dz = 0.0;
    lake->activenod = 1;
  }
  else {
    lake->runoff_out += ( lake->volume -lake->ice_water_eq) * 1000. / (soil_con.cell_area * lake_con.Cl[0] * lakefrac);
    lake->volume = lake->ice_water_eq;
    lake->surfdz = 0.0;
    lake->dz = 0.0;
    lake->activenod = 0;
    lake->ldepth = 0.0;
  } 	

  // lake_con.basin equals the surface area at specific depths as input by
  // the user in the lake parameter file or calculated in vic411_read_lakeparam(), 
  // lake->surface equals the area at the top of each dynamic solution layer 
  
  /* Re-calculate lake surface area  */
  for(k=0; k<=lake->activenod; k++) {
    if(k==0)
      ldepth = lake->ldepth;
    else
      ldepth = lake->dz*(lake->activenod - k);
    
    ErrorFlag = vic411_get_sarea(lake_con, ldepth, &(lake->surface[k]));
    if ( ErrorFlag == ERROR ) {
      fprintf(stderr, "Something went wrong in vic411_get_sarea; record = %d, depth = %f, sarea = %e\n",rec,ldepth,lake->surface[k]);
      return ( ErrorFlag );
    }
  }
  if (lake->new_ice_area > lake->surface[0])
    newfraction = lake->new_ice_area/lake_con.basin[0];
  else
    newfraction = lake->surface[0]/lake_con.basin[0];
  
  /*******************************************************************/  
  /* Adjust temperature distribution if number of nodes has changed. 
     Note:  This approach (and the lake model in general) does not preserve 
     the thermal energy of the water column. */
 
  index = 0.;
  if ( lake->activenod != isave_n ) {
    for(k=0; k< lake->activenod; k++) {
      Tnew[k] = 0.0;
      for(i=0; i< isave_n; i++) {
	index += (1./lake->activenod);
	Tnew[k] += lake->temp[(int)floor(index)];
      }
    }
    for(k=0; k< lake->activenod; k++) {
      if(isave_n > 0)
	lake->temp[k] = Tnew[k]/isave_n;
      else
	lake->temp[k] = prcp->energy[iveg][band].Tsurf;
    }	
  }
 
  if(lake->activenod == isave_n && isave_n == 0)
    lake->temp[k] = prcp->energy[iveg][band].Tsurf;	

  // Copy moisture fluxes into lake->soil structure
  lake->soil.vic411_runoff = lake->runoff_out;
  lake->soil.baseflow = lake->baseflow_out;
  lake->soil.inflow = lake->baseflow_out;
 
  /**********************************************************************
      5. Rescale the fluxes in the lake and the wetland by the change in lake area; advect the storages
   **********************************************************************/
  // Wetland
  if (newfraction < 1.0) { // wetland exists at end of time step
    for (dist=0; dist<Ndist; dist++) {
      vic411_advect_soil_veg_storage(lakefrac, max_newfraction, newfraction, delta_moist, &soil_con, &veg_con, &(cell[dist][iveg][band]), &(veg_var[dist][iveg][band])); 
      vic411_rescale_soil_veg_fluxes((1-lakefrac), (1-newfraction), &(cell[dist][iveg][band]), &(veg_var[dist][iveg][band])); 
    }
    vic411_advect_snow_storage(lakefrac, max_newfraction, newfraction, &(snow[iveg][band])); 
    vic411_rescale_snow_energy_fluxes((1-lakefrac), (1-newfraction), &(snow[iveg][band]), &(energy[iveg][band])); 
    for (j=0; j<vic411_options.Nlayer; j++) moist[j] = cell[0][iveg][band].layer[j].moist;
    ErrorFlag = vic411_distribute_node_moisture_properties(energy[iveg][band].moist, energy[iveg][band].ice,
                                                    energy[iveg][band].kappa_node, energy[iveg][band].Cs_node,
                                                    soil_con.Zsum_node, energy[iveg][band].T,
                                                    soil_con.max_moist_node,
#if QUICK_FS
                                                    soil_con.ufwc_table_node,
#else
                                                    soil_con.expt_node,
                                                    soil_con.bubble_node,
#endif // QUICK_FS
#if EXCESS_ICE
                                                    soil_con.porosity_node,
                                                    soil_con.effective_porosity_node,
#endif // EXCESS_ICE
                                                    moist, soil_con.depth,
                                                    soil_con.soil_density,
                                                    soil_con.bulk_density,
                                                    soil_con.quartz, vic411_options.Nnode,
                                                    vic411_options.Nlayer, soil_con.FS_ACTIVE);
    if ( ErrorFlag == ERROR ) return (ERROR);
  }

  // Lake
  if (newfraction > 0.0) { // lake exists at end of time step
    if (lakefrac > 0.0) { // lake existed at beginning of time step
      vic411_rescale_soil_veg_fluxes(lakefrac, newfraction, &(lake->soil), NULL); 
      vic411_rescale_snow_energy_fluxes(lakefrac, newfraction, &(lake->snow), &(lake->energy)); 
      vic411_rescale_lake_fluxes(lakefrac, newfraction, lake); 
    }
    else { // lake didn't exist at beginning of time step; create new lake
      depth_in_save = lake_con.depth_in;
      lake_con.depth_in = lake->ldepth;
      vic411_initialize_lake(lake, lake_con, &soil_con, energy[iveg][band].T[0]);
      lake_con.depth_in = depth_in_save;
    }
  }

  free((char*)delta_moist);
  free((char*)moist);

  return(0);

}

void vic411_advect_soil_veg_storage(double lakefrac,
                             double max_newfraction,
                             double newfraction,
                             double *delta_moist,
                             vic411_soil_con_struct  *soil_con,
                             vic411_veg_con_struct   *veg_con,
                             vic411_cell_data_struct *cell,
                             vic411_veg_var_struct   *veg_var)
{

  extern vic411_option_struct   vic411_options;
  int lidx;
  double new_moist;

  if ((1-newfraction) < SMALL) {
    newfraction = 1 - SMALL;
  }

  if (lakefrac < 1.0) { // wetland existed during this step

    for (lidx=0; lidx<vic411_options.Nlayer; lidx++) {
      if (lakefrac >= max_newfraction) {
        new_moist = cell->layer[lidx].moist*(1-lakefrac) + soil_con->max_moist[lidx]*(lakefrac-newfraction);
      }
      else {
        new_moist = cell->layer[lidx].moist*(1-max_newfraction);
        if (newfraction < max_newfraction) {
          if (newfraction > lakefrac) {
            new_moist += (cell->layer[lidx].moist+delta_moist[lidx])*(max_newfraction-newfraction);
          }
          else {
            new_moist += (cell->layer[lidx].moist+delta_moist[lidx])*(max_newfraction-lakefrac) + (soil_con->max_moist[lidx])*(lakefrac-newfraction);
          }
        }
      }
      new_moist /= (1-newfraction);
#if SPATIAL_FROST
      for (k=0; k<FROST_SUBAREAS; k++) {
        cell->layer[lidx].ice[k] *= new_moist/cell->layer[lidx].moist;
      }
#else
      cell->layer[lidx].ice     *= new_moist/cell->layer[lidx].moist;
#endif
      cell->layer[lidx].moist = new_moist;
    }

    if (newfraction < lakefrac) { // lake receded
      cell->asat = (cell->asat*(1-lakefrac) + lakefrac-newfraction) / (1-newfraction);
      if (veg_var != NULL) {
        veg_var->Wdew *= (1-lakefrac)/(1-newfraction);
      }
    }
    else {
      if (veg_var != NULL) {
        veg_var->Wdew *= (1-max_newfraction)/(1-newfraction);
      }
    }

  }
  else { // Wetland didn't exist until now; create new wetland

    for (lidx=0; lidx<vic411_options.Nlayer; lidx++) {
      cell->layer[lidx].moist = soil_con->max_moist[lidx];
#if SPATIAL_FROST
      for (k=0; k<FROST_SUBAREAS; k++) {
        cell->layer[lidx].ice[k]     = 0.0;
        cell->layer[lidx].min_liq[k] = soil_con->resid_moist[lidx];
      }
#else
      cell->layer[lidx].ice      = 0.0;
      cell->layer[lidx].min_liq  = soil_con->resid_moist[lidx];
#endif
    }
    cell->asat = 1.0;

    if (veg_var != NULL) {
      veg_var->Wdew = 0.0;
    }

  }

  // Compute rootmoist and wetness
  cell->rootmoist = 0;
  cell->wetness = 0;
  for(lidx=0;lidx<vic411_options.Nlayer;lidx++) {
    if (veg_con->root[lidx] > 0)
      cell->rootmoist += cell->layer[lidx].moist;
#if EXCESS_ICE
    cell->wetness += (cell->layer[lidx].moist - soil_con->Wpwp[lidx])/(soil_con->effective_porosity[lidx]*soil_con->depth[lidx]*1000 - soil_con->Wpwp[lidx]);
#else
    cell->wetness += (cell->layer[lidx].moist - soil_con->Wpwp[lidx])/(soil_con->porosity[lidx]*soil_con->depth[lidx]*1000 - soil_con->Wpwp[lidx]);
#endif
  }
  cell->wetness /= vic411_options.Nlayer;

}


void vic411_rescale_soil_veg_fluxes(double oldfrac,
                             double newfrac,
                             vic411_cell_data_struct *cell,
                             vic411_veg_var_struct   *veg_var)
{

  extern vic411_option_struct   vic411_options;
  int lidx;

  if (newfrac < SMALL) {
    newfrac = SMALL;
  }

  if (oldfrac > 0.0) { // existed at beginning of time step
    for (lidx=0; lidx<vic411_options.Nlayer; lidx++) {
      cell->layer[lidx].evap *= oldfrac/newfrac;
    }
    cell->baseflow *= oldfrac/newfrac;
    cell->inflow *= oldfrac/newfrac;
    cell->vic411_runoff *= oldfrac/newfrac;
    if (veg_var != NULL) {
      veg_var->canopyevap *= oldfrac/newfrac;
      veg_var->throughfall *= oldfrac/newfrac;
    }
  }
  else { // didn't exist at beginning of time step; set fluxes to 0
    for (lidx=0; lidx<vic411_options.Nlayer; lidx++) {
      cell->layer[lidx].evap = 0.0;
    }
    cell->baseflow = 0.0;
    cell->inflow = 0.0;
    cell->vic411_runoff = 0.0;
    if (veg_var != NULL) {
      veg_var->canopyevap = 0.0;
      veg_var->throughfall = 0.0;
    }
  }

}

void vic411_advect_snow_storage(double lakefrac,
                         double max_newfraction,
                         double newfraction,
                         vic411_snow_data_struct  *snow)
{

  if ((1-newfraction) < SMALL) {
    newfraction = 1 - SMALL;
  }

  if (lakefrac < 1.0) { // wetland existed during this step

    if (lakefrac >= max_newfraction) { // lake receded
      snow->depth *= (1-lakefrac)/(1-newfraction);
      snow->pack_water *= (1-lakefrac)/(1-newfraction);
      snow->snow_canopy *= (1-lakefrac)/(1-newfraction);
      snow->surf_water *= (1-lakefrac)/(1-newfraction);
      snow->swq *= (1-lakefrac)/(1-newfraction);
    }
    else {
      snow->depth *= (1-max_newfraction)/(1-newfraction);
      snow->pack_water *= (1-max_newfraction)/(1-newfraction);
      snow->snow_canopy *= (1-max_newfraction)/(1-newfraction);
      snow->surf_water *= (1-max_newfraction)/(1-newfraction);
      snow->swq *= (1-max_newfraction)/(1-newfraction);
    }

  }

  else { // Wetland didn't exist until now; create new wetland

    snow->depth = 0.0;
    snow->pack_water = 0.0;
    snow->snow_canopy = 0.0;
    snow->surf_water = 0.0;
    snow->swq = 0.0;

  }

}

void vic411_rescale_snow_energy_fluxes(double oldfrac,
                                double newfrac,
                                vic411_snow_data_struct  *snow,
                                vic411_energy_bal_struct *energy)
{

  int nidx;

  if (newfrac < SMALL) {
    newfrac = SMALL;
  }

  if (oldfrac > 0.0) { // existed at beginning of time step
    snow->blowing_flux *= oldfrac/newfrac;
    snow->melt *= oldfrac/newfrac;
    snow->surface_flux *= oldfrac/newfrac;
    snow->vapor_flux *= oldfrac/newfrac;
    energy->advected_sensible *= oldfrac/newfrac;
    energy->advection *= oldfrac/newfrac;
    energy->AtmosError *= oldfrac/newfrac;
    energy->AtmosLatent *= oldfrac/newfrac;
    energy->AtmosLatentSub *= oldfrac/newfrac;
    energy->AtmosSensible *= oldfrac/newfrac;
    energy->canopy_advection *= oldfrac/newfrac;
    energy->canopy_latent *= oldfrac/newfrac;
    energy->canopy_latent_sub *= oldfrac/newfrac;
    energy->canopy_refreeze *= oldfrac/newfrac;
    energy->canopy_sensible *= oldfrac/newfrac;
    energy->deltaCC *= oldfrac/newfrac;
    energy->deltaH *= oldfrac/newfrac;
    energy->error *= oldfrac/newfrac;
    energy->fusion *= oldfrac/newfrac;
    energy->grnd_flux *= oldfrac/newfrac;
    energy->latent *= oldfrac/newfrac;
    energy->latent_sub *= oldfrac/newfrac;
    energy->longwave *= oldfrac/newfrac;
    energy->LongOverIn *= oldfrac/newfrac;
    energy->LongUnderIn *= oldfrac/newfrac;
    energy->LongUnderOut *= oldfrac/newfrac;
    energy->melt_energy *= oldfrac/newfrac;
    energy->NetLongAtmos *= oldfrac/newfrac;
    energy->NetLongOver *= oldfrac/newfrac;
    energy->NetLongUnder *= oldfrac/newfrac;
    energy->NetShortAtmos *= oldfrac/newfrac;
    energy->NetShortGrnd *= oldfrac/newfrac;
    energy->NetShortOver *= oldfrac/newfrac;
    energy->NetShortUnder *= oldfrac/newfrac;
    energy->out_long_canopy *= oldfrac/newfrac;
    energy->out_long_surface *= oldfrac/newfrac;
    energy->refreeze_energy *= oldfrac/newfrac;
    energy->sensible *= oldfrac/newfrac;
    energy->shortwave *= oldfrac/newfrac;
    energy->ShortOverIn *= oldfrac/newfrac;
    energy->ShortUnderIn *= oldfrac/newfrac;
    energy->snow_flux *= oldfrac/newfrac;
  }
  else { // didn't exist at beginning of time step; set fluxes to 0
    snow->blowing_flux = 0.0;
    snow->melt = 0.0;
    snow->surface_flux = 0.0;
    snow->vapor_flux = 0.0;
    energy->advected_sensible = 0.0;
    energy->advection = 0.0;
    energy->AtmosError = 0.0;
    energy->AtmosLatent = 0.0;
    energy->AtmosLatentSub = 0.0;
    energy->AtmosSensible = 0.0;
    energy->canopy_advection = 0.0;
    energy->canopy_latent = 0.0;
    energy->canopy_latent_sub = 0.0;
    energy->canopy_refreeze = 0.0;
    energy->canopy_sensible = 0.0;
    energy->deltaCC = 0.0;
    energy->deltaH = 0.0;
    energy->error = 0.0;
    energy->fusion = 0.0;
    energy->grnd_flux = 0.0;
    energy->latent = 0.0;
    energy->latent_sub = 0.0;
    energy->longwave = 0.0;
    energy->LongOverIn = 0.0;
    energy->LongUnderIn = 0.0;
    energy->LongUnderOut = 0.0;
    energy->melt_energy = 0.0;
    energy->NetLongAtmos = 0.0;
    energy->NetLongOver = 0.0;
    energy->NetLongUnder = 0.0;
    energy->NetShortAtmos = 0.0;
    energy->NetShortGrnd = 0.0;
    energy->NetShortOver = 0.0;
    energy->NetShortUnder = 0.0;
    energy->out_long_canopy = 0.0;
    energy->out_long_surface = 0.0;
    energy->refreeze_energy = 0.0;
    energy->sensible = 0.0;
    energy->shortwave = 0.0;
    energy->ShortOverIn = 0.0;
    energy->ShortUnderIn = 0.0;
    energy->snow_flux = 0.0;
  }

}

void vic411_rescale_lake_fluxes(double oldfrac,
                         double newfrac,
                         vic411_lake_var_struct *lake)
{

  if (newfrac < SMALL) {
    newfrac = SMALL;
  }

  if (oldfrac > 0.0) { // existed at beginning of time step
    lake->baseflow_in *= oldfrac/newfrac;
    lake->baseflow_out *= oldfrac/newfrac;
    lake->evapw *= oldfrac/newfrac;
    lake->recharge *= oldfrac/newfrac;
    lake->runoff_in *= oldfrac/newfrac;
    lake->runoff_out *= oldfrac/newfrac;
  }
  else { // didn't exist at beginning of time step; set fluxes to 0
    lake->baseflow_in = 0.0;
    lake->baseflow_out = 0.0;
    lake->evapw = 0.0;
    lake->recharge = 0.0;
    lake->runoff_in = 0.0;
    lake->runoff_out = 0.0;
  }

}

