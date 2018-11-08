#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>
#include <math.h>

static char vcid[] = "$Id: store_moisture_for_debug.c,v 4.2 2000/05/16 21:57:54 vicadmin Exp $";

#if LINK_DEBUG
void store_moisture_for_debug(int                 iveg,
		              int                 Nveg,
			      double             *mu,
			      vic411_cell_data_struct ***cell,
			      vic411_veg_var_struct   ***veg_var,
			      vic411_snow_data_struct  **snow,
			      vic411_soil_con_struct    *soil_con) {
/****************************************************************
  This subroutine was written to save the current water storage
  terms for use in calculating the model water balance error
****************************************************************/

  extern vic411_option_struct vic411_options;
  extern vic411_debug_struct  debug;

  int               Ndist;
  int               i;
  int               band;
  int               dist;
  int               Nbands;

  if(vic411_options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
  Nbands = vic411_options.SNOW_BAND;

  for(band=0;band<Nbands;band++) {
    if(soil_con->AreaFract[band]>0) {
      for(dist=0;dist<Ndist;dist++) 
	for(i=0;i<vic411_options.Nlayer+3;i++)
	  debug.store_moist[dist][band][i] = 0.;
      if(iveg<Nveg) {
	for(dist=0;dist<Ndist;dist++)
	  debug.store_moist[dist][band][0] 
	    += veg_var[dist][iveg][band].Wdew;
	debug.store_moist[WET][band][0] += (snow[iveg][band].snow_canopy) 
	  * 1000.;
      }
      for(dist=0;dist<Ndist;dist++)
	debug.store_moist[dist][band][vic411_options.Nlayer+2] 
	  += debug.store_moist[dist][band][0];
      debug.store_moist[WET][band][1] += (snow[iveg][band].swq*1000.);
      for(dist=0;dist<Ndist;dist++)
	debug.store_moist[dist][band][vic411_options.Nlayer+2] 
	  += debug.store_moist[dist][band][1];
      for(i=0;i<vic411_options.Nlayer;i++) {
	for(dist=0;dist<Ndist;dist++) {
	  debug.store_moist[dist][band][i+2] 
	    = cell[dist][iveg][band].layer[i].moist;
	  debug.store_moist[dist][band][vic411_options.Nlayer+2] 
	    += debug.store_moist[dist][band][i+2];
	}
      }
    }
  }
}
#endif

