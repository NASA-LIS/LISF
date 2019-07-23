#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_compute_pot_evap.c,v 1.1.2.2 2009/08/22 00:02:23 vicadmin Exp $";

void vic411_compute_pot_evap(int veg_class, 
		      vic411_dmy_struct *dmy, 
		      int rec, 
		      int dt, 
		      double shortwave,
		      double net_longwave,
		      double tair, 
		      double vpd,
		      double elevation,
		      double **aero_resist,
		      double *pot_evap)
/****************************************************************************
                                                                           
  vic411_compute_pot_evap: computes potential evaporation for several different
                    reference land cover types (which are defined in
                    vic411_vicNl_def.h and vic411_global.h).

  modifications:
  2009-Aug-21 Fixed bug in assignment of albedo for natural vegetation.	TJB

****************************************************************************/
{
  extern vic411_veg_lib_struct *vic411_veg_lib;
  extern char vic411_ref_veg_ref_crop[];

  int NVegLibTypes;
  int i;
  double albedo;
  double net_short;
  double net_rad;
  double rs;
  double rarc;
  double RGL;
  double lai;
  double gsm_inv;
  char ref_crop;
  double rc;
  double ra;

  /************************************************
  Estimate and store potential evap estimates using vic411_penman equation
  ************************************************/

  NVegLibTypes = vic411_veg_lib[0].NVegLibTypes;
  for (i=0; i<N_PET_TYPES; i++) {
    if (i < N_PET_TYPES_NON_NAT) {
      rs = vic411_veg_lib[NVegLibTypes+i].rmin;
      rarc = vic411_veg_lib[NVegLibTypes+i].rarc;
      RGL = vic411_veg_lib[NVegLibTypes+i].RGL;
      lai = vic411_veg_lib[NVegLibTypes+i].LAI[dmy[rec].month-1];
      albedo = vic411_veg_lib[NVegLibTypes+i].albedo[dmy[rec].month-1];
    }
    else {
      rs = vic411_veg_lib[veg_class].rmin;
      if (i == PET_VEGNOCR) rs = 0;
      rarc = vic411_veg_lib[veg_class].rarc;
      RGL = vic411_veg_lib[veg_class].RGL;
      lai = vic411_veg_lib[veg_class].LAI[dmy[rec].month-1];
      albedo = vic411_veg_lib[veg_class].albedo[dmy[rec].month-1];
    }
    gsm_inv = 1.0;
    ref_crop = vic411_ref_veg_ref_crop[i];
    rc = vic411_calc_rc(rs, net_short, RGL, tair, vpd, lai, gsm_inv, ref_crop);
    if (i < N_PET_TYPES_NON_NAT || !vic411_veg_lib[veg_class].overstory)
      ra = aero_resist[i][0];
    else
      ra = aero_resist[i][1];
    net_short = (1.0 - albedo) * shortwave;
    net_rad = net_short + net_longwave;
    pot_evap[i] = vic411_penman(tair, elevation, net_rad, vpd, ra, rc, rarc) * dt/24.0;
  }

}
