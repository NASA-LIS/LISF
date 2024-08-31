//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <stdio.h>
#include <stdlib.h>
#include "vic411_vicNl.h"

//BOP
//
// !ROUTINE: write_prcp
// \label{write_prcp411}
// 
// !REVISION HISTORY: 
//  Aug 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
void vic411_write_prcp(vic411_dist_prcp_struct *prcp, int Nveg)
// !DESCRIPTION: 
//  This routine writes contents of the prcp data structure to stdout, 
//  used primarily for debugging.
//
//  The arguments are:
//  \begin{description}
//  \item[prcp] given prcp data structure
//  \item[Nveg] number of vegetation types
//  \end{description}
//
//EOP
{
   extern vic411_option_struct vic411_options;

   int dist, veg, band, lidx, i;
   int Ndist, Nbands;

   double tmpval;
   FILE * fp;

   if(vic411_options.DIST_PRCP) 
      Ndist = 2;
   else 
      Ndist = 1;
   Nbands = vic411_options.SNOW_BAND;

   fp=fopen("prcp.log","a");
   fprintf(fp,"Prcp Values:\n");

   for ( dist = 0; dist < Ndist; ++dist )
   {
      for ( veg = 0; veg <= Nveg; ++veg )
      {
         for ( band = 0; band < Nbands; ++band )
         {
            for ( i = 0; i < 2; ++i )
            {
               tmpval = prcp->cell[dist][veg][band].aero_resist[i];
               fprintf(fp,"\tcell[%d][%d][%d].aero_resist[%d]     = %f\n", dist, veg, band, i, tmpval);
            }

            tmpval = prcp->cell[dist][veg][band].asat;
            fprintf(fp,"\tcell[%d][%d][%d].asat                   = %f\n", dist, veg, band, tmpval);

            tmpval = prcp->cell[dist][veg][band].baseflow;
            fprintf(fp,"\tcell[%d][%d][%d].baseflow               = %f\n", dist, veg, band, tmpval);

            tmpval = prcp->cell[dist][veg][band].inflow;
            fprintf(fp,"\tcell[%d][%d][%d].inflow                 = %f\n", dist, veg, band, tmpval);

            for ( i = 0; i < N_PET_TYPES; ++i )
            {
               tmpval = prcp->cell[dist][veg][band].pot_evap[i];
               fprintf(fp,"\tcell[%d][%d][%d].pot_evap[%d]        = %f\n", dist, veg, band, i, tmpval);
            }

            for ( lidx = 0; lidx < vic411_options.Nlayer; ++lidx )
            {
               tmpval = prcp->cell[dist][veg][band].layer[lidx].Cs;
               fprintf(fp,"\tcell[%d][%d][%d].layer[%d].Cs                  = %f\n", dist, veg, band, lidx, tmpval);

               tmpval = prcp->cell[dist][veg][band].layer[lidx].T;
               fprintf(fp,"\tcell[%d][%d][%d].layer[%d].T                   = %f\n", dist, veg, band, lidx, tmpval);

               tmpval = prcp->cell[dist][veg][band].layer[lidx].evap;
               fprintf(fp,"\tcell[%d][%d][%d].layer[%d].evap                = %f\n", dist, veg, band, lidx, tmpval);

               tmpval = prcp->cell[dist][veg][band].layer[lidx].ice;
               fprintf(fp,"\tcell[%d][%d][%d].layer[%d].ice                 = %f\n", dist, veg, band, lidx, tmpval);

               tmpval = prcp->cell[dist][veg][band].layer[lidx].min_liq;
               fprintf(fp,"\tcell[%d][%d][%d].layer[%d].min_liq             = %f\n", dist, veg, band, lidx, tmpval);

               tmpval = prcp->cell[dist][veg][band].layer[lidx].kappa;
               fprintf(fp,"\tcell[%d][%d][%d].layer[%d].kappa               = %f\n", dist, veg, band, lidx, tmpval);

               tmpval = prcp->cell[dist][veg][band].layer[lidx].moist;
               fprintf(fp,"\tcell[%d][%d][%d].layer[%d].moist               = %f\n", dist, veg, band, lidx, tmpval);

               tmpval = prcp->cell[dist][veg][band].layer[lidx].phi;
               fprintf(fp,"\tcell[%d][%d][%d].layer[%d].phi                 = %f\n", dist, veg, band, lidx, tmpval);
            }

            tmpval = prcp->cell[dist][veg][band].rootmoist;
            fprintf(fp,"\tcell[%d][%d][%d].rootmoist              = %f\n", dist, veg, band, tmpval);

            tmpval = prcp->cell[dist][veg][band].wetness;
            fprintf(fp,"\tcell[%d][%d][%d].wetness                = %f\n", dist, veg, band, tmpval);
         }
      }
   }

   for ( veg = 0; veg <= Nveg; ++veg )
   {
      tmpval = prcp->mu[veg];
      fprintf(fp,"\tmu[%d] = %f\n",  veg, tmpval);
   }

   for ( veg = 0; veg <= Nveg; ++veg )
   {
      for ( band = 0; band < Nbands; ++band )
      {
         tmpval = prcp->energy[veg][band].AlbedoLake;
         fprintf(fp,"\tenergy[%d][%d].AlbedoLake                  = %f\n", veg, band, tmpval);

         tmpval = prcp->energy[veg][band].AlbedoOver;
         fprintf(fp,"\tenergy[%d][%d].AlbedoOver                  = %f\n", veg, band, tmpval);

         tmpval = prcp->energy[veg][band].AlbedoUnder;
         fprintf(fp,"\tenergy[%d][%d].AlbedoUnder;                = %f\n", veg, band, tmpval);

         for ( i = 0; i < 2; ++i )
         {
            tmpval = prcp->energy[veg][band].Cs[i];
            fprintf(fp,"\tenergy[%d][%d].Cs[%d]                   = %f\n", veg, band, i, tmpval);
         }

         for ( i = 0; i < MAX_NODES; ++i )
         {
            tmpval = prcp->energy[veg][band].Cs_node[i];
            fprintf(fp,"\tenergy[%d][%d].Cs_node[%d]              = %f\n", veg, band, i, tmpval);
         }

         for ( i = 0; i < MAX_FRONTS; ++i )
         {
            tmpval = prcp->energy[veg][band].fdepth[i];
            fprintf(fp,"\tenergy[%d][%d].fdepth[%d]               = %f\n", veg, band, i, tmpval);
         }

         tmpval = prcp->energy[veg][band].frozen;
         fprintf(fp,"\tenergy[%d][%d].frozen;                     = %f\n", veg, band, tmpval);

         for ( i = 0; i < MAX_NODES; ++i )
         {
            tmpval = prcp->energy[veg][band].ice[i];
            fprintf(fp,"\tenergy[%d][%d].ice[%d]                  = %f\n", veg, band, i, tmpval);
         }

         for ( i = 0; i < 2; ++i )
         {
            tmpval = prcp->energy[veg][band].kappa[i];
            fprintf(fp,"\tenergy[%d][%d].kappa[%d]                = %f\n", veg, band, i, tmpval);
         }

         for ( i = 0; i < MAX_NODES; ++i )
         {
            tmpval = prcp->energy[veg][band].kappa_node[i];
            fprintf(fp,"\tenergy[%d][%d].kappa_node[%d]           = %f\n", veg, band, i, tmpval);
         }

         for ( i = 0; i < MAX_NODES; ++i )
         {
            tmpval = prcp->energy[veg][band].moist[i];
            fprintf(fp,"\tenergy[%d][%d].moist[%d]                = %f\n", veg, band, i, tmpval);
         }

         tmpval = prcp->energy[veg][band].Nfrost;
         fprintf(fp,"\tenergy[%d][%d].Nfrost;                     = %f\n", veg, band, tmpval);

         tmpval = prcp->energy[veg][band].Nthaw;
         fprintf(fp,"\tenergy[%d][%d].Nthaw;                      = %f\n", veg, band, tmpval);

         for ( i = 0; i < MAX_NODES; ++i )
         {
            tmpval = prcp->energy[veg][band].T[i];
            fprintf(fp,"\tenergy[%d][%d].T[%d]                    = %f\n", veg, band, i, tmpval);
         }

         for ( i = 0; i < MAX_NODES; ++i )
         {
            tmpval = prcp->energy[veg][band].T_fbflag[i];
            fprintf(fp,"\tenergy[%d][%d].T_fbflag[%d]             = %f\n", veg, band, i, tmpval);
         }

         for ( i = 0; i < MAX_NODES; ++i )
         {
            tmpval = prcp->energy[veg][band].T_fbcount[i];
            fprintf(fp,"\tenergy[%d][%d].T_fbcount[%d]            = %f\n", veg, band, i, tmpval);
         }

         tmpval = prcp->energy[veg][band].T1_index;
         fprintf(fp,"\tenergy[%d][%d].T1_index;                   = %f\n", veg, band, tmpval);

         tmpval = prcp->energy[veg][band].Tcanopy;
         fprintf(fp,"\tenergy[%d][%d].Tcanopy;                    = %f\n", veg, band, tmpval);

         tmpval = prcp->energy[veg][band].Tcanopy_fbflag;
         fprintf(fp,"\tenergy[%d][%d].Tcanopy_fbflag              = %f\n", veg, band, tmpval);

         tmpval = prcp->energy[veg][band].Tcanopy_fbcount;
         fprintf(fp,"\tenergy[%d][%d].Tcanopy_fbcount             = %f\n", veg, band, tmpval);

         for ( i = 0; i < MAX_FRONTS; ++i )
         {
            tmpval = prcp->energy[veg][band].tdepth[i];
            fprintf(fp,"\tenergy[%d][%d].tdepth[%d]               = %f\n", veg, band, i, tmpval);
         }

         tmpval = prcp->energy[veg][band].Tfoliage;
         fprintf(fp,"\tenergy[%d][%d].Tfoliage                    = %f\n", veg, band, tmpval);

         tmpval = prcp->energy[veg][band].Tfoliage_fbflag;
         fprintf(fp,"\tenergy[%d][%d].Tfoliage_fbflag             = %f\n", veg, band, tmpval);

         tmpval = prcp->energy[veg][band].Tfoliage_fbcount;
         fprintf(fp,"\tenergy[%d][%d].Tfoliage_fbcount            = %f\n", veg, band, tmpval);

         tmpval = prcp->energy[veg][band].Tsurf;
         fprintf(fp,"\tenergy[%d][%d].Tfoliage                    = %f\n", veg, band, tmpval);

         tmpval = prcp->energy[veg][band].Tsurf_fbflag;
         fprintf(fp,"\tenergy[%d][%d].Tsurf_fbflag                = %f\n", veg, band, tmpval);

         tmpval = prcp->energy[veg][band].Tsurf_fbcount;
         fprintf(fp,"\tenergy[%d][%d].Tsurf_fbcount               = %f\n", veg, band, tmpval);

         tmpval = prcp->energy[veg][band].unfrozen;
         fprintf(fp,"\tenergy[%d][%d].unfrozen         = %f\n", veg, band, tmpval);              

         tmpval = prcp->energy[veg][band].advected_sensible;
         fprintf(fp,"\tenergy[%d][%d].advected_sensible         = %f\n", veg, band, tmpval);     

         tmpval = prcp->energy[veg][band].advection;
         fprintf(fp,"\tenergy[%d][%d].advection         = %f\n", veg, band, tmpval);             

         tmpval = prcp->energy[veg][band].AtmosError;
         fprintf(fp,"\tenergy[%d][%d].AtmosError         = %f\n", veg, band, tmpval);

         tmpval = prcp->energy[veg][band].AtmosLatent;
         fprintf(fp,"\tenergy[%d][%d].AtmosLatent         = %f\n", veg, band, tmpval);           

         tmpval = prcp->energy[veg][band].AtmosLatentSub;
         fprintf(fp,"\tenergy[%d][%d].AtmosLatentSub         = %f\n", veg, band, tmpval);        

         tmpval = prcp->energy[veg][band].AtmosSensible;
         fprintf(fp,"\tenergy[%d][%d].AtmosSensible         = %f\n", veg, band, tmpval);         

         tmpval = prcp->energy[veg][band].canopy_advection;
         fprintf(fp,"\tenergy[%d][%d].canopy_advection         = %f\n", veg, band, tmpval);      

         tmpval = prcp->energy[veg][band].canopy_latent;
         fprintf(fp,"\tenergy[%d][%d].canopy_latent         = %f\n", veg, band, tmpval);         

         tmpval = prcp->energy[veg][band].canopy_latent_sub;
         fprintf(fp,"\tenergy[%d][%d].canopy_latent_sub         = %f\n", veg, band, tmpval);     

         tmpval = prcp->energy[veg][band].canopy_refreeze;
         fprintf(fp,"\tenergy[%d][%d].canopy_refreeze         = %f\n", veg, band, tmpval);       

         tmpval = prcp->energy[veg][band].canopy_sensible;
         fprintf(fp,"\tenergy[%d][%d].canopy_sensible         = %f\n", veg, band, tmpval);       

         tmpval = prcp->energy[veg][band].deltaCC;
         fprintf(fp,"\tenergy[%d][%d].deltaCC         = %f\n", veg, band, tmpval);               

         tmpval = prcp->energy[veg][band].deltaH;
         fprintf(fp,"\tenergy[%d][%d].deltaH         = %f\n", veg, band, tmpval);                

         tmpval = prcp->energy[veg][band].error;
         fprintf(fp,"\tenergy[%d][%d].error         = %f\n", veg, band, tmpval);                 

         tmpval = prcp->energy[veg][band].fusion;
         fprintf(fp,"\tenergy[%d][%d].fusion         = %f\n", veg, band, tmpval);                

         tmpval = prcp->energy[veg][band].grnd_flux;
         fprintf(fp,"\tenergy[%d][%d].grnd_flux         = %f\n", veg, band, tmpval);             

         tmpval = prcp->energy[veg][band].latent;
         fprintf(fp,"\tenergy[%d][%d].latent         = %f\n", veg, band, tmpval);                

         tmpval = prcp->energy[veg][band].latent_sub;
         fprintf(fp,"\tenergy[%d][%d].latent_sub         = %f\n", veg, band, tmpval);            

         tmpval = prcp->energy[veg][band].longwave;
         fprintf(fp,"\tenergy[%d][%d].longwave         = %f\n", veg, band, tmpval);              

         tmpval = prcp->energy[veg][band].LongOverIn;
         fprintf(fp,"\tenergy[%d][%d].LongOverIn         = %f\n", veg, band, tmpval);            

         tmpval = prcp->energy[veg][band].LongUnderIn;
         fprintf(fp,"\tenergy[%d][%d].LongUnderIn         = %f\n", veg, band, tmpval);           

         tmpval = prcp->energy[veg][band].LongUnderOut;
         fprintf(fp,"\tenergy[%d][%d].LongUnderOut         = %f\n", veg, band, tmpval);          

         tmpval = prcp->energy[veg][band].melt_energy;
         fprintf(fp,"\tenergy[%d][%d].melt_energy         = %f\n", veg, band, tmpval);           

         tmpval = prcp->energy[veg][band].NetLongAtmos;
         fprintf(fp,"\tenergy[%d][%d].NetLongAtmos         = %f\n", veg, band, tmpval);          

         tmpval = prcp->energy[veg][band].NetLongOver;
         fprintf(fp,"\tenergy[%d][%d].NetLongOver         = %f\n", veg, band, tmpval);           

         tmpval = prcp->energy[veg][band].NetLongUnder;
         fprintf(fp,"\tenergy[%d][%d].NetLongUnder         = %f\n", veg, band, tmpval);          

         tmpval = prcp->energy[veg][band].NetShortAtmos;
         fprintf(fp,"\tenergy[%d][%d].NetShortAtmos         = %f\n", veg, band, tmpval);         

         tmpval = prcp->energy[veg][band].NetShortGrnd;
         fprintf(fp,"\tenergy[%d][%d].NetShortGrnd         = %f\n", veg, band, tmpval);          

         tmpval = prcp->energy[veg][band].NetShortOver;
         fprintf(fp,"\tenergy[%d][%d].NetShortOver         = %f\n", veg, band, tmpval);          

         tmpval = prcp->energy[veg][band].NetShortUnder;
         fprintf(fp,"\tenergy[%d][%d].NetShortUnder         = %f\n", veg, band, tmpval);         

         tmpval = prcp->energy[veg][band].out_long_canopy;
         fprintf(fp,"\tenergy[%d][%d].out_long_canopy         = %f\n", veg, band, tmpval);       

         tmpval = prcp->energy[veg][band].out_long_surface;
         fprintf(fp,"\tenergy[%d][%d].out_long_surface         = %f\n", veg, band, tmpval);      

         tmpval = prcp->energy[veg][band].refreeze_energy;
         fprintf(fp,"\tenergy[%d][%d].refreeze_energy         = %f\n", veg, band, tmpval);       

         tmpval = prcp->energy[veg][band].sensible;
         fprintf(fp,"\tenergy[%d][%d].sensible         = %f\n", veg, band, tmpval);              

         tmpval = prcp->energy[veg][band].shortwave;
         fprintf(fp,"\tenergy[%d][%d].shortwave         = %f\n", veg, band, tmpval);             

         tmpval = prcp->energy[veg][band].ShortOverIn;
         fprintf(fp,"\tenergy[%d][%d].ShortOverIn         = %f\n", veg, band, tmpval);           

         tmpval = prcp->energy[veg][band].ShortUnderIn;
         fprintf(fp,"\tenergy[%d][%d].ShortUnderIn         = %f\n", veg, band, tmpval);          

         tmpval = prcp->energy[veg][band].snow_flux;
         fprintf(fp,"\tenergy[%d][%d].snow_flux         = %f\n", veg, band, tmpval);             
      }
   }

   tmpval = prcp->lake_var.activenod;
   fprintf(fp,"\tlake_var.activenod = %f\n", tmpval);

   tmpval = prcp->lake_var.dz;
   fprintf(fp,"\tlake_var.dz = %f\n", tmpval);

   tmpval = prcp->lake_var.surfdz;
   fprintf(fp,"\tlake_var.surfdz = %f\n", tmpval);

   tmpval = prcp->lake_var.ldepth;
   fprintf(fp,"\tlake_var.ldepth = %f\n", tmpval);

   for ( i = 0; i < MAX_LAKE_NODES+1; ++i )
   {
      tmpval = prcp->lake_var.surface[i];
      fprintf(fp,"\tlake_var.surface[%d] = %f\n", i, tmpval);
   }

   tmpval = prcp->lake_var.sarea;
   fprintf(fp,"\tlake_var.sarea = %f\n", tmpval);

   tmpval = prcp->lake_var.volume;
   fprintf(fp,"\tlake_var.volume = %f\n", tmpval);

   for ( i = 0; i < MAX_LAKE_NODES; ++i )
   {
      tmpval = prcp->lake_var.temp[i];
      fprintf(fp,"\tlake_var.temp[%d] = %f\n", i, tmpval);
   }

   tmpval = prcp->lake_var.tempavg;
   fprintf(fp,"\tlake_var.tempavg = %f\n", tmpval);

   tmpval = prcp->lake_var.areai;
   fprintf(fp,"\tlake_var.areai = %f\n", tmpval);

   tmpval = prcp->lake_var.new_ice_area;
   fprintf(fp,"\tlake_var.new_ice_area = %f\n", tmpval);

   tmpval = prcp->lake_var.ice_water_eq;
   fprintf(fp,"\tlake_var.ice_water_eq = %f\n", tmpval);

   tmpval = prcp->lake_var.hice;
   fprintf(fp,"\tlake_var.hice = %f\n", tmpval);

   tmpval = prcp->lake_var.tempi;
   fprintf(fp,"\tlake_var.tempi = %f\n", tmpval);

   tmpval = prcp->lake_var.swe;
   fprintf(fp,"\tlake_var.swe = %f\n", tmpval);

   tmpval = prcp->lake_var.surf_temp;
   fprintf(fp,"\tlake_var.surf_temp = %f\n", tmpval);

   tmpval = prcp->lake_var.pack_temp;
   fprintf(fp,"\tlake_var.pack_temp = %f\n", tmpval);

   tmpval = prcp->lake_var.coldcontent;
   fprintf(fp,"\tlake_var.coldcontent = %f\n", tmpval);

   tmpval = prcp->lake_var.surf_water;
   fprintf(fp,"\tlake_var.surf_water = %f\n", tmpval);

   tmpval = prcp->lake_var.pack_water;
   fprintf(fp,"\tlake_var.pack_water = %f\n", tmpval);

   tmpval = prcp->lake_var.SAlbedo;
   fprintf(fp,"\tlake_var.SAlbedo = %f\n", tmpval);

   tmpval = prcp->lake_var.sdepth;
   fprintf(fp,"\tlake_var.sdepth = %f\n", tmpval);

   tmpval = prcp->lake_var.aero_resist;
   fprintf(fp,"\tlake_var.aero_resist = %f\n", tmpval);

   for ( i = 0; i < MAX_LAKE_NODES; ++i )
   {
      tmpval = prcp->lake_var.density[i];
      fprintf(fp,"\tlake_var.density[%d] = %f\n", i, tmpval);
   }

   tmpval = prcp->lake_var.baseflow_in;
   fprintf(fp,"\tlake_var.baseflow_in = %f\n", tmpval);

   tmpval = prcp->lake_var.baseflow_out;
   fprintf(fp,"\tlake_var.baseflow_out = %f\n", tmpval);

   tmpval = prcp->lake_var.evapw;
   fprintf(fp,"\tlake_var.evapw = %f\n", tmpval);

   tmpval = prcp->lake_var.recharge;
   fprintf(fp,"\tlake_var.recharge = %f\n", tmpval);

   tmpval = prcp->lake_var.runoff_in;
   fprintf(fp,"\tlake_var.runoff_in = %f\n", tmpval);

   tmpval = prcp->lake_var.runoff_out;
   fprintf(fp,"\tlake_var.runoff_out = %f\n", tmpval);

   tmpval = prcp->lake_var.snowmlt;
   fprintf(fp,"\tlake_var.snowmlt = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.albedo;
   fprintf(fp,"\tlake_var.snow.albedo = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.canopy_albedo;
   fprintf(fp,"\tlake_var.snow.canopy_albedo = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.coldcontent;
   fprintf(fp,"\tlake_var.snow.coldcontent = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.coverage;
   fprintf(fp,"\tlake_var.snow.coverage = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.density;
   fprintf(fp,"\tlake_var.snow.density = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.depth;
   fprintf(fp,"\tlake_var.snow.depth = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.last_snow;
   fprintf(fp,"\tlake_var.snow.last_snow = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.max_swq;
   fprintf(fp,"\tlake_var.snow.max_swq = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.MELTING;
   fprintf(fp,"\tlake_var.snow.MELTING = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.pack_temp;
   fprintf(fp,"\tlake_var.snow.pack_temp = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.pack_water;
   fprintf(fp,"\tlake_var.snow.pack_water = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.snow;
   fprintf(fp,"\tlake_var.snow.snow = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.snow_canopy;
   fprintf(fp,"\tlake_var.snow.snow_canopy = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.store_coverage;
   fprintf(fp,"\tlake_var.snow.store_coverage = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.store_snow;
   fprintf(fp,"\tlake_var.snow.store_snow = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.store_swq;
   fprintf(fp,"\tlake_var.snow.store_swq = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.surf_temp;
   fprintf(fp,"\tlake_var.snow.surf_temp = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.surf_temp_fbcount;
   fprintf(fp,"\tlake_var.snow.surf_temp_fbcount = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.surf_temp_fbflag;
   fprintf(fp,"\tlake_var.snow.surf_temp_fbflag = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.surf_water;
   fprintf(fp,"\tlake_var.snow.surf_water = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.swq;
   fprintf(fp,"\tlake_var.snow.swq = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.swq_slope;
   fprintf(fp,"\tlake_var.snow.swq_slope = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.tmp_int_storage;
   fprintf(fp,"\tlake_var.snow.tmp_int_storage = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.blowing_flux;
   fprintf(fp,"\tlake_var.snow.blowing_flux = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.canopy_vapor_flux;
   fprintf(fp,"\tlake_var.snow.canopy_vapor_flux = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.mass_error;
   fprintf(fp,"\tlake_var.snow.mass_error = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.melt;
   fprintf(fp,"\tlake_var.snow.melt = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.Qnet;
   fprintf(fp,"\tlake_var.snow.Qnet = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.surface_flux;
   fprintf(fp,"\tlake_var.snow.surface_flux = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.transport;
   fprintf(fp,"\tlake_var.snow.transport = %f\n", tmpval);

   tmpval = prcp->lake_var.snow.vapor_flux;
   fprintf(fp,"\tlake_var.snow.vapor_flux = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.AlbedoLake;
   fprintf(fp,"\tlake_var.energy.AlbedoLake = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.AlbedoOver;
   fprintf(fp,"\tlake_var.energy.AlbedoOver = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.AlbedoUnder;
   fprintf(fp,"\tlake_var.energy.AlbedoUnder = %f\n", tmpval);

   for ( i = 0; i < 2; ++i )
   {
      tmpval = prcp->lake_var.energy.Cs[i];
      fprintf(fp,"\tlake_var.energy.Cs[%d] = %f\n", i, tmpval);
   }

   for ( i = 0; i < MAX_NODES; ++i )
   {
      tmpval = prcp->lake_var.energy.Cs_node[i];
      fprintf(fp,"\tlake_var.energy.Cs_node[%d] = %f\n", i, tmpval);
   }

   for ( i = 0; i < MAX_FRONTS; ++i )
   {
      tmpval = prcp->lake_var.energy.fdepth[i];
      fprintf(fp,"\tlake_var.energy.fdepth[%d] = %f\n", i, tmpval);
   }

   tmpval = prcp->lake_var.energy.frozen;
   fprintf(fp,"\tlake_var.energy.frozen = %f\n", tmpval);

   for ( i = 0; i < MAX_NODES; ++i )
   {
      tmpval = prcp->lake_var.energy.ice[i];
      fprintf(fp,"\tlake_var.energy.ice[%d] = %f\n", i, tmpval);
   }

   for ( i = 0; i < 2; ++i )
   {
      tmpval = prcp->lake_var.energy.kappa[i];
      fprintf(fp,"\tlake_var.energy.kappa[%d] = %f\n", i, tmpval);
   }

   for ( i = 0; i < MAX_NODES; ++i )
   {
      tmpval = prcp->lake_var.energy.kappa_node[i];
      fprintf(fp,"\tlake_var.energy.kappa_node[%d] = %f\n", i, tmpval);
   }

   for ( i = 0; i < MAX_NODES; ++i )
   {
      tmpval = prcp->lake_var.energy.moist[i];
      fprintf(fp,"\tlake_var.energy.moist[%d] = %f\n", i, tmpval);
   }

   tmpval = prcp->lake_var.energy.Nfrost;
   fprintf(fp,"\tlake_var.energy.Nfrost = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.Nthaw;
   fprintf(fp,"\tlake_var.energy.Nthaw = %f\n", tmpval);

   for ( i = 0; i < MAX_NODES; ++i )
   {
      tmpval = prcp->lake_var.energy.T[i];
      fprintf(fp,"\tlake_var.energy.T[%d] = %f\n", i, tmpval);
   }

   for ( i = 0; i < MAX_NODES; ++i )
   {
      tmpval = prcp->lake_var.energy.T_fbflag[i];
      fprintf(fp,"\tlake_var.energy.T_fbflag[%d] = %f\n", i, tmpval);
   }

   for ( i = 0; i < MAX_NODES; ++i )
   {
      tmpval = prcp->lake_var.energy.T_fbcount[i];
      fprintf(fp,"\tlake_var.energy.T_fbcount[%d] = %f\n", i, tmpval);
   }

   tmpval = prcp->lake_var.energy.T1_index;
   fprintf(fp,"\tlake_var.energy.T1_index = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.Tcanopy;
   fprintf(fp,"\tlake_var.energy.Tcanopy = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.Tcanopy_fbflag;
   fprintf(fp,"\tlake_var.energy.Tcanopy_fbflag = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.Tcanopy_fbcount;
   fprintf(fp,"\tlake_var.energy.Tcanopy_fbcount = %f\n", tmpval);

   for ( i = 0; i < MAX_FRONTS; ++i )
   {
      tmpval = prcp->lake_var.energy.tdepth[i];
      fprintf(fp,"\tlake_var.energy.tdepth[%d] = %f\n", i, tmpval);
   }

   tmpval = prcp->lake_var.energy.Tfoliage;
   fprintf(fp,"\tlake_var.energy.Tfoliage = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.Tfoliage_fbflag;
   fprintf(fp,"\tlake_var.energy.Tfoliage_fbflag = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.Tfoliage_fbcount;
   fprintf(fp,"\tlake_var.energy.Tfoliage_fbcount = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.Tsurf;
   fprintf(fp,"\tlake_var.energy.Tsurf = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.Tsurf_fbflag;
   fprintf(fp,"\tlake_var.energy.Tsurf_fbflag = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.Tsurf_fbcount;
   fprintf(fp,"\tlake_var.energy.Tsurf_fbcount = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.unfrozen;
   fprintf(fp,"\tlake_var.energy.unfrozen = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.advected_sensible;
   fprintf(fp,"\tlake_var.energy.advected_sensible = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.advection;
   fprintf(fp,"\tlake_var.energy.advection = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.AtmosError;
   fprintf(fp,"\tlake_var.energy.AtmosError = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.AtmosLatent;
   fprintf(fp,"\tlake_var.energy.AtmosLatent = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.AtmosLatentSub;
   fprintf(fp,"\tlake_var.energy.AtmosLatentSub = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.AtmosSensible;
   fprintf(fp,"\tlake_var.energy.AtmosSensible = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.canopy_advection;
   fprintf(fp,"\tlake_var.energy.canopy_advection = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.canopy_latent;
   fprintf(fp,"\tlake_var.energy.canopy_latent = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.canopy_latent_sub;
   fprintf(fp,"\tlake_var.energy.canopy_latent_sub = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.canopy_refreeze;
   fprintf(fp,"\tlake_var.energy.canopy_refreeze = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.canopy_sensible;
   fprintf(fp,"\tlake_var.energy.canopy_sensible = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.deltaCC;
   fprintf(fp,"\tlake_var.energy.deltaCC = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.deltaH;
   fprintf(fp,"\tlake_var.energy.deltaH = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.error;
   fprintf(fp,"\tlake_var.energy.error = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.fusion;
   fprintf(fp,"\tlake_var.energy.fusion = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.grnd_flux;
   fprintf(fp,"\tlake_var.energy.grnd_flux = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.latent;
   fprintf(fp,"\tlake_var.energy.latent = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.latent_sub;
   fprintf(fp,"\tlake_var.energy.latent_sub = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.longwave;
   fprintf(fp,"\tlake_var.energy.longwave = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.LongOverIn;
   fprintf(fp,"\tlake_var.energy.LongOverIn = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.LongUnderIn;
   fprintf(fp,"\tlake_var.energy.LongUnderIn = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.LongUnderOut;
   fprintf(fp,"\tlake_var.energy.LongUnderOut = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.melt_energy;
   fprintf(fp,"\tlake_var.energy.melt_energy = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.NetLongAtmos;
   fprintf(fp,"\tlake_var.energy.NetLongAtmos = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.NetLongOver;
   fprintf(fp,"\tlake_var.energy.NetLongOver = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.NetLongUnder;
   fprintf(fp,"\tlake_var.energy.NetLongUnder = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.NetShortAtmos;
   fprintf(fp,"\tlake_var.energy.NetShortAtmos = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.NetShortGrnd;
   fprintf(fp,"\tlake_var.energy.NetShortGrnd = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.NetShortOver;
   fprintf(fp,"\tlake_var.energy.NetShortOver = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.NetShortUnder;
   fprintf(fp,"\tlake_var.energy.NetShortUnder = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.out_long_canopy;
   fprintf(fp,"\tlake_var.energy.out_long_canopy = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.out_long_surface;
   fprintf(fp,"\tlake_var.energy.out_long_surface = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.refreeze_energy;
   fprintf(fp,"\tlake_var.energy.refreeze_energy = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.sensible;
   fprintf(fp,"\tlake_var.energy.sensible = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.shortwave;
   fprintf(fp,"\tlake_var.energy.shortwave = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.ShortOverIn;
   fprintf(fp,"\tlake_var.energy.ShortOverIn = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.ShortUnderIn;
   fprintf(fp,"\tlake_var.energy.ShortUnderIn = %f\n", tmpval);

   tmpval = prcp->lake_var.energy.snow_flux;
   fprintf(fp,"\tlake_var.energy.snow_flux = %f\n", tmpval);

   for ( i = 0; i < 2; ++i )
   {
      tmpval = prcp->lake_var.soil.aero_resist[i];
      fprintf(fp,"\tlake_var.soil.aero_resist[%d] = %f\n", i, tmpval);
   }

   tmpval = prcp->lake_var.soil.asat;
   fprintf(fp,"\tlake_var.soil.asat = %f\n", tmpval);

   tmpval = prcp->lake_var.soil.baseflow;
   fprintf(fp,"\tlake_var.soil.baseflow = %f\n", tmpval);

   tmpval = prcp->lake_var.soil.inflow;
   fprintf(fp,"\tlake_var.soil.inflow = %f\n", tmpval);

   for ( i = 0; i < N_PET_TYPES; ++i )
   {
      tmpval = prcp->lake_var.soil.pot_evap[i];
      fprintf(fp,"\tlake_var.soil.pot_evap[%d] = %f\n", i, tmpval);
   }

   tmpval = prcp->lake_var.soil.vic411_runoff;
   fprintf(fp,"\tlake_var.soil.vic411_runoff = %f\n", tmpval);

   for ( lidx = 0; lidx < vic411_options.Nlayer; ++lidx )
   {
      tmpval = prcp->lake_var.soil.layer[lidx].Cs;
      fprintf(fp,"\tlake_var.soil.layer[%d].Cs = %f\n", lidx, tmpval);

      tmpval = prcp->lake_var.soil.layer[lidx].T;
      fprintf(fp,"\tlake_var.soil.layer[%d].T = %f\n", lidx, tmpval);

      tmpval = prcp->lake_var.soil.layer[lidx].evap;
      fprintf(fp,"\tlake_var.soil.layer[%d].evap = %f\n", lidx, tmpval);

      tmpval = prcp->lake_var.soil.layer[lidx].ice;
      fprintf(fp,"\tlake_var.soil.layer[%d].ice = %f\n", lidx, tmpval);

      tmpval = prcp->lake_var.soil.layer[lidx].min_liq;
      fprintf(fp,"\tlake_var.soil.layer[%d].min_liq = %f\n", lidx, tmpval);

      tmpval = prcp->lake_var.soil.layer[lidx].kappa;
      fprintf(fp,"\tlake_var.soil.layer[%d].kappa = %f\n", lidx, tmpval);

      tmpval = prcp->lake_var.soil.layer[lidx].moist;
      fprintf(fp,"\tlake_var.soil.layer[%d].moist = %f\n", lidx, tmpval);

      tmpval = prcp->lake_var.soil.layer[lidx].phi;
      fprintf(fp,"\tlake_var.soil.layer[%d].phi = %f\n", lidx, tmpval);
   }

   tmpval = prcp->lake_var.soil.rootmoist;
   fprintf(fp,"\tlake_var.soil.rootmoist = %f\n", tmpval);

   tmpval = prcp->lake_var.soil.wetness;
   fprintf(fp,"\tlake_var.soil.wetness = %f\n", tmpval);


   for ( veg = 0; veg <= Nveg; ++veg )
   {
      for ( band = 0; band < Nbands; ++band )
      {
         tmpval = prcp->snow[veg][band].albedo;
         fprintf(fp,"\tsnow[%d][%d].albedo = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].canopy_albedo;
         fprintf(fp,"\tsnow[%d][%d].canopy_albedo = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].coldcontent;
         fprintf(fp,"\tsnow[%d][%d].coldcontent = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].coverage;
         fprintf(fp,"\tsnow[%d][%d].coverage = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].density;
         fprintf(fp,"\tsnow[%d][%d].density = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].depth;
         fprintf(fp,"\tsnow[%d][%d].depth = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].last_snow;
         fprintf(fp,"\tsnow[%d][%d].last_snow = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].max_swq;
         fprintf(fp,"\tsnow[%d][%d].max_swq = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].MELTING;
         fprintf(fp,"\tsnow[%d][%d].MELTING = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].pack_temp;
         fprintf(fp,"\tsnow[%d][%d].pack_temp = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].pack_water;
         fprintf(fp,"\tsnow[%d][%d].pack_water = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].snow;
         fprintf(fp,"\tsnow[%d][%d].snow = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].snow_canopy;
         fprintf(fp,"\tsnow[%d][%d].snow_canopy = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].store_coverage;
         fprintf(fp,"\tsnow[%d][%d].store_coverage = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].store_snow;
         fprintf(fp,"\tsnow[%d][%d].store_snow = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].store_swq;
         fprintf(fp,"\tsnow[%d][%d].store_swq = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].surf_temp;
         fprintf(fp,"\tsnow[%d][%d].surf_temp = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].surf_temp_fbcount;
         fprintf(fp,"\tsnow[%d][%d].surf_temp_fbcount = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].surf_temp_fbflag;
         fprintf(fp,"\tsnow[%d][%d].surf_temp_fbflag = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].surf_water;
         fprintf(fp,"\tsnow[%d][%d].surf_water = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].swq;
         fprintf(fp,"\tsnow[%d][%d].swq = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].swq_slope;
         fprintf(fp,"\tsnow[%d][%d].swq_slope = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].tmp_int_storage;
         fprintf(fp,"\tsnow[%d][%d].tmp_int_storage = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].blowing_flux;
         fprintf(fp,"\tsnow[%d][%d].blowing_flux = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].canopy_vapor_flux;
         fprintf(fp,"\tsnow[%d][%d].canopy_vapor_flux = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].mass_error;
         fprintf(fp,"\tsnow[%d][%d].mass_error = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].melt;
         fprintf(fp,"\tsnow[%d][%d].melt = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].Qnet;
         fprintf(fp,"\tsnow[%d][%d].Qnet = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].surface_flux;
         fprintf(fp,"\tsnow[%d][%d].surface_flux = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].transport;
         fprintf(fp,"\tsnow[%d][%d].transport = %f\n", veg, band, tmpval);

         tmpval = prcp->snow[veg][band].vapor_flux;
         fprintf(fp,"\tsnow[%d][%d].vapor_flux = %f\n", veg, band, tmpval);
      }
   }

   for ( dist = 0; dist < Ndist; ++dist )
   {
      for ( veg = 0; veg <= Nveg; ++veg )
      {
         for ( i = 0; i < 2; ++i )
         {
            tmpval = prcp->veg_var[dist][veg][i].canopyevap;
            fprintf(fp,"\tveg_var[%d][%d][%d].canopyevap = %f\n", dist, veg, i, tmpval);

            tmpval = prcp->veg_var[dist][veg][i].throughfall;
            fprintf(fp,"\tveg_var[%d][%d][%d].throughfall = %f\n", dist, veg, i, tmpval);

            tmpval = prcp->veg_var[dist][veg][i].Wdew;
            fprintf(fp,"\tveg_var[%d][%d][%d].Wdew = %f\n", dist, veg, i, tmpval);
         }
      }
   }

   fclose(fp);
}
