//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.4
//
// Copyright (c) 2022 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <stdio.h>
#include <stdlib.h>
#include "vic411_vicNl.h" 
#include <string.h>

#define BARETYPE ( 12 ) // UMD -- HARDCODED!

extern int vic411_read_vegparam_lis(FILE *, int, int, int, vic411_veg_con_struct *);

int vic411_real_2_vic(int real_type, int Nveg_type)
{
   extern vic411_veg_lib_struct *vic411_veg_lib;
   int vic_type;
   int k, vic411_flag;

   vic411_flag = 0;
   for(k=0; k<Nveg_type; k++)
    { 
        if(real_type == vic411_veg_lib[k].veg_class)
        {
            vic_type = k;
            vic411_flag = 1;
            break;
        }
    }

    if(vic411_flag==0)
    {
        //printf("LIS-VIC error (vic411_set_vegparam.c - vic411_real_2_vic): undefined land cover type in vic411_veg_lib!\n");
        //printf("Check the vegetation libary of VIC.\n");
        //exit(1);
        vic_type = -1;        
    }

    return(vic_type);
}

int vic411_is_understory_type(int real_type, int Nveg_type)
{
    extern vic411_veg_lib_struct *vic411_veg_lib;
    int vic_type;
    if(real_type == BARETYPE)
    {
        return(0);
    }
    vic_type=vic411_real_2_vic(real_type, Nveg_type);
    if(vic_type==-1)
    {
        return(0);
    }
    if( vic411_veg_lib[vic_type].overstory==1)
    {
        return(0);
    }
    else
    {
        return(1);
    }
}

int vic411_is_overstory_type(int real_type, int Nveg_type)
{
    extern vic411_veg_lib_struct *vic411_veg_lib;
    int vic_type;
    if(real_type == BARETYPE)
    {
        return(0);
    } 
    vic_type=vic411_real_2_vic(real_type, Nveg_type);
    if(vic_type == -1)
    {
        return(0);
    }
    if(vic411_veg_lib[vic_type].overstory==1)
    {
        return(1);
    }
    else
    {
        return(0);
    }
}

int vic411_number_of_understory_types(int Nveg_type, int Nlc_type, float *veg_fracs)
{
    int num_understory = 0;
    int n, real_type;

    // loop all land cover types, n=1:Nlc_type
    for(n=0; n<Nlc_type; n++)
    {
        if(veg_fracs[n]>0)
        {
            real_type = n + 1; 
            if(vic411_is_understory_type(real_type, Nveg_type)==1)
            {
               num_understory++;
            }
        }
    }
    return(num_understory);
}


int vic411_number_of_overstory_types(int Nveg_type, int Nlc_type, float *veg_fracs)
{
    int num_overstory = 0;
    int n, real_type;
    
    // loop all land cover types, n=1:Nlc_type
    for(n=0; n<Nlc_type; n++)
    {
        if(veg_fracs[n]>0)
        {
            real_type = n + 1; 
            if(vic411_is_overstory_type(real_type, Nveg_type)==1)
            {
               num_overstory++;
            }
        }
    }
    return(num_overstory);
}

double vic411_normalize_cv(int real_type000, int Nveg_type, int Nlc_type, float *fracs)
{
    double understory_sum_cv = 0.0;
    int n, real_type;
    
    // loop all lanc cover types n=1:Nlc_type
    for(n=0; n<Nlc_type; n++)
    {
        real_type = n + 1;
        if(vic411_is_understory_type(real_type, Nveg_type))
        {   
            understory_sum_cv += fracs[n];
        }
    }
    // fracs indexes from 0
    return(fracs[real_type000-1]/understory_sum_cv);
}


int vic411_has_above_treeline_band(int tile_idx)
{
    extern vic411_soil_con_struct *vic411_lis_soil_con;
    extern vic411_option_struct   vic411_options;
    
    int band, vic411_flag;
    
    vic411_flag = 0; 
    for(band=0; band <vic411_options.SNOW_BAND; band++)
    {
        if(vic411_lis_soil_con[tile_idx].AboveTreeLine[band] == TRUE)
        {
            vic411_flag = 1;
            break;        
        }
    }
    
    return(vic411_flag);
}

//BOP
//
// !ROUTINE: vic411_set_vegparam
// \label{set_vegparam411}
// 
// !REVISION HISTORY: 
//  22 Feb 2012  James Geiger; Initial adaptation from VIC 4.1.1 (vic411_read_vegparam.c).
//  01 Jun 2012  Shugong Wang; Add support to COMPUTE_TREELINE option 
// !INTERFACE:
vic411_veg_con_struct *vic411_set_vegparam(int tile_idx, 
                             int vegclass, 
                             int Nveg_type,
                             int gridcel, 
                             FILE *fp_vegparam, 
                             float *veg_fracs, 
                             int Nlc_type)
// !DESCRIPTION: 
//  This routine initializes the vegetation parameters veg\_con data structure based on the given
//  vegclass value.
//
//  The arguments are:
//  \begin{description}
//  \item[tile\_idx] index of current tile 
//  \item[vegclass] vegetation classification, as determined by LIS, of the current tile,
//                  considered intent(in)
//  \item[Nveg\_type]    number of vegetation types in the veg\_lib parameter file, considered intent(in)
//  \item[gridcel]       grid id in VIC soil file
//  \item[fp\_veg\_param]  file pointer to vegetation parameter file
//  \item[veg\_fracs]     fractions of all vegetation classes 
//  \item[Nlc\_type]      number of vegetation classes in landcover dataset
//  \end{description}
//
//EOP
{

   extern vic411_veg_lib_struct *vic411_veg_lib;
   extern vic411_option_struct   vic411_options;
#if LINK_DEBUG
   extern vic411_debug_struct    debug;
#endif

   vic411_veg_con_struct *temp;

   int vegetat_type_num;
   int MaxVeg;
   int tempclass;
   float depth_sum;
   float sum;
   int i, j, n, n0, k;
   char str[500];
   char ErrStr[MAXSTRING];
   int AllOverstory;
   int vic_type;
   float overstory_sum;
   int n_under; 
   int real_type;
   int vic411_flag;
   int n_root_zone;
   double tmp1; 
   /*default root zone configuration for grass*/
   float grass_root_zone_depth[3] = {0.10, 1.00,  0.50}; 
   float grass_root_zone_fract[3] = {0.10, 0.70,  0.20};

   // by default, set the number of vegetation in a tile to be 1
   vegetat_type_num = 1;
   overstory_sum    = 0.0;
   
   /* Code added by Shugong Wang on 06/01/2012 to support COMPUTE_TREELINE option. 
      If the option is turned on, there should be multiple elevation/snow bands in the grid */
    // (1) If the vegetation class of current tile is an overstory vegetation class, then we need
    //     to consider the COMPUTE_TREELINE option
   if(vic411_options.SNOW_BAND>1                       &&  // there are multiple snow bands
      vic411_options.COMPUTE_TREELINE > 0              &&  // COMPUTE_TREELINE optin is turned on
      vegclass != BARETYPE                      &&  // not bare soil - BARETYPE may be replaced with LIS option later
      vic411_is_overstory_type(vegclass, Nveg_type)==1 &&  // current tile has overstory vegetation
      vic411_has_above_treeline_band(tile_idx))          // current grid has above-treeline band(s)
   {
        // (2) Determine whether there are understory vegetation class in the grid or not
        if(vic411_number_of_understory_types(Nveg_type, Nlc_type, veg_fracs)>=1)
        {
            AllOverstory = 0;
        }
        else
        {
            AllOverstory = 1;
        }
        
        // (3) If all vegetation classes in the grid are overstory
        //     replace the vegetation class of current tile with prescribed vegetation class or bare soil.
        //     In this case, there is only one vegetation class in the tile. Otherwise, put all understory vegetation 
        //     classes to the tile 
        if(AllOverstory==1)
        {
            // there are no understory classes and a prescribed understory type is specified in VIC global control file
            if(vic411_options.AboveTreelineVeg >=0)
            {
                vegetat_type_num = 2;   
            }
            // no prescribed understory vegation class type and will use bare soil to fill the space of overstory vegation class
            // there is no need to vegation class
            else
            {
                vegetat_type_num = 1;
            }
        }
        else // there are understory vegetation classes in the tile, put all understory vegetation classes in the tile.  
        {
            n_under = vic411_number_of_understory_types(Nveg_type, Nlc_type, veg_fracs);         
            vegetat_type_num = 1 + n_under;
        }
   }
      
   // Make sure to allocate extra memory for bare soil tile
   // and optionally an above-treeline veg tile
   MaxVeg = vegetat_type_num+1;
   
   /** Allocate memory for vegetation grid cell parameters **/
   temp = (vic411_veg_con_struct*) calloc( MaxVeg, sizeof(vic411_veg_con_struct));
   temp[0].Cv_sum = 1.0;            // by default, there is only one vegetation class when using LIS-based tiling 
   temp[0].Cv     = 1.0;            // So we can set Cv_sum = 1.0 and CV=1.0 for temp[0]
   temp[0].veg_class = vegclass;    // At this moment, veg_class is still the real vegetation class, VIC will
                                    // convert it into the index of corresponding vegetation class in vic411_veg_lib
                                    // in vic411_read_vegparam.c
   /* vic411_compute_treeline option, continue */
   // (4) adjust fractions  if vegetation class of current tile is overstory type and entire grid covered by 
   //     overstory types 
   if(vic411_options.SNOW_BAND>1                       &&  // there are multiple snow bands
      vic411_options.COMPUTE_TREELINE > 0              &&  // COMPUTE_TREELINE optin is turned on
      vegclass != BARETYPE                      &&  // not bare soil - BARETYPE may be replaced with LIS option later
      vic411_is_overstory_type(vegclass, Nveg_type)==1 &&  // current tile has overstory vegetation
      vic411_has_above_treeline_band(tile_idx))            // current grid has above-treeline band(s)
   {
      
      // substract 0.001 from the fraction of current tile if its vegeation class is an overstory type
      temp[0].Cv = 0.999;

      // if there are understory classes in the grid, in other words, not all vegetation classes are overstory types.
      if(AllOverstory==0)
      {
        // (4.1) determine the number of understory types in the grid
        n_under = vic411_number_of_understory_types(Nveg_type, Nlc_type, veg_fracs);
        n0=0;
        tmp1 = temp[0].Cv;
        for(n=1; n<vegetat_type_num; n++)
        {
            //printf("n = %d total_cv = %lf\n", n, tmp1);
            for(k=n0; k<Nveg_type; k++)
            {
                real_type = k + 1;
                // if current type is an understory type, add to vegcon
                if(veg_fracs[k]>0 && vic411_is_understory_type(real_type, Nveg_type)==1)
                {
                    //printf("vegetat_type_num = %d veg_class = %d\n",vegetat_type_num, real_type);
                    temp[n].veg_class = real_type;
                    temp[n].Cv        = 0.001 * vic411_normalize_cv(real_type, Nveg_type, Nlc_type, veg_fracs);
                    tmp1 += temp[n].Cv;
                    n0 = k + 1;
                    break;
                }
            }
        }
      }
      // there is no understory classes in the grid
      else
      {
        // if an understory type is prescribed in global control file
        if(vic411_options.AboveTreelineVeg >=0) 
        {
            temp[1].veg_class = vic411_options.AboveTreelineVeg;
            temp[1].Cv        = 0.001;
        }
        else // leave a very small portion to bare soil    
        {
            temp[0].Cv_sum = 0.999;
        }
      }
   }
   
   if ( vegclass != BARETYPE )
   {
      for ( i = 0; i < vegetat_type_num; i++ )
      {
         temp[i].zone_depth = calloc(vic411_options.ROOT_ZONES,sizeof(float));
         temp[i].zone_fract = calloc(vic411_options.ROOT_ZONES,sizeof(float));
         temp[i].vegetat_type_num = vegetat_type_num;
    
         temp[i].LAKE = 0;

         depth_sum = 0;

         // added by Shugong Wang to read in root zone information and update vic411_lis_veg_lib
         vic411_flag = vic411_read_vegparam_lis(fp_vegparam, gridcel , Nveg_type, vegclass, &(temp[i]));
         if(vic411_flag==0) // cannot find veg parameters
         {
            // Since root zones are not defined they are copied from the last
            // vegetation type.
            if(vic411_options.ROOT_ZONES>3)
            {
                n_root_zone = 3;
            }
            else
            {
                n_root_zone = vic411_options.ROOT_ZONES;
            }
            for ( j = 0; j < n_root_zone; j++ )
            {
              temp[i].zone_depth[j] = grass_root_zone_depth[j];
              temp[i].zone_fract[j] = grass_root_zone_fract[j];
            }
         }

         for ( j = 0; j < vic411_options.ROOT_ZONES; ++j)
         {
            depth_sum += temp[i].zone_depth[j];
         }

         sum = 0.;
         for ( j = 0; j < vic411_options.ROOT_ZONES; ++j)
         {
            sum += temp[i].zone_fract[j];
         }

         if ( depth_sum <= 0)
         {
            sprintf(str,"Root zone depths must sum to a value greater than 0.");
            vic411_nrerror(str);
         }

         if ( sum != 1.)
         {
            fprintf(stderr,"WARNING: Root zone fractions sum to more than 1 ( = %f), normalizing fractions.  If the sum is large, check that your vegetation parameter file is in the form - <zone 1 depth> <zone 1 fract> <zone 2 depth> <zone 2 fract> ...\n", sum);
            for ( j = 0;j < vic411_options.ROOT_ZONES; j++)
            {
               temp[i].zone_fract[j] /= sum;
            }
         }

         // check missing
         tempclass = MISSING;
         for ( j = 0; j < Nveg_type; j++ )
         {
            if(temp[i].veg_class == vic411_veg_lib[j].veg_class)
            {
               tempclass = j;
            }
         }

         if ( tempclass == MISSING )
         {
            sprintf(ErrStr,"The vegetation class id %i in vegetation tile %i from cell %i is not defined in the vegetation library file.", temp[i].veg_class, i, gridcel);
            vic411_nrerror(ErrStr);
         }
         else // convert real type into vic type
         {
            temp[i].veg_class = tempclass;
         }

         //temp[0].Cv_sum += temp[i].Cv;

      }
   }
   else
   {
      // For baresoil set the 0th index to some sane values.
      temp[0].zone_depth = calloc(vic411_options.ROOT_ZONES,sizeof(float));
      temp[0].zone_fract = calloc(vic411_options.ROOT_ZONES,sizeof(float));
      temp[0].vegetat_type_num = vegetat_type_num;
      temp[0].LAKE = 0;
      temp[0].veg_class = 1;
      temp[0].Cv = 0.0;
      // added by Shugong Wang, is it necessary? 
      temp[0].Cv_sum = 0.0;
       
      for ( j = 0; j < vic411_options.ROOT_ZONES; ++j)
      {
         temp[0].zone_depth[j] = 1.0 / vic411_options.ROOT_ZONES;
         temp[0].zone_fract[j] = 1.0 / vic411_options.ROOT_ZONES;
      }
   }


   // Bare soil tile
   if ( vegclass == BARETYPE )
   {
      j = vegetat_type_num;
      temp[j].veg_class = Nveg_type; // Create a veg_class ID for bare soil, which is not mentioned in the veg library
      temp[j].Cv = 1.0 - temp[0].Cv_sum;
      // Don't allocate any root-zone-related arrays
      if(vic411_options.BLOWING)
      {
         if (vegetat_type_num > 0) 
         {
            temp[j].sigma_slope = temp[0].sigma_slope;
            temp[j].lag_one = temp[0].lag_one;
            temp[j].fetch = temp[0].fetch;
         }
         else
         {
            temp[j].sigma_slope = 0.005;
            temp[j].lag_one = 0.95;
            temp[j].fetch = 2000;
         }
      }
   }

  return temp;
} 
