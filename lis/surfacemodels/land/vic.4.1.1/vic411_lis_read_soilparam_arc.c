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
#include <vic411_vicNl.h>
#include <string.h>
//BOP
//
// !ROUTINE: vic411_lis_read_soilparam_arc
// \label{vic411_lis_read_soilparam_arc}
// 
// !REVISION HISTORY: 
//  Sept 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
vic411_soil_con_struct vic411_lis_read_soilparam_arc(FILE  *soilparam, 
				                           char  *soilparamdir, 
                                       float *lat,
                                       float *lng)
//
// !DESCRIPTION: 
//  This routine reads soil parameters for each grid cell from an ASCII
//  ARC/INFO output grid.
//
//  This routine was modified from the VIC 4.1.1 read\_soilparam\_arc.c file.
//
//  The arguments are:
//  \begin{description}
//  \item[soilparam] filename of the soil parameter file
//  \item[soilparamdir] directory containing the soil parameter file
//  \item[lat] latitude of the desired grid-cell
//  \item[lng] longitude of the desired grid-cell
//  \end{description}
//
//EOP
{
  extern vic411_option_struct vic411_options;
  extern vic411_global_param_struct vic411_global_param;
  extern vic411_veg_lib_struct *vic411_veg_lib;
#if LINK_DEBUG
  extern vic411_debug_struct debug;
#endif

  int             layer;
  int             i,j;
  int             Nbands,band;
  vic411_soil_con_struct temp; 
  double          Wcr_FRACT[MAX_LAYERS];
  double          Wpwp_FRACT[MAX_LAYERS];
  double          off_gmt;
  double          clay[MAX_LAYERS];
  double          sand[MAX_LAYERS];
  double          sum_depth;
  double          extra_depth;
  char            ErrStr[MAXSTRING];
  char            namestr[MAXSTRING];
  char            tmpstr[MAXSTRING];
  double          tmp_lat;
  double          tmp_lng;
  double          start_lat;
  double          right_lng;
  double          left_lng;
  double          delta;
  double          dist;
  double          tmp_bubble;
#if EXCESS_ICE
  double          init_ice_fract[MAX_LAYERS];
#endif

    tmp_bubble = 0;
  
    rewind(soilparam);
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
  
    temp.lat = *lat;
    temp.lng = *lng;
    temp.gridcel = (int)vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);

    fscanf(soilparam,"%s",tmpstr);

#if VERBOSE
    /* add print statements for grid cell number */
    fprintf(stderr,"\ncell: %d,  lat: %.4f, long: %.4f\n",temp.gridcel,temp.lat,temp.lng);
#endif

    /** Get Average Grid Cell Elevation **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.elevation = (float)vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
  
    /** Get Grid Cell Infiltration Parameter **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.b_infilt = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
  
    /** Get Maximum Baseflow Fraction **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.Ds = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
  
    /** Get Maximum Baseflow Velocity **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.Dsmax = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
  
    /** Get Maximum Soil Moisture Fraction **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.Ws = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
  
    /** Get Exponential **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.c = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);

    /** Get Average Soil Temperature **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.avg_temp = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
    if(vic411_options.FULL_ENERGY && (temp.avg_temp>100. || temp.avg_temp<-50)) {
      fprintf(stderr,"Need valid average soil temperature in degrees C to run");
      fprintf(stderr," Full Energy model, %f is not acceptable.\n",
	      temp.avg_temp);
      exit(0);
    }

    /** Get Soil Thermal Damping Depth **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.dp = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);

    /** Get Data Time Zone Offset from GMT **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    off_gmt = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);

    /** Get Critical Soil Moisture Fraction for each layer **/
    for(layer = 0; layer < vic411_options.Nlayer; layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,"/");
      strcat(namestr,tmpstr);
      Wcr_FRACT[layer] = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Wilting Point Soil Moisture Fraction for each layer **/
    for(layer = 0; layer < vic411_options.Nlayer; layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,"/");
      strcat(namestr,tmpstr);
      Wpwp_FRACT[layer] = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Bare Soil Roughness **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.rough = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
#if !OUTPUT_FORCE
    /* Overwrite default bare soil aerodynamic resistance parameters
       with the values taken from the soil parameter file */
    for (j=0; j<12; j++) {
      vic411_veg_lib[vic411_veg_lib[0].NVegLibTypes].roughness[j] = temp.rough;
      vic411_veg_lib[vic411_veg_lib[0].NVegLibTypes].displacement[j] = temp.rough*0.667/0.123;
    }
#endif // !OUTPUT_FORCE

    /** Get Snow Surface Roughness **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.snow_rough = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);

    /** Get Average Annual Precipitation **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.annual_prec = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);

    /** Get Layer Percent Sand **/
    for(layer=0;layer<vic411_options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,"/");
      strcat(namestr,tmpstr);
      sand[layer] = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Layer Percent Clay **/
    for(layer=0;layer<vic411_options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,"/");
      strcat(namestr,tmpstr);
      clay[layer] = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Layer Saturated Hydraulic Conductivity **/
    for(layer=0;layer<vic411_options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,"/");
      strcat(namestr,tmpstr);
      temp.Ksat[layer] = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Layer Soil Moisture Diffusion Parameter **/
    for(layer=0;layer<vic411_options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,"/");
      strcat(namestr,tmpstr);
      temp.phi_s[layer] = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Initial Layer Moisture **/
    for(layer=0;layer<vic411_options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,"/");
      strcat(namestr,tmpstr);
      temp.init_moist[layer] = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Layer Thickness **/
    for(layer=0;layer<vic411_options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,"/");
      strcat(namestr,tmpstr);
      temp.depth[layer] = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
#if !EXCESS_ICE
      temp.depth[layer] = (float)(int)(temp.depth[layer] * 1000 + 0.5) / 1000;
#endif
    }

    /** Get Layer Bulk Density **/
    for(layer=0;layer<vic411_options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,"/");
      strcat(namestr,tmpstr);
      temp.bulk_density[layer] = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Layer Porosity **/
    for(layer=0;layer<vic411_options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,"/");
      strcat(namestr,tmpstr);
      temp.porosity[layer] = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Activate Frozen Soil for Grid Cell **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.FS_ACTIVE = (char)vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);

    /** Minimum Snow Depth of Full Coverage **/
    fscanf(soilparam,"%s",tmpstr);
#if SPATIAL_SNOW
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.depth_full_snow_cover = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
#endif // SPATIAL_SNOW

    /** Slope of Frozen Soil Distribution **/
    fscanf(soilparam,"%s",tmpstr);
#if SPATIAL_FROST
    strcpy(namestr,soilparamdir);
    strcat(namestr,"/");
    strcat(namestr,tmpstr);
    temp.frost_slope = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
#endif // SPATIAL_FROST

    /** Layer Initial Volumetric Ice Fraction **/
    for(layer=0;layer<vic411_options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
#if EXCESS_ICE
      strcpy(namestr,soilparamdir);
      strcat(namestr,"/");
      strcat(namestr,tmpstr);
      init_ice_fract[layer] = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
#endif // EXCESS_ICE
    }

    if (vic411_options.COMPUTE_TREELINE && vic411_options.JULY_TAVG_SUPPLIED && (fscanf(soilparam,"%s",tmpstr)) != EOF) {
      /** Get Avg July Air Temperature **/
      strcpy(namestr,soilparamdir);
      strcat(namestr,"/");
      strcat(namestr,tmpstr);
      temp.avgJulyAirTemp = vic411_read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

#if !OUTPUT_FORCE
    /*******************************************
      Compute Soil Layer Properties
      *******************************************/
    sum_depth   = 0.;
    for(layer = 0; layer < vic411_options.Nlayer; layer++) {
      sum_depth += temp.depth[layer];
      temp.bulk_density[layer] *= 1000.;
      temp.soil_density[layer] = temp.bulk_density[layer] 
	/ (1.0 - temp.porosity[layer]);
      temp.quartz[layer] = sand[layer] / 100.;
#if !EXCESS_ICE
      temp.max_moist[layer] = temp.depth[layer] * temp.porosity[layer] * 1000.;
#endif
      temp.bubble[layer] = exp(5.3396738 + 0.1845038*clay[layer] 
			 - 2.48394546*temp.porosity[layer] 
			 - 0.00213853*pow(clay[layer],2.)
			 - 0.04356349*sand[layer]*temp.porosity[layer]
			 - 0.61745089*clay[layer]*temp.porosity[layer]
			 + 0.00143598*pow(sand[layer],2.)
			 * pow(temp.porosity[layer],2.)
			 - 0.00855375*pow(clay[layer],2.)
			 * pow(temp.porosity[layer],2.)
			 - 0.00001282*pow(sand[layer],2.)*clay[layer]
			 + 0.00895359*pow(clay[layer],2.)*temp.porosity[layer]
			 - 0.00072472*pow(sand[layer],2.)*temp.porosity[layer]
			 + 0.00000540*pow(clay[layer],2.)*sand[layer]
			 + 0.50028060*pow(temp.porosity[layer],2.)*clay[layer]);
      temp.expt[layer] = exp(-0.7842831 + 0.0177544*sand[layer] 
			     - 1.062498*temp.porosity[layer] 
			     - 0.00005304*pow(sand[layer],2.)
			     - 0.00273493*pow(clay[layer],2.)
			     + 1.11134946*pow(temp.porosity[layer],2.)
			     - 0.03088295*sand[layer]*temp.porosity[layer]
			     + 0.00026587*pow(sand[layer],2.)
			     * pow(temp.porosity[layer],2.)
			     - 0.00610522*pow(clay[layer],2.)
			     * pow(temp.porosity[layer],2.)
			     - 0.00000235*pow(sand[layer],2.)*clay[layer]
			     + 0.00798746*pow(clay[layer],2.)*temp.porosity[layer]
			     - 0.00674491*pow(temp.porosity[layer],2.)*clay[layer]);
      temp.expt[layer] = 2. / temp.expt[layer] + 3.;
      temp.resid_moist[layer] = - 0.0182482 + 0.00087269 * sand[layer]
	                      + 0.00513488 * clay[layer] 
	                      + 0.02939286 * temp.porosity[layer] 
                              - 0.00015395 * pow(clay[layer],2.) 
                              - 0.00108270 * sand[layer] * temp.porosity[layer] 
                              - 0.00018233 * pow(clay[layer],2.) 
                                           * pow(temp.porosity[layer],2.) 
                              + 0.00030703 * pow(clay[layer],2.0) 
                                           * temp.porosity[layer] 
                              - 0.00235840 * pow(temp.porosity[layer],2.) 
                                           * clay[layer];

      /** Check for valid values of generated parameters **/
      if(temp.bubble[layer]<1.36) {
	fprintf(stderr,"WARNING: estimated bubbling pressure too low (%f), resetting to minimum value (%f).\n",temp.bubble[layer],1.36);
	temp.bubble[layer] = 1.36;
      }
      if(temp.bubble[layer]>187.2) {
	fprintf(stderr,"WARNING: estimated bubbling pressure too high (%f), resetting to maximum value (%f).\n",temp.bubble[layer],187.2);
	temp.bubble[layer] = 187.2;
      }
      if(temp.expt[layer] < 2. / 1.090 + 3.) {
	fprintf(stderr,"WARNING: estimated exponential (expt) too low (%f), resetting to minimum value (%f).\n", temp.expt[layer], 2. / 1.090 + 3.);
	temp.expt[layer] = 2. / 1.090 + 3.;
      }
      if(temp.expt[layer] > 2. / 0.037 + 3.) {
	fprintf(stderr,"WARNING: estimated exponential (expt) too high (%f), resetting to maximum value (%f).\n",temp.expt[layer], 2. / 0.037 + 3.);
	temp.expt[layer] = 2. / 0.037 + 3.;
      }
      if(temp.resid_moist[layer] < -0.038) {
	fprintf(stderr,"WARNING: estimated residual soil moisture too low (%f), resetting to minimum value (%f).\n",temp.resid_moist[layer],-0.038);
	temp.resid_moist[layer] = -0.038;
      }
      if(temp.resid_moist[layer] > 0.205) {
	fprintf(stderr,"WARNING: estimated residual soil moisture too high (%f), resetting to maximum value (%f).\n",temp.resid_moist[layer],0.205);
	temp.resid_moist[layer] = 0.205;
      }
      tmp_bubble += temp.bubble[layer];
    }
    for(layer=0;layer<vic411_options.Nlayer;layer++) temp.bubble[layer] = tmp_bubble/3.;

#if !EXCESS_ICE
    /*******************************************
      Validate Initial Soil Layer Moisture Content for !EXCESS_ICE option
    *******************************************/
    if (!vic411_options.INIT_STATE) {
      for(layer = 0; layer < vic411_options.Nlayer; layer++) {
	if(temp.init_moist[layer] > temp.max_moist[layer]) {
	  fprintf(stderr,"Initial soil moisture (%f mm) is greater than the maximum moisture (%f mm) for layer %i.\n\tResetting soil moisture to maximum.\n",
		temp.init_moist[layer], temp.max_moist[layer], layer);
	  temp.init_moist[layer] = temp.max_moist[layer];
	}
	if(temp.init_moist[layer] < temp.resid_moist[layer] * temp.depth[layer] * 1000.) {
	  fprintf(stderr,"Initial soil moisture (%f mm) is less than calculated residual moisture (%f mm) for layer %i.\n\tResetting soil moisture to residual moisture.\n",
		temp.init_moist[layer], temp.resid_moist[layer] * temp.depth[layer] * 1000., layer);
	  temp.init_moist[layer] = temp.resid_moist[layer] * temp.depth[layer] * 1000.;
	}
      }
    }
#endif

#if EXCESS_ICE
    /*******************************************
      Compute Soil Layer Properties for EXCESS_ICE option
    *******************************************/
    extra_depth=0;
    for(layer = 0; layer < vic411_options.Nlayer; layer++) {
      temp.min_depth[layer]=temp.depth[layer];
      if(init_ice_fract[layer]>MAX_ICE_INIT){ // validate amount based on physical constraints
	fprintf(stderr,"Initial ice fraction (%f) is greater than maximum ice content for layer %d.\n\tResetting to maximum of %f\n",init_ice_fract[layer],layer,MAX_ICE_INIT);
	init_ice_fract[layer]=MAX_ICE_INIT;
      }
      if(init_ice_fract[layer]>=temp.porosity[layer]){//excess ground ice present
	fprintf(stderr,"Excess ground ice present in layer %d:\n",layer+1);
	fprintf(stderr,"\t\tSubsidence will occur when the average soil layer\n\t\t  temperature exceeds %.2f degrees Celsius.\n",powf((1.-ICE_AT_SUBSIDENCE),(3.-temp.expt[layer])/2.)*273.16*9.81*temp.bubble[layer]/(-Lf*100.));	
	temp.effective_porosity[layer]=init_ice_fract[layer];
	fprintf(stderr,"\t\tEffective porosity increased from %.2f to %.2f.\n",temp.porosity[layer],temp.effective_porosity[layer]);
	temp.depth[layer] = temp.min_depth[layer]*(1.0 - temp.porosity[layer])/(1.0 - temp.effective_porosity[layer]); //adjust soil layer depth
	extra_depth += temp.depth[layer]-temp.min_depth[layer]; //net increase in depth due to excess ice
	fprintf(stderr,"\t\tDepth of soil layer adjusted for excess ground ice: from %.2f m to %.2f m.\n",temp.min_depth[layer],temp.depth[layer]);
	fprintf(stderr,"\t\tBulk density adjusted for excess ground ice: from %.2f kg/m^3 to %.2f kg/m^3.\n",temp.bulk_density[layer],(1.0-temp.effective_porosity[layer])*temp.soil_density[layer]);
	temp.bulk_density[layer] = (1.0-temp.effective_porosity[layer])*temp.soil_density[layer]; //adjust bulk density
      }
      else //excess ground ice not present
	temp.effective_porosity[layer]=temp.porosity[layer];
    }
    if(extra_depth>0) {
      fprintf(stderr,"Damping depth adjusted for excess ground ice: from %.2f m to %.2f m.\n",temp.dp,temp.dp+extra_depth);
      temp.dp = temp.dp + extra_depth;  //adjust damping depth
    }
    
    /* final soil layer thicknesses for EXCESS_ICE option */
    for(layer = 0; layer < vic411_options.Nlayer; layer++) 
      temp.depth[layer] = (float)(int)(temp.depth[layer] * 1000 + 0.5) / 1000;
    
    /* Calculate and Validate Maximum Initial Soil Layer Moisture Content for EXCESS_ICE option */
    for(layer = 0; layer < vic411_options.Nlayer; layer++) {
      temp.max_moist[layer] = temp.depth[layer] * temp.effective_porosity[layer] * 1000.;
      if(temp.effective_porosity[layer]>temp.porosity[layer])//excess ground ice present
	temp.init_moist[layer] = temp.max_moist[layer]; 
      else {//excess ground ice not present
	if(temp.depth[layer] * init_ice_fract[layer] * 1000. > temp.init_moist[layer])
	  temp.init_moist[layer] = temp.depth[layer] * init_ice_fract[layer] * 1000.;
      }
    }
    for(layer = 0; layer < vic411_options.Nlayer; layer++) {
      if(temp.init_moist[layer] > temp.max_moist[layer]) {
	fprintf(stderr,"Initial soil moisture (%f mm) is greater than the maximum moisture (%f mm) for layer %d.\n\tResetting soil moisture to maximum.\n",
		temp.init_moist[layer], temp.max_moist[layer], layer);
	temp.init_moist[layer] = temp.max_moist[layer];
      }
      if(temp.init_moist[layer] < temp.resid_moist[layer] * temp.depth[layer] * 1000.) { 
	fprintf(stderr,"Initial soil moisture (%f mm) is less than calculated residual moisture (%f mm) for layer %d.\n\tResetting soil moisture to residual moisture.\n",
		temp.init_moist[layer], temp.resid_moist[layer] * temp.depth[layer] * 1000., layer);
	temp.init_moist[layer] = temp.resid_moist[layer] * temp.depth[layer] * 1000.;
      }
    }
    /* print final values for each layer */
    //for(layer = 0; layer < vic411_options.Nlayer; layer++)
    //fprintf(stderr,"final soil values: %d %.2f %.2f %.2f %.2f %.2f\n",layer,temp.effective_porosity[layer],temp.depth[layer],temp.bulk_density[layer],temp.max_moist[layer],temp.init_moist[layer]);
#endif // EXCESS_ICE   

    /*******************************************
      Validate Soil Layer Thicknesses                    
    *******************************************/
    for(layer=0;layer<vic411_options.Nlayer;layer++) {
      if(temp.depth[layer] < MINSOILDEPTH) {
	fprintf(stderr,"Model will not function with layer depth %f < %f m.\n",
		temp.depth[layer],MINSOILDEPTH);
	exit(0);
      }
    }
    if(temp.depth[0] > temp.depth[1]) {
      sprintf(ErrStr,"ERROR: Model will not function with layer %i depth (%f m) < layer %i depth (%f m).\n",
	      0,temp.depth[0],1,temp.depth[1]);
      vic411_nrerror(ErrStr);
    }
#if EXCESS_ICE
    for(layer = 0; layer < vic411_options.Nlayer; layer++) {
      if(temp.min_depth[layer] < MINSOILDEPTH) {
	sprintf(ErrStr,"ERROR: Model will not function with layer %d depth %f < %f m.\n",
		layer,temp.min_depth[layer],MINSOILDEPTH);
	vic411_nrerror(ErrStr);
      }
    }
    if(temp.min_depth[0] > temp.min_depth[1]) {
      sprintf(ErrStr,"ERROR: Model will not function with layer %d depth (%f m) < layer %d depth (%f m).\n",
	      0,temp.min_depth[0],1,temp.min_depth[1]);
      vic411_nrerror(ErrStr);
    }
#endif /* EXCESS_ICE */

    /**********************************************
      Compute Maximum Infiltration for Upper Layers
      **********************************************/
    if(vic411_options.Nlayer==2)
      temp.max_infil = (1.0+temp.b_infilt)*temp.max_moist[0];
    else
      temp.max_infil = (1.0+temp.b_infilt)*(temp.max_moist[0]
					    + temp.max_moist[1]);

    /****************************************************************
      Compute Soil Layer Critical and Wilting Point Moisture Contents
      ****************************************************************/
    for(layer=0;layer<vic411_options.Nlayer;layer++) {
      temp.Wcr[layer]  = Wcr_FRACT[layer] * temp.max_moist[layer];
      temp.Wpwp[layer] = Wpwp_FRACT[layer] * temp.max_moist[layer];
#if EXCESS_ICE
      temp.Wcr_FRACT[layer]  = Wcr_FRACT[layer];
      temp.Wpwp_FRACT[layer] = Wpwp_FRACT[layer]; 
#endif
      if(temp.Wpwp[layer] > temp.Wcr[layer]) {
	sprintf(ErrStr,"Calculated wilting point moisture (%f mm) is greater than calculated critical point moisture (%f mm) for layer %i.\n\tIn the soil parameter file, Wpwp_FRACT MUST be <= Wcr_FRACT.\n",
		temp.Wpwp[layer], temp.Wcr[layer], layer);
	vic411_nrerror(ErrStr);
      }
      if(temp.Wpwp[layer] < temp.resid_moist[layer] * temp.depth[layer] * 1000.) {
	sprintf(ErrStr,"Calculated wilting point moisture (%f mm) is less than calculated residual moisture (%f mm) for layer %i.\n\tIn the soil parameter file, Wpwp_FRACT MUST be >= resid_moist / (1.0 - bulk_density/soil_density).\n",
		temp.Wpwp[layer], temp.resid_moist[layer] * temp.depth[layer] * 1000., layer);
	vic411_nrerror(ErrStr);
      }
    }

    /*************************************************
    if BASEFLOW == NIJSSEN2001 then convert the baseflow
    parameters d1, d2, d3, d4 to Ds, Dsmax, Ws, and c.  JA
    *************************************************/
#if EXCESS_ICE
    temp.Dsmax_orig = temp.Dsmax;
    temp.Ds_orig = temp.Ds;
    temp.Ws_orig = temp.Ws;
#endif
    if(vic411_options.BASEFLOW == NIJSSEN2001) {
      layer = vic411_options.Nlayer-1;
      temp.Dsmax = temp.Dsmax * 
	pow((double)(1./(temp.max_moist[layer]-temp.Ws)), -temp.c) +
	temp.Ds * temp.max_moist[layer];
      temp.Ds = temp.Ds * temp.Ws / temp.Dsmax;
      temp.Ws = temp.Ws/temp.max_moist[layer];
    }

    /*******************************************************************
      Calculate grid cell area.
    ******************************************************************/

    if (vic411_options.EQUAL_AREA) {

      temp.cell_area = vic411_global_param.resolution * 1000. * 1000.; /* Grid cell area in m^2. */

    }
    else {

      tmp_lat = fabs(*lat);
      tmp_lng = fabs(*lng);

      start_lat = tmp_lat - vic411_global_param.resolution / 2;
      right_lng = tmp_lng + vic411_global_param.resolution / 2;
      left_lng  = tmp_lng - vic411_global_param.resolution / 2;

      delta = vic411_get_dist(tmp_lat,tmp_lng,tmp_lat+vic411_global_param.resolution/10.,tmp_lng);

      dist = 0.;

      for ( i = 0; i < 10; i++ ) {
        dist += vic411_get_dist(start_lat,left_lng,start_lat,right_lng) * delta;
        start_lat += vic411_global_param.resolution/10;
      }

      temp.cell_area = dist * 1000. * 1000.; /* Grid cell area in m^2. */

    }

#endif /* !OUTPUT_FORCE */
    
    /*************************************************
      Determine Central Longitude of Current Time Zone 
    *************************************************/
    temp.time_zone_lng = off_gmt * 360./24.;
    
    /*************************************************
      Allocate and Initialize Snow Band Parameters
    *************************************************/
    Nbands         = vic411_options.SNOW_BAND;
    temp.AreaFract     = (double *)calloc(Nbands,sizeof(double));
    temp.BandElev      = (float *)calloc(Nbands,sizeof(float));
    temp.Tfactor       = (double *)calloc(Nbands,sizeof(double));
    temp.Pfactor       = (double *)calloc(Nbands,sizeof(double));
    temp.AboveTreeLine = (char *)calloc(Nbands,sizeof(char));

    if (temp.Tfactor == NULL || temp.Pfactor == NULL || temp.AreaFract == NULL) 
      vic411_nrerror("Memory allocation failure in vic411_read_snowband");

    if ( Nbands <= 0 ) {
      sprintf(ErrStr,"Number of snow bands must be > 0 (%d)",Nbands);
      vic411_nrerror(ErrStr);
    }

    /** Set default values for factors to use unmodified forcing data **/
    for (band = 0; band < Nbands; band++) {
      temp.AreaFract[band] = 0.;
      temp.BandElev[band]  = temp.elevation;
      temp.Tfactor[band]   = 0.;
      temp.Pfactor[band]   = 1.;
    }
    temp.AreaFract[0] = 1.;

  return temp;
} 
