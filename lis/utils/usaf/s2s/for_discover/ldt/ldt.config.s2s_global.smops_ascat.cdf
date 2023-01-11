# Overall driver options
# - Based on USAF-Global SM CDF LDT options
# - For a ~25KM Global domain
# ___________________________________________________________

LDT running mode:                       "DA preprocessing"
Processed LSM parameter filename:       ./LDT_Params/lis_input.s2s_global.noahmp401_hymap2.25km.nc
LIS number of nests:                    1
Number of surface model types:          2
Surface model types:                    "LSM"  "Openwater"
Land surface model:                     "Noah-MP.4.0.1"
Routing model:                          none
Lake model:                             none
Water fraction cutoff value:            0.5
Incorporate crop information:           .false.

Number of met forcing sources:          0
Met forcing sources:                    none

LDT diagnostic file:                    ldtlog.s2s_global.ascat_cdf_201401-202201
Mask-parameter fill diagnostic file:    MaskParamFill.s2s_global.ascat_cdf_2014-202201.log
LDT output directory:                   OUTPUT
Undefined value:                        -9999.0

#Global 25KM domain:
Map projection of the LIS domain:       latlon
Run domain lower left lat:             -89.875
Run domain upper right lat:             89.875
Run domain lower left lon:            -179.875
Run domain upper right lon:            179.875
Run domain resolution (dx):              0.25
Run domain resolution (dy):              0.25

# ----------------------
# CDF generation options:

DA preprocessing method:                    "CDF generation"
DA observation source:                      "SMOPS soil moisture"
Name of the preprocessed DA file:           ./smops_ascat_cdf_100obs_201401-202201

Apply anomaly correction to obs:            0
Temporal resolution of CDFs:               "monthly"
Number of bins to use in the CDF:           100
Observation count threshold:                100
Temporal averaging interval:                1da
Apply external mask:                        0
External mask directory:                    none

#SMOPS soil moisture observation directory:  ./RS_DATA/SMOPS_v3
SMOPS soil moisture observation directory:  ../MET_FORCING/SMOPS
SMOPS soil moisture use ASCAT data:         1
SMOPS soil moisture use WindSat data:       0
SMOPS soil moisture use SMOS data:          0
SMOPS soil moisture use AMSR2 data:         0
SMOPS soil moisture use SMAP data:          0
SMOPS search radius for openwater proximity detection: 3
#Search radius to remove island points:      1

Starting year:                                    2014
Starting month:                                      1
Starting day:                                        1 
Starting hour:                                       0
Starting minute:                                     0
Starting second:                                     0
Ending year:                                      2022
Ending month:                                        1
Ending day:                                          1
Ending hour:                                         0
Ending minute:                                       0
Ending second:                                       0
LIS output timestep:                               1da

Number of ensembles per tile:           1

#The following options are used for subgrid tiling based on vegetation
Maximum number of surface type tiles per grid:    1
Minimum cutoff percentage (surface type tiles):   0.05
Maximum number of soil texture tiles per grid:    1
Minimum cutoff percentage (soil texture tiles):   0.05
Maximum number of soil fraction tiles per grid:   1
Minimum cutoff percentage (soil fraction tiles):  0.05
Maximum number of elevation bands per grid:       1
Minimum cutoff percentage (elevation bands):      0.05
Maximum number of slope bands per grid:           1
Minimum cutoff percentage (slope bands):          0.05
Maximum number of aspect bands per grid:          1
Minimum cutoff percentage (aspect bands):         0.05

# ----------------------

#Landcover parameter inputs
Landcover data source:                  "MODIS_Native"
Landcover classification:               "IGBPNCEP"
Landcover file:                         ./LDT_Params/input/LS_PARAMETERS/noah_2dparms/igbp.bin
Landcover spatial transform:            tile
Landcover map projection:               latlon
Landcover fill option:                  neighbor
Landcover fill radius:                  5
Landcover fill value:                   10

#Landmask parameter inputs
Create or readin landmask:              "readin"
Landmask data source:                   "UKMO_CAP_Netcdf"
Landmask file:                          ./LDT_Params/input/data/cap2ldt_ps41.nc
#Landmask spatial transform:             none   # none | mode | neighbor
Landmask spatial transform:             mode    # Go from ~10 KM to 25KM domain ...
Landmask map projection:                latlon
Landmask lower left lat:               -89.9531250
Landmask lower left lon:              -179.9296875
Landmask upper right lat:               89.9531250
Landmask upper right lon:              179.9296875
Landmask resolution (dx):                0.1406250
Landmask resolution (dy):                0.0937500

#Soil parameter inputs
Soil fraction data source:              none
Soils spatial transform:                none
Soils map projection:                   latlon
Soils fill option:                      none
Porosity data source:                   none
Porosity map:                           none

#Soil texture map:
Soil texture data source:               "STATSGOFAO_Native"
Soil texture map:                       ./LDT_Params/input/LS_PARAMETERS/noah_2dparms/topsoil30snew
Soil texture spatial transform:         mode
Soil texture map projection:            latlon
Soil texture fill option:               neighbor
Soil texture fill radius:               5
Soil texture fill value:                6
Soil texture fill value for Antarctica: 16
Soil texture force exclusion of water points during fill: true

# ISRIC Soils span the latitude range of:  -62S to 87N
# ** Issues for full global domain ... 
#Soil texture data source:    ISRIC
#Soil texture map:        ./LDT_Params/input/LS_PARAMETERS/soil_parms/ISRIC/v2017/TEXMHT_M_sl1_250m.tif     # v2017 file
#Soil texture map projection:     latlon
#Soil texture spatial transform:   mode                  # none | mode | neighbor | tile
#Soil texture fill option:       neighbor                # none | neighbor
#Soil texture fill radius:         5                     # Number of pixels to search for neighbor
#Soil texture fill value:          6                     # Static value to fill where missing
#Soil texture fill value for Antarctica:  4              # 4 -- USDA value for silty-loam
#Soil texture force exclusion of water points during fill: true

#Topography parameter inputs
Elevation data source:                  "MERIT_1K"
Elevation number of bands:              1
Elevation map:                          ./LDT_Params/input/LS_PARAMETERS/topo_parms/MERIT
Elevation fill option:                  average
Elevation fill radius:                  5
Elevation fill value:                   300.

Slope data source:                      "MERIT_1K"
Slope number of bands:                  1
Slope map:                              ./LDT_Params/input/LS_PARAMETERS/topo_parms/MERIT
Slope fill option:                      average
Slope fill radius:                      5
Slope fill value:                       0

Aspect data source:                     "MERIT_1K"
Aspect number of bands:                 1
Aspect map:                             ./LDT_Params/input/LS_PARAMETERS/topo_parms/MERIT
Aspect fill option:                     average
Aspect fill radius:                     5
Aspect fill value:                      3.14159

Topography spatial transform:           average
Topography map projection:              latlon

#Albedo inputs
Albedo data source:                     "NCEP_Native"
Albedo map:                             ./LDT_Params/input/LS_PARAMETERS/noah_2dparms/albedo
Albedo climatology interval:            monthly
Albedo spatial transform:               "budget-bilinear"
Albedo map projection:                  latlon
Albedo fill option:                     neighbor
Albedo fill radius:                     5
Albedo fill value:                      0.15

#Maximum snow albedo inputs
Max snow albedo data source:            "Barlage_Native"
Max snow albedo map:                    ./LDT_Params/input/LS_PARAMETERS/noah_2dparms/maximum_snow_albedo.hdf
Max snow albedo spatial transform:      average
Max snow albedo map projection:         latlon
Max snow albedo fill option:            neighbor
Max snow albedo fill radius:            5
Max snow albedo fill value:             0.3

#Greenness inputs
Greenness data source:                  "NCEP_Native"
Greenness fraction map:                 ./LDT_Params/input/LS_PARAMETERS/noah_2dparms/gfrac
Greenness climatology interval:         monthly
Calculate min-max greenness fraction:   .false.
Greenness maximum map:                  ./LDT_Params/input/LS_PARAMETERS/noah_2dparms/gfrac_max.asc
Greenness minimum map:                  ./LDT_Params/input/LS_PARAMETERS/noah_2dparms/gfrac_min.asc
Greenness spatial transform:            "budget-bilinear"
Greenness map projection:               latlon
Greenness fill option:                  neighbor
Greenness fill radius:                  5
Greenness fill value:                   0.3
Greenness maximum fill value:           1.0
Greenness minimum fill value:           0.0

#Slope type inputs
Slope type data source:         "none"

#Bottom temperature inputs
Bottom temperature data source:         "ISLSCP1"
Bottom temperature map:                 ./LDT_Params/input/LS_PARAMETERS/noah_2dparms/SOILTEMP.60
Bottom temperature spatial transform:   "budget-bilinear"
Bottom temperature map projection:      latlon
Bottom temperature fill option:         average
Bottom temperature fill radius:         5
Bottom temperature fill value:          287.0
Bottom temperature topographic downscaling:  "lapse-rate"

#Noah-MP LSM inputs
Noah-MP PBL Height Value:               900.

