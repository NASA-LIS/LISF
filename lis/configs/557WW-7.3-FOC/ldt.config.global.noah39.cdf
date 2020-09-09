#Overall driver options
LDT running mode:                       "DA preprocessing"
Processed LSM parameter filename:       ./output/lis_input.global.noah39.nc
LIS number of nests:                    1
Number of surface model types:          2
Surface model types:                    "LSM"  "Openwater"
Land surface model:                     "Noah.3.9"
Routing model:                          none
Lake model:                             none
Water fraction cutoff value:            0.5
Incorporate crop information:           .false.
Number of met forcing sources:          0
Met forcing sources:                    none
LDT diagnostic file:                    ldtlog.global.noah39.cdf
Mask-parameter fill diagnostic file:    MaskParamFill.global.noah39.cdf.log
LDT output directory:                   OUTPUT
Undefined value:                        -9999.0

#Rotated GALWEM 10-km domain
Map projection of the LIS domain:       latlon
Run domain lower left lat:               -89.9531250
Run domain lower left lon:              -179.9296875
Run domain upper right lat:               89.9531250
Run domain upper right lon:              179.9296875
Run domain resolution (dx):                0.1406250
Run domain resolution (dy):                0.0937500

#Landcover parameter inputs
Landcover data source:                  "MODIS_Native"
Landcover classification:               "IGBPNCEP"
Landcover file:                         ./input/LS_PARAMETERS/noah_2dparms/igbp.bin
Landcover spatial transform:            tile
Landcover map projection:               latlon
Landcover fill option:                  neighbor
Landcover fill radius:                  5
Landcover fill value:                   10

#Landmask parameter inputs
Create or readin landmask:              "readin"
Landmask data source:                   "UKMO_CAP_Netcdf"
Landmask file:                          ./data/cap2ldt_ps41.nc
Landmask spatial transform:             none
Landmask map projection:                latlon
Landmask fill option:                   none
Landmask lower left lat:                 -89.9531250
Landmask lower left lon:                -179.9296875
Landmask upper right lat:                 89.9531250
Landmask upper right lon:                179.9296875
Landmask resolution (dx):                  0.1406250
Landmask resolution (dy):                  0.0937500

#Soil parameter inputs
Soil fraction data source:              none
Soils spatial transform:                none
Soils map projection:                   latlon
Soils fill option:                      none
Porosity data source:                   none
Porosity map:                           none

#Soil texture inputs
Soil texture data source:               "STATSGOFAO_Native"
Soil texture map:                       ./input/LS_PARAMETERS/noah_2dparms/topsoil30snew
Soil texture spatial transform:         mode
Soil texture map projection:            latlon
Soil texture fill option:               neighbor
Soil texture fill radius:               5
Soil texture fill value:                6
Soil texture fill value for Antarctica: 16
Soil texture force exclusion of water points during fill: true

#Topography parameter inputs
Elevation data source:                  "SRTM_Native"
Elevation number of bands:              1
Elevation map:                          ./input/LS_PARAMETERS/topo_parms/SRTM/SRTM30/raw_wgtopo30antarc
Elevation fill option:                  none
Elevation fill radius:                  5
Elevation fill value:                   0

Slope data source:                      "SRTM_Native"
Slope number of bands:                  1
Slope map:                              ./input/LS_PARAMETERS/topo_parms/SRTM/SRTM30/raw_wgtopo30antarc
Slope fill option:                      none
Slope fill radius:                      5
Slope fill value:                       0

Aspect data source:                     "SRTM_Native"
Aspect number of bands:                 1
Aspect map:                             ./input/LS_PARAMETERS/topo_parms/SRTM/SRTM30/raw_wgtopo30antarc
Aspect fill option:                     none
Aspect fill radius:                     5
Aspect fill value:                      3.14159

Topography spatial transform:           average
Topography map projection:              latlon

#Albedo inputs
Albedo data source:                     "NCEP_Native"
Albedo map:                             ./input/LS_PARAMETERS/noah_2dparms/albedo
Albedo climatology interval:            monthly
Albedo spatial transform:               "budget-bilinear"
Albedo map projection:                  latlon
Albedo fill option:                     neighbor
Albedo fill radius:                     5
Albedo fill value:                      0.15

#Maximum snow albedo inputs
Max snow albedo data source:            "Barlage_Native"
Max snow albedo map:                    ./input/LS_PARAMETERS/noah_2dparms/maximum_snow_albedo.hdf
Max snow albedo spatial transform:      average
Max snow albedo map projection:         latlon
Max snow albedo fill option:            neighbor
Max snow albedo fill radius:            5
Max snow albedo fill value:             0.3

#Greenness inputs
Greenness data source:                  "NCEP_Native"
Greenness fraction map:                 ./input/LS_PARAMETERS/noah_2dparms/gfrac
Greenness climatology interval:         monthly
Calculate min-max greenness fraction:   .false.
Greenness maximum map:                  ./input/LS_PARAMETERS/noah_2dparms/gfrac_max.asc
Greenness minimum map:                  ./input/LS_PARAMETERS/noah_2dparms/gfrac_min.asc
Greenness spatial transform:            "budget-bilinear"
Greenness map projection:               latlon
Greenness fill option:                  neighbor
Greenness fill radius:                  5
Greenness fill value:                   0.3
Greenness maximum fill value:           1.0
Greenness minimum fill value:           0.0

#Slope type inputs
Slope type data source:                 none

#Bottom temperature inputs
Bottom temperature data source:         "ISLSCP1"
Bottom temperature map:                 ./input/LS_PARAMETERS/noah_2dparms/SOILTEMP.60
Bottom temperature spatial transform:   "budget-bilinear"
Bottom temperature map projection:      latlon
Bottom temperature fill option:         average
Bottom temperature fill radius:         5
Bottom temperature fill value:          287.0
Bottom temperature topographic downscaling:  "lapse-rate"

#Runtime options
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

# CDF generation options
DA preprocessing method:                          "CDF generation"
DA observation source:                            "LIS LSM soil moisture"
Name of the preprocessed DA file:                 ./noah39_cdf_200obs.nc
Apply anomaly correction to obs:                  0
Temporal resolution of CDFs:                      "monthly"
Number of bins to use in the CDF:                 100
Observation count threshold:                      200
Temporal averaging interval:                      1da
Apply external mask:                              0
External mask directory:                          none

LIS soil moisture output model name:              "Noah.3.9"
LIS soil moisture output directory:               ./output_noah39_6hr
LIS soil moisture output format:                  "netcdf"
LIS soil moisture output map projection:          "latlon"
LIS soil moisture output methodology:             "2d gridspace"
LIS soil moisture output naming style:            "3 level hierarchy"
LIS soil moisture output nest index:              1

LIS soil moisture domain lower left lat:       -89.9531250
LIS soil moisture domain lower left lon:      -179.9296875
LIS soil moisture domain upper right lat:       89.9531250
LIS soil moisture domain upper right lon:      179.9296875
LIS soil moisture domain resolution (dx):        0.1406250
LIS soil moisture domain resolution (dy):        0.0937500

Starting year:                                    2008
Starting month:                                      1
Starting day:                                        2
Starting hour:                                       0
Starting minute:                                     0
Starting second:                                     0
Ending year:                                      2019
Ending month:                                        6
Ending day:                                          1
Ending hour:                                         0
Ending minute:                                       0
Ending second:                                       0
LIS output timestep:                              6hr

