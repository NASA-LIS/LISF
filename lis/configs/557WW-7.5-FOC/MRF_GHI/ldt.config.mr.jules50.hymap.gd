# Overall driver options
LDT running mode:                       "Ensemble restart processing"
Processed LSM parameter filename:       ./output/lis_input.mr.jules50.nc

LIS number of nests:                    1
Number of surface model types:          2
Surface model types:                    "LSM"  "Openwater"
Land surface model:                     "JULES.5.0"
Routing model:                          "HYMAP2"
Lake model:                             none
Water fraction cutoff value:            0.5
Incorporate crop information:           .false.
Number of met forcing sources:          0
Met forcing sources:                    none
LDT diagnostic file:                    output/ldtlog.mr.jules50.hymap
Mask-parameter fill diagnostic file:    output/MaskParamFill.mr.jules50.hymap.log
LDT output directory:                   output
Undefined value:                        -9999.0

#Rotated GALWEM 10-km domain
Map projection of the LIS domain:       latlon
Run domain lower left lat:               -89.9531250
Run domain lower left lon:              -179.9296875
Run domain upper right lat:               89.9531250
Run domain upper right lon:              179.9296875
Run domain resolution (dx):                0.1406250
Run domain resolution (dy):                0.0937500

#Runtime options
Number of ensembles per tile:           12

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

# Ensemble restart options
LIS restart source:                               "LSM"
Ensemble restart generation mode:                 "downscale"
Number of ensembles per tile (input restart):     12
Number of ensembles per tile (output restart):    1
Ensemble restart generation sampling strategy:    "random sampling"

Input restart filename:                           ./input/rstfiles/LIS_RST_JULES50_202306010000.d01.nc
Output restart filename: ./output/LIS_RST_JULES50_202306010000_EN01.d01.nc

# ________________

#Landcover parameter inputs
Landcover data source:                  "UM_Native_Ancillary"
Landcover classification:               "JULES_PFT"
Landcover file:                         ./input/cap2ldt_ps41.nc
Landcover spatial transform:            none
Landcover map projection:               latlon
Landcover fill option:                  none

#Landmask parameter inputs
Create or readin landmask:              "create"
Landmask data source:                   "UM_Ancillary"

#Soil parameter inputs
Soil fraction data source:              none
Soils spatial transform:                none
Soils map projection:                   latlon
Soils fill option:                      none
Porosity data source:                   none
Porosity map:                           none

#Topography parameter inputs
Elevation data source:                  "MERIT_1K"
Elevation number of bands:              1
Elevation map:                          ./input/LS_PARAMETERS/topo_parms/MERIT
Elevation fill option:                  average
Elevation fill radius:                  5
Elevation fill value:                   300.

Slope data source:                      "MERIT_1K"
Slope number of bands:                  1
Slope map:                              ./input/LS_PARAMETERS/topo_parms/MERIT
Slope fill option:                      average
Slope fill radius:                      5
Slope fill value:                       0

Aspect data source:                     "MERIT_1K"
Aspect number of bands:                 1
Aspect map:                             ./input/LS_PARAMETERS/topo_parms/MERIT
Aspect fill option:                     average
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

#JULES LSM inputs
JULES soil parameter mode:              readin
JULES ancillary file:                   ./input/cap2ldt_ps41.nc

# HYMAP River routing parameters:
HYMAP river width map:                   ./input/hymap_params/lis_rivwth.bin
HYMAP river height map:                  ./input/hymap_params/lis_rivhgt.bin
HYMAP river roughness map:               ./input/hymap_params/lis_rivman.bin
HYMAP floodplain roughness map:          ./input/hymap_params/lis_fldman.bin
HYMAP river length map:                  ./input/hymap_params/lis_rivlen.bin
HYMAP floodplain height map:             ./input/hymap_params/lis_fldhgt.bin
HYMAP floodplain height levels:          10
HYMAP flow direction x map:              ./input/hymap_params/lis_nextx.bin
HYMAP flow direction y map:              ./input/hymap_params/lis_nexty.bin
HYMAP grid elevation map:                ./input/hymap_params/lis_elevtn.bin
HYMAP grid distance map:                 ./input/hymap_params/lis_nxtdst.bin
HYMAP grid area map:                     ./input/hymap_params/lis_grarea.bin
HYMAP runoff time delay map:             ./input/hymap_params/lis_getirana_paiva.bin
HYMAP runoff time delay multiplier map:  ./input/hymap_params/lis_trunoff.bin
HYMAP baseflow time delay map:           ./input/hymap_params/lis_tbasflw.bin
HYMAP reference discharge map:           ./input/hymap_params/lis_qrefer.bin
HYMAP basin mask map:                    ./input/hymap_params/lis_mask.bin
HYMAP drainage area map:                 ./input/hymap_params/lis_uparea.bin
HYMAP basin map:                         ./input/hymap_params/lis_basin.bin
HYMAP river flow type map:               ./input/hymap_params/lis_trunoff.bin
HYMAP baseflow dwi ratio map:            ./input/hymap_params/lis_trunoff.bin
HYMAP runoff dwi ratio map:              ./input/hymap_params/lis_trunoff.bin
HYMAP params spatial transform:          none
HYMAP params map projection:             latlon
HYMAP params lower left lat:               -89.9531250
HYMAP params lower left lon:              -179.9296875
HYMAP params upper right lat:               89.9531250
HYMAP params upper right lon:              179.9296875
HYMAP params resolution (dx):                0.1406250
HYMAP params resolution (dy):                0.0937500
