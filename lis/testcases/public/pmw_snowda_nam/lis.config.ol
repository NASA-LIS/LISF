#Overall driver options
Running mode: 		         "retrospective"
Map projection of the LIS domain: "latlon"
Number of nests:                 1
Number of surface model types:    1
Surface model types:            "LSM"
Surface model output interval:  "1da" #"1hr"
Land surface model:             "Noah.3.3"
Number of met forcing sources:   1
Blending method for forcings:    "overlay"
Met forcing sources:             "NAM242"
Topographic correction method (met forcing):  "lapse-rate"
Enable spatial downscaling of precipitation:  0
Spatial interpolation method (met forcing):   "bilinear"
Spatial upscaling method (met forcing):       "average"
Temporal interpolation method (met forcing):  "linear"

#Runtime options
Forcing variables list file:               ./forcing_variables.txt
Output forcing:                            1   #1-yes
Output parameters:                         0   #0- no
Output methodology:                        "2d gridspace"
Output model restart files:                1
Output data format:                        "netcdf"
Output naming style:                       "3 level hierarchy"
Start mode:                                "restart"
Starting year:                             2011
Starting month:                            3
Starting day:                              1
Starting hour:                             6
Starting minute:                           0
Starting second:                           0
Ending year:                               2011
Ending month:                              5
Ending day:                                31
Ending hour:                               0
Ending minute:                             0
Ending second:                             0
Undefined value:                          -9999
Output directory:                         'OUTPUT/ol' 
Diagnostic output file:                   'diag/lislog.ol'
Number of ensembles per tile:              20

#The following options are used for subgrid tiling based on vegetation
Maximum number of surface type tiles per grid:     1
Minimum cutoff percentage (surface type tiles):    0.10 
Maximum number of soil texture tiles per grid:     1
Minimum cutoff percentage (soil texture tiles):    0.10
Maximum number of soil fraction tiles per grid:    1
Minimum cutoff percentage (soil fraction tiles):   0.10
Maximum number of elevation bands per grid:        1
Minimum cutoff percentage (elevation bands):       0.10
Maximum number of slope bands per grid:            1
Minimum cutoff percentage (slope bands):           0.10
Maximum number of aspect bands per grid:           1
Minimum cutoff percentage (aspect bands):          0.10

#Processor Layout	
#Should match the total number of processors used

Number of processors along x:    16
Number of processors along y:    12
Halo size along x: 0
Halo size along y: 0 

#------------------------ ROUTING -------------------------------------

Routing model:                    "none"

#------------------------RADIATIVE TRANSFER MODELS--------------------------

Radiative transfer model:   "none"

#------------------------APPLICATION MODELS---------------------------------

Number of application models: 0

#---------------------DATA ASSIMILATION ----------------------------------
#Data Assimilation Options

Number of data assimilation instances:               0
Data assimilation algorithm:                         "Direct insertion"
Data assimilation set:                               "SNODEP" 
Number of state variables:                           2
Data assimilation exclude analysis increments:       1
Data assimilation output interval for diagnostics:   "1da"  
Data assimilation number of observation types:       1 
Data assimilation output ensemble members:           0
Data assimilation output processed observations:     0
Data assimilation output innovations:                0

Apply perturbation bias correction:        0
Bias estimation algorithm:                "none"
Bias estimation attributes file:          "none"
Bias estimation restart output frequency:
Bias estimation start mode:
Bias estimation restart file:

#Perturbation options
Perturbations start mode:                 "coldstart"
Perturbations restart output interval:    "1mo"
Perturbations restart filename:           "none"

Forcing perturbation algorithm:           "GMAO scheme" 
Forcing perturbation frequency:           "1hr"
Forcing attributes file:                  ./forcing_attribs.txt
Forcing perturbation attributes file:     ./forcing_pert_attribs.txt

State perturbation algorithm:             "none"
State perturbation frequency:             "3hr"
State attributes file:                    "none"
State perturbation attributes file:       "none"

Observation perturbation algorithm:       "none"
Observation perturbation frequency:       "6hr"
Observation attributes file:              "none"
Observation perturbation attributes file: "none"


#------------------------DOMAIN SPECIFICATION--------------------------
#Definition of Running Domain
#Specify the domain extremes in latitude and longitude

Run domain lower left lat:                  50.025 #   30.125
Run domain lower left lon:                -172.975 
Run domain upper right lat:                 75.725 #   83.625
Run domain upper right lon:               -130.025 #  -63.125
Run domain resolution (dx):                  0.05
Run domain resolution (dy):                  0.05

#The following options list the choice of parameter maps to be 
#used

Landmask data source:            "LDT"
Landcover data source:           "LDT"
Soil texture data source:        "LDT"
Soil fraction data source:       "none"
Soil color data source:          "none"
Elevation data source:           "LDT"
Slope data source:               "none"
Aspect data source:              "none"
Curvature data source:           "none"
LAI data source:                 "none"
SAI data source:                 "none"
Albedo data source:              "LDT"
Max snow albedo data source:     "none"
Greenness data source:           "LDT"  
Roughness data source:           "none"  
Porosity data source:            "none"
Ksat data source:                "none"
B parameter data source:         "none"
Quartz data source:              "none"
Emissivity data source:          "none"

LIS domain and parameter data file: ./lis_input.d01.nc
Use greenness fraction climatology: 1
Use albedo climatology: 1
Albedo climatology interval type: "monthly"

#--------------------------------FORCINGS----------------------------------
NAM242 forcing directory: ./input/MET_FORCING/NAM242

#-----------------------LAND SURFACE MODELS--------------------------
TEMPLATE model timestep:              30mn
TEMPLATE restart output interval:     1da
TEMPLATE restart file:                none

Noah.3.3 model timestep:                  "1hr"
Noah.3.3 restart output interval:         "1mo"
Noah.3.3 restart file:                    ./input/LIS_RST_NOAH33_201102282300.ens20.d01.nc
Noah.3.3 vegetation parameter table:      ./input/PARAMETERS/noah33_parms/VEGPARM.TBL
Noah.3.3 soil parameter table:            ./input/PARAMETERS/noah33_parms/SOILPARM.TBL
Noah.3.3 general parameter table:         ./input/PARAMETERS/noah33_parms/GENPARM.TBL
Noah.3.3 use PTF for mapping soil properties: 0
Noah.3.3 soils scheme:                    2
Noah.3.3 number of soil layers:           4
Noah.3.3 layer thicknesses:               0.1  0.3  0.6  1.0
Noah.3.3 initial skin temperature:        288.0
Noah.3.3 initial soil temperatures:       288.0  288.0  288.0  288.0
Noah.3.3 initial total soil moistures:    0.080 0.080 0.080 0.080
Noah.3.3 initial liquid soil moistures:   0.080 0.080 0.080 0.080
Noah.3.3 initial canopy water:            0.0
Noah.3.3 initial snow depth:              0.0
Noah.3.3 initial snow equivalent:         0.0
Noah.3.3 fixed max snow albedo:           0.0
Noah.3.3 fixed deep soil temperature:     0.0
Noah.3.3 fixed vegetation type:           0
Noah.3.3 fixed soil type:                 0
Noah.3.3 fixed slope type:                0
Noah.3.3 sfcdif option:                   1
Noah.3.3 z0 veg-type dependence option:   0
Noah.3.3 greenness fraction:  0.01  0.02  0.07  0.17  0.27  0.58  0.93  0.96  0.65  0.24  0.11  0.02
Noah.3.3 background albedo:   0.18  0.17  0.16  0.15  0.15  0.15  0.15  0.16  0.16  0.17  0.17  0.18
Noah.3.3 background roughness length: 0.020 0.020 0.025 0.030 0.035 0.036 0.035 0.030 0.027 0.025 0.020 0.020
Noah.3.3 reference height for forcing T and q:   6.0
Noah.3.3 reference height for forcing u and v:   6.0

#---------------------------MODEL OUTPUT CONFIGURATION-----------------------
#Specify the list of ALMA variables that need to be featured in the 
#LSM model output

Output start year:
Output start month:
Output start day:
Output start hour:
Output start minutes:
Output start seconds:

Output GRIB Table Version: 128
Output GRIB Center Id:     57
Output GRIB Subcenter Id:  2
Output GRIB Process Id:    88
Output GRIB Grid Id:       255

Model output attributes file: './MODEL_OUTPUT_LIST.TBL'      
