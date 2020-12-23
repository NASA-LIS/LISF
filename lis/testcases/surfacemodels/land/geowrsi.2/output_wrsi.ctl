* DSET wrsi time steps folder
DSET ^OUTPUT/SURFACEMODEL/%ch
options template
DTYPE netcdf
TITLE WRSI output
UNDEF -9999.0
XDEF  294  LINEAR    22.05  0.1
YDEF  348  LINEAR   -11.65  0.1
ZDEF 1 LINEAR 1 1
* TDEF wrsi number of time steps
TDEF 36 LINEAR 06:00Z10jun2009 10dy
* CHSUB wrsi time steps
chsub  1   1  200906/LIS_HIST_200906100000.d01.nc
chsub  2   2  200906/LIS_HIST_200906200000.d01.nc
chsub  3   3  200906/LIS_HIST_200906300000.d01.nc
chsub  4   4  200907/LIS_HIST_200907100000.d01.nc
chsub  5   5  200907/LIS_HIST_200907200000.d01.nc
chsub  6   6  200907/LIS_HIST_200907310000.d01.nc
chsub  7   7  200908/LIS_HIST_200908100000.d01.nc
chsub  8   8  200908/LIS_HIST_200908200000.d01.nc
chsub  9   9  200908/LIS_HIST_200908310000.d01.nc
chsub  10  10 200909/LIS_HIST_200909100000.d01.nc
chsub  11  11 200909/LIS_HIST_200909200000.d01.nc
chsub  12  12 200909/LIS_HIST_200909300000.d01.nc
chsub  13  13 200910/LIS_HIST_200910100000.d01.nc
chsub  14  14 200910/LIS_HIST_200910200000.d01.nc
chsub  15  15 200910/LIS_HIST_200910310000.d01.nc
chsub  16  16 200911/LIS_HIST_200911100000.d01.nc
chsub  17  17 200911/LIS_HIST_200911200000.d01.nc
chsub  18  18 200911/LIS_HIST_200911300000.d01.nc
chsub  19  19 200912/LIS_HIST_200912100000.d01.nc
chsub  20  20 200912/LIS_HIST_200912200000.d01.nc
chsub  21  21 200912/LIS_HIST_200912310000.d01.nc
chsub  22  22 201001/LIS_HIST_201001100000.d01.nc
chsub  23  23 201001/LIS_HIST_201001200000.d01.nc
chsub  24  24 201001/LIS_HIST_201001310000.d01.nc
chsub  25  25 201002/LIS_HIST_201002100000.d01.nc
chsub  26  26 201002/LIS_HIST_201002200000.d01.nc
chsub  27  27 201002/LIS_HIST_201002280000.d01.nc
chsub  28  28 201003/LIS_HIST_201003100000.d01.nc
chsub  29  29 201003/LIS_HIST_201003200000.d01.nc
chsub  30  30 201003/LIS_HIST_201003310000.d01.nc
chsub  31  31 201004/LIS_HIST_201004100000.d01.nc
chsub  32  32 201004/LIS_HIST_201004200000.d01.nc
chsub  33  33 201004/LIS_HIST_201004300000.d01.nc
chsub  34  34 201005/LIS_HIST_201005100000.d01.nc
chsub  35  35 201005/LIS_HIST_201005200000.d01.nc
chsub  36  36 201005/LIS_HIST_201005310000.d01.nc
VARS 42
lat=>lat                     1 y,x    Latitude  [-]
lon=>lon                     1 y,x    Longitude [-]
Rainf_inst=>Rainf            1 y,x  ** Total rainfall precip [kg/m2]
PotEvap_inst=>PotEvap        1 y,x  ** Total potential evaporation [kg/m2]
SOS_inst=>SOS                1 y,x  ** Start-of-season [in dekads]
WRSI_inst=>WRSI              1 y,x  ** Water requirement satisfaction index [ratio]
KF2_inst=>KF2                1 y,x  ** Percent of Season [%]
SumWR_inst=>SumWR            1 y,x  ** Sum of Water Requirement [mm]
SumET_inst=>SumET            1 y,x  ** Sum of Evapotranspiration [mm]
SWI_inst=>SWI                1 y,x  ** Soil Water Index [%]
SOSa_inst=>SOSa              1 y,x  ** Start-of-season Anomaly [in dekads]
TotalSurplusWater_inst=>TotalSurplusWat      1 y,x  ** Total surplus water ~ different stages of crop growth [mm]
MaxSurplusWater_inst=>MaxSurplusWater        1 y,x  ** Max surplus water experienced in 1 dekad [mm]
TotalWaterDeficit_inst=>TotalWaterDefic      1 y,x  ** Total water deficit ~ different stages of crop growth [mm]
MaxWaterDeficit_inst=>MaxWaterDeficit        1 y,x  ** Max water deficit experienced in 1 dekad [mm]
TotalAETInitial_inst=>TotalAETInitial        1 y,x  ** Actual evapotranspiration ~ Initial stage [mm]
TotalWRInitial_inst=>TotalWRInitial          1 y,x  ** Water requirement         ~ Initial stage [mm]
TotalSurplusWaterInitial_inst=>InitTotSurWat 1 y,x  ** Surplus water      ~ Initial stage [mm]
TotalWaterDeficitInitial_inst=>InitTotWatDef 1 y,x  ** Water deficit      ~ Initial stage [mm]
TotalAETVeg_inst=>TotalAETVeg                1 y,x  ** Actual evapotranspiration ~ Vegetative stage [mm]
TotalWRVeg_inst=>TotalWRVeg                  1 y,x  ** Water requirement         ~ Vegetative stage [mm]
TotalSurplusWaterVeg_inst=>VegtTotSurWat     1 y,x  ** Surplus water      ~ Vegetative stage [mm]
TotalWaterDeficitVeg_inst=>VegtTotWatDef     1 y,x  ** Water deficit      ~ Vegetative stage [mm]
TotalAETFlower_inst=>TotalAETFlower          1 y,x  ** Actual evapotranspiration ~ Flowering stage [mm]
TotalWRFlower_inst=>TotalWRFlower            1 y,x  ** Water requirement         ~ Flowering stage [mm]
TotalSurplusWaterFlower_inst=>FlwrTotSurWat  1 y,x  ** Surplus water       ~ Flowering stage [mm]
TotalWaterDeficitFlower_inst=>FlwrTotWatDef  1 y,x  ** Water deficit       ~ Flowering stage [mm]
TotalAETRipe_inst=>TotalAETRipe              1 y,x  ** Actual evapotranspiration ~ Ripening stage [mm]
TotalWRRipe_inst=>TotalWRRipe                1 y,x  ** Water requirement         ~ Ripening stage [mm]
TotalSurplusWaterRipe_inst=>RipeTotSurWat    1 y,x  ** Surplus water       ~ Ripening stage [mm]
TotalWaterDeficitRipe_inst=>RipeTotWatDef    1 y,x  ** Water deficit       ~ Ripening stage [mm]
PermWiltDate_inst=>PermWiltDate              1 y,x  ** Permanent wilting date [dekad]
Wilting1_inst=>Wilting1                      1 y,x  ** First wilting date [dekad]
Wilting2_inst=>Wilting2                      1 y,x  ** Second wilting date [dekad]
WRSIa_inst=>WRSIa                            1 y,x  ** WRSI anomaly [-]
growing_season_inst=>growing_season          1 y,x  ** Growing season [year]
WHC_inst=>WHC                                1 y,x  ** Water holding capacity; parameter [mm]
LGP_inst=>LGP                                1 y,x  ** Length of growing period; parameter [dekad]
WR_TimeStep_inst=>WR_TS                      1 y,x  ** Water resource; TS - timestep []; added by SY
AET_TimeStep_inst=>AET_TS                    1 y,x  ** Actual ET timestep [ ]
WRSI_TimeStep_inst=>WRSI_TS                  1 y,x  ** WRSI timestep [ ]
SurplusWater_TimeStep_inst=>SurplWat_TS      1 y,x  ** Surplus water timestep [ ]
ENDVARS
