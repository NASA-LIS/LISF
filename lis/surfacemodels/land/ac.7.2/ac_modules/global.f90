module ac_global

use ac_kinds, only: int8, &
                    int16, &
                    int32, &
                    intEnum, &
                    sp
use ac_project_input, only: GetNumberSimulationRuns, &
                            ProjectInput
use ac_utils, only: roundc, &
                    GetReleaseDate, &
                    GetVersionString, &
                    trunc
use iso_fortran_env, only: iostat_end
implicit none


real(sp), parameter :: equiv = 0.64_sp
    !! conversion factor: 1 dS/m = 0.64 g/l
integer(int32), parameter :: max_SoilLayers = 5
integer(int32), parameter :: max_No_compartments = 12
real(sp), parameter :: undef_double = -9.9_sp
    !! value for 'undefined' real(sp) variables
integer(int32), parameter :: undef_int = -9
    !! value for 'undefined' int32 variables
real(sp), parameter :: PI = 3.1415926535_sp
real(sp), parameter :: CO2Ref = 369.41_sp;
    !! reference CO2 in ppm by volume for year 2000 for Mauna Loa
    !! (Hawaii,USA)
real(sp), parameter :: EvapZmin = 15._sp
    !! cm  minimum soil depth for water extraction by evaporation
real(sp), parameter :: eps =10E-08
real(sp), parameter :: ac_zero_threshold = 0.000001_sp
real(sp), dimension(12), parameter :: ElapsedDays = [0._sp, 31._sp, 59.25_sp, &
                                                    90.25_sp, 120.25_sp, 151.25_sp, &
                                                    181.25_sp, 212.25_sp, 243.25_sp, &
                                                    273.25_sp, 304.25_sp, 334.25_sp]
integer(int32), dimension (12), parameter :: DaysInMonth = [31,28,31,30,31,30, &
                                                            31,31,30,31,30,31]

character(len=10), dimension(12), parameter :: NameMonth = &
    [character(len=10) :: 'January','February','March','April','May','June',&
                'July','August','September','October','November','December']


integer(intEnum), parameter :: modeCycle_GDDays = 0
    !! index of GDDays in modeCycle enumerated type
integer(intEnum), parameter :: modeCycle_CalendarDays = 1
    !! index of CalendarDays in modeCycle enumerated type

integer(intEnum), parameter :: pMethod_NoCorrection = 0
    !! index of NoCorrection in pMethod enumerated type
integer(intEnum), parameter :: pMethod_FAOCorrection = 1
    !! index of FAOCorrection in pMethod enumerated type

integer(intEnum), parameter :: subkind_Vegetative = 0
    !! index of Vegetative in subkind enumerated type
integer(intEnum), parameter :: subkind_Grain = 1
    !! index of Grain in subkind enumerated type
integer(intEnum), parameter :: subkind_Tuber = 2
    !! index of Tuber in subkind enumerated type
integer(intEnum), parameter :: subkind_Forage = 3
    !! index of Forage in subkind enumerated type

integer(intEnum), parameter :: plant_seed = 0
    !! index of seed in planting enumerated type
integer(intEnum), parameter :: plant_transplant = 1
    !! index of transplant in planting enumerated type
integer(intEnum), parameter :: plant_regrowth= 2
    !! index of regrowth in planting enumerated type

integer(intEnum), parameter :: Method_full = 0
    !! index of full in Method enumerated type
integer(intEnum), parameter :: Method_usda = 1
    !! index of usda in Method enumerated type
integer(intEnum), parameter :: Method_percentage= 2
    !! index of percentage in Method enumerated type

integer(intEnum), parameter :: EffectiveRainMethod_full = 0
    !! index of full in EffectiveRainMethod enumerated type
integer(intEnum), parameter :: EffectiveRainMethod_usda = 1
    !! index of usda in EffectiveRainMethod enumerated type
integer(intEnum), parameter :: EffectiveRainMethod_percentage= 2
    !! index of percentage in EffectiveRainMethod enumerated type

integer(intEnum), parameter :: TimeCuttings_NA = 0
    !! index of NA in TimeCuttings enumerated type
integer(intEnum), parameter :: TimeCuttings_IntDay = 1
    !! index of IntDay in TimeCuttings enumerated type
integer(intEnum), parameter :: TimeCuttings_IntGDD = 2
    !! index of IntGDD in TimeCuttings enumerated type
integer(intEnum), parameter :: TimeCuttings_DryB = 3
    !! index of DryB in TimeCuttings enumerated type
integer(intEnum), parameter :: TimeCuttings_DryY = 4
    !! index of DryY in TimeCuttings enumerated type
integer(intEnum), parameter :: TimeCuttings_FreshY= 5
    !! index of FreshY in TimeCuttings enumerated type

integer(intEnum), parameter :: Criterion_CumulRain = 0
    !! index of CumulRain in Criterion enumerated type
integer(intEnum), parameter :: Criterion_RainPeriod = 1
    !! index of RainPeriod in Criterion enumerated type
integer(intEnum), parameter :: Criterion_RainDecade = 2
    !! index of RainDecade in Criterion enumerated type
integer(intEnum), parameter :: Criterion_RainVsETo = 3
    !! index of RainVsETo in Criterion enumerated type

integer(intEnum), parameter :: AirTCriterion_TminPeriod = 0
    !! index of TminPeriod in AirTCriterion enumerated type
integer(intEnum), parameter :: AirTCriterion_TmeanPeriod = 1
    !! index of TmeanPeriod in AirTCriterion enumerated type
integer(intEnum), parameter :: AirTCriterion_GDDPeriod = 2
    !! index of GDDPeriod in AirTCriterion enumerated type
integer(intEnum), parameter :: AirTCriterion_CumulGDD = 3
    !! index of CumulGDD in AirTCriterion enumerated type

integer(intEnum), parameter :: GenerateTimeMode_FixInt = 0
    !! index of FixInt in GenerateTimeMode enumerated type
integer(intEnum), parameter :: GenerateTimeMode_AllDepl = 1
    !! index of AllDepl in GenerateTimeMode enumerated type
integer(intEnum), parameter :: GenerateTimeMode_AllRAW = 2
    !! index of AllRAW in GenerateTimeMode enumerated type
integer(intEnum), parameter :: GenerateTimeMode_WaterBetweenBunds = 3
    !! index of WaterBetweenBunds in GenerateTimeMode enumerated type

integer(intEnum), parameter :: GenerateDepthMode_ToFC = 0
    !! index of ToFC in GenerateDepthMode enumerated type
integer(intEnum), parameter :: GenerateDepthMode_FixDepth = 1
    !! index of FixDepth in GenerateDepthMode enumerated type

integer(intEnum), parameter :: IrriMode_NoIrri = 0
    !! index of NoIrri in IrriMode enumerated type
integer(intEnum), parameter :: IrriMode_Manual = 1
    !! index of Manual in IrriMode enumerated type
integer(intEnum), parameter :: IrriMode_Generate = 2
    !! index of Generate in IrriMode enumerated type
integer(intEnum), parameter :: IrriMode_Inet = 3
    !! index of inet in IrriMode enumerated type

integer(intEnum), parameter :: IrriMethod_MBasin = 0
    !! index of MBasin in IrriMethod enumerated type
integer(intEnum), parameter :: IrriMethod_MBorder = 1
    !! index of MBorder in IrriMethod enumerated type
integer(intEnum), parameter :: IrriMethod_MDrip = 2
    !! index of MDrip in IrriMethod enumerated type
integer(intEnum), parameter :: IrriMethod_MFurrow = 3
    !! index of MFurrow in IrriMethod enumerated type
integer(intEnum), parameter :: IrriMethod_MSprinkler = 4
    !! index of MSprinkler in IrriMethod enumerated type

integer(intEnum), parameter :: datatype_daily = 0
    !! index of daily in datatype enumerated type
integer(intEnum), parameter :: datatype_decadely = 1
    !! index of decadely in datatype enumerated type
integer(intEnum), parameter :: datatype_monthly= 2
    !! index of monthly in datatype enumerated type

integer(intEnum), parameter :: typeproject_typepro = 0
    !! index of TypePRO in typeproject enumerated type
integer(intEnum), parameter :: typeproject_typeprm = 1
    !! index of TypePRM in typeproject enumerated type
integer(intEnum), parameter :: typeproject_typenone = 2
    !! index of TypeNone in typeproject enumerated type

integer(intEnum), parameter :: typeObsSim_ObsSimCC = 0
    !! index of ObsSimCC in typeObsSim enumerated type
integer(intEnum), parameter :: typeObsSim_ObsSimB = 1
    !! index of ObsSimB in typeObsSim enumerated type
integer(intEnum), parameter :: typeObsSim_ObsSimSWC = 2
    !! index of ObsSimSWC in typeObsSim enumerated type

type rep_DayEventInt
    integer(int32) :: DayNr
        !! Undocumented
    integer(int32) :: param
        !! Undocumented
end type rep_DayEventInt

type CompartmentIndividual
    real(sp) :: Thickness
        !! meter
    real(sp) :: theta
        !! m3/m3
    real(sp) :: fluxout
        !! mm/day
    integer(int32) :: Layer
        !! Undocumented
    real(sp) :: Smax
        !! Maximum root extraction m3/m3.day
    real(sp) :: FCadj
        !! Vol % at Field Capacity adjusted to Aquifer
    integer(int32) :: DayAnaero
        !! number of days under anaerobic conditions
    real(sp) :: WFactor
        !! weighting factor 0 ... 1
        !! Importance of compartment in calculation of
        !! - relative wetness (RUNOFF)
        !! - evaporation process
        !! - transpiration process *)
    !! salinity factors
    real(sp), dimension(11) :: Salt
        !! salt content in solution in cells (g/m2)
    real(sp), dimension(11) :: Depo
        !! salt deposit in cells (g/m2)
end type CompartmentIndividual

type SoilLayerIndividual
    character(len=25) :: Description
        !! Undocumented
    real(sp) :: Thickness
        !! meter
    real(sp) :: SAT
        !! Vol % at Saturation
    real(sp) :: FC
        !! Vol % at Field Capacity
    real(sp) :: WP
        !! Vol % at Wilting Point
    real(sp) :: tau
        !! drainage factor 0 ... 1
    real(sp) :: InfRate
        !! Infiltration rate at saturation mm/day
    integer(int8) :: Penetrability
        !! root zone expansion rate in percentage
    integer(int8) :: GravelMass
        !! mass percentage of gravel
    real(sp) :: GravelVol
        !! volume percentage of gravel
    real(sp) :: WaterContent
        !! mm
    ! salinity parameters (cells)
    integer(int8) :: Macro
        !! Macropores : from Saturation to Macro [vol%]
    real(sp), dimension(11) :: SaltMobility
        !! Mobility of salt in the various salt cellS
    integer(int8) :: SC
        !! number of Saltcels between 0 and SC/(SC+2)*SAT vol%
    integer(int8) :: SCP1
        !! SC + 1   (1 extra Saltcel between SC/(SC+2)*SAT vol% and SAT)
        !! THis last celL is twice as large as the other cels *)
    real(sp) :: UL
        !! Upper Limit of SC salt cells = SC/(SC+2) * (SAT/100) in m3/m3
    real(sp) :: Dx
        !! Size of SC salt cells [m3/m3] = UL/SC
    ! capilary rise parameters
    integer(int8) :: SoilClass
        !! 1 = sandy, 2 = loamy, 3 = sandy clayey, 4 - silty clayey soils
    real(sp) :: CRa, CRb
        !! coefficients for Capillary Rise
end type SoilLayerIndividual

type rep_Shapes
    integer(int8) :: Stress
        !! Percentage soil fertility stress for calibration
    real(sp) :: ShapeCGC
        !! Shape factor for the response of Canopy Growth Coefficient to soil
        !! fertility stress
    real(sp) :: ShapeCCX
        !! Shape factor for the response of Maximum Canopy Cover to soil
        !! fertility stress
    real(sp) :: ShapeWP
        !! Shape factor for the response of Crop Water Producitity to soil
        !! fertility stress
    real(sp) :: ShapeCDecline
        !! Shape factor for the response of Decline of Canopy Cover to soil
        !! fertility stress
    logical :: Calibrated
        !! Undocumented
end type rep_Shapes

type rep_soil
    integer(int8) :: REW
        !! (* Readily evaporable water mm *)
    integer(int8) :: NrSoilLayers
    integer(int8) :: CNvalue
    real(sp) :: RootMax
        !! maximum rooting depth in soil profile for selected crop
end type rep_soil

type rep_Assimilates
    logical :: On
        !! Undocumented
    integer(int32) :: Period
        !! Number of days at end of season during which assimilates are stored in root system
    integer(int8) :: Stored
        !! Percentage of assimilates, transferred to root system at last day of season
    integer(int8) :: Mobilized
        !! Percentage of stored assimilates, transferred to above ground parts in next season
end type rep_Assimilates

type rep_Onset
    logical :: GenerateOn
        !! by rainfall or temperature criterion
    logical :: GenerateTempOn
        !! by temperature criterion
    integer(intEnum) :: Criterion
        !! Undocumented
    integer(intEnum) :: AirTCriterion
        !! Undocumented
    integer(int32) :: StartSearchDayNr
        !! daynumber
    integer(int32) :: StopSearchDayNr
        !! daynumber
    integer(int32) :: LengthSearchPeriod
        !! days
end type rep_Onset

type rep_EndSeason
    integer(int32) :: ExtraYears
        !! to add to YearStartCropCycle
    logical :: GenerateTempOn
        !! by temperature criterion
    integer(intEnum) :: AirTCriterion
        !! Undocumented
    integer(int32) :: StartSearchDayNr
        !! daynumber
    integer(int32) :: StopSearchDayNr
        !! daynumber
    integer(int32) :: LengthSearchPeriod
        !! days
end type rep_EndSeason

type rep_Content
    real(sp) :: BeginDay
        !! at the beginning of the day
    real(sp) :: EndDay
        !! at the end of the day
    real(sp) :: ErrorDay
        !! error on WaterContent or SaltContent over the day
end type rep_Content

type rep_EffectStress
    integer(int8) :: RedCGC
        !! Reduction of CGC (%)
    integer(int8) :: RedCCX
        !! Reduction of CCx (%)
    integer(int8) :: RedWP
        !! Reduction of WP (%)
    real(sp) :: CDecline
        !! Average decrease of CCx in mid season (%/day)
    integer(int8) :: RedKsSto
        !! Reduction of KsSto (%)
end type rep_EffectStress

type rep_EffectiveRain
    !!  for 10-day or monthly rainfall data
    integer(intEnum) :: Method
        !! Undocumented
    integer(int8) :: PercentEffRain
        !! IF Method = Percentage
    integer(int8) :: ShowersInDecade
        !! adjustment of surface run-off
    integer(int8) :: RootNrEvap
        !! Root for reduction in soil evaporation
end type rep_EffectiveRain

type rep_RootZoneWC
    real(sp) :: Actual
        !! actual soil water content in rootzone [mm]
    real(sp) :: FC
        !! soil water content [mm] in rootzone at FC
    real(sp) :: WP
        !! soil water content [mm] in rootzone at WP
    real(sp) :: SAT
        !! soil water content [mm] in rootzone at Sat
    real(sp) :: Leaf
        !! soil water content [mm] in rootzone at upper Threshold for leaf
        !! expansion
    real(sp) :: Thresh
        !! soil water content [mm] in rootzone at Threshold for stomatal
        !! closure
    real(sp) :: Sen
        !! soil water content [mm] in rootzone at Threshold for canopy
        !! senescence
    real(sp) :: ZtopAct
        !! actual soil water content [mm] in top soil (= top compartment)
    real(sp) :: ZtopFC
        !! soil water content [mm] at FC in top soil (= top compartment)
    real(sp) :: ZtopWP
        !! soil water content [mm] at WP in top soil (= top compartment)
    real(sp) :: ZtopThresh
        !! soil water content [mm] at Threshold for stomatal closure in top
        !! soil
end type rep_RootZoneWC

type rep_IrriECw
    real(sp) :: PreSeason
        !! Undocumented
    real(sp) :: PostSeason
        !! Undocumented
end type rep_IrriECw

type rep_clim
    integer(intEnum) :: DataType
        !! Undocumented
    integer(int32) :: FromD, FromM, FromY
        !! D = day or decade, Y=1901 is not linked to specific year
    integer(int32) :: ToD, ToM, ToY
        !! Undocumented
    integer(int32) :: FromDayNr, ToDayNr
        !! daynumber
    character(len=:), allocatable :: FromString, ToString
        !! Undocumented
    integer(int32) :: NrObs
        !! number of observations
end type rep_clim

type rep_CropFileSet
    integer(int32) :: DaysFromSenescenceToEnd
        !! Undocumented
    integer(int32) :: DaysToHarvest
        !! given or calculated from GDD
    integer(int32) :: GDDaysFromSenescenceToEnd
        !! Undocumented
    integer(int32) :: GDDaysToHarvest
        !! given or calculated from Calendar Days
end type rep_CropFileSet

type rep_Cuttings
    logical :: Considered
        !! Undocumented
    integer(int32) :: CCcut
        !! Canopy cover (%) after cutting
    integer(int32) :: Day1
        !! first day after time window for generating cuttings (1 = start crop cycle)
    integer(int32) :: NrDays
        !! number of days of time window for generate cuttings (-9 is whole crop cycle)
    logical :: Generate
        !! ture: generate cuttings; false : schedule for cuttings
    integer(intEnum) :: Criterion
        !! time criterion for generating cuttings
    logical :: HarvestEnd
        !! final harvest at crop maturity
    integer(int32) :: FirstDayNr
        !! first dayNr of list of specified cutting events (-9 = onset growing cycle)
end type rep_Cuttings

type rep_Manag
    integer(int8) :: Mulch
        !! percent soil cover by mulch in growing period
    integer(int8) :: SoilCoverBefore
        !! percent soil cover by mulch before growing period
    integer(int8) :: SoilCoverAfter
        !! percent soil cover by mulch after growing period
    integer(int8) :: EffectMulchOffS
        !! effect Mulch on evaporation before and after growing period
    integer(int8) :: EffectMulchInS
        !! effect Mulch on evaporation in growing period
    integer(int32) :: FertilityStress
        !! Undocumented
    real(sp) :: BundHeight
        !! meter;
    logical :: RunoffOn
        !! surface runoff
    integer(int32) :: CNcorrection
        !! percent increase/decrease of CN
    integer(int8) :: WeedRC
        !! Relative weed cover in percentage at canopy closure
    integer(int32) :: WeedDeltaRC
        !! Increase/Decrease of Relative weed cover in percentage during mid season
    real(sp) :: WeedShape
        !! Shape factor for crop canopy suppression
    integer(int8) :: WeedAdj
        !! replacement (%) by weeds of the self-thinned part of the Canopy Cover - only for perennials
    type(rep_Cuttings) :: Cuttings
        !! Multiple cuttings
end type rep_Manag

type rep_param
    !!  DEFAULT.PAR
    !! crop parameters IN CROP.PAR - with Reset option
    integer(int8) :: EvapDeclineFactor
        !! exponential decline with relative soil water [1 = small ... 8 = sharp]
    real(sp) :: KcWetBare
        !! Soil evaporation coefficients from wet bare soil
    integer(int8) :: PercCCxHIfinal
        !! CC threshold below which HI no longer increase (% of 100)
    integer(int32) :: RootPercentZmin
        !! starting depth of root sine function in % of Zmin (sowing depth)
    real(sp) :: MaxRootZoneExpansion
        !! maximum root zone expansion in cm/day - fixed at 5 cm/day
    integer(int8) :: KsShapeFactorRoot
        !! shape factro for the effect of water stress on root zone expansion
    integer(int8) :: TAWGermination
        !! Soil water content (% TAW) required at sowing depth for germination
    real(sp) :: pAdjFAO
        !! Adjustment factor for FAO-adjustment of soil water depletion (p) for various ET
    integer(int32) :: DelayLowOxygen
        !! delay [days] for full effect of anaeroby
    real(sp) :: ExpFsen
        !! exponent of senescence factor adjusting drop in photosynthetic activity of dying crop
    integer(int8) :: Beta
        !! Percentage decrease of p(senescence) once early canopy senescence is triggered
    integer(int8) :: ThicknessTopSWC
        !! Thickness of top soil for determination of its Soil Water Content (cm)
    !! Field parameter IN FIELD.PAR  - with Reset option
    integer(int8) :: EvapZmax
        !! cm  maximum soil depth for water extraction by evaporation
    !! Runoff parameters IN RUNOFF.PAR  - with Reset option
    real(sp) :: RunoffDepth
        !! considered depth (m) of soil profile for calculation of mean soil water content for CN adjustment
    logical :: CNcorrection
        !! correction Antecedent Moisture Class (On/Off)
    !! Temperature parameters IN TEMPERATURE.PAR  - with Reset option
    real(sp) :: Tmin
        !! Default Minimum and maximum air temperature (degC) if no temperature file
    real(sp) :: Tmax
        !! Default Minimum and maximum air temperature (degC) if no temperature file
    integer(int8) :: GDDMethod
        !! 1 for Method 1, 2 for Method 2, 3 for Method 3
    !! General parameters IN GENERAL.PAR
    integer(int32) :: PercRAW
        !! allowable percent RAW depletion for determination Inet
    real(sp) :: CompDefThick
        !! Default thickness of soil compartments [m]
    integer(int32) :: CropDay1
        !! First day after sowing/transplanting (DAP = 1)
    real(sp) :: Tbase
        !! Default base and upper temperature (degC) assigned to crop
    real(sp) :: Tupper
        !! Default base and upper temperature (degC) assigned to crop
    integer(int8) :: IrriFwInSeason
        !! Percentage of soil surface wetted by irrigation in crop season
    integer(int8) :: IrriFwOffSeason
        !! Percentage of soil surface wetted by irrigation off-season
    !! Showers parameters (10-day or monthly rainfall) IN SHOWERS.PAR
    integer(int32), dimension(12) :: ShowersInDecade !! 10-day or Monthly rainfall --> Runoff estimate
    type(rep_EffectiveRain) :: EffectiveRain
        !! 10-day or Monthly rainfall --> Effective rainfall
    !! Salinity
    integer(int8) :: SaltDiff
        !! salt diffusion factor (capacity for salt diffusion in micro pores) [%]
    integer(int8) :: SaltSolub
        !! salt solubility [g/liter]
    !! Groundwater table
    logical :: ConstGwt
        !! groundwater table is constant (or absent) during the simulation period
    !! Capillary rise
    integer(int8) :: RootNrDF
        !! Undocumented
    !! Initial abstraction for surface runoff
    integer(int8) :: IniAbstract
        !! Undocumented
end type rep_param

type rep_sum
    real(sp) :: Epot, Tpot, Rain, Irrigation, Infiltrated
        !! Undocumented
    real(sp) :: Runoff, Drain, Eact, Tact, TrW, ECropCycle, CRwater
        !! mm
    real(sp) :: Biomass, YieldPart, BiomassPot, BiomassUnlim, BiomassTot
        !! ton/ha
    real(sp) :: SaltIn, SaltOut, CRsalt
        !! ton/ha
end type rep_sum

type rep_RootZoneSalt
    real(sp) :: ECe
        !! Electrical conductivity of the saturated soil-paste extract (dS/m)
    real(sp) :: ECsw
        !! Electrical conductivity of the soil water (dS/m)
    real(sp) :: ECswFC
        !! Electrical conductivity of the soil water at Field Capacity(dS/m)
    real(sp) :: KsSalt
        !! stress coefficient for salinity
end type rep_RootZoneSalt

type rep_IniSWC
    logical :: AtDepths
        !! at specific depths or for specific layers
    integer(int8) :: NrLoc
        !! number of depths or layers considered
    real(sp), dimension(max_No_compartments) :: Loc
        !! depth or layer thickness [m]
    real(sp), dimension(max_No_compartments) :: VolProc
        !! soil water content (vol%)
    real(sp), dimension(max_No_compartments) :: SaltECe
        !! ECe in dS/m
    logical :: AtFC
        !! If iniSWC is at FC
end type rep_IniSWC

type rep_storage
    real(sp) :: Btotal
        !! assimilates (ton/ha) stored in root systemn by CropString in Storage-Season
    character(len=:), allocatable :: CropString
        !! full name of crop file which stores Btotal during Storage-Season
    integer(int8) :: Season
        !! season in which Btotal is stored
end type rep_storage

type rep_sim
    integer(int32) :: FromDayNr
        !! daynumber
    integer(int32) :: ToDayNr
        !! daynumber
    type(rep_IniSWC) :: IniSWC
        !! Undocumented
    real(sp), dimension(max_No_compartments) :: ThetaIni
        !! dS/m
    real(sp), dimension(max_No_compartments) :: ECeIni
        !! dS/m
    real(sp) :: SurfaceStorageIni
        !! Undocumented
    real(sp) :: ECStorageIni
        !! Undocumented
    real(sp) :: CCini
        !! Undocumented
    real(sp) :: Bini
        !! Undocumented
    real(sp) :: Zrini
        !! Undocumented
    logical :: LinkCropToSimPeriod
        !! Undocumented
    logical :: ResetIniSWC
        !! soil water and salts
    integer(int32) :: InitialStep
        !! Undocumented
    logical :: EvapLimitON
        !! soil evap is before late season stage limited due to sheltering effect of (partly) withered canopy cover
    real(sp) :: EvapWCsurf
        !! remaining water (mm) in surface soil layer for stage 1 evaporation [REW .. 0]
    integer(int8) :: EvapStartStg2
        !! % extra to define upper limit of soil water content at start of stage 2 [100 .. 0]
    real(sp) :: EvapZ
        !! actual soil depth (m) for water extraction by evaporation  [EvapZmin/100 .. EvapZmax/100]
    integer(int32) :: HIfinal
        !! final Harvest Index might be smaller than HImax due to early canopy decline
    integer(int32) :: DelayedDays
        !! delayed days since sowing/planting due to water stress (crop cannot germinate)
    logical :: Germinate
        !! germinate is false when crop cannot germinate due to water stress
    real(sp) :: SumEToStress
        !! Sum ETo during stress period to delay canopy senescence
    real(sp) :: SumGDD
        !! Sum of Growing Degree-days
    real(sp) :: SumGDDfromDay1
        !! Sum of Growing Degree-days since Crop.Day1
    real(sp) :: SCor
        !! correction factor for Crop.SmaxBot if restrictive soil layer inhibit root development
    logical :: MultipleRun
        !! Project with a sequence of simulation runs
    integer(int32) :: NrRuns
        !! Undocumented
    logical :: MultipleRunWithKeepSWC
        !! Project with a sequence of simulation runs and initial SWC is once or more KeepSWC
    real(sp) :: MultipleRunConstZrx
        !! Maximum rooting depth for multiple projects with KeepSWC
    real(sp) :: IrriECw
        !! quality of irrigation water (dS/m)
    integer(int8) :: DayAnaero
        !! number of days under anaerobic conditions
    type(rep_EffectStress) :: EffectStress
        !! effect of soil fertility and salinity stress on CC, WP and KsSto
    logical :: SalinityConsidered
        !! Undocumented
    logical :: ProtectedSeedling
        !! IF protected (before CC = 1.25 CC0), seedling triggering of early senescence is switched off
    logical :: SWCtopSoilConsidered
        !! Top soil is relative wetter than root zone and determines water stresses
    integer(int32) :: LengthCuttingInterval
        !! Default length of cutting interval (days)
    integer(int8) :: YearSeason
        !! year number for perennials (1 = 1st year, 2, 3, 4, max = 127)
    integer(int8) :: RCadj
        !! adjusted relative cover of weeds with self thinning for perennials
    type(rep_storage) :: Storage
        !! Undocumented
    integer(int32) :: YearStartCropCycle
        !! calendar year in which crop cycle starts
    integer(int32) :: CropDay1Previous
        !! previous daynumber at the start of teh crop cycle
end type rep_sim

type rep_DayEventDbl
    integer(int32) :: DayNr
        !! Undocumented
    real(sp) :: Param
        !! Undocumented
end type rep_DayEventDbl

type rep_Crop
    integer(intEnum) :: subkind
        !! Undocumented
    integer(intEnum) :: ModeCycle
        !! Undocumented
    integer(intEnum) :: Planting
        !! 1 = sown, 0 = transplanted, -9 = regrowth
    integer(intEnum) :: pMethod
        !! Undocumented
    real(sp) :: pdef
        !! soil water depletion fraction for no stomatal stress as defined (ETo = 5 mm/day)
    real(sp) :: pActStom
        !! actual p for no stomatal stress for ETo of the day
    real(sp) :: KsShapeFactorLeaf
        !! Undocumented
    real(sp) :: KsShapeFactorStomata
        !! Undocumented
    real(sp) :: KsShapeFactorSenescence
        !! Undocumented
    real(sp) :: pLeafDefUL
        !! soil water depletion fraction for leaf expansion (ETo = 5 mm/day)
    real(sp) :: pLeafDefLL
        !! soil water depletion fraction for leaf expansion (ETo = 5 mm/day)
    real(sp) :: pLeafAct
        !! actual p for upper limit leaf expansion for ETo of the day
    real(sp) :: pSenescence
        !! soil water depletion fraction for canopys senescence (ETo = 5 mm/day)
    real(sp) :: pSenAct
        !! actual p for canopy senescence for ETo of the day
    real(sp) :: pPollination
        !! soil water depletion fraction for failure of pollination
    integer(int32) :: SumEToDelaySenescence
        !! Undocumented
    integer(int32) :: AnaeroPoint
        !! (SAT - [vol%]) at which deficient aeration
    type(rep_Shapes) :: StressResponse
        !! is reponse to soil fertility stress
    integer(int8) :: ECemin
        !! lower threshold for salinity stress (dS/m)
    integer(int8) :: ECemax
        !! upper threshold for salinity stress (dS/m)
    integer(int8) :: CCsaltDistortion
        !! distortion canopy cover for calibration for simulation of effect of salinity stress (%)
    integer(int32) :: ResponseECsw
        !! Response of Ks stomata to ECsw for calibration: From 0 (none) to +200 (very strong)
    real(sp) :: SmaxTopQuarter
        !! Smax Top 1/4 root zone HOOGLAND
    real(sp) :: SmaxBotQuarter
        !! Smax Bottom 1/4 root zone HOOGLAND
    real(sp) :: SmaxTop
        !! Smax Top root zone HOOGLAND
    real(sp) :: SmaxBot
        !! Smax Bottom root zone HOOGLAND
    real(sp) :: KcTop
        !! Undocumented
    real(sp) :: KcDecline
        !! Reduction Kc (%CCx/day) as result of ageing effects, nitrogen defficiency, etc.
    integer(int32) :: CCEffectEvapLate
        !! %
    integer(int32) :: Day1
        !! Daynummer: first day of croping period starting from sowing/transplanting
    integer(int32) :: DayN
        !! Daynummer: last day = harvest day
    integer(int32), dimension(4) :: Length
        !! Length : rep_int_array; @ 1 .. 4 :  = the four growth stages
    real(sp) :: RootMin
        !! rooting depth in meter
    real(sp) :: RootMax
        !! rooting depth in meter
    integer(int8) :: RootShape
        !! 10 times the root of the root function
    real(sp) :: Tbase
        !! Base Temperature (degC)
    real(sp) :: Tupper
        !! Upper temperature threshold (degC)
    integer(int8) :: Tcold
        !! Minimum air temperature below which pollination starts to fail (cold stress) (degC)
    integer(int8) :: Theat
        !! Maximum air temperature above which pollination starts to fail (heat stress) (degC)
    real(sp) :: GDtranspLow
        !! Minimum growing degrees required for full crop transpiration (degC - day)
    real(sp) :: SizeSeedling
        !! Canopy cover per seedling (cm2)
    real(sp) :: SizePlant
        !! Canopy cover of plant on 1st day (cm2) when regrowth
    integer(int32) :: PlantingDens
        !! number of plants per hectare
    real(sp) :: CCo
        !! starting canopy size  (fraction canopy cover)
    real(sp) :: CCini
        !! starting canopy size for regrowth (fraction canopy cover)
    real(sp) :: CGC
        !! Canopy growth coefficient (increase of CC in fraction per day)
    real(sp) :: GDDCGC
        !! Canopy growth coefficient (increase of CC in fraction per growing-degree day)
    real(sp) :: CCx
        !! expected maximum canopy cover  (fraction canopy cover)
    real(sp) :: CDC
        !! Canopy Decline Coefficient (decrease of CC in fraction per day)
    real(sp) :: GDDCDC
        !! Canopy Decline Coefficient (decrease of CC in fraction per growing-degree day)
    real(sp) :: CCxAdjusted
        !! maximum canopy cover given water stress
    real(sp) :: CCxWithered
        !! maximum existed CC during season (for correction Evap for withered canopy)
    real(sp) :: CCoAdjusted
        !! initial canopy size after soil water stress
    integer(int32) :: DaysToCCini
        !! required for regrowth (if CCini > CCo)
    integer(int32) :: DaysToGermination
        !! given or calculated from GDD
    integer(int32) :: DaysToFullCanopy
        !! given or calculated from GDD
    integer(int32) :: DaysToFullCanopySF
        !! adjusted to soil fertility
    integer(int32) :: DaysToFlowering
        !! given or calculated from GDD
    integer(int32) :: LengthFlowering
        !! given or calculated from GDD
    integer(int32) :: DaysToSenescence
        !! given or calculated from GDD
    integer(int32) :: DaysToHarvest
        !! given or calculated from GDD
    integer(int32) :: DaysToMaxRooting
        !! given or calculated from GDD
    integer(int32) :: DaysToHIo
        !! given or calculated from GDD
    integer(int32) :: GDDaysToCCini
        !! required for regrowth (if CCini > CCo)
    integer(int32) :: GDDaysToGermination
        !! given or calculated from Calendar Days
    integer(int32) :: GDDaysToFullCanopy
        !! given or calculated from Calendar Days
    integer(int32) :: GDDaysToFullCanopySF
        !! adjusted to soil fertility
    integer(int32) :: GDDaysToFlowering
        !! given or calculated from Calendar Days
    integer(int32) :: GDDLengthFlowering
        !! given or calculated from Calendar Days
    integer(int32) :: GDDaysToSenescence
        !! given or calculated from Calendar Days
    integer(int32) :: GDDaysToHarvest
        !! given or calculated from Calendar Days
    integer(int32) :: GDDaysToMaxRooting
        !! given or calculated from Calendar Days
    integer(int32) :: GDDaysToHIo
        !! given or calculated from Calendar Days
    real(sp) :: WP
        !! (normalized) water productivity (gram/m2)
    integer(int32) :: WPy
        !! (normalized) water productivity during yield formation (Percent WP)
    integer(int8) :: AdaptedToCO2
        !! Crop performance under elevated atmospheric CO2 concentration (%)
    integer(int32) :: HI
        !! HI harvest index (percentage)
    real(sp) :: dHIdt
        !! average rate of change in harvest index (% increase per calendar day)
    integer(int8) :: HIincrease
        !! possible increase (%) of HI due to water stress before flowering
    real(sp) :: aCoeff
        !! coefficient describing impact of restricted vegetative growth at flowering on HI
    real(sp) :: bCoeff
        !! coefficient describing impact of stomatal closure at flowering on HI
    integer(int8) :: DHImax
        !! allowable maximum increase (%) of specified HI
    logical :: DeterminancyLinked
        !! linkage of determinancy with flowering
    integer(int16) :: fExcess
        !! potential excess of fruits (%) ranging form
    integer(int8) :: DryMatter
        !! dry matter content (%) of fresh yield
    real(sp) :: RootMinYear1
        !! minimum rooting depth in first year in meter (for perennial crops)
    logical :: SownYear1
        !! True = Sown, False = transplanted (for perennial crops)
    integer(int8) :: YearCCx
        !! number of years at which CCx declines to 90 % of its value due to self-thinning - Perennials
    real(sp) :: CCxRoot
        !! shape factor of the decline of CCx over the years due to self-thinning - Perennials
    type(rep_Assimilates) :: Assimilates
        !! Undocumented
end type rep_Crop

type rep_PerennialPeriod
    logical :: GenerateOnset
        !! onset is generated by air temperature criterion
    integer(intEnum) :: OnsetCriterion
        !! another doscstring
    integer(int32) :: OnsetFirstDay
        !! Undocumented
    integer(int32) :: OnsetFirstMonth
        !! Undocumented
    integer(int32) :: OnsetStartSearchDayNr
        !! daynumber
    integer(int32) :: OnsetStopSearchDayNr
        !! daynumber
    integer(int32) :: OnsetLengthSearchPeriod
        !! days
    real(sp) :: OnsetThresholdValue
        !! degC or degree-days
    integer(int32) :: OnsetPeriodValue
        !! number of successive days
    integer(int8) :: OnsetOccurrence
        !! number of occurrences (1,2 or 3)
    logical :: GenerateEnd
        !! end is generate by air temperature criterion
    integer(intEnum) :: EndCriterion
        !! Undocumented
    integer(int32) :: EndLastDay
        !! Undocumented
    integer(int32) :: EndLastMonth
        !! Undocumented
    integer(int32) :: ExtraYears
        !! number of years to add to the onset year
    integer(int32) :: EndStartSearchDayNr
        !! daynumber
    integer(int32) :: EndStopSearchDayNr
        !! daynumber
    integer(int32) :: EndLengthSearchPeriod
        !! days
    real(sp) :: EndThresholdValue
        !! degC or degree-days
    integer(int32) :: EndPeriodValue
        !! number of successive days
    integer(int8) :: EndOccurrence
        !! number of occurrences (1,2 or 3)
    integer(int32) :: GeneratedDayNrOnset
        !! Undocumented
    integer(int32) :: GeneratedDayNrEnd
        !! Undocumented
end type rep_PerennialPeriod

type rep_FileOK
    logical :: Climate_Filename
    logical :: Temperature_Filename
    logical :: ETo_Filename
    logical :: Rain_Filename
    logical :: CO2_Filename
    logical :: Calendar_Filename
    logical :: Crop_Filename
    logical :: Irrigation_Filename
    logical :: Management_Filename
    logical :: GroundWater_Filename
    logical :: Soil_Filename
    logical :: SWCIni_Filename
    logical :: OffSeason_Filename
    logical :: Observations_Filename
end type rep_FileOK
    
character(len=:), allocatable :: RainFile
character(len=:), allocatable :: RainFileFull
character(len=:), allocatable :: RainDescription
character(len=:), allocatable :: EToFile
character(len=:), allocatable :: EToFileFull
character(len=:), allocatable :: EToDescription
character(len=:), allocatable :: CalendarFile
character(len=:), allocatable :: CalendarFileFull
character(len=:), allocatable :: CalendarDescription
character(len=:), allocatable :: CO2File
character(len=:), allocatable :: CO2FileFull
character(len=:), allocatable :: CO2Description
character(len=:), allocatable :: IrriFile
character(len=:), allocatable :: IrriFileFull
character(len=:), allocatable :: CropFile
character(len=:), allocatable :: CropFileFull
character(len=:), allocatable :: CropDescription
character(len=:), allocatable :: PathNameProg
character(len=:), allocatable :: PathNameOutp
character(len=:), allocatable :: PathNameSimul
character(len=:), allocatable :: ProfFile
character(len=:), allocatable :: ProfFilefull
character(len=:), allocatable :: ProfDescription
character(len=:), allocatable :: ManFile
character(len=:), allocatable :: ManFilefull
character(len=:), allocatable :: ObservationsFile
character(len=:), allocatable :: ObservationsFilefull
character(len=:), allocatable :: ObservationsDescription
character(len=:), allocatable :: OffSeasonFile
character(len=:), allocatable :: OffSeasonFilefull
character(len=:), allocatable :: OutputName
character(len=:), allocatable :: GroundWaterFile
character(len=:), allocatable :: GroundWaterFilefull
character(len=:), allocatable :: ClimateFile
character(len=:), allocatable :: ClimateFileFull
character(len=:), allocatable :: ClimateDescription
character(len=:), allocatable :: IrriDescription
character(len=:), allocatable :: ClimFile
character(len=:), allocatable :: SWCiniFile
character(len=:), allocatable :: SWCiniFileFull
character(len=:), allocatable :: SWCiniDescription
character(len=:), allocatable :: ProjectDescription
character(len=:), allocatable :: ProjectFile
character(len=:), allocatable :: ProjectFileFull
character(len=:), allocatable :: MultipleProjectDescription
character(len=:), allocatable :: MultipleProjectFile
character(len=:), allocatable :: TemperatureFile
character(len=:), allocatable :: TemperatureFileFull
character(len=:), allocatable :: TemperatureDescription
character(len=:), allocatable :: MultipleProjectFileFull
character(len=:), allocatable :: FullFileNameProgramParameters
character(len=:), allocatable :: ManDescription
character(len=:), allocatable :: ClimDescription
character(len=:), allocatable :: OffSeasonDescription
character(len=:), allocatable :: GroundwaterDescription
character(len=:), allocatable :: TnxReferenceFile
character(len=:), allocatable :: TnxReferenceFileFull
character(len=:), allocatable :: TnxReference365DaysFile
character(len=:), allocatable :: TnxReference365DaysFileFull

type(rep_IrriECw) :: IrriECw
type(rep_Manag) :: Management
type(rep_PerennialPeriod) :: perennialperiod
type(rep_param) :: simulparam
type(rep_Cuttings) :: Cuttings
type(rep_Onset) :: onset
type(rep_EndSeason) :: endseason
type(rep_Crop) :: crop
type(rep_Content) :: TotalSaltContent
type(rep_Content) :: TotalWaterContent
type(rep_EffectiveRain) :: effectiverain
type(rep_soil) :: Soil
type(rep_RootZoneWC) :: RootZoneWC
type(rep_CropFileSet) :: CropFileSet
type(rep_sum) :: SumWaBal
type(rep_RootZoneSalt) :: RootZoneSalt
type(rep_clim)  :: TemperatureRecord, ClimRecord, RainRecord, EToRecord
type(rep_sim) :: Simulation


integer(intEnum) :: GenerateTimeMode
integer(intEnum) :: GenerateDepthMode
integer(intEnum) :: IrriMode
integer(intEnum) :: IrriMethod

integer(int32) :: TnxReferenceYear
integer(int32) :: DaySubmerged
integer(int32) :: MaxPlotNew
integer(int32) :: NrCompartments
integer(int32) :: IrriFirstDayNr
integer(int32) :: ZiAqua ! Depth of Groundwater table below
                         ! soil surface in centimeter

integer(int8) :: IniPercTAW ! Default Value for Percentage TAW for Initial
                            ! Soil Water Content Menu
integer(int8) :: MaxPlotTr
integer(int8) :: OutputAggregate

integer :: fTnxReference  ! file handle
integer :: fTnxReference_iostat  ! IO status
integer :: fTnxReference365Days  ! file handle
integer :: fTnxReference365Days_iostat  ! IO status

real(sp) :: CCiActual
real(sp) :: CCiprev
real(sp) :: CCiTopEarlySen
real(sp) :: CRsalt ! gram/m2
real(sp) :: CRwater ! mm/day
real(sp) :: ECdrain ! EC drain water dS/m
real(sp) :: ECiAqua ! EC of the groundwater table in dS/m
real(sp) :: ECstorage !EC surface storage dS/m
real(sp) :: Eact ! mm/day
real(sp) :: Epot ! mm/day
real(sp) :: ETo ! mm/day
real(sp) :: Drain  ! mm/day
real(sp) :: Infiltrated ! mm/day
real(sp) :: Irrigation ! mm/day
real(sp) :: Rain  ! mm/day
real(sp) :: RootingDepth
real(sp) :: Runoff  ! mm/day
real(sp) :: SaltInfiltr ! salt infiltrated in soil profile Mg/ha
real(sp) :: Surf0 ! surface water [mm] begin day
real(sp) :: SurfaceStorage !mm/day
real(sp) :: Tact ! mm/day
real(sp) :: Tpot ! mm/day
real(sp) :: TactWeedInfested !mm/day
real(sp) :: Tmax ! degC
real(sp) :: Tmin ! degC
real(sp) :: TmaxCropReference ! degC
real(sp) :: TminCropReference ! degC
real(sp) :: TmaxTnxReference365Days ! degC
real(sp) :: TminTnxReference365Days ! degC
real(sp), dimension(1:366) :: TmaxRun, TminRun
real(sp), dimension(1:12) ::  TmaxTnxReference12MonthsRun, TminTnxReference12MonthsRun
real(sp), dimension(1:365) :: TmaxCropReferenceRun, TminCropReferenceRun
real(sp), dimension(1:365) :: TmaxTnxReference365DaysRun, TminTnxReference365DaysRun

logical :: EvapoEntireSoilSurface ! True of soil wetted by RAIN (false = IRRIGATION and fw < 1)
logical :: PreDay, OutDaily, Out8Irri
logical :: Out1Wabal
logical :: Out2Crop
logical :: Out3Prof
logical :: Out4Salt
logical :: Out5CompWC
logical :: Out6CompEC
logical :: Out7Clim
logical :: Part1Mult,Part2Eval

character(len=:), allocatable :: PathNameList,PathNameParam

type(CompartmentIndividual), dimension(max_No_compartments) :: Compartment
type(SoilLayerIndividual), dimension(max_SoilLayers) :: soillayer


type(rep_DayEventInt), dimension(5) :: IrriBeforeSeason
type(rep_DayEventInt), dimension(5) :: IrriAfterSeason


contains


real(sp) function DeduceAquaCropVersion(FullNameXXFile)
    character(len=*), intent(in) :: FullNameXXFile

    integer :: fhandle
    real(sp) :: VersionNr

    open(newunit=fhandle, file=trim(FullNameXXFile), status='old', &
         action='read')

    read(fhandle, *)  ! Description
    read(fhandle, *) VersionNr  ! AquaCrop version

    close(fhandle)

    DeduceAquaCropVersion = VersionNr
end function DeduceAquaCropVersion


real(sp) function RootMaxInSoilProfile(ZmaxCrop, TheNrSoilLayers, TheSoilLayer)
    real(sp), intent(in) :: ZmaxCrop
    integer(int8), intent(in) :: TheNrSoilLayers
    type(SoilLayerIndividual), dimension(max_SoilLayers), intent(in) :: &
                                                                TheSoilLayer

    real(sp) :: Zmax
    real(sp) :: Zsoil
    integer :: layi

    Zmax = ZmaxCrop
    Zsoil = 0._sp

    layi = 0
    do while ((layi < TheNrSoilLayers) .and. (Zmax > 0._sp))
        layi = layi + 1

        if ((TheSoilLayer(layi)%Penetrability < 100) .and. &
            (roundc(Zsoil*1000, mold=1_int32) < roundc(ZmaxCrop*1000, mold=1_int32))) then
            Zmax = real(undef_int, kind=sp)
        end if

        Zsoil = Zsoil + TheSoilLayer(layi)%Thickness
    end do

    if (Zmax < 0._sp) then
        call ZrAdjustedToRestrictiveLayers(ZmaxCrop, TheNrSoilLayers, &
                                           TheSoilLayer, Zmax)
    end if

    RootMaxInSoilProfile = real(Zmax, kind=sp)
end function RootMaxInSoilProfile


subroutine ZrAdjustedToRestrictiveLayers(ZrIN, TheNrSoilLayers, TheLayer, ZrOUT)
    real(sp), intent(in) :: ZrIN
    integer(int8), intent(in) :: TheNrSoilLayers
    type(SoilLayerIndividual), dimension(max_SoilLayers), intent(in) :: TheLayer
    real(sp), intent(inout) :: ZrOUT

    integer :: Layi
    real(sp) :: Zsoil, ZrAdj, ZrRemain, DeltaZ, ZrTest
    logical :: TheEnd

    ZrOUT = ZrIn

    ! initialize (layer 1)
    layi = 1
    Zsoil = TheLayer(layi)%Thickness
    ZrAdj = 0
    ZrRemain = ZrIN
    DeltaZ = Zsoil
    TheEnd = .false.

    ! check succesive layers
    do while (.not. TheEnd)
        ZrTest = ZrAdj + ZrRemain * (TheLayer(layi)%Penetrability/100._sp)

        if ((layi == TheNrSoilLayers) &
            .or. (TheLayer(layi)%Penetrability == 0) &
            .or. (roundc(ZrTest*10000, mold=1_int32) &
                 <= roundc(Zsoil*10000, mold=1_int32))) then
            ! no root expansion in layer
            TheEnd = .true.
            ZrOUT = ZrTest
        else
            ZrAdj = Zsoil
            ZrRemain = ZrRemain - DeltaZ/(TheLayer(layi)%Penetrability/100._sp)
            layi = layi + 1
            Zsoil = Zsoil + TheLayer(layi)%Thickness
            DeltaZ = TheLayer(layi)%Thickness
        end if
    end do
end subroutine ZrAdjustedToRestrictiveLayers


subroutine DeclareInitialCondAtFCandNoSalt()

    integer(int32) :: layeri, compi, celli, ind

    call SetSWCiniFile('(None)')
    call SetSWCiniFileFull(GetSWCiniFile()) ! no file
    call SetSWCiniDescription('Soil water profile at Field Capacity')
    call SetSimulation_IniSWC_AtDepths(.false.)
    call SetSimulation_IniSWC_NrLoc(GetSoil_NrSoilLayers())
    do layeri = 1, GetSoil_NrSoilLayers()
        call SetSimulation_IniSWC_Loc_i(layeri, GetSoilLayer_Thickness(layeri))
        call SetSimulation_IniSWC_VolProc_i(layeri, GetSoilLayer_FC(layeri))
        call SetSimulation_IniSWC_SaltECe_i(layeri, 0._sp)
    end do
    call SetSimulation_IniSWC_AtFC(.true.)
    do layeri = (GetSoil_NrSoilLayers()+1), max_No_compartments
        call SetSimulation_IniSWC_Loc_i(layeri, undef_double)
        call SetSimulation_IniSWC_VolProc_i(layeri, undef_double)
        call SetSimulation_IniSWC_SaltECe_i(layeri, undef_double)
    end do
    do compi = 1, GetNrCompartments()
        if (GetCompartment_Layer(compi) == 0) then
            ind = 1 ! LB: added an if statement to avoid having index=0
        else
            ind = GetCompartment_Layer(compi)
        end if
        do celli = 1, GetSoilLayer_SCP1(ind)
            ! salinity in cells
            call SetCompartment_Salt(compi, celli, 0.0_sp)
            call SetCompartment_Depo(compi, celli, 0.0_sp)
        end do
    end do
end subroutine DeclareInitialCondAtFCandNoSalt


subroutine set_layer_undef(LayerData)
    type(SoilLayerIndividual), intent(inout) :: LayerData

    integer(int32) :: i

    LayerData%Description = ''
    LayerData%Thickness = undef_double
    LayerData%SAT = undef_double
    LayerData%FC = undef_double
    LayerData%WP = undef_double
    LayerData%tau = undef_double
    LayerData%InfRate = undef_double
    LayerData%Penetrability = undef_int
    LayerData%GravelMass = undef_int
    LayerData%GravelVol = undef_int
    LayerData%Macro = undef_int
    LayerData%UL = undef_double
    LayerData%Dx = undef_double
    do i = 1, 11
        LayerData%SaltMobility(i) = undef_double  ! maximum 11 salt cells
    end do
    LayerData%SoilClass = undef_int
    LayerData%CRa = undef_int
    LayerData%CRb = undef_int
    LayerData%WaterContent = undef_double
end subroutine set_layer_undef


subroutine CropStressParametersSoilFertility(CropSResp, &
                    StressLevel, StressOUT)
    type(rep_Shapes), intent(in) :: CropSResp
    integer(int32), intent(in)    :: StressLevel
    type(rep_EffectStress), intent(inout) :: StressOUT

    real(sp) :: Ksi, pULActual, pLLActual

    pLLActual = 1._sp

    ! decline canopy growth coefficient (CGC)
    pULActual = 0._sp
    Ksi = KsAny(StressLevel/100._sp, pULActual, pLLActual, CropSResp%ShapeCGC)
    StressOUT%RedCGC = roundc((1._sp-Ksi)*100._sp, mold=1_int8)
    ! decline maximum canopy cover (CCx)
    pULActual = 0._sp
    Ksi = KsAny(StressLevel/100._sp, pULActual, pLLActual, CropSResp%ShapeCCX)
    StressOUT%RedCCX = roundc((1._sp-Ksi)*100._sp, mold=1_int8)
    ! decline crop water productivity (WP)
    pULActual = 0._sp
    Ksi = KsAny(StressLevel/100._sp, pULActual, pLLActual, CropSResp%ShapeWP)
    StressOUT%RedWP = roundc((1._sp-Ksi)*100._sp, mold=1_int8)
    ! decline Canopy Cover (CDecline)
    pULActual = 0._sp
    Ksi = KsAny(StressLevel/100._sp, pULActual, pLLActual, CropSResp%ShapeCDecline)
    StressOUT%CDecline = 1._sp - Ksi
    ! inducing stomatal closure (KsSto) not applicable
    Ksi = 1._sp
    StressOUT%RedKsSto = roundc((1._sp-Ksi)*100._sp, mold=1_int8)
end subroutine CropStressParametersSoilFertility


real(sp) function TimeRootFunction(t, ShapeFactor, tmax, t0)
    real(sp), intent(in) :: t
    integer(int8), intent(in) :: ShapeFactor
    real(sp), intent(in) :: tmax
    real(sp), intent(in) :: t0

    TimeRootFunction = exp((10._sp / ShapeFactor) * log((t-t0) / (tmax-t0)))
end function TimeRootFunction


real(sp) function TimeToReachZroot(Zi, Zo, Zx, ShapeRootDeepening, Lo, LZxAdj)
    real(sp), intent(in) :: Zi
    real(sp), intent(in) :: Zo
    real(sp), intent(in) :: Zx
    integer(int8), intent(in) :: ShapeRootDeepening
    integer(int32), intent(in) :: Lo
    integer(int32), intent(in) :: LZxAdj

    real(sp) :: ti, T1

    ti = real(undef_int, kind=sp)

    if (roundc(Zi*100, mold=1_int32) >= roundc(Zx*100, mold=1_int32)) then
        ti = real(LZxAdj, kind=sp)
    else
        if (((Zo+0.0001_sp) < Zx) .and. (LZxAdj > Lo/2._sp) .and. (LZxAdj > 0) &
            .and. (ShapeRootDeepening > 0)) then
            T1 = exp((ShapeRootDeepening/10._sp) * log((Zi-Zo) / (Zx-Zo)))
            ti = T1 * (LZxAdj - Lo/2._sp) + Lo/2._sp
        end if
    end if

    TimeToReachZroot = ti
end function TimeToReachZroot


real(sp) function CanopyCoverNoStressSF(DAP, L0, L123, &
       LMaturity, GDDL0, GDDL123, GDDLMaturity, CCo, CCx,&
       CGC, CDC, GDDCGC, GDDCDC, SumGDD, TypeDays, SFRedCGC, SFRedCCx)
    integer(int32), intent(in) :: DAP
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: L123
    integer(int32), intent(in) :: LMaturity
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDL123
    integer(int32), intent(in) :: GDDLMaturity
    real(sp), intent(in) :: CCo
    real(sp), intent(in) :: CCx
    real(sp), intent(in) :: CGC
    real(sp), intent(in) :: CDC
    real(sp), intent(in) :: GDDCGC
    real(sp), intent(in) :: GDDCDC
    real(sp), intent(in) :: SumGDD
    integer(intEnum), intent(in) :: TypeDays
    integer(int8), intent(in) :: SFRedCGC
    integer(int8), intent(in) :: SFRedCCx

    select case (TypeDays)
    case (modeCycle_GDDays)
        CanopyCoverNoStressSF = CanopyCoverNoStressGDDaysSF(GDDL0, GDDL123,&
            GDDLMaturity, SumGDD, CCo, CCx, GDDCGC, GDDCDC, SFRedCGC, SFRedCCx)
    case default
        CanopyCoverNoStressSF = CanopyCoverNoStressDaysSF(DAP, L0, L123,&
            LMaturity, CCo, CCx, CGC, CDC, SFRedCGC, SFRedCCx)
    end select


    contains


    real(sp) function CanopyCoverNoStressDaysSF(DAP, L0, L123,&
           LMaturity, CCo, CCx, CGC, CDC, SFRedCGC, SFRedCCx)
        integer(int32), intent(in) :: DAP
        integer(int32), intent(in) :: L0
        integer(int32), intent(in) :: L123
        integer(int32), intent(in) :: LMaturity
        real(sp), intent(in) :: CCo
        real(sp), intent(in) :: CCx
        real(sp), intent(in) :: CGC
        real(sp), intent(in) :: CDC
        integer(int8), intent(in) :: SFRedCGC
        integer(int8), intent(in) :: SFRedCCx

        real(sp) :: CC, CCxAdj, CDCadj
        integer(int32) :: t

        ! CanopyCoverNoStressDaysSF
        CC = 0.0_sp
        t = DAP - GetSimulation_DelayedDays()
        ! CC refers to canopy cover at the end of the day

        if ((t >= 1) .and. (t <= LMaturity) .and. (CCo > epsilon(1._sp))) then
            if (t <= L0) then ! before germination or recovering of transplant
                CC = 0._sp
            else
                if (t < L123) then ! Canopy development and Mid-season stage
                    CC = CCatTime((t-L0), CCo, ((1._sp-SFRedCGC/100._sp)*CGC),&
                                  ((1._sp-SFRedCCx/100._sp)*CCx))
                else
                   ! Late-season stage  (t <= LMaturity)
                    if (CCx < 0.001) then
                        CC = 0._sp
                    else
                        CCxAdj = CCatTime((L123-L0), CCo, &
                                  ((1._sp-SFRedCGC/100._sp)*CGC),&
                                  ((1._sp-SFRedCCx/100._sp)*CCx))
                        CDCadj = CDC*(CCxAdj+2.29_sp)/(CCx+2.29_sp)
                        if (CCxAdj < 0.001) then
                            CC = 0._sp
                        else
                            CC = CCxAdj * (1._sp - 0.05_sp *&
                                 (exp((t-L123)*3.33_sp*CDCAdj/(CCxAdj+2.29_sp))-1._sp))
                        end if
                    end if
                end if
            end if
        end if
        if (CC > 1._sp) then
            CC = 1._sp
        end if
        if (CC < epsilon(1._sp)) then
            CC = 0._sp
        end if
        CanopyCoverNoStressDaysSF = CC
    end function CanopyCoverNoStressDaysSF
end function CanopyCoverNoStressSF


real(sp) function CCiNoWaterStressSF(Dayi, L0, L12SF, L123, L1234, GDDL0,&
    GDDL12SF, GDDL123, GDDL1234, CCo, CCx, CGC, GDDCGC, CDC, GDDCDC, SumGDD,&
    RatDGDD, SFRedCGC, SFRedCCx, SFCDecline, TheModeCycle)

    integer(int32), intent(in) :: Dayi
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: L12SF
    integer(int32), intent(in) :: L123
    integer(int32), intent(in) :: L1234
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDL12SF
    integer(int32), intent(in) :: GDDL123
    integer(int32), intent(in) :: GDDL1234
    real(sp), intent(in) :: CCo
    real(sp), intent(in) :: CCx
    real(sp), intent(in) :: CGC
    real(sp), intent(in) :: GDDCGC
    real(sp), intent(in) :: CDC
    real(sp), intent(in) :: GDDCDC
    real(sp), intent(in) :: SumGDD
    real(sp), intent(in) :: RatDGDD
    integer(int8), intent(in) :: SFRedCGC
    integer(int8), intent(in) :: SFRedCCx
    real(sp), intent(in) :: SFCDecline
    integer(intEnum), intent(in) :: TheModeCycle

    real(sp) :: CCi, CCibis, CCxAdj, CDCadj, GDDCDCadj

    ! Calculate CCi
    CCi = CanopyCoverNoStressSF(Dayi, L0, L123, L1234, GDDL0, GDDL123,&
                                GDDL1234, CCo, CCx, CGC, CDC, GDDCGC,&
                                GDDCDC, SumGDD, TheModeCycle, SFRedCGC,&
                                SFRedCCX)

    ! Consider CDecline for limited soil fertiltiy
    if ((Dayi > L12SF) .and. (SFCDecline > ac_zero_threshold) .and. (L12SF < L123)) then
        if (Dayi < L123) then
            if (TheModeCycle == modeCycle_CalendarDays) then
                CCi = CCi - (SFCDecline/100.0_sp)&
                            * exp(2.0_sp*log(real(Dayi-L12SF, kind=sp)))&
                            / real(L123-L12SF, kind=sp)
            else
                if ((SumGDD > GDDL12SF) .and. (GDDL123 > GDDL12SF)) then
                    CCi = CCi - (RatDGDD*SFCDecline/100.0_sp)&
                                * exp(2.0_sp*log(SumGDD-GDDL12SF))&
                                / real(GDDL123-GDDL12SF, kind=sp)
                end if
            end if
            if (CCi < 0.0_sp) then
                CCi = 0.0_sp
            end if
        else
            if (TheModeCycle == modeCycle_CalendarDays) then
                CCi = CCatTime((L123-L0), CCo, (CGC*(1.0_sp-SFRedCGC/100.0_sp)),&
                               ((1.0_sp-SFRedCCX/100.0_sp)*CCx))
                ! CCibis is CC in late season when Canopy decline continues
                CCibis = CCi  - (SFCDecline/100.0_sp)&
                                * (exp(2.0_sp*log(real(Dayi-L12SF, kind=sp)))&
                                   / real(L123-L12SF, kind=sp))
                if (CCibis < 0.0_sp) then
                    CCi = 0.0_sp
                else
                    CCi = CCi  - ((SFCDecline/100.0_sp) * (L123-L12SF))
                end if
                if (CCi < 0.001_sp) then
                    CCi = 0.0_sp
                else
                    ! is CCx at start of late season, adjusted for canopy
                    ! decline with soil fertility stress
                    CCxAdj = CCi
                    CDCadj = CDC * (CCxAdj + 2.29_sp)/(CCx + 2.29_sp)
                    if (Dayi < (L123 + LengthCanopyDecline(CCxAdj, CDCadj))) then
                        CCi = CCxAdj * (1.0_sp &
                                        - 0.05_sp*(exp((Dayi-L123)*3.33_sp*CDCadj&
                                                       /(CCxAdj+2.29_sp))&
                                                   -1.0_sp))
                        if (CCibis < CCi) then
                            CCi = CCibis ! accept smallest Canopy Cover
                        end if
                    else
                        CCi = 0.0_sp
                    end if
                end if
            else
                CCi = CCatTime((GDDL123-GDDL0), CCo,&
                               (GDDCGC*(1.0_sp-SFRedCGC/100.0_sp)),&
                               ((1.0_sp-SFRedCCX/100.0_sp)*CCx))
                ! CCibis is CC in late season when Canopy decline continues
                if ((SumGDD > GDDL12SF) .and. (GDDL123 > GDDL12SF)) then
                    CCibis = CCi  - (RatDGDD*SFCDecline/100.0_sp)&
                                    * (exp(2.0_sp*log(SumGDD-GDDL12SF))&
                                      /real(GDDL123-GDDL12SF, kind=sp))
                else
                    CCibis = CCi
                end if
                if (CCibis < 0.0_sp) then
                    CCi = 0.0_sp
                else
                    CCi = CCi - ((RatDGDD*SFCDecline/100.0_sp) * (GDDL123-GDDL12SF))
                end if
                if (CCi < 0.001_sp) then
                    CCi = 0.0_sp
                else
                    ! is CCx at start of late season, adjusted for canopy
                    ! decline with soil fertility stress
                    CCxAdj = CCi
                    GDDCDCadj = GDDCDC * (CCxAdj + 2.29_sp)/(CCx + 2.29_sp)
                    if (SumGDD < (GDDL123 + LengthCanopyDecline(CCxAdj, GDDCDCadj))) then
                        CCi = CCxAdj * (1.0_sp&
                                        - 0.05_sp*(exp((SumGDD-GDDL123)*3.33_sp&
                                                        *GDDCDCadj/(CCxAdj+2.29_sp))&
                                                    -1.0_sp))
                        if (CCibis < CCi) then
                            CCi = CCibis ! accept smallest Canopy Cover
                        end if
                    else
                        CCi = 0.0_sp
                    end if
                end if
            end if
            if (CCi < 0.0_sp) then
                CCi = 0.0_sp
            end if
        end if
    end if
    CCiNoWaterStressSF = CCi
end function CCiNoWaterStressSF


real(sp) function FromGravelMassToGravelVolume(PorosityPercent,&
                                               GravelMassPercent)
    real(sp), intent(in)      :: PorosityPercent
    integer(int8), intent(in) :: GravelMassPercent

    real(sp), parameter ::  MineralBD = 2.65 !! Mg/m3
    real(sp) :: MatrixBD, SoilBD

    if (GravelMassPercent > 0) then
        MatrixBD = MineralBD * (1._sp - PorosityPercent/100._sp)
        SoilBD = 100._sp/(GravelMassPercent/MineralBD + &
                          (100._sp-GravelMassPercent)/MatrixBD)
        FromGravelMassToGravelVolume = GravelMassPercent * (SoilBD/MineralBD)
   else
       FromGravelMassToGravelVolume = 0.0_sp
   end if
end function FromGravelMassToGravelVolume


subroutine CheckForWaterTableInProfile(DepthGWTmeter, ProfileComp, WaterTableInProfile)
    real(sp), intent(in) :: DepthGWTmeter
    type(CompartmentIndividual), dimension(max_No_compartments), &
                                                    intent(in) :: ProfileComp
    logical, intent(inout) :: WaterTableInProfile

    real(sp) :: Ztot, Zi
    integer(int32) :: compi

    WaterTableInProfile = .false.
    Ztot = 0._sp
    compi = 0

    if (DepthGWTmeter >= epsilon(0._sp)) then
        ! groundwater table is present
        do while ((.not. WaterTableInProfile) .and. (compi < getNrCompartments()))
            compi = compi + 1
            Ztot = Ztot + ProfileComp(compi)%Thickness
            Zi = Ztot - ProfileComp(compi)%Thickness/2._sp
            if (Zi >= DepthGWTmeter) then
                WaterTableInProfile = .true.
            end if
        end do
    end if
end subroutine CheckForWaterTableInProfile


real(sp) function GetWeedRC(TheDay, GDDayi, fCCx, TempWeedRCinput, TempWeedAdj,&
                            TempWeedDeltaRC, L12SF, TempL123, GDDL12SF, &
                            TempGDDL123, TheModeCycle)
    integer(int32), intent(in) :: TheDay
    real(sp), intent(in) :: GDDayi
    real(sp), intent(in) :: fCCx
    integer(int8), intent(in) :: TempWeedRCinput
    integer(int8), intent(in) :: TempWeedAdj
    integer(int32), intent(inout) :: TempWeedDeltaRC
    integer(int32), intent(in) :: L12SF
    integer(int32), intent(in) :: TempL123
    integer(int32), intent(in) :: GDDL12SF
    integer(int32), intent(in) :: TempGDDL123
    integer(intEnum), intent(in) :: TheModeCycle

    real(sp) :: WeedRCDayCalc

    WeedRCDayCalc = TempWeedRCinput

    if ((TempWeedRCinput > 0) .and. (TempWeedDeltaRC /= 0)) then
        ! daily RC when increase/decline of RC in season (i.e. TempWeedDeltaRC <> 0)
        ! adjust the slope of increase/decline of RC in case of self-thinning (i.e. fCCx < 1)
        if ((TempWeedDeltaRC /= 0) .and. (fCCx < 0.999_sp)) then
            ! only when self-thinning and there is increase/decline of RC
            if (fCCx < 0.005_sp) then
                TempWeedDeltaRC = 0
            else
                TempWeedDeltaRC = roundc(TempWeedDeltaRC *&
                     exp(log(fCCx) * (1+TempWeedAdj/100._sp)), mold=1_int32)
            end if
        end if

        ! calculate WeedRCDay by considering (adjusted) decline/increase of RC
        if (TheModeCycle == modeCycle_CalendarDays) then
            if (TheDay > L12SF) then
                if (TheDay >= TempL123) then
                    WeedRCDayCalc = TempWeedRCinput * (1 + &
                                        TempWeedDeltaRC/100._sp)
                else
                    WeedRCDayCalc = TempWeedRCinput * (1 + &
                                        (TempWeedDeltaRC/100._sp) &
                                         * (TheDay-L12SF) / real(TempL123-L12SF,kind=sp))
                end if
            end if
        else
            if (GDDayi > GDDL12SF) then
                if (GDDayi > TempGDDL123) then
                    WeedRCDayCalc = TempWeedRCinput * (1 + &
                                        TempWeedDeltaRC/100._sp)
                else
                    WeedRCDayCalc = TempWeedRCinput * (1 + &
                                        (TempWeedDeltaRC/100._sp) &
                                         * (GDDayi-GDDL12SF) &
                                         / real(TempGDDL123-GDDL12SF,kind=sp))
                end if
            end if
        end if

        ! fine-tuning for over- or undershooting in case of self-thinning
        if (fCCx < 0.999_sp) then
            ! only for self-thinning
            if ((fCCx < 1) .and. (fCCx > 0) .and. (WeedRCDayCalc > 98)) then
                WeedRCDayCalc = 98._sp
            end if
            if (WeedRCDayCalc < 0) then
                WeedRCDayCalc = 0._sp
            end if
            if (fCCx <= 0) then
                WeedRCDayCalc = 100._sp
            end if
        end if
    end if

    GetWeedRC = WeedRCDayCalc
end function GetWeedRC


subroutine DetermineLengthGrowthStages(CCoVal, CCxVal, CDCVal, L0, TotalLength, &
                                        CGCgiven, TheDaysToCCini, ThePlanting, &
                                         Length123, StLength, Length12, CGCVal)
    real(sp), intent(in) :: CCoVal
    real(sp), intent(in) :: CCxVal
    real(sp), intent(in) :: CDCVal
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: TotalLength
    logical, intent(in) :: CGCgiven
    integer(int32), intent(in) :: TheDaysToCCini
    integer(intEnum), intent(in) :: ThePlanting
    integer(int32), intent(inout) :: Length123
    integer(int32), dimension(4), intent(inout) :: StLength
    integer(int32), intent(inout) :: Length12
    real(sp), intent(inout) :: CGCVal

    real(sp) :: CCxVal_scaled
    real(sp) :: CCToReach
    integer(int32) :: L12Adj

    if (Length123 < Length12) then
        Length123 = Length12
    end if

    ! 1. Initial and 2. Crop Development stage
    ! CGC is given and Length12 is already adjusted to it
    ! OR Length12 is given and CGC has to be determined
    if ((CCoVal >= CCxVal) .or. (Length12 <= L0)) then
        Length12 = 0
        StLength(1) = 0
        StLength(2) = 0
        CGCVal = Undef_int
    else
        if (.not. CGCgiven) then ! Length12 is given and CGC has to be determined
            CGCVal = real((log((0.25_sp*CCxVal/CCoVal)/(1._sp-0.98_sp))&
                           /real(Length12-L0, kind=sp)), kind=sp)
            ! Check if CGC < maximum value (0.40) and adjust Length12 if required
            if (CGCVal > 0.40_sp) then
                CGCVal = 0.40_sp
                CCxVal_scaled = 0.98_sp*CCxVal
                Length12 = DaysToReachCCwithGivenCGC(CCxVal_scaled , CCoVal, &
                                                             CCxVal, CGCVal, L0)
                if (Length123 < Length12) then
                    Length123 = Length12
                end if
            end if
        end if
        ! find StLength[1]
        CCToReach = 0.10_sp
        StLength(1) = DaysToReachCCwithGivenCGC(CCToReach, CCoVal, CCxVal, &
                                                                    CGCVal, L0)
        ! find StLength[2]
        StLength(2) = Length12 - StLength(1)
    end if
    L12Adj = Length12

    ! adjust Initial and Crop Development stage, in case crop starts as regrowth
    if (ThePlanting == plant_regrowth) then
        if (TheDaystoCCini == undef_int) then
            ! maximum canopy cover is already reached at start season
            L12Adj = 0
            StLength(1) = 0
            StLength(2) = 0
        else
            if (TheDaystoCCini == 0) then
                ! start at germination
                L12Adj = Length12 - L0
                StLength(1) = StLength(1) - L0
            else
                ! start after germination
                L12Adj = Length12 - (L0 + TheDaysToCCini)
                StLength(1) = StLength(1) - (L0 + TheDaysToCCini)
            end if
            if (StLength(1) < 0) then
                StLength(1) = 0
            end if
            StLength(2) = L12Adj - StLength(1)
        end if
    end if

    ! 3. Mid season stage
    StLength(3) = Length123 - L12Adj

    ! 4. Late season stage
    StLength(4) = LengthCanopyDecline(CCxVal, CDCVal)

    ! final adjustment
    if (StLength(1) > TotalLength) then
        StLength(1) = TotalLength
        StLength(2) = 0
        StLength(3) = 0
        StLength(4) = 0
    else
        if ((StLength(1)+StLength(2)) > TotalLength) then
            StLength(2) = TotalLength - StLength(1)
            StLength(3) = 0
            StLength(4) = 0
        else
            if ((StLength(1)+StLength(2)+StLength(3)) > TotalLength) then
                StLength(3) = TotalLength - StLength(1) - StLength(2)
                StLength(4) = 0
            elseif ((StLength(1)+StLength(2)+StLength(3)+StLength(4)) > &
                                                              TotalLength) then
                StLength(4) = TotalLength - StLength(1) - StLength(2) - StLength(3)
            end if
        end if
    end if
end subroutine DetermineLengthGrowthStages


integer(int32) function TimeToCCini(ThePlantingType, TheCropPlantingDens, &
                          TheSizeSeedling, TheSizePlant, TheCropCCx, TheCropCGC)
    integer(intEnum), intent(in) :: ThePlantingType
    integer(int32), intent(in) :: TheCropPlantingDens
    real(sp), intent(in) :: TheSizeSeedling
    real(sp), intent(in) :: TheSizePlant
    real(sp), intent(in) :: TheCropCCx
    real(sp), intent(in) :: TheCropCGC

    integer(int32) :: ElapsedTime
    real(sp) :: TheCropCCo
    real(sp) :: TheCropCCini

    if ((ThePlantingType == plant_seed) .or. (ThePlantingType == plant_transplant) &
                                   .or. (TheSizeSeedling >= TheSizePlant)) then
        ElapsedTime = 0
    else
        TheCropCCo = (TheCropPlantingDens/10000._sp) * (TheSizeSeedling/10000._sp)
        TheCropCCini = (TheCropPlantingDens/10000._sp) * (TheSizePlant/10000._sp)
        if (TheCropCCini >= (0.98_sp*TheCropCCx)) then
            ElapsedTime = undef_int
        else
            ElapsedTime = DaysToReachCCwithGivenCGC(TheCropCCini, TheCropCCo, &
                                                    TheCropCCx, TheCropCGC, 0)
        end if
    end if
    TimeToCCini = ElapsedTime
end function TimeToCCini


real(sp) function MultiplierCCxSelfThinning(Yeari, Yearx, ShapeFactor)
    integer(int32), intent(in) :: Yeari
    integer(int32), intent(in) :: Yearx
    real(sp), intent(in) :: ShapeFactor

    real(sp) :: fCCx, Year0

    fCCx = 1
    if ((Yeari >= 2) .and. (Yearx >= 2) .and. &
        (roundc(100._sp*ShapeFactor, mold=1_int32) /= 0)) then
        Year0 = 1._sp + (Yearx-1._sp) * exp(ShapeFactor*log(10._sp))
        if (Yeari >= Year0) then
            fCCx = 0
        else
            fCCx = 0.9_sp + 0.1_sp * (1._sp - exp((1._sp/ShapeFactor) &
                        *log((Yeari-1._sp)/(Yearx-1._sp))))
        end if
        if (fCCx < 0) then
            fCCx = 0
        end if
    end if
    MultiplierCCxSelfThinning = fCCx
end function MultiplierCCxSelfThinning


integer(int32) function DaysToReachCCwithGivenCGC(CCToReach, CCoVal, &
                                                        CCxVal, CGCVal, L0)
    real(sp), intent(in) :: CCToReach
    real(sp), intent(in) :: CCoVal
    real(sp), intent(in) :: CCxVal
    real(sp), intent(in) :: CGCVal
    integer(int32), intent(in) :: L0

    real(sp) :: L
    real(sp) :: CCToReach_local

    CCToReach_local = CCToReach
    if ((CCoVal > CCToReach_local) .or. (CCoVal >= CCxVal)) then
        L = 0
    else
        if (CCToReach_local > (0.98_sp*CCxVal)) then
            CCToReach_local = 0.98_sp*CCxVal
        end if
        if (CCToReach_local <= CCxVal/2._sp) then
            L = log(CCToReach_local/CCoVal)/CGCVal
        else
            L = log((0.25_sp*CCxVal*CCxVal/CCoVal)/(CCxVal-CCToReach_local))&
                                                                        /CGCVal
        end if

    end if
    DaysToReachCCwithGivenCGC = L0 + roundc(L, mold=1_int32)
end function DaysToReachCCwithGivenCGC


integer(int32) function LengthCanopyDecline(CCx, CDC)
    real(sp), intent(in) :: CCx
    real(sp), intent(in) :: CDC

    integer(int32) :: ND

    ND = 0
    if (CCx > 0) then
        if (CDC <= epsilon(1._sp)) then
            ND = undef_int
        else
            ND = roundc((((CCx+2.29_sp)/(CDC*3.33_sp))* &
                         log(1._sp + 1._sp/0.05_sp) + 0.50_sp), mold=1_int32)
                         ! + 0.50 to guarantee that CC is zero
        end if

    end if
    LengthCanopyDecline = ND
end function LengthCanopyDecline


real(sp) function HarvestIndexGrowthCoefficient(HImax, dHIdt)
    real(sp), intent(in) :: HImax
    real(sp), intent(in) :: dHIdt

    real(sp) :: HIo, HIvar, HIGC, t

    HIo = 1

    if (HImax > HIo) then
        t = HImax/dHIdt
        HIGC = 0.001_sp
        HIGC = HIGC + 0.001_sp
        HIvar = (HIo*HImax)/(HIo+(HImax-HIo)*exp(-HIGC*t))
        do while (HIvar <= (0.98_sp*HImax))
            HIGC = HIGC + 0.001_sp
            HIvar = (HIo*HImax)/(HIo+(HImax-HIo)*exp(-HIGC*t))
        end do

        if (HIvar >= HImax) then
            HIGC = HIGC - 0.001_sp
        end if
    else
        HIGC = undef_int
    end if
    HarvestIndexGrowthCoefficient = HIGC

end function HarvestIndexGrowthCoefficient


real(sp) function TauFromKsat(Ksat)
    real(sp), intent(in) :: Ksat

    integer(int32) :: TauTemp

    if (abs(Ksat) < epsilon(1._sp)) then
        TauFromKsat = 0
    else
        TauTemp = roundc(100.0_sp*0.0866_sp*exp(0.35_sp*log(Ksat)), mold=1_int32)
        if (TauTemp < 0) then
            TauTemp = 0
        end if
        if (TauTemp > 100) then
            TauTemp = 100
        end if
        TauFromKsat = TauTemp/100.0_sp
    end if
end function TauFromKsat


integer(int8) function NumberSoilClass(SatvolPro, FCvolPro, PWPvolPro, Ksatmm)
    real(sp), intent(in) :: SatvolPro
    real(sp), intent(in) :: FCvolPro
    real(sp), intent(in) :: PWPvolPro
    real(sp), intent(in) :: Ksatmm

    if (SATvolPro <= 55.0_sp) then
        if (PWPvolPro >= 20.0_sp) then
            if ((SATvolPro >= 49.0_sp) .and. (FCvolPro >= 40.0_sp)) then
                NumberSoilClass = 4  ! silty clayey soils
            else
                NumberSoilClass = 3  ! sandy clayey soils
            end if
        else
            if (FCvolPro < 23.0_sp) then
                NumberSoilClass = 1 ! sandy soils
            else
                if ((PWPvolPro > 16.0_sp) .and. (Ksatmm < 100.0_sp)) then
                    NumberSoilClass = 3 ! sandy clayey soils
                else
                    if ((PWPvolPro < 6.0_sp) .and. (FCvolPro < 28.0_sp) &
                        .and. (Ksatmm >750.0_sp)) then
                        NumberSoilClass = 1 ! sandy soils
                    else
                        NumberSoilClass = 2  ! loamy soils
                    end if
                end if
            end if
        end if
    else
        NumberSoilClass = 4 ! silty clayey soils
    end if
end function NumberSoilClass


subroutine DeriveSmaxTopBottom(SxTopQ, SxBotQ, SxTop, SxBot)
    real(sp), intent(in) :: SxTopQ
    real(sp), intent(in) :: SxBotQ
    real(sp), intent(inout) :: SxTop
    real(sp), intent(inout) :: SxBot

    real(sp) :: x, V1, V2, V11, V22

    V1 = SxTopQ
    V2 = SxBotQ
    if (abs(V1 - V2) < 1e-12_sp) then
        SxTop = V1
        SxBot = V2
    else
        if (SxTopQ < SxBotQ) then
            V1 = SxBotQ
            V2 = SxTopQ
        end if
        x = 3.0_sp * V2/(V1-V2)
        if (x < 0.5_sp) then
            V11 = (4.0_sp/3.5_sp) * V1
            V22 = 0.0_sp
        else
            V11 = (x + 3.5_sp) * V1/(x+3.0_sp)
            V22 = (x - 0.5_sp) * V2/x
        end if
        if (SxTopQ > SxBotQ) then
            SxTop = V11
            SxBot = V22
        else
            SxTop = V22
            SxBot = V11
        end if
    end if
end subroutine DeriveSmaxTopBottom


real(sp) function KsTemperature(T0, T1, Tin)
    real(sp), intent(in) :: T0
    real(sp), intent(in) :: T1
    real(sp), intent(in) :: Tin

    real(sp) :: M
    integer(int8) :: a

    M = 1._sp ! no correction applied (TO and/or T1 is undefined, or T0=T1)
    if (((roundc(T0, mold=1_int32) /= undef_int) .and. &
         (roundc(T1, mold=1_int32) /= undef_int)) .and. abs(T0-T1)> eps) then
        if (T0 < T1) then
            a =  1  ! cold stress
        else
            a = -1 ! heat stress
        end if
        if ((a*Tin > a*T0) .and. (a*Tin < a*T1)) then
           ! within range for correction
            M = GetKs(T0, T1, Tin)
            if (M < 0) then
                M = 0._sp
            end if
            if (M > 1) then
                M = 1._sp
            end if
        else
            if (a*Tin <= a*T0) then
                M = 0._sp
            end if
            if (a*Tin >= a*T1) then
                M = 1._sp
            end if
        end if
    end if
    KsTemperature = M


    contains


    real(sp) function GetKs(T0, T1, Tin)
        real(sp), intent(in) :: T0
        real(sp), intent(in) :: T1
        real(sp), intent(in) :: Tin

        real(sp), parameter  :: Mo = 0.02_sp
        real(sp), parameter  :: Mx = 1.0_sp

        real(sp) :: MRate, Ksi, Trel

        Trel = (Tin-T0)/(T1-T0)
        ! derive rate of increase (MRate)
        MRate = (-1._sp)*(log((Mo*Mx-0.98_sp*Mo)/(0.98_sp*(Mx-Mo))))
        ! get Ks from logistic equation
        Ksi = (Mo*Mx)/(Mo+(Mx-Mo)*exp(-MRate*Trel))
        ! adjust for Mo
        Ksi = Ksi - Mo * (1._sp - Trel)
        GetKs = Ksi
     end function GetKs
end function KsTemperature


real(sp) function KsSalinity(SalinityResponsConsidered, &
                ECeN, ECeX, ECeVAR, KsShapeSalinity)
    logical, intent(in) :: SalinityResponsConsidered
    integer(int8), intent(in) :: ECeN
    integer(int8), intent(in) :: ECeX
    real(sp), intent(in) :: ECeVAR
    real(sp), intent(in) :: KsShapeSalinity

    real(sp) :: M, tmp_var

    M = 1._sp ! no correction applied
    if (SalinityResponsConsidered) then
        if ((ECeVAR > ECeN) .and. (ECeVar < ECeX)) then
            ! within range for correction
            if ((roundc(KsShapeSalinity*10._sp, mold=1_int32) /= 0) .and. &
                (roundc(KsShapeSalinity*10._sp, mold=1_int32) /= 990)) then
                tmp_var = real(ECeN, kind=sp)
                M = KsAny(ECeVar, tmp_var, real(ECeX, kind=sp), KsShapeSalinity)
                ! convex or concave
            else
                if (roundc(KsShapeSalinity*10._sp, mold=1_int32) == 0) then
                    M = 1._sp - (ECeVAR-ECeN)/(ECeX-ECeN)
                    ! linear (KsShapeSalinity = 0)
                else
                    M = KsTemperature(real(ECeX, kind=sp), real(ECeN, kind=sp), ECeVAR)
                    ! logistic equation (KsShapeSalinity = 99)
                end if
            end if
        else
            if (ECeVAR <= ECeN) then
                M = 1._sp  ! no salinity stress
            end if
            if (ECeVar >= ECeX) then
                M = 0._sp  ! full salinity stress
            end if
        end if
    end if
    if (M > 1) then
        M = 1._sp
    end if
    if (M < 0) then
        M = 0._sp
    end if
    KsSalinity = M
end function KsSalinity


subroutine TimeToMaxCanopySF(CCo, CGC, CCx, L0, L12, L123, LToFlor, LFlor, DeterminantCrop, L12SF, RedCGC, RedCCx, ClassSF)
    real(sp), intent(in) :: CCo
    real(sp), intent(in) :: CGC
    real(sp), intent(in) :: CCx
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: L12
    integer(int32), intent(in) :: L123
    integer(int32), intent(in) :: LToFlor
    integer(int32), intent(in) :: LFlor
    logical, intent(in) :: DeterminantCrop
    integer(int32), intent(inout) :: L12SF
    integer(int8), intent(inout) :: RedCGC
    integer(int8), intent(inout) :: RedCCx
    integer(int32), intent(inout) :: ClassSF

    real(sp) :: CCToReach
    integer(int32) :: L12SFmax

    if ((ClassSF == 0) .or. ((RedCCx == 0) .and. (RedCGC == 0))) then
        L12SF = L12
    else
        CCToReach = 0.98_sp*(1-RedCCX/100._sp)*CCx
        L12SF = DaysToReachCCwithGivenCGC(CCToReach, CCo, ((1-RedCCX/100._sp)*CCx), (CGC*(1-(RedCGC)/100._sp)), L0)
        ! determine L12SFmax
        if (DeterminantCrop) then
            L12SFmax = LToFlor + roundc(LFlor/2._sp, mold=1_int32)
        else
            L12SFmax = L123
        end if
        ! check for L12SFmax
        if (L12SF > L12SFmax) then
            ! full canopy cannot be reached in potential period for vegetative growth
            ! ClassSF := undef_int; ! switch to user defined soil fertility
            ! 1. increase CGC(soil fertility)
            do while ((L12SF > L12SFmax) .and. (RedCGC > 0))
                RedCGC = RedCGC - 1_int8
                L12SF = DaysToReachCCwithGivenCGC(CCToReach, CCo, ((1-RedCCX/100._sp)*CCx), (CGC*(1-(RedCGC)/100._sp)), L0)
            end do
            ! 2. if not sufficient decrease CCx(soil fertility)
            do while ((L12SF > L12SFmax) .and. ( ((1-RedCCX/100._sp)*CCx) > 0.10_sp) .and. (RedCCx <= 50))
                RedCCx = RedCCx + 1_int8
                CCToReach = 0.98_sp*(1-RedCCX/100._sp)*CCx
                L12SF = DaysToReachCCwithGivenCGC(CCToReach, CCo, ((1-RedCCX/100._sp)*CCx), (CGC*(1-(RedCGC)/100._sp)), L0)
            end do
        end if
    end if
end subroutine TimeToMaxCanopySF


real(sp) function SoilEvaporationReductionCoefficient(Wrel, Edecline)
    real(sp), intent(in) :: Wrel
    real(sp), intent(in) :: Edecline

    if (Wrel <= 0.00001_sp) then
        SoilEvaporationReductionCoefficient = 0.0_sp
    else
        if (Wrel >= 0.99999_sp) then
            SoilEvaporationReductionCoefficient = 1.0_sp
        else
            SoilEvaporationReductionCoefficient =&
                (exp(Edecline*Wrel) - 1.0_sp)/(exp(Edecline) - 1.0_sp)
        end if
    end if
end function SoilEvaporationReductionCoefficient


real(sp) function MaxCRatDepth(ParamCRa, ParamCRb, Ksat, Zi, DepthGWT)
    real(sp), intent(in) :: ParamCRa
    real(sp), intent(in) :: ParamCRb
    real(sp), intent(in) :: Ksat
    real(sp), intent(in) :: Zi
    real(sp), intent(in) :: DepthGWT

    real(sp) :: CRmax

    CRmax = 0._sp
    if ((Ksat > 0._sp) .and. (DepthGWT > 0._sp) .and. ((DepthGWT-Zi) < 4._sp)) then
        if (Zi >= DepthGWT) then
            CRmax = 99._sp
        else
            CRmax = exp((log(DepthGWT - Zi) - ParamCRb)/ParamCRa)
            if (CRmax > 99._sp) then
                CRmax = 99._sp
            end if
        end if
    end if
    MaxCRatDepth = CRmax
end function MaxCRatDepth


real(sp) function CCmultiplierWeed(ProcentWeedCover, CCxCrop, FshapeWeed)
    integer(int8), intent(in) :: ProcentWeedCover
    real(sp), intent(in) :: CCxCrop
    real(sp), intent(in) :: FshapeWeed

    real(sp) :: fWeed

    if ((ProcentWeedCover > 0) .and. (CCxCrop < 0.9999_sp) .and. (CCxCrop > 0.001_sp)) then
        if (ProcentWeedCover == 100) then
            fWeed = 1._sp/CCxCrop
        else
            fWeed = 1._sp - (1._sp - 1._sp/CCxCrop) * &
              (exp(FshapeWeed*ProcentWeedCover/100._sp) - 1._sp)/(exp(FshapeWeed) - 1._sp)
            if (fWeed > (1._sp/CCxCrop)) then
                fWeed = 1._sp/CCxCrop
            end if
        end if
    else
        fWeed = 1._sp
    end if
    CCmultiplierWeed = fWeed
end function CCmultiplierWeed


real(sp) function CCmultiplierWeedAdjusted(ProcentWeedCover, CCxCrop, FshapeWeed, fCCx, Yeari, MWeedAdj, RCadj)
    integer(int8), intent(in) :: ProcentWeedCover
    real(sp), intent(in) :: CCxCrop
    real(sp) :: FshapeWeed
    real(sp), intent(in) :: fCCx
    integer(int8), intent(in) :: Yeari
    integer(int8), intent(in) :: MWeedAdj
    integer(int8), intent(inout) :: RCadj

    real(sp) :: fWeedi, CCxTot100, CCxTot0, CCxTotM, fweedMax, RCadjD, FshapeMinimum

    fWeedi = 1._sp
    RCadj = ProcentWeedCover
    if (ProcentWeedCover > 0) then
        fweedi = CCmultiplierWeed(ProcentWeedCover, CCxCrop, FshapeWeed)
        ! FOR perennials when self-thinning
        if ((GetCrop_subkind() == subkind_Forage) .and. (Yeari > 1) .and. (fCCx < 0.995)) then
            ! need for adjustment
            ! step 1 - adjusment of shape factor to degree of crop replacement by weeds
            FshapeMinimum = 10 - 20*( (exp(fCCx*3._sp)-1)/(exp(3._sp)-1) + sqrt(MWeedAdj/100._sp))
            if (roundc(FshapeMinimum*10,mold=1_int32) == 0) then
                FshapeMinimum = 0.1
            end if
            FshapeWeed = FshapeWeed;
            if (FshapeWeed < FshapeMinimum) then
                FshapeWeed = FshapeMinimum
            end if

            ! step 2 - Estimate of CCxTot
            ! A. Total CC (crop and weeds) when self-thinning and 100% weed take over
            fweedi = CCmultiplierWeed(ProcentWeedCover, CCxCrop, FshapeWeed)
            CCxTot100 = fweedi * CCxCrop
            ! B. Total CC (crop and weeds) when self-thinning and 0% weed take over
            if (fCCx > 0.005) then
                fweedi = CCmultiplierWeed(roundc(fCCx*ProcentWeedCover,mold=1_int8),&
                    (fCCx*CCxCrop), FshapeWeed)
            else
                fweedi = 1
            end if
            CCxTot0 = fweedi * (fCCx*CCxCrop)
            ! C. total CC (crop and weeds) with specified weed take over (MWeedAdj)
            CCxTotM = CCxTot0 + (CCxTot100 - CCxTot0)* MWeedAdj/100
            if (CCxTotM < (fCCx*CCxCrop*(1-ProcentWeedCover/100._sp))) then
                CCxTotM = fCCx*CCxCrop*(1-ProcentWeedCover/100._sp)
            end if
            if (fCCx > 0.005) then
                fweedi = CCxTotM/(fCCx*CCxCrop)
                fweedMax = 1._sp/(fCCx*CCxCrop)
                if (roundc(fweedi*1000,mold=1_int32) > roundc(fWeedMax*1000,mold=1_int32)) then
                    fweedi = fweedMax
                end if
            end if

            ! step 3 - Estimate of adjusted weed cover
            RCadjD = ProcentWeedCover + (1-fCCx)*CCxCrop*MWeedAdj
            if (fCCx > 0.005) then
                if (RCadjD < (100*(CCxTotM - fCCx*CCxCrop)/CCxTotM)) then
                    RCadjD = 100*(CCxTotM - fCCx*CCxCrop)/CCxTotM
                end if
                if (RCadjD > (100 * (1- (fCCx*CCxCrop*(1-ProcentWeedCover/100._sp)/CCxTotM)))) then
                    RCadjD = 100*(1- fCCx*CCxCrop*(1-ProcentWeedCover/100._sp)/CCxTotM)
                end if
            end if
            RCadj = roundc(RCadjD,mold=1_int8)
            if (RCadj > 100) then
                RCadj = 100
            end if
        end if
    end if
    CCmultiplierWeedAdjusted = fWeedi
    ! CCmultiplierWeedAdjusted
end function CCmultiplierWeedAdjusted


real(sp) function BMRange(HIadj)
    integer(int32), intent(in) :: HIadj

    real(sp) :: BMR

    if (HIadj <= 0) then
        BMR = 0.0_sp
    else
        BMR = (log(real(HIadj, kind=sp))/0.0562_sp)/100.0_sp
    end if
    if (BMR > 1.0_sp) then
        BMR = 1.0_sp
    end if
    BMRange = BMR
end function BMRange


real(sp) function HImultiplier(RatioBM, RangeBM, HIadj)
    real(sp), intent(in) :: RatioBM
    real(sp), intent(in) :: RangeBM
    integer(int8), intent(in) :: HIadj

    real(sp) :: Rini, Rmax, Rend

    Rini = 1.0_sp - RangeBM
    REnd = 1.0_sp
    RMax = Rini + (2.0_sp/3.0_sp) * (REnd-Rini)
    if (RatioBM <= RIni) then
        HImultiplier = 1.0_sp
    elseif (RatioBM <= RMax) then
        HImultiplier = 1.0_sp&
            + (1.0_sp + sin(PI*(1.5_sp-(RatioBM-RIni)/(RMax-RIni))))&
            *(HIadj/200.0_sp)
    elseif (RatioBM <= REnd) then
        HImultiplier = 1.0_sp&
            + (1.0_sp + sin(PI*(0.5_sp+(RatioBM-RMax)/(REnd-RMax))))&
            *(HIadj/200.0_sp)
    else
        HImultiplier = 1.0_sp
    end if
end function HImultiplier


real(sp) function AdjustedKsStoToECsw(ECeMin, ECeMax, ResponseECsw, ECei, &
            ECswi, ECswFCi, Wrel, Coeffb0Salt, Coeffb1Salt, Coeffb2Salt, KsStoIN)
    integer(int8), intent(in) :: ECeMin
    integer(int8), intent(in) :: ECeMax
    integer(int32), intent(in) :: ResponseECsw
    real(sp), intent(in) :: ECei
    real(sp), intent(in) :: ECswi
    real(sp), intent(in) :: ECswFCi
    real(sp), intent(in) :: Wrel
    real(sp), intent(in) :: Coeffb0Salt
    real(sp), intent(in) :: Coeffb1Salt
    real(sp), intent(in) :: Coeffb2Salt
    real(sp), intent(in) :: KsStoIN

    real(sp) :: ECswRel, LocalKsShapeFactorSalt, KsSalti, SaltStressi, StoClosure, KsStoOut

    if ((ResponseECsw > 0) .and. (Wrel > epsilon(1._sp)) .and. &
                            (GetSimulation_SalinityConsidered() .eqv. .true.)) then
        ! adjustment to ECsw considered
        ECswRel = ECswi - (ECswFCi - ECei) + (ResponseECsw-100._sp)*Wrel
        if ((ECswRel > ECeMin) .and. (ECswRel < ECeMax)) then
            ! stomatal closure at ECsw relative
            LocalKsShapeFactorSalt = 3._sp ! CONVEX give best ECsw response
            KsSalti = KsSalinity(GetSimulation_SalinityConsidered(), ECeMin, &
                                        ECeMax, ECswRel, LocalKsShapeFactorSalt)
            SaltStressi = (1._sp-KsSalti)*100._sp
            StoClosure = Coeffb0Salt + Coeffb1Salt * SaltStressi + Coeffb2Salt &
                                                    * SaltStressi * SaltStressi
            ! adjusted KsSto
            KsStoOut = (1._sp - StoClosure/100._sp)
            if (KsStoOut < 0.0_sp) then
                KsStoOut = 0._sp
            end if
            if (KsStoOut > KsStoIN) then
                KsStoOut = KsStoIN
            end if
        else
            if (ECswRel >= ECeMax) then
                KsStoOut = 0._sp ! full stress
            else
                KsStoOut = KsStoIN ! no extra stress
            end if
        end if
    else
        KsStoOut = KsStoIN  ! no adjustment to ECsw
    end if
    AdjustedKsStoToECsw = KsStoOut
end function AdjustedKsStoToECsw


real(sp) function CCatTime(Dayi, CCoIN, CGCIN, CCxIN)
    integer(int32), intent(in) :: Dayi
    real(sp), intent(in) :: CCoIN
    real(sp), intent(in) :: CGCIN
    real(sp), intent(in) :: CCxIN

    real(sp) :: CCi

    CCi = CCoIN * exp(CGCIN * Dayi)
    if (CCi > CCxIN/2._sp) then
        CCi = CCxIN - 0.25_sp * (CCxIN/CCoIN) * CCxIN * exp(-CGCIN*Dayi)
    end if
    CCatTime = CCi
end function CCatTime


subroutine DetermineDayNr(Dayi, Monthi, Yeari, DayNr)
    integer(int32), intent(in) :: Dayi
    integer(int32), intent(in) :: Monthi
    integer(int32), intent(in) :: Yeari
    integer(int32), intent(inout) :: DayNr

    DayNr = trunc((Yeari - 1901)*365.25_sp + ElapsedDays(Monthi) + Dayi + 0.05_sp)
end subroutine DetermineDayNr


subroutine DetermineDate(DayNr, Dayi, Monthi, Yeari)
    integer(int32), intent(in) :: DayNr
    integer(int32), intent(inout) :: Dayi
    integer(int32), intent(inout) :: Monthi
    integer(int32), intent(inout) :: Yeari

    real(sp) :: SumDayMonth

    Yeari = trunc((DayNr-0.05_sp)/365.25_sp)
    SumDayMonth = (DayNr - Yeari*365.25_sp)
    Yeari = 1901 + Yeari
    Monthi = 1

    do while (Monthi < 12)
        if (SumDayMonth <= ElapsedDays(Monthi+1)) exit
        Monthi = Monthi + 1
    end do
    Dayi = roundc(SumDayMonth - ElapsedDays(Monthi) + 0.25_sp + 0.06_sp, &
                  mold=1_int32)
end subroutine DetermineDate


real(sp) function DegreesDay(Tbase, Tupper, TDayMin, TDayMax, GDDSelectedMethod)
    real(sp), intent(in) :: Tbase
    real(sp), intent(in) :: Tupper
    real(sp), intent(in) :: TDayMin
    real(sp), intent(in) :: TDayMax
    integer(int8), intent(in) :: GDDSelectedMethod

    real(sp) :: TstarMax, TstarMin
    real(sp) :: Tavg, DgrD

    select case (GDDSelectedMethod)
    case (1)
        ! Method 1. - No adjustemnt of Tmax, Tmin before calculation of Taverage
        Tavg = (TDayMax+TDayMin)/2._sp
        if (Tavg > Tupper) then
            Tavg = Tupper
        end if
        if (Tavg < Tbase) then
            Tavg = Tbase
        end if

    case (2)
        ! Method 2. -  Adjustment for Tbase before calculation of Taverage
        TstarMax = TDayMax
        if (TDayMax < Tbase) then
            TstarMax = Tbase
        end if
        if (TDayMax > Tupper) then
            TstarMax = Tupper
        end if
        TstarMin = TDayMin
        if (TDayMin < Tbase) then
            TstarMin = Tbase
        end if
        if (TDayMin > Tupper) then
            TstarMin = Tupper
        end if
        Tavg = (TstarMax+TstarMin)/2._sp

    case default
       ! Method 3.
        TstarMax = TDayMax
        if (TDayMax < Tbase) then
             TstarMax = Tbase
        end if
        if (TDayMax > Tupper) then
            TstarMax = Tupper
        end if
        TstarMin = TDayMin
        if (TDayMin > Tupper) then
            TstarMin = Tupper
        end if
        Tavg = (TstarMax+TstarMin)/2._sp
        if (Tavg < Tbase) then
            Tavg = Tbase
        end if
    end select
    DgrD =  Tavg - Tbase
    DegreesDay =  DgrD
end function DegreesDay


subroutine DetermineCNIandIII(CN2, CN1, CN3)
    integer(int8), intent(in) :: CN2
    integer(int8), intent(inout) :: CN1
    integer(int8), intent(inout) :: CN3

    CN1 = roundc(1.4_sp*(exp(-14*log(10._sp))) + 0.507_sp*CN2 &
                  - 0.00374_sp*CN2*CN2 + 0.0000867_sp*CN2*CN2*CN2, mold=1_int8)
    CN3 = roundc(5.6_sp*(exp(-14*log(10._sp))) + 2.33_sp*CN2 &
                 - 0.0209_sp*CN2*CN2 + 0.000076_sp*CN2*CN2*CN2, mold=1_int8)

    if (CN1 <= 0) then
        CN1 = 1
    elseif (CN1 > 100) then
        CN1 = 100
    end if
    if (CN3 <= 0) then
        CN3 = 1
    elseif (CN3 > 100) then
        CN3 = 100
    end if
    if (CN3 < CN2) then
        CN3 = CN2
    end if
end subroutine DetermineCNIandIII


subroutine DetermineCN_default(Infiltr, CN2)
    real(sp), intent(in) :: Infiltr
    integer(int8), intent(inout) :: CN2

    if (Infiltr > 864) then
        CN2 = 46
    elseif (Infiltr >= 347) then
        CN2 = 61
    elseif (Infiltr >= 36) then
        CN2 = 72
    else
        CN2 = 77
    end if
end subroutine DetermineCN_default


real(sp) function ECeComp(Comp)
    type(CompartmentIndividual), intent(in) :: Comp

    real(sp) :: volSat, TotSalt, denominator
    integer(int32) :: i

    volSAT = GetSoilLayer_SAT(Comp%Layer)
    TotSalt = 0._sp
    do i = 1, GetSoilLayer_SCP1(Comp%Layer)
        TotSalt = TotSalt + Comp%Salt(i) + Comp%Depo(i) ! g/m2
    end do

    denominator = volSAT*10._sp * Comp%Thickness * &
                  (1._sp - GetSoilLayer_GravelVol(Comp%Layer)/100._sp)
    TotSalt = TotSalt / denominator  ! g/l

    if (TotSalt > GetSimulParam_SaltSolub()) then
        TotSalt = GetSimulParam_SaltSolub()
    end if

    ECeComp = TotSalt / Equiv ! dS/m
end function ECeComp


real(sp) function ECswComp(Comp, atFC)
    type(CompartmentIndividual), intent(in) :: Comp
    logical, intent(in) :: atFC

    real(sp) :: TotSalt
    integer(int32) :: i

    TotSalt = 0
    do i = 1, GetSoilLayer_SCP1(Comp%Layer)
        TotSalt = TotSalt + Comp%Salt(i) + Comp%Depo(i) ! g/m2
    end do
    if (atFC .eqv. .true.) then
        TotSalt = TotSalt/ (GetSoilLayer_FC(Comp%Layer)*10._sp*Comp%Thickness &
                    *(1._sp-GetSoilLayer_GravelVol(Comp%Layer)/100._sp)) ! g/l
    else
        TotSalt = TotSalt/(Comp%theta*1000._sp*Comp%Thickness* &
                     (1._sp-GetSoilLayer_GravelVol(Comp%Layer)/100._sp)) ! g/l
    end if
    if (TotSalt > GetSimulParam_SaltSolub()) then
        TotSalt = GetSimulParam_SaltSolub()
    end if
    ECswComp = TotSalt/Equiv
end function ECswComp


subroutine SaltSolutionDeposit(mm, SaltSolution, SaltDeposit) ! mm = l/m2, SaltSol/Saltdepo = g/m2
    real(sp), intent(in) :: mm
    real(sp), intent(inout) :: SaltSolution
    real(sp), intent(inout) :: SaltDeposit


    SaltSolution = SaltSolution + SaltDeposit
    if (SaltSolution > GetSimulParam_SaltSolub() * mm) then
        SaltDeposit = SaltSolution - GetSimulParam_SaltSolub() * mm
        SaltSolution = GetSimulParam_SaltSolub() * mm
    else
        SaltDeposit = 0._sp
    end if
end subroutine SaltSolutionDeposit


real(sp) function MultiplierCCoSelfThinning(Yeari, Yearx, ShapeFactor)
    integer(int32), intent(in) :: Yeari
    integer(int32), intent(in) :: Yearx
    real(sp), intent(in) :: ShapeFactor

    real(sp) :: fCCo, Year0

    fCCo = 1._sp
    if ((Yeari >= 1) .and. (Yearx >= 2) .and. (roundc(100*ShapeFactor, mold=1_int32) /= 0)) then
        Year0 = 1._sp + (Yearx-1) * exp(ShapeFactor*log(10._sp))
        if ((Yeari >= Year0) .or. (Year0 <= 1)) then
            fCCo = 0._sp
        else
            fCCo = 1._sp - (Yeari-1)/(Year0-1._sp)
        end if
        if (fCCo < 0._sp) then
            fCCo = 0._sp
        end if
    end if
    MultiplierCCoSelfThinning = fCCo
end function MultiplierCCoSelfThinning


real(sp) function KsAny(Wrel, pULActual, pLLActual, ShapeFactor)
    real(sp), intent(in) :: Wrel
    real(sp), intent(in) :: pULActual
    real(sp), intent(in) :: pLLActual
    real(sp), intent(in) :: ShapeFactor

    real(sp) :: pRelativeLLUL, KsVal
    real(sp) :: pULActual_local
    ! Wrel : WC in rootzone (negative .... 0=FC ..... 1=WP .... > 1)
    !        FC .. UpperLimit ... LowerLimit .. WP
    ! p relative (negative .... O=UpperLimit ...... 1=LowerLimit .....>1)

    pULActual_local = pULActual

    if ((pLLActual - pULActual_local) < 0.0001_sp) then
        pULActual_local = pLLActual - 0.0001_sp
    end if

    pRelativeLLUL = (Wrel - pULActual_local)/(pLLActual - pULActual_local)

    if (pRelativeLLUL <= epsilon(0._sp)) then
        KsVal = 1._sp
    elseif (pRelativeLLUL >= 1._sp) then
        KsVal = 0._sp
    else
        if (roundc(10*ShapeFactor, mold=1_int32) == 0) then ! straight line
            KsVal = 1._sp - &
                    (exp(pRelativeLLUL*0.01_sp)-1._sp)/(exp(0.01_sp)-1._sp)
        else
            KsVal = 1._sp - &
                    (exp(pRelativeLLUL*ShapeFactor)-1._sp)/(exp(ShapeFactor)-1._sp)
        end if
        if (KsVal > 1._sp) then
            KsVal = 1._sp
        end if
        if (KsVal < 0._sp) then
            KsVal = 0._sp
        end if
    end if
    KsAny = KsVal
end function KsAny


real(sp) function CCatGDD(GDDi, CCoIN, GDDCGCIN, CCxIN)
    real(sp), intent(in) :: GDDi
    real(sp), intent(in) :: CCoIN
    real(sp), intent(in) :: GDDCGCIN
    real(sp), intent(in) :: CCxIN

    real(sp) :: CCi

    CCi = CCoIN * exp(GDDCGCIN * GDDi)
    if (CCi > CCxIN/2._sp) then
        CCi = CCxIN - 0.25_sp * (CCxIN/CCoIN) * CCxIN * exp(-GDDCGCIN*GDDi)
    end if
    CCatGDD = CCi
end function CCatGDD


real(sp) function CanopyCoverNoStressGDDaysSF(GDDL0, GDDL123, GDDLMaturity, SumGDD, &
        CCo, CCx, GDDCGC, GDDCDC, SFRedCGC, SFRedCCx)
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDL123
    integer(int32), intent(in) :: GDDLMaturity
    real(sp), intent(in) :: SumGDD
    real(sp), intent(in) :: CCo
    real(sp), intent(in) :: CCx
    real(sp), intent(in) :: GDDCGC
    real(sp), intent(in) :: GDDCDC
    integer(int8), intent(in) :: SFRedCGC
    integer(int8), intent(in) :: SFRedCCx

    real(sp) :: CC, CCxAdj, GDDCDCadj

    ! SumGDD refers to the end of the day and Delayed days are not considered
    CC = 0._sp
    if ((SumGDD > 0._sp) .and. (roundc(SumGDD, mold=1_int32) <= GDDLMaturity) .and. (CCo > 0._sp)) then
        if (SumGDD <= GDDL0) then ! before germination or recovering of transplant
            CC = 0._sp
        else
            if (SumGDD < GDDL123) then ! Canopy development and Mid-season stage
                CC = CCatGDD(real((SumGDD-GDDL0),sp), CCo, ((1._sp-SFRedCGC/100._sp)*GDDCGC), &
                  ((1._sp-SFRedCCx/100._sp)*CCx))
            else
                ! Late-season stage  (SumGDD <= GDDLMaturity)
                if (CCx < 0.001_sp) then
                    CC = 0._sp
                else
                    CCxAdj = CCatGDD(real(GDDL123-GDDL0,sp), CCo, ((1_sp-SFRedCGC/100._sp)*GDDCGC), &
                      ((1._sp-SFRedCCx/100._sp)*CCx))
                    GDDCDCadj = GDDCDC*(CCxadj+2.29_sp)/(CCx+2.29_sp)
                    if (CCxAdj < 0.001_sp) then
                        CC = 0._sp
                    else
                        CC = CCxAdj * (1._sp - 0.05_sp*(exp((SumGDD-GDDL123)*3.33_sp*GDDCDCadj/&
                          (CCxAdj+2.29_sp))-1._sp))
                    end if
                end if
            end if
        end if
    end if
    if (CC > 1._sp) then
        CC = 1._sp
    end if
    if (CC < 0._sp) then
        CC = 0._sp
    end if
    CanopyCoverNoStressGDDaysSF = CC
end function CanopyCoverNoStressGDDaysSF


real(sp) function HIadjWStressAtFlowering(KsVeg, KsSto, a, b)
    real(sp), intent(in) :: KsVeg
    real(sp), intent(in) :: KsSto
    integer(int8), intent(in) :: a
    real(sp), intent(in) :: b

    if (a == undef_int) then
        if (roundc(b, mold=1_int32) == undef_int) then
            HIadjWStressAtFlowering = 1._sp
        elseif (KsSto > 0.001_sp) then
            HIadjWStressAtFlowering = (exp(0.10_sp*log(KsSto))) * (1._sp-(1._sp-KsSto)/b)
        else
            HIadjWStressAtFlowering = 0.0_sp
        end if
    else
        if (roundc(b, mold=1_int32) == undef_int) then
            HIadjWStressAtFlowering = (1._sp + (1._sp-KsVeg)/a)
        elseif (KsSto > 0.001_sp) then
            HIadjWStressAtFlowering = (1._sp + (1._sp-KsVeg)/a) * (exp(0.10_sp*log(KsSto))) &
              * (1._sp-(1._sp-KsSto)/b)
        else
            HIadjWStressAtFlowering = 0._sp
        end if
    end if
end function HIadjWStressAtFlowering


real(sp) function fAdjustedForCO2(CO2i, WPi, PercentA)
    real(sp), intent(in) :: CO2i
    real(sp), intent(in) :: WPi
    integer(int8), intent(in) :: PercentA

    real(sp) :: fW, fCO2Old, fType, fSink, fShape, CO2rel, fCO2adj, fCO2

    ! 1. Correction for crop type: fType
    if (WPi >= 40._sp) then
        fType = 0._sp ! no correction for C4 crops
    else
        if (WPi <= 20._sp) then
            fType = 1._sp ! full correction for C3 crops
        else
            fType = (40._sp-WPi)/(40._sp-20._sp)
        end if
    end if

    ! 2. crop sink strength coefficient: fSink
    fSink = PercentA/100._sp
    if (fSink < 0._sp) then
        fSink = 0._sp ! based on FACE expirements
    end if
    if (fSink > 1._sp) then
        fSink = 1._sp ! theoretical adjustment
    end if

    ! 3. Correction coefficient for CO2: fCO2Old
    fCO2Old = undef_int
    if (CO2i <= 550._sp) then
        ! 3.1 weighing factor for CO2
        if (CO2i <= CO2Ref) then
            fw = 0._sp
        else
            if (CO2i >= 550._sp) then
                fW = 1._sp
            else
                fw = 1._sp - (550._sp - CO2i)/(550._sp - CO2Ref)
            end if
        end if

        ! 3.2 adjustment for CO2
        fCO2Old = (CO2i/CO2Ref)/(1._sp+(CO2i-CO2Ref)*((1._sp-fW)*0.000138_sp&
            + fW*(0.000138_sp*fSink + 0.001165_sp*(1._sp-fSink))))
    end if

    ! 4. Adjusted correction coefficient for CO2: fCO2adj
    fCO2adj = undef_int
    if (CO2i > CO2Ref) then
        ! 4.1 Shape factor
        fShape = -4.61824_sp - 3.43831_sp*fSink - 5.32587_sp*fSink*fSink

        ! 4.2 adjustment for CO2
        if (CO2i >= 2000._sp) then
            fCO2adj = 1.58_sp  ! maximum is reached
        else
            CO2rel = (CO2i-CO2Ref)/(2000._sp-CO2Ref)
            fCO2adj = 1._sp + 0.58_sp * ((exp(CO2rel*fShape) - 1._sp)/&
                (exp(fShape) - 1._sp))
        end if
    end if

    ! 5. Selected adjusted coefficient for CO2: fCO2
    if (CO2i <= CO2Ref) then
        fCO2 = fCO2Old
    else
        fCO2 = fCO2adj
        if ((CO2i <= 550._sp) .and. (fCO2Old < fCO2adj)) then
            fCO2 = fCO2Old
        end if
    end if

    ! 6. final adjustment
    fAdjustedForCO2 = 1._sp + fType*(fCO2-1._sp)
end function fAdjustedForCO2


logical function FullUndefinedRecord(FromY, FromD, FromM, ToD, ToM)
    integer(int32), intent(in) :: FromY
    integer(int32), intent(in) :: FromD
    integer(int32), intent(in) :: FromM
    integer(int32), intent(in) :: ToD
    integer(int32), intent(in) :: ToM

    FullUndefinedRecord = ((FromY == 1901) .and. (FromD == 1)&
        .and. (FromM == 1) .and. (ToD == 31) .and. (ToM == 12))
end function FullUndefinedRecord


subroutine NoIrrigation()

    integer(int32) :: Nri

    call SetIrriMode(IrriMode_NoIrri)
    call SetIrriDescription('Rainfed cropping')
    call SetIrriMethod(IrriMethod_MSprinkler)
    call SetSimulation_IrriECw(0.0_sp) ! dS/m
    call SetGenerateTimeMode(GenerateTimeMode_AllRAW)
    call SetGenerateDepthMode(GenerateDepthMode_ToFC)
    IrriFirstDayNr = undef_int
    do Nri = 1, 5
        call SetIrriBeforeSeason_DayNr(Nri, 0)
        call SetIrriBeforeSeason_Param(Nri, 0)
        call SetIrriAfterSeason_DayNr(Nri, 0)
        call SetIrriAfterSeason_Param(Nri, 0)
    end do
    call SetIrriECw_PreSeason(0.0_sp) ! dS/m
    call SetIrriECw_PostSeason(0.0_sp) ! dS/m
end subroutine NoIrrigation


subroutine LoadIrriScheduleInfo(FullName)
    character(len=*), intent(in) :: FullName

    integer(int32) :: fhandle
    integer(int32) :: i, rc
    real(sp) :: VersionNr
    integer(int8) :: simul_irri_in
    integer(int32) :: simul_percraw

    open(newunit=fhandle, file=trim(FullName), status='old', action='read')
    read(fhandle, '(a)', iostat=rc) IrriDescription
    read(fhandle, *, iostat=rc) VersionNr  ! AquaCrop version

    ! irrigation method
    read(fhandle, *, iostat=rc) i
    select case (i)
    case(1)
        call SetIrriMethod(IrriMethod_MSprinkler)
    case(2)
        call SetIrriMethod(IrriMethod_MBasin)
    case(3)
        call SetIrriMethod(IrriMethod_MBorder)
    case(4)
        call SetIrriMethod(IrriMethod_MFurrow)
    case default
        call SetIrriMethod(IrriMethod_MDrip)
    end select
    ! fraction of soil surface wetted
    read(fhandle, *, iostat=rc) simul_irri_in
    call SetSimulParam_IrriFwInSeason(simul_irri_in)

    ! irrigation mode and parameters
    read(fhandle, *, iostat=rc) i
    select case (i)
    case(0)
        call SetIrriMode(IrriMode_NoIrri) ! rainfed
    case(1)
        call SetIrriMode(IrriMode_Manual)
    case(2)
        call SetIrriMode(IrriMode_Generate)
    case default
        call SetIrriMode(IrriMode_Inet)
    end select

    ! 1. Irrigation schedule
    if ((i == 1) .and. (roundc(VersionNr*10,mold=1) >= 70)) then
        read(fhandle, *, iostat=rc) IrriFirstDayNr ! line 6
    else
        IrriFirstDayNr = undef_int ! start of growing period
    end if


    ! 2. Generate
    if (GetIrriMode() == IrriMode_Generate) then
        read(fhandle, *, iostat=rc) i ! time criterion
        select case (i)
        case(1)
            call SetGenerateTimeMode(GenerateTimeMode_FixInt)
        case(2)
            call SetGenerateTimeMode(GenerateTimeMode_AllDepl)
        case(3)
            call SetGenerateTimeMode(GenerateTimeMode_AllRAW)
        case(4)
            call SetGenerateTimeMode(GenerateTimeMode_WaterBetweenBunds)
        case default
            call SetGenerateTimeMode(GenerateTimeMode_AllRAW)
        end select
        read(fhandle, *, iostat=rc) i ! depth criterion
        select case (i)
        case(1)
            call SetGenerateDepthMode(GenerateDepthMode_ToFc)
        case default
            call SetGenerateDepthMode(GenerateDepthMode_FixDepth)
        end select
        IrriFirstDayNr = undef_int ! start of growing period
    end if

    ! 3. Net irrigation requirement
    if (GetIrriMode() == IrriMode_Inet) then
        read(fhandle, *, iostat=rc) simul_percraw
        call SetSimulParam_PercRAW(simul_percraw)
        IrriFirstDayNr = undef_int  ! start of growing period
    end if
    close(fhandle)
    ! LoadIrriScheduleInfo
end subroutine LoadIrriScheduleInfo


subroutine GenerateCO2Description(CO2FileFull, CO2Description)
    character(len=*), intent(in) :: CO2FileFull
    character(len=*), intent(inout) :: CO2Description

    integer :: fhandle

    open(newunit=fhandle, file=trim(CO2FileFull), status='old', &
         action='read')
    read(fhandle, *) CO2Description
    close(fhandle)

    if (trim(GetCO2File()) == 'MaunaLoa.CO2') then
        ! since this is an AquaCrop file, the Description is determined by AquaCrop
        CO2Description = 'Default atmospheric CO2 concentration from 1902 to 2099'
    end if
end subroutine GenerateCO2Description


subroutine GetIrriDescription(IrriFileFull, IrriDescription)
    character(len=*), intent(in) :: IrriFileFull
    character(len=*), intent(inout) :: IrriDescription

    integer :: fhandle

    open(newunit=fhandle, file=trim(IrriFileFull), status='old', &
         action='read')
    read(fhandle, *) IrriDescription
    close(fhandle)
end subroutine GetIrriDescription


subroutine SetIrriDescription(str)
    !! Setter for the "IrriDescription" global variable.
    character(len=*), intent(in) :: str

    IrriDescription = str
end subroutine SetIrriDescription


subroutine GetDaySwitchToLinear(HImax, dHIdt, HIGC, tSwitch, HIGClinear)
    integer(int32), intent(in) :: HImax
    real(sp), intent(in) :: dHIdt
    real(sp), intent(in) :: HIGC
    integer(int32), intent(inout) :: tSwitch
    real(sp), intent(inout) :: HIGClinear

    real(sp) :: HIi, HiM1, HIfinal
    integer(int32) :: tmax, ti
    integer(int32), parameter :: HIo = 1

    tmax = roundc(HImax/dHIdt, mold=1_int32)
    ti = 0
    HiM1 = HIo
    if (tmax > 0) then
        loop: do
            ti = ti + 1
            HIi = (HIo*HImax)/ (HIo+(HImax-HIo)*exp(-HIGC*ti))
            HIfinal = HIi + (tmax - ti)*(HIi-HIM1)
            HIM1 = HIi
            if ((HIfinal > HImax) .or. (ti >= tmax)) exit loop
        end do loop
        tSwitch = ti - 1
    else
        tSwitch = 0
    end if
    if (tSwitch > 0) then
        HIi = (HIo*HImax)/ (HIo+(HImax-HIo)*exp(-HIGC*tSwitch))
    else
        HIi = 0
    end if
    HIGClinear = (HImax-HIi)/(tmax-tSwitch)
end subroutine GetDaySwitchToLinear


logical function FileExists(full_name)
    character(len=*), intent(in) :: full_name

    inquire(file=trim(full_name), exist=FileExists)
end function FileExists


subroutine SplitStringInTwoParams(StringIN, Par1, Par2)
    character(len=*), intent(in) :: StringIN
    real(sp), intent(inout) :: Par1
    real(sp), intent(inout) :: Par2

    integer(int32) :: LengthS, i, Parami
    character :: CharA
    character(len=255) :: StringNumber

    LengthS = len(StringIN)
    i = 0
    Parami = 0
    ! divide the line in parameters
    do while ((i < LengthS) .and. (Parami < 2))
        i = i + 1
        CharA = StringIN(i:i)
        if (ichar(CharA) > 32) then
            ! next Parameter
            Parami = Parami + 1
            StringNumber = ''
            do while ((ichar(CharA) > 32) .and. (i <= LengthS))
                StringNumber = trim(StringNumber) // CharA
                i = i + 1
                if (i <= LengthS) then
                    CharA = StringIN(i:i)
                end if
            end do
            if (Parami == 1) then
                read(StringNumber, *) Par1
            end if
            if (Parami == 2) then
                read(StringNumber, *) Par2
            end if
            ! next Parameter
        end if
        ! end of line
    end do
end subroutine SplitStringInTwoParams


subroutine SplitStringInThreeParams(StringIN, Par1, Par2, Par3)
    character(len=*), intent(in) :: StringIN
    real(sp), intent(inout) :: Par1
    real(sp), intent(inout) :: Par2
    real(sp), intent(inout) :: Par3

    integer(int32) :: LengthS, i, Parami
    character :: CharA
    character(len=255) :: StringNumber

    LengthS = len(StringIN)
    i = 0
    Parami = 0
    ! divide the line in parameters
    do while ((i < LengthS) .and. (Parami < 3))
        i = i + 1
        CharA = StringIN(i:i)
        if (ichar(CharA) > 32) then
            ! next Parameter
            Parami = Parami + 1
            StringNumber = ''
            do while ((ichar(CharA) > 32) .and. (i <= LengthS))
                StringNumber = trim(StringNumber) // CharA
                i = i + 1
                if (i <= LengthS) then
                    CharA = StringIN(i:i)
                end if
            end do
            if (Parami == 1) then
                read(StringNumber, *) Par1
            end if
            if (Parami == 2) then
                read(StringNumber, *) Par2
            end if
            if (Parami == 3) then
                read(StringNumber, *) Par3
            end if
            ! next Parameter
        end if
        ! end of line
    end do
end subroutine SplitStringInThreeParams


real(sp) function CO2ForSimulationPeriod(FromDayNr, ToDayNr)
    integer(int32), intent(in) :: FromDayNr
    integer(int32), intent(in) :: ToDayNr

    integer(int32) :: i, Dayi, Monthi, FromYi, ToYi, rc
    real(sp) :: CO2From, CO2To, CO2a, CO2b, YearA, YearB
    integer :: fhandle
    character(len=255) :: TempString

    call DetermineDate(FromDayNr, Dayi, Monthi, FromYi)
    call DetermineDate(ToDayNr, Dayi, Monthi, ToYi)

    if ((FromYi == 1901) .or. (ToYi == 1901)) then
        CO2ForSimulationPeriod = CO2Ref
    else
        open(newunit=fhandle, file=trim(GetCO2FileFull()), status='old', &
                                                    action='read',iostat=rc)
        do i= 1, 3
            read(fhandle, *, iostat=rc) ! Description and Title
        end do
        ! from year
        read(fhandle, '(a)', iostat=rc) TempString
        call SplitStringInTwoParams(trim(TempString), YearB, CO2b)
        if (roundc(YearB, mold=1) >= FromYi) then
            CO2From = CO2b
            YearA = YearB
            CO2a = CO2b
        else
            loop: do
                YearA = YearB
                CO2a = Co2b
                read(fhandle, '(a)', iostat=rc) TempString
                call SplitStringInTwoParams(trim(TempString), YearB, CO2b)
                if ((roundc(YearB, mold=1) >= FromYi) .or. (rc == iostat_end)) exit loop
            end do loop
            if (FromYi > roundc(YearB, mold=1)) then
                CO2From = CO2b
            else
                CO2From = CO2a + (CO2b-CO2a)* (FromYi - &
                roundc(YearA, mold=1))/real(roundc(YearB, mold=1)-roundc(YearA, mold=1),kind=sp)
            end if
        end if
        ! to year
        CO2To = CO2From
        if ((ToYi > FromYi) .and. (ToYi > roundc(YearA, mold=1))) then
            if (roundc(YearB, mold=1) >= ToYi) then
                CO2To = CO2a + (CO2b-CO2a)* (ToYi - &
                roundc(YearA, mold=1))/real(roundc(YearB, mold=1)-roundc(YearA, mold=1), kind=sp)
            elseif (.not. (rc == iostat_end)) then
                loop_2: do
                    YearA = YearB
                    CO2a = Co2b
                    read(fhandle, '(a)', iostat=rc) TempString
                    call SplitStringInTwoParams(trim(TempString), YearB, CO2b)
                    if (((roundc(YearB, mold=1) >= ToYi) .or. (rc == iostat_end))) &
                                                                    exit loop_2
                end do loop_2
                if (ToYi > roundc(YearB, mold=1)) then
                    CO2To = CO2b
                else
                    CO2To = CO2a + (CO2b-CO2a)* (ToYi - &
                    roundc(YearA, mold=1))/real(roundc(YearB, mold=1)-roundc(YearA, mold=1), kind=sp)
                end if
            end if
        end if
        Close(fhandle)
        CO2ForSimulationPeriod = (CO2From+CO2To)/2._sp
    end if
end function CO2ForSimulationPeriod


subroutine ReadRainfallSettings()

    integer :: fhandle
    character :: fullname
    integer(int8) :: NrM, effrainperc,effrainshow,effrainrootE

    fullName = trim(GetPathNameSimul()) // 'Rainfall.PAR'

    open(newunit=fhandle, file=trim(fullname), status='old', action='read')
    read(fhandle, *)! Settings for processing 10-day or monthly rainfall data
    read(fhandle, *) NrM
    select case (NrM)
        case (0)
            call SetSimulParam_EffectiveRain_Method(EffectiveRainMethod_full)
        case (1)
            call SetSimulParam_EffectiveRain_Method(EffectiveRainMethod_usda)
        case (2)
            call SetSimulParam_EffectiveRain_Method(EffectiveRainMethod_percentage)
    end select
    read(fhandle, *) effrainperc ! IF Method is Percentage
    call SetSimulParam_EffectiveRain_PercentEffRain(effrainperc)
    read(fhandle, *) effrainshow  ! For estimation of surface run-off
    call SetSimulParam_EffectiveRain_ShowersInDecade(effrainshow)
    read(fhandle, *) effrainrootE ! For reduction of soil evaporation
    call SetSimulParam_EffectiveRain_RootNrEvap(effrainrootE)
    close(fhandle)
end subroutine ReadRainfallSettings


subroutine ReadSoilSettings()

    integer :: fhandle
    character(len=:), allocatable :: fullName
    integer(int8) :: i, simul_saltdiff, simul_saltsolub, simul_root, simul_iniab
    real(sp) :: simul_rod

    fullName = trim(GetPathNameSimul()) // 'Soil.PAR'

    open(newunit=fhandle, file=trim(fullname), status='old', action='read')
    read(fhandle,*) simul_rod ! considered depth (m) of soil profile for calculation of mean soil water content
    call SetSimulParam_RunoffDepth(simul_rod)
    read(fhandle, *) i   ! correction CN for Antecedent Moisture Class
    if (i == 1) then
        call SetSimulParam_CNcorrection(.true.)
    else
        call SetSimulParam_CNcorrection(.false.)
    end if
    read(fhandle, *) simul_saltdiff ! salt diffusion factor (%)
    read(fhandle, *) simul_saltsolub ! salt solubility (g/liter)
    read(fhandle, *) simul_root ! shape factor capillary rise factor
    call SetSimulParam_SaltDiff(simul_saltdiff)
    call SetSimulParam_SaltSolub(simul_saltsolub)
    call SetSimulParam_RootNrDF(simul_root)
    ! new Version 4.1
    read(fhandle, *) simul_iniab ! Percentage of S for initial abstraction for surface runoff
    call SetSimulParam_IniAbstract(simul_iniab)
    call SetSimulParam_IniAbstract(5_int8) ! fixed in Version 5.0 cannot be changed since linked with equations for CN AMCII and CN converions
    close(fhandle)
end subroutine ReadSoilSettings


subroutine LoadClimate(FullName, ClimateDescription, TempFile, EToFile, RainFile, CO2File)
    character(len=*), intent(in) :: FullName
    character(len=*), intent(inout) :: ClimateDescription
    character(len=*), intent(inout) :: TempFile
    character(len=*), intent(inout) :: EToFile
    character(len=*), intent(inout) :: RainFile
    character(len=*), intent(inout) :: CO2File

    integer :: fhandle

    open(newunit=fhandle, file=trim(FullName), status='old', &
         action='read')
    read(fhandle, '(a)') ClimateDescription
    read(fhandle) ! AquaCrop Version
    read(fhandle, *) TempFile
    read(fhandle, *) EToFile
    read(fhandle, *) RainFile
    read(fhandle, *) CO2File
    close(fhandle)
end subroutine LoadClimate


subroutine LoadCropCalendar(FullName, GetOnset, GetOnsetTemp, DayNrStart, YearStart)
    character(len=*), intent(in) :: FullName
    logical, intent(inout) :: GetOnset
    logical, intent(inout) :: GetOnsetTemp
    integer(int32), intent(inout) :: DayNrStart
    integer(int32), intent(in) :: YearStart

    integer :: fhandle
    integer(int8) :: Onseti
    integer(int32) :: Dayi, Monthi, Yeari, CriterionNr
    integer(int32) :: DayNr
    GetOnset = .false.
    GetOnsetTemp = .false.

    open(newunit=fhandle, file=trim(FullName), status='old', action='read')
    read(fhandle, '(a)') CalendarDescription
    read(fhandle, *) ! AquaCrop Version

    ! Specification of Onset and End growing season
    read(fhandle, *) Onseti ! specification of the onset

    ! Onset growing season
    if (Onseti == 0) then
        ! onset on a specific day
        read(fhandle, *) ! start search period - not applicable
        read(fhandle, *) ! length search period - not applicable
        read(fhandle, *) DayNr ! day-number
        call DetermineDate(DayNr, Dayi, Monthi, Yeari)
        call DetermineDayNr(Dayi, Monthi, YearStart, DayNrStart)
    else
        ! onset is generated
        GetOnset = .true.
        read(fhandle, *) ! start search period
        read(fhandle, *) ! length search period
        read(fhandle, *) CriterionNr ! criterion number to decide if based on rainfall or air temperature
        if (CriterionNr > 10) then
            GetOnsetTemp = .true.
        end if
    end if
    close(fhandle)
end subroutine LoadCropCalendar


subroutine NoManagement()

    type(rep_EffectStress) :: EffectStress_temp

    call SetManDescription('No specific field management')
    ! mulches
    call SetManagement_Mulch(0_int8)
    call SetManagement_EffectMulchInS(50_int8)
    ! soil fertility
    call SetManagement_FertilityStress(0_int32)
    EffectStress_temp = GetSimulation_EffectStress()
    call CropStressParametersSoilFertility(GetCrop_StressResponse(), &
                                      GetManagement_FertilityStress(), &
                                      EffectStress_temp)
    call SetSimulation_EffectStress(EffectStress_temp)
    ! soil bunds
    call SetManagement_BundHeight(0._sp)
    call SetSimulation_SurfaceStorageIni(0.0_sp)
    call SetSimulation_ECStorageIni(0.0_sp)
    ! surface run-off
    call SetManagement_RunoffOn(.true.)
    call SetManagement_CNcorrection(0)
    ! weed infestation
    call SetManagement_WeedRC(0_int8)
    call SetManagement_WeedDeltaRC(0)
    call SetManagement_WeedShape(-0.01_sp)
    call SetManagement_WeedAdj(100_int8)
    ! multiple cuttings
    call SetManagement_Cuttings_Considered(.false.)
    call SetManagement_Cuttings_CCcut(30)
    call SetManagement_Cuttings_Day1(1)
    call SetManagement_Cuttings_NrDays(undef_int)
    call SetManagement_Cuttings_Generate(.false.)
    call SetManagement_Cuttings_Criterion(TimeCuttings_NA)
    call SetManagement_Cuttings_HarvestEnd(.false.)
    call SetManagement_Cuttings_FirstDayNr(undef_int)
end subroutine NoManagement


subroutine LoadManagement(FullName)
    character(len=*), intent(in) :: FullName

    integer :: fhandle
    integer(int8) :: i
    real(sp) :: VersionNr
    integer(int8) :: TempShortInt
    integer(int32) :: TempInt
    real(sp) :: TempDouble
    type(rep_EffectStress) :: EffectStress_temp
    character(len=1025) :: mandescription_temp

    open(newunit=fhandle, file=trim(FullName), status='old', action='read')
    read(fhandle, '(a)') mandescription_temp
    call SetManDescription(trim(mandescription_temp))
    read(fhandle, *) VersionNr ! AquaCrop Version
    ! mulches
    read(fhandle, *) TempShortInt
    call SetManagement_Mulch(TempShortInt)
    read(fhandle, *) TempShortInt
    call SetManagement_EffectMulchInS(TempShortInt)
    ! soil fertility
    read(fhandle, *) TempInt ! effect is crop specific
    call SetManagement_FertilityStress(TempInt)
    EffectStress_temp = GetSimulation_EffectStress()
    call CropStressParametersSoilFertility(GetCrop_StressResponse(), &
                                      GetManagement_FertilityStress(), &
                                      EffectStress_temp)
    call SetSimulation_EffectStress(EffectStress_temp)
    ! soil bunds
    read(fhandle, *) TempDouble
    call SetManagement_BundHeight(TempDouble)
    call SetSimulation_SurfaceStorageIni(0.0_sp)
    call SetSimulation_ECStorageIni(0.0_sp)
    ! surface run-off
    read(fhandle, *) i
    if (i == 1) then
        call SetManagement_RunoffON(.false.)   ! prevention of surface runoff
    else
        call SetManagement_RunoffON(.true.)   ! surface runoff is not prevented
    end if
    if (roundc(VersionNr*10, mold=1) < 50) then
        ! UPDATE required for CN adjustment
        call SetManagement_CNcorrection(0)
    else
        read(fhandle, *) TempInt ! project increase/decrease of CN
        call SetManagement_CNcorrection(TempInt)
    end if
    ! weed infestation
    if (roundc(VersionNr*10, mold=1) < 50) then
        ! UPDATE required for Version 3.0, 3.1 and 4.0
        call SetManagement_WeedRC(0_int8) ! relative cover of weeds (%)
        call SetManagement_WeedDeltaRC(0)
        call SetManagement_WeedShape(-0.01_sp) ! shape factor of the CC expansion
                                          ! function in a weed infested field
    else
        read(fhandle, *) TempShortInt ! relative cover of weeds (%)
        call SetManagement_WeedRC(TempShortInt)
        if (roundc(VersionNr*10, mold=1) < 51) then
            call SetManagement_WeedDeltaRC(0)
        else
            read(fhandle, *) TempInt
            call SetManagement_WeedDeltaRC(TempInt)
        end if
        read(fhandle, *) TempDouble ! shape factor of the CC expansion
                                    ! function in a weed infested field
        call SetManagement_WeedShape(TempDouble)
    end if
    if (roundc(VersionNr*10, mold=1) < 70) then
        ! UPDATE required for versions below 7
        call SetManagement_WeedAdj(100_int8) ! replacement (%) by weeds of the
                                        ! self-thinned part of the Canopy Cover
                                        ! - only for perennials
    else
        read(fhandle, *) TempShortInt
        call SetManagement_WeedAdj(TempShortInt)
    end if
    ! multiple cuttings
    if (roundc(VersionNr*10, mold=1) >= 70) then
        ! UPDATE required for multiple cuttings
        read(fhandle, *) i  ! Consider multiple cuttings: True or False
        if (i == 0) then
            call SetManagement_Cuttings_Considered(.false.)
        else
            call SetManagement_Cuttings_Considered(.true.)
        end if
        read(fhandle, *) TempInt  ! Canopy cover (%) after cutting
        call SetManagement_Cuttings_CCcut(TempInt)
        ! Next line is expected to be present in the input file, however
        ! A PARAMETER THAT IS NO LONGER USED since AquaCrop version 7.1
        read(fhandle, *) TempInt ! Increase (percentage) of CGC after cutting
        read(fhandle, *) TempInt ! Considered first day when generating cuttings
                                 ! (1 = start of growth cycle)
        call SetManagement_Cuttings_Day1(TempInt)
        read(fhandle, *) TempInt  ! Considered number owhen generating cuttings
                                  ! (-9 = total growth cycle)
        call SetManagement_Cuttings_NrDays(TempInt)
        read(fhandle, *) i  ! Generate multiple cuttings: True or False
        if (i == 1) then
            call SetManagement_Cuttings_Generate(.true.)
        else
            call SetManagement_Cuttings_Generate(.false.)
        end if
        read(fhandle, *) i  ! Time criterion for generating cuttings
        select case (i)
            case(0)
                call SetManagement_Cuttings_Criterion(TimeCuttings_NA)
                ! not applicable
            case(1)
                call SetManagement_Cuttings_Criterion(TimeCuttings_IntDay)
                ! interval in days
            case(2)
                call SetManagement_Cuttings_Criterion(TimeCuttings_IntGDD)
                ! interval in Growing Degree Days
            case(3)
                call SetManagement_Cuttings_Criterion(TimeCuttings_DryB)
                ! produced dry above ground biomass (ton/ha)
            case(4)
                call SetManagement_Cuttings_Criterion(TimeCuttings_DryY)
                ! produced dry yield (ton/ha)
            case(5)
                call SetManagement_Cuttings_Criterion(TimeCuttings_FreshY)
                ! produced fresh yield (ton/ha)
        end select
        read(fhandle, *) i  ! final harvest at crop maturity:
                            ! True or False (When generating cuttings)
        if (i == 1) then
            call SetManagement_Cuttings_HarvestEnd(.true.)
        else
            call SetManagement_Cuttings_HarvestEnd(.false.)
        end if
        read(fhandle, *) TempInt ! dayNr for Day 1 of list of cuttings
                                 ! (-9 = Day1 is start growing cycle)
        call SetManagement_Cuttings_FirstDayNr(TempInt)
    else
        call SetManagement_Cuttings_Considered(.false.)
        call SetManagement_Cuttings_CCcut(30)
        call SetManagement_Cuttings_Day1(1)
        call SetManagement_Cuttings_NrDays(undef_int)
        call SetManagement_Cuttings_Generate(.false.)
        call SetManagement_Cuttings_Criterion(TimeCuttings_NA)
        call SetManagement_Cuttings_HarvestEnd(.false.)
        call SetManagement_Cuttings_FirstDayNr(undef_int)
    end if
    close(fhandle)
end subroutine LoadManagement


subroutine SaveCrop(totalname)
    character(len=*), intent(in) :: totalname

    integer :: fhandle
    integer(int32) :: i, j
    character(len=:), allocatable :: TempString

    open(newunit=fhandle, file=trim(totalname), status='replace', action='write')
    write(fhandle, '(a)') GetCropDescription()
    ! AquaCrop version
    write(fhandle, '("     ", a, "       : AquaCrop Version (", a, ")")') &
          GetVersionString(), GetReleaseDate()
    write(fhandle, '(a)') '     1         : File not protected'

    ! SubKind
    i = 2
    select case (GetCrop_subkind())
        case(subkind_Vegetative)
            i = 1
            TempString = '         : leafy vegetable crop'
        case(subkind_Grain)
            i = 2
            TempString = '         : fruit/grain producing crop'
        case(subkind_Tuber)
            i = 3
            TempString = '         : root/tuber crop'
        case(subkind_Forage)
            i = 4
            TempString = '         : forage crop'
    end select
    write(fhandle, '(i6,a)') i, TempString

    ! Sown, transplanting or regrowth
    if (GetCrop_Planting() == plant_Seed) then
        i = 1
        if (GetCrop_subkind() == subkind_Forage) then
            write(fhandle, '(i6,a)') i, '         : Crop is sown in 1st year'
        else
            write(fhandle, '(i6,a)') i, '         : Crop is sown'
        end if
    else
        if (GetCrop_Planting() == plant_Transplant) then
            i = 0
            if (GetCrop_subkind() == subkind_Forage) then
                write(fhandle, '(i6,a)') i, &
                                    '         : Crop is transplanted in 1st year'
            else
                write(fhandle, '(i6,a)') i, &
                                    '         : Crop is transplanted'
            end if
        else
            i = -9
            write(fhandle, '(i6,a)') i, &
                                    '         : Crop is regrowth'
        end if
    end if

    ! Mode (description crop cycle)
    i = 1
    TempString = '         : Determination of crop cycle : by calendar days'
    if (GetCrop_ModeCycle() == ModeCycle_GDDays) then
        i = 0
        TempString = '         : Determination of crop cycle : by growing degree-days'
    end if
    write(fhandle, '(i6,a)') i, TempString

    ! p correction for ET
    if (GetCrop_pMethod() == pMethod_NoCorrection) then
        j = 0
        write(fhandle, '(i6,a)') j, &
            '         : No adjustment by ETo of soil water depletion factors (p)'
    else
        j = 1
        write(fhandle, '(i6,a)') j, &
            '         : Soil water depletion factors (p) are adjusted by ETo'
    end if

    ! temperatures controlling crop development
    write(fhandle, '(f8.1,a)') GetCrop_Tbase(), &
    '       : Base temperature (degC) below which crop development does not progress'
    write(fhandle, '(f8.1,a)') GetCrop_Tupper(), &
    '       : Upper temperature (degC) above which crop development no longer increases with an increase in temperature'

    ! required growing degree days to complete the crop cycle (is identical as to maturity)
    write(fhandle, '(i6,a)') GetCrop_GDDaysToHarvest(), &
    '         : Total length of crop cycle in growing degree-days'

    ! water stress
    write(fhandle, '(f9.2,a)') GetCrop_pLeafDefUL(), &
    '      : Soil water depletion factor for canopy expansion (p-exp) - Upper threshold'
    write(fhandle, '(f9.2,a)') GetCrop_pLeafDefLL(), &
    '      : Soil water depletion factor for canopy expansion (p-exp) - Lower threshold'
    write(fhandle, '(f8.1,a)') GetCrop_KsShapeFactorLeaf(), &
    '       : Shape factor for water stress coefficient for canopy expansion (0.0 = straight line)'
    write(fhandle, '(f9.2,a)') GetCrop_pdef(), &
    '      : Soil water depletion fraction for stomatal control (p - sto) - Upper threshold'
    write(fhandle, '(f8.1,a)') GetCrop_KsShapeFactorStomata(), &
    '       : Shape factor for water stress coefficient for stomatal control (0.0 = straight line)'
    write(fhandle, '(f9.2,a)') GetCrop_pSenescence(), &
    '      : Soil water depletion factor for canopy senescence (p - sen) - Upper threshold'
    write(fhandle, '(f8.1,a)') GetCrop_KsShapeFactorSenescence(), &
    '       : Shape factor for water stress coefficient for canopy senescence (0.0 = straight line)'
    write(fhandle, '(i6,a)') GetCrop_SumEToDelaySenescence(), &
    '         : Sum(ETo) during dormant period to be exceeded before crop is permanently wilted'
    if (abs(GetCrop_pPollination() - undef_int) < epsilon(0._sp)) then
        write(fhandle, '(f9.2,a)') GetCrop_pPollination(), &
        '      : Soil water depletion factor for pollination - Not Applicable'
    else
        write(fhandle, '(f9.2,a)') GetCrop_pPollination(), &
        '      : Soil water depletion factor for pollination (p - pol) - Upper threshold'
    end if
    write(fhandle, '(i6,a)') GetCrop_AnaeroPoint(), &
    '         : Vol% for Anaerobiotic point (* (SAT - [vol%]) at which deficient aeration occurs *)'

    ! stress response
    write(fhandle, '(i6,a)') GetCrop_StressResponse_Stress(), &
    '         : Considered soil fertility stress for calibration of stress response (%)'
    if (GetCrop_StressResponse_ShapeCGC() > 24.9_sp) then
        write(fhandle, '(f9.2,a)') GetCrop_StressResponse_ShapeCGC(), &
        '      : Response of canopy expansion is not considered'
    else
        write(fhandle, '(f9.2,a)') GetCrop_StressResponse_ShapeCGC(), &
        '      : Shape factor for the response of canopy expansion to soil fertility stress'
    end if
    if (GetCrop_StressResponse_ShapeCCX() > 24.9_sp) then
        write(fhandle, '(f9.2,a)') GetCrop_StressResponse_ShapeCCX(), &
        '      : Response of maximum canopy cover is not considered'
    else
        write(fhandle, '(f9.2,a)') GetCrop_StressResponse_ShapeCCX(), &
        '      : Shape factor for the response of maximum canopy cover to soil fertility stress'
    end if
    if (GetCrop_StressResponse_ShapeWP() > 24.9_sp) then
        write(fhandle, '(f9.2,a)') GetCrop_StressResponse_ShapeWP(), &
        '      : Response of crop Water Productivity is not considered'
    else
        write(fhandle, '(f9.2,a)') GetCrop_StressResponse_ShapeWP(), &
        '      : Shape factor for the response of crop Water Productivity to soil fertility stress'
    end if
    if (GetCrop_StressResponse_ShapeCDecline() > 24.9_sp) then
        write(fhandle, '(f9.2,a)') GetCrop_StressResponse_ShapeCDecline(), &
        '      : Response of decline of canopy cover is not considered'
    else
        write(fhandle, '(f9.2,a)') GetCrop_StressResponse_ShapeCDecline(), &
        '      : Shape factor for the response of decline of canopy cover to soil fertility stress'
    end if
    write(fhandle, '(a)') '    -9         : dummy - Parameter no Longer required'

    ! temperature stress
    if (int(GetCrop_Tcold(), int32) == undef_int) then
        write(fhandle, '(i6,a)') GetCrop_Tcold(), &
        '         : Cold (air temperature) stress affecting pollination - not considered'
    else
        write(fhandle, '(i6,a)') GetCrop_Tcold(), &
        '         : Minimum air temperature below which pollination starts to fail (cold stress) (degC)'
    end if
    if (int(GetCrop_Theat(), int32) == undef_int) then
        write(fhandle, '(i6,a)') GetCrop_Theat(), &
        '         : Heat (air temperature) stress affecting pollination - not considered'
    else
        write(fhandle, '(i6,a)') GetCrop_Theat(), &
        '         : Maximum air temperature above which pollination starts to fail (heat stress) (degC)'
    end if
    if (roundc(GetCrop_GDtranspLow(), mold=1) == undef_int) then
        write(fhandle, '(f8.1,a)') GetCrop_GDtranspLow(), &
        '       : Cold (air temperature) stress on crop transpiration not considered'
    else
        write(fhandle, '(f8.1,a)') GetCrop_GDtranspLow(), &
        '       : Minimum growing degrees required for full crop transpiration (degC - day)'
    end if

    ! salinity stress
    write(fhandle, '(i6,a)') GetCrop_ECemin(), &
    '         : Electrical Conductivity of soil saturation extract at which crop starts to be affected by soil salinity (dS/m)'
    write(fhandle, '(i6,a)') GetCrop_ECemax(), &
    '         : Electrical Conductivity of soil saturation extract at which crop can no longer grow (dS/m)'
    write(fhandle, '(a)') '    -9         : Dummy - no longer applicable' ! shape factor Ks(salt)-ECe
    write(fhandle, '(i6,a)') GetCrop_CCsaltDistortion(), &
    '         : Calibrated distortion (%) of CC due to salinity stress (Range: 0 (none) to +100 (very strong))'
    write(fhandle, '(i6,a)') GetCrop_ResponseECsw(), &
    '         : Calibrated response (%) of stomata stress to ECsw (Range: 0 (none) to +200 (extreme))'

    ! evapotranspiration
    write(fhandle, '(f9.2,a)') GetCrop_KcTop(),  &
    '      : Crop coefficient when canopy is complete but prior to senescence (KcTr,x)'
    write(fhandle, '(f10.3,a)') GetCrop_KcDecline(), &
    '     : Decline of crop coefficient (%/day) as a result of ageing, nitrogen deficiency, etc.'
    write(fhandle, '(f9.2,a)') GetCrop_RootMin(), &
    '      : Minimum effective rooting depth (m)'
    write(fhandle, '(f9.2,a)') GetCrop_RootMax(), &
    '      : Maximum effective rooting depth (m)'
    write(fhandle, '(i6,a)') GetCrop_RootShape(), &
    '         : Shape factor describing root zone expansion'
    write(fhandle, '(f10.3,a)') GetCrop_SmaxTopQuarter(), &
    '     : Maximum root water extraction (m3water/m3soil.day) in top quarter of root zone'
    write(fhandle, '(f10.3,a)') GetCrop_SmaxBotQuarter(), &
    '     : Maximum root water extraction (m3water/m3soil.day) in bottom quarter of root zone'
    write(fhandle, '(i6,a)') GetCrop_CCEffectEvapLate(), &
    '         : Effect of canopy cover in reducing soil evaporation in late season stage'

    ! canopy development
    write(fhandle, '(f9.2,a)') GetCrop_SizeSeedling(), &
    '      : Soil surface covered by an individual seedling at 90 % emergence (cm2)'
    write(fhandle, '(f9.2,a)') GetCrop_SizePlant(), &
    '      : Canopy size of individual plant (re-growth) at 1st day (cm2)'
    write(fhandle, '(i9,a)') GetCrop_PlantingDens(), &
    '      : Number of plants per hectare'
    write(fhandle, '(f12.5,a)') GetCrop_CGC(), &
    '   : Canopy growth coefficient (CGC): Increase in canopy cover (fraction soil cover per day)'
    if (GetCrop_YearCCx() == undef_int) then
        write(fhandle, '(i6,a)') GetCrop_YearCCx(), &
        '         : Number of years at which CCx declines to 90 % of its value due to self-thinning - Not Applicable'
    else
        write(fhandle, '(i6,a)') GetCrop_YearCCx(), &
        '         : Number of years at which CCx declines to 90 % of its value due to self-thinning - for Perennials'
    end if
    if (roundc(GetCrop_CCxRoot(), mold=1) == undef_int) then
        write(fhandle, '(f9.2,a)') GetCrop_CCxRoot(), &
        '      : Shape factor of the decline of CCx over the years due to self-thinning - Not Applicable'
    else
        write(fhandle, '(f9.2,a)') GetCrop_CCxRoot(), &
        '      : Shape factor of the decline of CCx over the years due to self-thinning - for Perennials'
    end if
    write(fhandle, '(a)') '    -9         : dummy - Parameter no Longer required'

    write(fhandle, '(f9.2,a)') GetCrop_CCx(), &
    '      : Maximum canopy cover (CCx) in fraction soil cover'
    write(fhandle, '(f12.5,a)') GetCrop_CDC(), &
    '   : Canopy decline coefficient (CDC): Decrease in canopy cover (in fraction per day)'
    if (GetCrop_Planting() == plant_Seed) then
        write(fhandle, '(i6,a)') GetCrop_DaysToGermination(), &
        '         : Calendar Days: from sowing to emergence'
        write(fhandle, '(i6,a)') GetCrop_DaysToMaxRooting(), &
        '         : Calendar Days: from sowing to maximum rooting depth'
        write(fhandle, '(i6,a)') GetCrop_DaysToSenescence(), &
        '         : Calendar Days: from sowing to start senescence'
        write(fhandle, '(i6,a)') GetCrop_DaysToHarvest(), &
        '         : Calendar Days: from sowing to maturity (length of crop cycle)'
        if (GetCrop_subkind() == subkind_Tuber) then
            write(fhandle, '(i6,a)') GetCrop_DaysToFlowering(), &
            '         : Calendar Days: from sowing to start of yield formation'
        else
            write(fhandle, '(i6,a)') GetCrop_DaysToFlowering(), &
            '         : Calendar Days: from sowing to flowering'
        end if
    else
        if (GetCrop_Planting() == plant_Transplant) then
            write(fhandle, '(i6,a)') GetCrop_DaysToGermination(), &
            '         : Calendar Days: from transplanting to recovered transplant'
            write(fhandle, '(i6,a)') GetCrop_DaysToMaxRooting(), &
            '         : Calendar Days: from transplanting to maximum rooting depth'
            write(fhandle, '(i6,a)') GetCrop_DaysToSenescence(), &
            '         : Calendar Days: from transplanting to start senescence'
            write(fhandle, '(i6,a)') GetCrop_DaysToHarvest(), &
            '         : Calendar Days: from transplanting to maturity'
            if (GetCrop_subkind() == subkind_Tuber) then
                write(fhandle, '(i6,a)') GetCrop_DaysToFlowering(), &
                '         : Calendar Days: from transplanting to start of yield formation'
            else
                write(fhandle, '(i6,a)') GetCrop_DaysToFlowering(), &
                '         : Calendar Days: from transplanting to flowering'
            end if
        else
            ! planting = regrowth
            write(fhandle, '(i6,a)') GetCrop_DaysToGermination(), &
            '         : Calendar Days: from regrowth to recovering'
            write(fhandle, '(i6,a)') GetCrop_DaysToMaxRooting(), &
            '         : Calendar Days: from regrowth to maximum rooting depth'
            write(fhandle, '(i6,a)') GetCrop_DaysToSenescence(), &
            '         : Calendar Days: from regrowth to start senescence'
            write(fhandle, '(i6,a)') GetCrop_DaysToHarvest(), &
            '         : Calendar Days: from regrowth to maturity'
            if (GetCrop_subkind() == subkind_Tuber) then
                write(fhandle, '(i6,a)') GetCrop_DaysToFlowering(), &
                '         : Calendar Days: from regrowth to start of yield formation'
            else
                write(fhandle, '(i6,a)') GetCrop_DaysToFlowering(), &
                '         : Calendar Days: from regrowth to flowering'
            end if
        end if
    end if
    write(fhandle, '(i6,a)') GetCrop_LengthFlowering(), &
    '         : Length of the flowering stage (days)'

    ! Crop.DeterminancyLinked
    if (GetCrop_DeterminancyLinked() .eqv. .true.) then
        i = 1
        TempString = '         : Crop determinancy linked with flowering'
    else
        i = 0
        TempString = '         : Crop determinancy unlinked with flowering'
    end if
    write(fhandle, '(i6,a)') i, TempString

    ! Potential excess of fruits (%)
    if ((GetCrop_subkind() == subkind_Vegetative) &
                        .or. (GetCrop_subkind() == subkind_Forage)) then
        write(fhandle, '(i6,a)') undef_int, &
        '         : parameter NO LONGER required' ! Building up of Harvest Index (% of growing cycle)')
    else
        if (GetCrop_fExcess() == undef_int) then
            TempString = '         : Excess of potential fruits - Not Applicable'
        else
            TempString = '         : Excess of potential fruits (%)'
        end if
        write(fhandle, '(i6, a)') GetCrop_fExcess(), TempString
    end if

    ! Building-up of Harvest Index
    if (GetCrop_DaysToHIo() == undef_int) then
        TempString = '         : Building up of Harvest Index - Not Applicable'
    else
        select case (GetCrop_subkind())
            case(subkind_Vegetative, subkind_Forage)
                TempString = '         : Building up of Harvest Index starting at sowing/transplanting (days)'
            case(subkind_Grain)
                TempString = '         : Building up of Harvest Index starting at flowering (days)'
            case(subkind_Tuber)
                TempString = '         : Building up of Harvest Index starting at root/tuber enlargement (days)'
            case default
                TempString = '         : Building up of Harvest Index during yield formation (days)'
        end select
    end if
    write(fhandle, '(i6, a)') GetCrop_DaysToHIo(), TempString

    ! yield response to water
    write(fhandle, '(f8.1,a)') GetCrop_WP(), &
    '       : Water Productivity normalized for ETo and CO2 (WP*) (gram/m2)'
    write(fhandle, '(i6,a)') GetCrop_WPy(), &
    '         : Water Productivity normalized for ETo and CO2 during yield formation (as % WP*)'
    write(fhandle, '(i6,a)') GetCrop_AdaptedToCO2(), &
    '         : Crop performance under elevated atmospheric CO2 concentration (%)'
    write(fhandle, '(i6,a)') GetCrop_HI(), &
    '         : Reference Harvest Index (HIo) (%)'
    if (GetCrop_subkind() == subkind_Tuber) then
        write(fhandle, '(i6,a)') GetCrop_HIincrease(), &
        '         : Possible increase (%) of HI due to water stress before start of yield formation'
    else
        write(fhandle, '(i6,a)') GetCrop_HIincrease(), &
        '         : Possible increase (%) of HI due to water stress before flowering'
    end if
    if (roundc(GetCrop_aCoeff(), mold=1) == undef_int) then
        write(fhandle, '(f8.1,a)') GetCrop_aCoeff(), &
            '       : No impact on HI of restricted vegetative growth during yield formation '
    else
        write(fhandle, '(f8.1,a)') GetCrop_aCoeff(), &
        '       : Coefficient describing positive impact on HI of restricted vegetative growth during yield formation'
    end if
    if (roundc(GetCrop_bCoeff(), mold=1) == undef_int) then
        write(fhandle, '(f8.1,a)') GetCrop_bCoeff(), &
        '       : No effect on HI of stomatal closure during yield formation'
    else
        write(fhandle, '(f8.1,a)') GetCrop_bCoeff(), &
        '       : Coefficient describing negative impact on HI of stomatal closure during yield formation'
    end if
    write(fhandle, '(i6,a)') GetCrop_DHImax(), &
    '         : Allowable maximum increase (%) of specified HI'

    ! growing degree days
    if (GetCrop_Planting() == plant_Seed) then
        write(fhandle, '(i6,a)') GetCrop_GDDaysToGermination(), &
        '         : GDDays: from sowing to emergence'
        write(fhandle, '(i6,a)') GetCrop_GDDaysToMaxRooting(), &
        '         : GDDays: from sowing to maximum rooting depth'
        write(fhandle, '(i6,a)') GetCrop_GDDaysToSenescence(), &
        '         : GDDays: from sowing to start senescence'
        write(fhandle, '(i6,a)') GetCrop_GDDaysToHarvest(), &
        '         : GDDays: from sowing to maturity (length of crop cycle)'
        if (GetCrop_subkind() == subkind_Tuber) then
            write(fhandle, '(i6,a)') GetCrop_GDDaysToFlowering(), &
            '         : GDDays: from sowing to start tuber formation'
        else
            write(fhandle, '(i6,a)') GetCrop_GDDaysToFlowering(), &
            '         : GDDays: from sowing to flowering'
        end if
    else
        if (GetCrop_Planting() == plant_Transplant) then
            write(fhandle, '(i6,a)') GetCrop_GDDaysToGermination(), &
            '         : GDDays: from transplanting to recovered transplant'
            write(fhandle, '(i6,a)') GetCrop_GDDaysToMaxRooting(), &
            '         : GDDays: from transplanting to maximum rooting depth'
            write(fhandle, '(i6,a)') GetCrop_GDDaysToSenescence(), &
            '         : GDDays: from transplanting to start senescence'
            write(fhandle, '(i6,a)') GetCrop_GDDaysToHarvest(), &
            '         : GDDays: from transplanting to maturity'
            if (GetCrop_subkind() == subkind_Tuber) then
                write(fhandle, '(i6,a)') GetCrop_GDDaysToFlowering(), &
                '         : GDDays: from transplanting to start yield formation'
            else
                write(fhandle, '(i6,a)') GetCrop_GDDaysToFlowering(), &
                '         : GDDays: from transplanting to flowering'
            end if
        else
            ! Crop.Planting = regrowth
            write(fhandle, '(i6,a)') GetCrop_GDDaysToGermination(), &
            '         : GDDays: from regrowth to recovering'
            write(fhandle, '(i6,a)') GetCrop_GDDaysToMaxRooting(), &
            '         : GDDays: from regrowth to maximum rooting depth'
            write(fhandle, '(i6,a)') GetCrop_GDDaysToSenescence(), &
            '         : GDDays: from regrowth to start senescence'
            write(fhandle, '(i6,a)') GetCrop_GDDaysToHarvest(), &
            '         : GDDays: from regrowth to maturity'
            if (GetCrop_subkind() == subkind_Tuber) then
                write(fhandle, '(i6,a)') GetCrop_GDDaysToFlowering(), &
                '         : GDDays: from regrowth to start yield formation'
            else
                write(fhandle, '(i6,a)') GetCrop_GDDaysToFlowering(), &
                '         : GDDays: from regrowth to flowering'
            end if
        end if
    end if
    write(fhandle, '(i6,a)') GetCrop_GDDLengthFlowering(), &
    '         : Length of the flowering stage (growing degree days)'
    write(fhandle, '(f13.6,a)') GetCrop_GDDCGC(), &
    '  : CGC for GGDays: Increase in canopy cover (in fraction soil cover per growing-degree day)'
    write(fhandle, '(f13.6,a)') GetCrop_GDDCDC(), &
    '  : CDC for GGDays: Decrease in canopy cover (in fraction per growing-degree day)'
    write(fhandle, '(i6,a)') GetCrop_GDDaysToHIo(), &
    '         : GDDays: building-up of Harvest Index during yield formation'

    ! added to 6.2
    write(fhandle, '(i6,a)') GetCrop_DryMatter(), &
    '         : dry matter content (%) of fresh yield'

    ! added to 7.0 - Perennial crops
    if (GetCrop_subkind() == subkind_Forage) then
        write(fhandle, '(f9.2,a)') GetCrop_RootMinYear1(), &
        '      : Minimum effective rooting depth (m) in first year (for perennials)'
    else
        write(fhandle, '(f9.2,a)') GetCrop_RootMinYear1(), &
        '      : Minimum effective rooting depth (m) in first year - required only in case of regrowth'
    end if
    if (GetCrop_SownYear1() .eqv. .true.) then
        i = 1
        if (GetCrop_subkind() == subkind_Forage) then
            write(fhandle, '(i6,a)') i, &
            '         : Crop is sown in 1st year (for perennials)'
        else
            write(fhandle, '(i6,a)') i, &
            '         : Crop is sown in 1st year - required only in case of regrowth'
        end if
    else
        i = 0
        if (GetCrop_subkind() == subkind_Forage) then
            write(fhandle, '(i6,a)') i, &
            '         : Crop is transplanted in 1st year (for perennials)'
        else
            write(fhandle, '(i6,a)') i, &
            '         : Crop is transplanted in 1st year - required only in case of regrowth'
        end if
    end if

    ! added to 7.0 - Assimilates
    if (GetCrop_Assimilates_On() .eqv. .false.) then
        i = 0
        write(fhandle, '(i6,a)') i, &
        '         : Transfer of assimilates from above ground parts to root system is NOT considered'
        write(fhandle, '(i6,a)') i, &
        '         : Number of days at end of season during which assimilates are stored in root system'
        write(fhandle, '(i6,a)') i, &
        '         : Percentage of assimilates transferred to root system at last day of season'
        write(fhandle, '(i6,a)') i, &
        '         : Percentage of stored assimilates transferred to above ground parts in next season'
    else
        i = 1
        write(fhandle, '(i6,a)') i, &
        '         : Transfer of assimilates from above ground parts to root system is considered'
        write(fhandle, '(i6,a)') GetCrop_Assimilates_Period(), &
        '         : Number of days at end of season during which assimilates are stored in root system'
        write(fhandle, '(i6,a)') GetCrop_Assimilates_Stored(), &
        '         : Percentage of assimilates transferred to root system at last day of season'
        write(fhandle, '(i6,a)') GetCrop_Assimilates_Mobilized(), &
        '         : Percentage of stored assimilates transferred to above ground parts in next season'
    end if
    close(fhandle)

    ! maximum rooting depth in given soil profile
    call SetSoil_RootMax(RootMaxInSoilProfile(GetCrop_RootMax(), &
                                    GetSoil_NrSoilLayers(), GetSoilLayer()))

    ! copy to CropFileSet
    call SetCropFileSet_DaysFromSenescenceToEnd(GetCrop_DaysToHarvest() &
                                            - GetCrop_DaysToSenescence())
    call SetCropFileSet_DaysToHarvest(GetCrop_DaysToHarvest())

    if (GetCrop_ModeCycle() == ModeCycle_GDDays) then
        call SetCropFileSet_GDDaysFromSenescenceToEnd( &
                GetCrop_GDDaysToHarvest() - GetCrop_GDDaysToSenescence())
        call SetCropFileSet_GDDaysToHarvest(GetCrop_GDDaysToHarvest())
    else
        call SetCropFileSet_GDDaysFromSenescenceToEnd(undef_int)
        call SetCropFileSet_GDDaysToHarvest(undef_int)
    end if
end subroutine SaveCrop


subroutine SaveProfile(totalname)

    character(len=*), intent(in) :: totalname

    integer :: fhandle
    integer(int32) :: i

    open(newunit=fhandle, file=trim(totalname), status='replace', action='write')
    write(fhandle, '(a)') GetProfDescription()

    write(fhandle, &
          '("        ", a, "                 : AquaCrop Version (", a, ")")') &
          GetVersionString(), GetReleaseDate()    ! AquaCrop version
    write(fhandle, '(i9, a)') GetSoil_CNvalue(), &
    '                   : CN (Curve Number)'
    write(fhandle, '(i9, a)') GetSoil_REW(), &
    '                   : Readily evaporable water from top layer (mm)'
    write(fhandle, '(i9, a)') GetSoil_NrSoilLayers(), &
    '                   : number of soil horizons'
    write(fhandle, '(i9, a)') undef_int, &
    '                   : variable no longer applicable'
    write(fhandle, '(a)') &
    '  Thickness  Sat   FC    WP     Ksat   Penetrability  Gravel  CRa       CRb           description'
    write(fhandle, '(a)') &
    '  ---(m)-   ----(vol %)-----  (mm/day)      (%)        (%)    -----------------------------------------'
    do i = 1, GetSoil_NrSoilLayers()
        write(fhandle, '(f8.2, f8.1, f6.1, f6.1, f8.1, i11, i10, f14.6, f10.6, a, a)') &
                            GetSoilLayer_Thickness(i), GetSoilLayer_SAT(i), &
                            GetSoilLayer_FC(i), GetSoilLayer_WP(i), &
                            GetSoilLayer_InfRate(i), GetSoilLayer_Penetrability(i), &
                            GetSoilLayer_GravelMass(i), GetSoilLayer_CRa(i), &
                            GetSoilLayer_CRb(i), '   ', trim(GetSoilLayer_Description(i))
    end do
    close(fhandle)

    ! maximum rooting depth in  soil profile for given crop
    call SetSoil_RootMax(RootMaxInSoilProfile(GetCrop_RootMax(), &
                            GetSoil_NrSoilLayers(), GetSoilLayer()))
end subroutine SaveProfile


subroutine DetermineParametersCR(SoilClass, KsatMM, aParam, bParam)
    integer(int8), intent(in) :: SoilClass
    real(sp), intent(in) :: KsatMM
    real(sp), intent(inout) :: aParam
    real(sp), intent(inout) :: bParam

    ! determine parameters
    if (roundc(KsatMM*1000, mold=1) <= 0) then
        aParam = undef_int
        bParam = undef_int
    else
        select case (SoilClass)
            case(1) ! sandy soils
                aParam = -0.3112_sp - KsatMM/100000._sp
                bParam = -1.4936_sp + 0.2416_sp*log(KsatMM)
            case(2) ! loamy soils
                aParam = -0.4986_sp + 9._sp*KsatMM/100000._sp
                bParam = -2.1320_sp + 0.4778_sp*log(KsatMM)
            case(3) ! sandy clayey soils
                aParam = -0.5677_sp - 4._sp*KsatMM/100000._sp
                bParam = -3.7189_sp + 0.5922_sp*log(KsatMM)
            case default ! silty clayey soils
            aParam = -0.6366_sp + 8._sp*KsatMM/10000._sp
            bParam = -1.9165_sp + 0.7063_sp*log(KsatMM)
        end select
    end if
end subroutine DetermineParametersCR


subroutine DetermineNrandThicknessCompartments()

    real(sp) :: TotalDepthL, TotalDepthC, DeltaZ
    integer(int32) :: i

    TotalDepthL = 0._sp
    do i = 1, GetSoil_NrSoilLayers()
        TotalDepthL = TotalDepthL + GetSoilLayer_Thickness(i)
    end do
    TotalDepthC = 0._sp
    call SetNrCompartments(0)
    loop: do
        DeltaZ = (TotalDepthL - TotalDepthC)
        call SetNrCompartments(GetNrCompartments() + 1)
        if (DeltaZ > GetSimulParam_CompDefThick()) then
            call SetCompartment_Thickness(GetNrCompartments(), GetSimulParam_CompDefThick())
        else
            call SetCompartment_Thickness(GetNrCompartments(), DeltaZ)
        end if
        TotalDepthC = TotalDepthC + GetCompartment_Thickness(GetNrCompartments())
        if ((GetNrCompartments() == max_No_compartments) &
                .or. (abs(TotalDepthC - TotalDepthL) < 0.0001_sp)) exit loop
        end do loop
end subroutine DetermineNrandThicknessCompartments


subroutine DetermineRootZoneSaltContent(RootingDepth, ZrECe, ZrECsw, ZrECswFC, ZrKsSalt)
    real(sp), intent(in) :: RootingDepth
    real(sp), intent(inout) :: ZrECe
    real(sp), intent(inout) :: ZrECsw
    real(sp), intent(inout) :: ZrECswFC
    real(sp), intent(inout) :: ZrKsSalt


    real(sp) :: CumDepth, Factor, frac_value
    integer(int32) :: compi

    CumDepth = 0._sp
    compi = 0._sp
    ZrECe = 0._sp
    ZrECsw = 0._sp
    ZrECswFC = 0._sp
    ZrKsSalt = 1._sp
    if (RootingDepth >= GetCrop_RootMin()) then
        loop: do
            compi = compi + 1
            CumDepth = CumDepth + GetCompartment_Thickness(compi)
            if (CumDepth <= RootingDepth) then
                Factor = 1._sp
            else
                frac_value = RootingDepth - (CumDepth - GetCompartment_Thickness(compi))
                if (frac_value > 0._sp) then
                    Factor = frac_value/GetCompartment_Thickness(compi)
                else
                    Factor = 0._sp
                end if
            end if
            Factor = Factor * (GetCompartment_Thickness(compi))/RootingDepth ! weighting factor
            ZrECe = ZrECe + Factor * ECeComp(GetCompartment_i(compi))
            ZrECsw = ZrECsw + Factor * ECswComp(GetCompartment_i(compi), (.false.)) ! not at FC
            ZrECswFC = ZrECswFC + Factor * ECswComp(GetCompartment_i(compi), (.true.)) ! at FC
            if ((CumDepth >= RootingDepth) .or. (compi == NrCompartments)) exit loop
        end do loop
        if (((GetCrop_ECemin() /= undef_int) .and. (GetCrop_ECemax() /= undef_int)) .and. &
                                    (GetCrop_ECemin() < GetCrop_ECemax())) then
            ZrKsSalt = KsSalinity((.true.), GetCrop_ECemin(), GetCrop_ECemax(), ZrECe, (0.0_sp))
        else
            ZrKsSalt = KsSalinity((.false.), GetCrop_ECemin(), GetCrop_ECemax(), ZrECe, (0.0_sp))
        end if
    else
        ZrECe = undef_int
        ZrECsw = undef_int
        ZrECswFC = undef_int
        ZrKsSalt = undef_int
    end if
end subroutine DetermineRootZoneSaltContent


subroutine CalculateAdjustedFC(DepthAquifer, CompartAdj)
    real(sp), intent(in) :: DepthAquifer
    type(CompartmentIndividual), dimension(max_No_compartments), &
                                    intent(inout) :: CompartAdj

    integer(int32) :: compi, ic
    real(sp) :: Zi, Depth, DeltaV, DeltaFC, Xmax

    Depth = 0._sp
    do compi = 1, GetNrCompartments()
        Depth = Depth + CompartAdj(compi)%Thickness
    end do

    compi = GetNrCompartments()

    loop: do
        Zi = Depth - CompartAdj(compi)%Thickness/2._sp
        Xmax = NoAdjustment(GetSoilLayer_FC(CompartAdj(compi)%Layer))

        if ((DepthAquifer < 0._sp) &
            .or. ((DepthAquifer - Zi) >= Xmax)) then
            do ic = 1, compi
                CompartAdj(ic)%FCadj = GetSoilLayer_FC(CompartAdj(ic)%Layer)
            end do

            compi = 0
        else
            if (GetSoilLayer_FC(CompartAdj(compi)%Layer) &
                    >= GetSoilLayer_SAT(CompartAdj(compi)%Layer)) then
                CompartAdj(compi)%FCadj = GetSoilLayer_FC(CompartAdj(compi)%Layer)
            else
                if (Zi >= DepthAquifer) then
                    CompartAdj(compi)%FCadj = GetSoilLayer_SAT(CompartAdj(compi)%Layer)
                else
                    DeltaV = GetSoilLayer_SAT(CompartAdj(compi)%Layer) &
                             - GetSoilLayer_FC(CompartAdj(compi)%Layer)
                    DeltaFC = (DeltaV/(Xmax**2)) &
                              * (Zi - (DepthAquifer - Xmax))**2
                    CompartAdj(compi)%FCadj = GetSoilLayer_FC(CompartAdj(compi)%Layer) &
                                              + DeltaFC
                end if
            end if

            Depth = Depth - CompartAdj(compi)%Thickness
            compi = compi - 1
        end if

        if (compi < 1) exit loop
    end do loop


    contains


    real(sp) function NoAdjustment(FCvolPr)
        real(sp), intent(in) :: FCvolPr

        real(sp) :: pF

        if (FCvolPr <= 10._sp) then
            NoAdjustment = 1._sp
        else
            if (FCvolPr >= 30._sp) then
                NoAdjustment = 2._sp
            else
                pF = 2._sp + 0.3_sp * (FCvolPr-10._sp)/20._sp
                NoAdjustment = (exp(pF*log(10._sp)))/100._sp
            end if
        end if
    end function NoAdjustment
end subroutine CalculateAdjustedFC


subroutine AdjustOnsetSearchPeriod()

    integer(int32) :: temp_Integer

    if (GetClimFile() == '(None)') then
        call SetOnset_StartSearchDayNr(1)
        call SetOnset_StopSearchDayNr(GetOnset_StartSearchDayNr() &
               + GetOnset_LengthSearchPeriod() - 1)
    else
        temp_Integer = GetOnset_StartSearchDayNr()
        call DetermineDayNr((1), (1), GetSimulation_YearStartCropCycle(), &
                              temp_Integer) ! 1 January
        call SetOnset_StartSearchDayNr(temp_Integer)
        if (GetOnset_StartSearchDayNr() < GetClimRecord_FromDayNr()) then
            call SetOnset_StartSearchDayNr(GetClimRecord_FromDayNr())
        end if
        call SetOnset_StopSearchDayNr(GetOnset_StartSearchDayNr() &
                        + GetOnset_LengthSearchPeriod() - 1)
        if (GetOnset_StopSearchDayNr() > GetClimRecord_ToDayNr()) then
            call SetOnset_StopSearchDayNr(GetClimRecord_ToDayNr())
            call SetOnset_LengthSearchPeriod(GetOnset_StopSearchDayNr() &
                     - GetOnset_StartSearchDayNr() + 1)
        end if
    end if
end subroutine AdjustOnsetSearchPeriod


integer(int32) function ActiveCells(Comp)
    type(CompartmentIndividual), intent(in) :: Comp

    integer(int32) ::  celi

    if (Comp%theta <= GetSoilLayer_UL(Comp%Layer)) then
        celi = 1
        do while (Comp%theta > (GetSoilLayer_Dx(Comp%Layer)) * celi)
            celi = celi + 1
        end do
    else
        celi = GetSoilLayer_SCP1(Comp%Layer)
    end if
    ActiveCells = celi
end function ActiveCells


subroutine DetermineSaltContent(ECe, Comp)
    real(sp), intent(in) :: ECe
    type(CompartmentIndividual), intent(inout) :: Comp

    real(sp) :: TotSalt, SumDF, SAT, UL, Dx, mm, mm1, mmN
    integer(int32) :: celn, i

    TotSalt = ECe*Equiv*(GetSoilLayer_SAT(Comp%Layer))*10._sp*Comp%Thickness
    celn = ActiveCells(Comp)
    SAT = (GetSoilLayer_SAT(Comp%Layer))/100._sp  ! m3/m3
    UL = GetSoilLayer_UL(Comp%Layer) ! m3/m3   ! Upper limit of SC salt cel
    Dx = GetSoilLayer_Dx(Comp%Layer)  ! m3/m3  ! Size of salts cel (expect last one)
    mm1 = Dx*1000._sp*Comp%Thickness &
                    * (1._sp - GetSoilLayer_GravelVol(Comp%Layer)/100._sp)
                     ! g/l ! volume [mm]=[l/m2] of cells
    mmN = (SAT-UL)*1000._sp*Comp%Thickness &
                    * (1._sp - GetSoilLayer_GravelVol(Comp%Layer)/100._sp)
                     ! g/l ! volume [mm]=[l/m2] of last cell
    SumDF = 0._sp
    do i = 1, GetSoilLayer_SCP1(Comp%Layer)
        Comp%Salt(i) = 0._sp
        Comp%Depo(i) = 0._sp
    end do
    do i = 1, celn
        SumDF = SumDF + GetSoilLayer_SaltMobility_i(Comp%Layer, i)
    end do
    do i = 1, celn
        Comp%Salt(i) = TotSalt * GetSoilLayer_SaltMobility_i(Comp%Layer, i)/SumDF
        mm = mm1
        if (i == GetSoilLayer_SCP1(Comp%Layer)) then
            mm = mmN
        end if
        call SaltSolutionDeposit(mm, Comp%Salt(i), Comp%Depo(i))
    end do
end subroutine DetermineSaltContent


subroutine SetClimData()

    type(rep_clim) :: SetARecord, SetBRecord
    integer(int32) :: tmptoD, tmpToM, tmpToY
    integer(int32) :: tmpFromD, tmpFromM, tmpFromY
    character(len=:), allocatable :: tmpstr

    call SetClimRecord_NrObs(999) ! (heeft geen belang)

    ! Part A - ETo and Rain files --> ClimFile
    if ((GetEToFile() == '(None)') .and. (GetRainFile() == '(None)')) then
        call SetClimFile('(None)')
        call SetClimDescription('Specify Climatic data when Running AquaCrop')
        call SetClimRecord_DataType(datatype_daily)
        call SetClimRecord_FromString('any date')
        call SetClimRecord_ToString('any date')
        call SetClimRecord_FromY(1901)
    else
        call SetClimFile('EToRainTempFile')
        call SetClimDescription('Read ETo/RAIN/TEMP data set')
        if (GetEToFile() == '(None)') then
            call SetClimRecord_FromY(GetRainRecord_FromY())
            call SetClimRecord_FromDayNr(GetRainRecord_FromDayNr())
            call SetClimRecord_ToDayNr(GetRainRecord_ToDayNr())
            call SetClimRecord_FromString(GetRainRecord_FromString())
            call SetClimRecord_ToString(GetRainRecord_ToString())
            if (FullUndefinedRecord(GetRainRecord_FromY(), &
                    GetRainRecord_FromD(), GetRainRecord_FromM(),&
                    GetRainRecord_ToD(), GetRainRecord_ToM())) then
                call SetClimRecord_NrObs(365)
            end if
        end if
        if (GetRainFile() == '(None)') then
            call SetClimRecord_FromY(GetEToRecord_FromY())
            call SetClimRecord_FromDayNr(GetEToRecord_FromDayNr())
            call SetClimRecord_ToDayNr(GetEToRecord_ToDayNr())
            call SetClimRecord_FromString(GetEToRecord_FromString())
            call SetClimRecord_ToString(GetEToRecord_ToString())
            if (FullUndefinedRecord(GetEToRecord_FromY(), &
                    GetEToRecord_FromD(), GetEToRecord_FromM(), &
                    GetEToRecord_ToD(), GetEToRecord_ToM())) then
                call SetClimRecord_NrObs(365)
            end if
        end if

        if ((GetEToFile() /= '(None)') .and. (GetRainFile() /= '(None)')) then
            SetARecord = GetEToRecord()
            SetBRecord = GetRainRecord()
            if (((GetEToRecord_FromY() == 1901) &
                 .and. FullUndefinedRecord(GetEToRecord_FromY(),&
                         GetEToRecord_FromD(), GetEToRecord_FromM(), &
                         GetEToRecord_ToD(), GetEToRecord_ToM())) &
                .and. ((GetRainRecord_FromY() == 1901) &
                  .and. FullUndefinedRecord(GetRainRecord_FromY(),&
                          GetRainRecord_FromD(), GetRainRecord_FromM(), &
                          GetRainRecord_ToD(), GetRainRecord_ToM()))) then
                call SetClimRecord_NrObs(365)
            end if

           if ((GetEToRecord_FromY() == 1901) .and. &
               (GetRainRecord_FromY() /= 1901)) then
                ! Jaartal van RainRecord ---> SetARecord (= EToRecord)
                ! FromY + adjust FromDayNr and FromString
                SetARecord%FromY = GetRainRecord_FromY()
                call DetermineDayNr(GetEToRecord_FromD(), GetEToRecord_FromM(), &
                               SetARecord%FromY, SetARecord%FromDayNr)
                if (((SetARecord%FromDayNr < GetRainRecord_FromDayNr())) &
                   .and. (GetRainRecord_FromY() < GetRainRecord_ToY())) then
                    SetARecord%FromY = GetRainRecord_FromY() + 1
                    call DetermineDayNr(GetEToRecord_FromD(), GetEToRecord_FromM(),&
                          SetARecord%FromY, SetARecord%FromDayNr)
                end if
                call SetClimRecord_FromY(SetARecord%FromY)
                ! nodig voor DayString (werkt met ClimRecord)
                SetARecord%FromString = DayString(SetARecord%FromDayNr)
                ! ToY + adjust ToDayNr and ToString
                if (FullUndefinedRecord(GetEToRecord_FromY(),&
                        GetEToRecord_FromD(), &
                        GetEToRecord_FromM(), GetEToRecord_ToD(),&
                        GetEToRecord_ToM())) then
                    SetARecord%ToY = GetRainRecord_ToY()
                else
                    SetARecord%ToY = SetARecord%FromY
                end if
                call DetermineDayNr(GetEToRecord_ToD(), GetEToRecord_ToM(), &
                               SetARecord%ToY, SetARecord%ToDayNr)
                SetARecord%ToString = DayString(SetARecord%ToDayNr)
            end if

            if ((GetEToRecord_FromY() /= 1901) .and. &
                (GetRainRecord_FromY() == 1901)) then
                ! Jaartal van EToRecord ---> SetBRecord (= RainRecord)
                ! FromY + adjust FromDayNr and FromString
                SetBRecord%FromY = GetEToRecord_FromY()
                call DetermineDayNr(GetRainRecord_FromD(), GetRainRecord_FromM(),&
                               SetBRecord%FromY, SetBRecord%FromDayNr)
                if (((SetBRecord%FromDayNr < GetEToRecord_FromDayNr())) &
                    .and. (GetEToRecord_FromY() < GetEToRecord_ToY())) then
                    SetBRecord%FromY = GetEToRecord_FromY() + 1
                    call DetermineDayNr(GetRainRecord_FromD(),&
                                   GetRainRecord_FromM(),&
                                   SetBRecord%FromY, SetBRecord%FromDayNr)
                end if
                call SetClimRecord_FromY(SetBRecord%FromY)
                ! nodig voor DayString (werkt met ClimRecord)
                SetBRecord%FromString = DayString(SetBRecord%FromDayNr)
                ! ToY + adjust ToDayNr and ToString
                if (FullUndefinedRecord(GetRainRecord_FromY(), &
                      GetRainRecord_FromD(), &
                      GetRainRecord_FromM(), GetRainRecord_ToD(),&
                      GetRainRecord_ToM())) then
                    SetBRecord%ToY = GetEToRecord_ToY()
                else
                    SetBRecord%ToY = SetBRecord%FromY
                end if
                call DetermineDayNr(GetRainRecord_ToD(), GetRainRecord_ToM(), &
                               SetBRecord%ToY, SetBRecord%ToDayNr)
                SetBRecord%ToString = DayString(SetBRecord%ToDayNr)
            end if
            ! bepaal characteristieken van ClimRecord
            call SetClimRecord_FromY(SetARecord%FromY)
            call SetClimRecord_FromDayNr(SetARecord%FromDayNr)
            tmpstr = SetARecord%FromString
            call SetClimRecord_FromString(tmpstr)
            if (GetClimRecord_FromDayNr() < SetBRecord%FromDayNr) then
                call SetClimRecord_FromY(SetBRecord%FromY)
                call SetClimRecord_FromDayNr(SetBRecord%FromDayNr)
                call SetClimRecord_FromString(SetBRecord%FromString)
            end if
            call SetClimRecord_ToDayNr(SetARecord%ToDayNr)
            tmpstr = SetARecord%ToString
            call SetClimRecord_ToString(tmpstr)
            if (GetClimRecord_ToDayNr() > SetBRecord%ToDayNr) then
                call SetClimRecord_ToDayNr(SetBRecord%ToDayNr)
                call SetClimRecord_ToString(SetBRecord%ToString)
            end if
            if (GetClimRecord_ToDayNr() < GetClimRecord_FromDayNr()) then
                call SetClimFile('(None)')
                call SetClimDescription(&
                       'ETo data set <--NO OVERLAP--> RAIN data set')
                call SetClimRecord_NrObs(0)
                call SetClimRecord_FromY(1901)
            end if
        end if
    end if

    ! Part B - ClimFile and Temperature files --> ClimFile
    if ((GetTemperatureFile() == '(None)') .or.&
        (GetTemperatureFile() == '(External)')) then
        ! no adjustments are required
    else
        if (GetClimFile() == '(None)') then
            call SetClimFile('EToRainTempFile')
            call SetClimDescription('Read ETo/RAIN/TEMP data set')
            call SetClimRecord_FromY(GetTemperatureRecord_FromY())
            call SetClimRecord_FromDayNr(GetTemperatureRecord_FromDayNr())
            call SetClimRecord_ToDayNr(GetTemperatureRecord_ToDayNr())
            call SetClimRecord_FromString(GetTemperatureRecord_FromString())
            call SetClimRecord_ToString(GetTemperatureRecord_ToString())
            if ((GetTemperatureRecord_FromY() == 1901) .and.&
                 FullUndefinedRecord(GetTemperatureRecord_FromY(),&
                     GetTemperatureRecord_FromD(),&
                     GetTemperatureRecord_FromM(),&
                     GetTemperatureRecord_ToD(), GetTemperatureRecord_ToM())) then
                call SetClimRecord_NrObs(365)
            else
                call SetClimRecord_NrObs(GetTemperatureRecord_ToDayNr() -&
                       GetTemperatureRecord_FromDayNr() + 1)
            end if
        else
            call DetermineDate(GetClimRecord_FromDayNr(), &
                       tmpFromD, tmpFromM, tmpFromY)
            call SetClimRecord_FromD(tmpFromD)
            call SetClimRecord_FromM(tmpFromM)
            call SetClimRecord_FromY(tmpFromY)
            call DetermineDate(GetClimRecord_ToDayNr(), tmpToD, tmpToM, tmpToY)
            call SetClimRecord_ToD(tmpToD)
            call SetClimRecord_ToM(tmpToM)
            call SetClimRecord_ToY(tmpToY)
            SetARecord = GetClimRecord()
            SetBRecord = GetTemperatureRecord()

            if ((GetClimRecord_FromY() == 1901) .and.&
                (GetTemperatureRecord_FromY() == 1901) &
                .and. (GetClimRecord_NrObs() == 365) &
                .and. FullUndefinedRecord(GetTemperatureRecord_FromY(), &
                           GetTemperatureRecord_FromD(),&
                           GetTemperatureRecord_FromM(), &
                           GetTemperatureRecord_ToD(),&
                           GetTemperatureRecord_ToM())) then
                call SetClimRecord_NrObs(365)
            else
                call SetClimRecord_NrObs(GetTemperatureRecord_ToDayNr() - &
                       GetTemperatureRecord_FromDayNr() + 1)
            end if
            if ((GetClimRecord_FromY() == 1901) .and.&
                (GetTemperatureRecord_FromY() /= 1901)) then
                ! Jaartal van TemperatureRecord ---> SetARecord (= ClimRecord)
                ! FromY + adjust FromDayNr and FromString
                SetARecord%FromY = GetTemperatureRecord_FromY()
                call DetermineDayNr(GetClimRecord_FromD(),&
                                    GetClimRecord_FromM(), &
                                    SetARecord%FromY, SetARecord%FromDayNr)
                if ((SetARecord%FromDayNr < GetTemperatureRecord_FromDayNr())&
                   .and. (GetTemperatureRecord_FromY() < GetTemperatureRecord_ToY())) then
                    SetARecord%FromY = GetTemperatureRecord_FromY() + 1
                    call DetermineDayNr(GetClimRecord_FromD(), GetClimRecord_FromM(),&
                                   SetARecord%FromY, SetARecord%FromDayNr)
                end if
                SetARecord%FromString = DayString(SetARecord%FromDayNr)
                ! ToY + adjust ToDayNr and ToString
                if (FullUndefinedRecord(GetClimRecord_FromY(),&
                     GetClimRecord_FromD(), &
                     GetClimRecord_FromM(), GetClimRecord_ToD(),&
                     GetClimRecord_ToM())) then
                    SetARecord%ToY = GetTemperatureRecord_ToY()
                else
                    SetARecord%ToY = SetARecord%FromY
                end if
                call DetermineDayNr(GetClimRecord_ToD(), GetClimRecord_ToM(), &
                          SetARecord%ToY, SetARecord%ToDayNr)
                SetARecord%ToString = DayString(SetARecord%ToDayNr)
            end if

            if ((GetClimRecord_FromY() /= 1901) .and.&
                (GetTemperatureRecord_FromY() == 1901)) then
                ! Jaartal van ClimRecord ---> SetBRecord (=
                ! GetTemperatureRecord())
                ! FromY + adjust FromDayNr and FromString
                SetBRecord%FromY = GetClimRecord_FromY()
                call DetermineDayNr(GetTemperatureRecord_FromD(), &
                       GetTemperatureRecord_FromM(), SetBRecord%FromY,&
                       SetBRecord%FromDayNr)
                if (((SetBRecord%FromDayNr < GetClimRecord_FromDayNr())) .and.&
                     (GetClimRecord_FromY() < GetClimRecord_ToY())) then
                    SetBRecord%FromY = GetClimRecord_FromY() + 1
                    call DetermineDayNr(GetTemperatureRecord_FromD(), &
                           GetTemperatureRecord_FromM(), SetBRecord%FromY,&
                           SetBRecord%FromDayNr)
                end if
                ! SetClimRecord_FromY(SetBRecord.FromY); ! nodig voor DayString
                ! (werkt met ClimRecord)
                SetBRecord%FromString = DayString(SetBRecord%FromDayNr)
                ! ToY + adjust ToDayNr and ToString
                if (FullUndefinedRecord(GetTemperatureRecord_FromY(),&
                      GetTemperatureRecord_FromD(),&
                      GetTemperatureRecord_FromM(), &
                      GetTemperatureRecord_ToD(), &
                      GetTemperatureRecord_ToM())) then
                    SetBRecord%ToY = GetClimRecord_ToY()
                else
                    SetBRecord%ToY = SetBRecord%FromY
                end if
                call DetermineDayNr(GetTemperatureRecord_ToD(), &
                       GetTemperatureRecord_ToM(), SetBRecord%ToY, &
                       SetBRecord%ToDayNr)
                SetBRecord%ToString = DayString(SetBRecord%ToDayNr)
            end if

            ! bepaal nieuwe characteristieken van ClimRecord
            call SetClimRecord_FromY(SetARecord%FromY)
            call SetClimRecord_FromDayNr(SetARecord%FromDayNr)
            call SetClimRecord_FromString(SetARecord%FromString)
            if (GetClimRecord_FromDayNr() < SetBRecord%FromDayNr) then
                call SetClimRecord_FromY(SetBRecord%FromY)
                call SetClimRecord_FromDayNr(SetBRecord%FromDayNr)
                call SetClimRecord_FromString(SetBRecord%FromString)
            end if
            call SetClimRecord_ToDayNr(SetARecord%ToDayNr)
            call SetClimRecord_ToString(SetARecord%ToString)
            if (GetClimRecord_ToDayNr() > SetBRecord%ToDayNr) then
                call SetClimRecord_ToDayNr(SetBRecord%ToDayNr)
                call SetClimRecord_ToString(SetBRecord%ToString)
            end if
            if (GetClimRecord_ToDayNr() < GetClimRecord_FromDayNr()) then
                call SetClimFile('(None)')
                call SetClimDescription(&
                        'Clim data <--NO OVERLAP--> TEMPERATURE data')
                call SetClimRecord_NrObs(0)
                call SetClimRecord_FromY(1901)
            end if
        end if
    end if
end subroutine SetClimData


character(len=17) function DayString(DNr)
    integer(int32), intent(in) :: DNr

    integer(int32) :: dayi, monthi, yeari
    character(2) :: strA
    character(len=:), allocatable :: strB
    integer(int32) :: DNr_t

    DNr_t = DNr
    if (GetClimFile() == '(None)') then
        do while (DNr_t > 365)
            DNr_t = DNr_t - 365
        end do
    end if
    call DetermineDate(DNr_t, dayi, monthi, yeari)
    write(strA, '(i2)') dayi
    if (GetClimRecord_FromY() == 1901) then
        strB = ''
    else
        write(strB, '(i4)') yeari
    end if
    strB = trim(strA)//' '//trim(NameMonth(monthi))//' '//trim(strB)
    do while (len(strB) < 17)
        strB = strB//' '
    end do
    DayString = strB
end function DayString


subroutine AdjustYearPerennials(TheYearSeason, Sown1stYear, TheCycleMode, &
          Zmax, ZminYear1, TheCCo, TheSizeSeedling, TheCGC, TheCCx, TheGDDCGC, &
          ThePlantingDens, TypeOfPlanting, Zmin, TheSizePlant, TheCCini,&
          TheDaysToCCini, TheGDDaysToCCini)
    integer(int8), intent(in) :: TheYearSeason
    logical, intent(in) :: Sown1stYear
    integer(intEnum), intent(in) :: TheCycleMode
    real(sp), intent(in) :: Zmax
    real(sp), intent(in) :: ZminYear1
    real(sp), intent(in) :: TheCCo
    real(sp), intent(in) :: TheSizeSeedling
    real(sp), intent(in) :: TheCGC
    real(sp), intent(in) :: TheCCx
    real(sp), intent(in) :: TheGDDCGC
    integer(int32), intent(in) :: ThePlantingDens
    integer(intEnum), intent(inout) :: TypeOfPlanting
    real(sp), intent(inout) :: Zmin
    real(sp), intent(inout) :: TheSizePlant
    real(sp), intent(inout) :: TheCCini
    integer(int32), intent(inout) :: TheDaysToCCini
    integer(int32), intent(inout) :: TheGDDaysToCCini

    if (TheYearSeason == 1) then
        if (Sown1stYear .eqv. .true.) then ! planting
            TypeOfPlanting = plant_seed
        else
            TypeOfPlanting = plant_transplant
        end if
        Zmin = ZminYear1  ! rooting depth
    else
        TypeOfPlanting = plant_regrowth ! planting
        Zmin = Zmax  ! rooting depth
        ! plant size by regrowth
        if (roundc(100._sp*TheSizePlant,mold=1_int32) < &
            roundc(100._sp*TheSizeSeedling,mold=1_int32)) then
            TheSizePlant = 10._sp * TheSizeSeedling
        end if
        if (roundc(100._sp*TheSizePlant,mold=1_int32) > &
            roundc((100._sp*TheCCx*10000._sp)/&
                   (ThePlantingDens/10000._sp),mold=1_int32)) then
            TheSizePlant = (TheCCx*10000._sp)/(ThePlantingDens/10000._sp)
            ! adjust size plant to maximum possible
        end if
    end if
    TheCCini = (ThePlantingDens/10000._sp) * (TheSizePlant/10000._sp)
    TheDaysToCCini = TimeToCCini(TypeOfPlanting, ThePlantingDens, &
                       TheSizeSeedling, TheSizePlant, TheCCx, TheCGC)
    if (TheCycleMode == modeCycle_GDDays) then
        TheGDDaysToCCini = TimeToCCini(TypeOfPlanting, ThePlantingDens, &
                       TheSizeSeedling, TheSizePlant, TheCCx, TheGDDCGC)
    else
        TheGDDaysToCCini = undef_int
    end if

end subroutine AdjustYearPerennials


subroutine NoCropCalendar()
    call SetCalendarFile('(None)')
    call SetCalendarFileFull(GetCalendarFile())  ! no file
    call SetCalendarDescription('')
    call SetOnset_GenerateOn(.false.)
    call SetOnset_GenerateTempOn(.false.)
    call SetEndSeason_GenerateTempOn(.false.)
    call SetCalendarDescription('No calendar for the Seeding/Planting year')
end subroutine NoCropCalendar


subroutine DetermineLinkedSimDay1(CropDay1, SimDay1)
    integer(int32), intent(in) :: CropDay1
    integer(int32), intent(inout) :: SimDay1

    SimDay1 = CropDay1
    if (GetClimFile() /= '(None)') then
        if ((SimDay1 < GetClimRecord_FromDayNr()) .or. &
            (SimDay1 > GetClimRecord_ToDayNr())) then
            call SetSimulation_LinkCropToSimPeriod(.false.)
            SimDay1 = GetClimRecord_FromDayNr()
        end if
    end if
end subroutine DetermineLinkedSimDay1


subroutine AdjustSimPeriod()
    integer(int32) :: IniSimFromDayNr
    character(len=:), allocatable :: FullFileName
    integer(int32) :: FromDayNr_temp
    type(CompartmentIndividual), dimension(max_No_compartments) :: Compartment_temp

    integer(int32) :: ZiAqua_tmp
    real(sp) :: ECiAqua_tmp

    IniSimFromDayNr = GetSimulation_FromDayNr()
    select case(GetSimulation_LinkCropToSimPeriod())
    case(.true.)
        FromDayNr_temp = GetSimulation_FromDayNr()
        call DetermineLinkedSimDay1(GetCrop_Day1(), FromDayNr_temp)
        call SetSimulation_FromDayNr(FromDayNr_temp)
        if (GetCrop_Day1() == GetSimulation_FromDayNr()) then
            call SetSimulation_ToDayNr(GetCrop_DayN())
        else
            call SetSimulation_ToDayNr(GetSimulation_FromDayNr() + 30) ! 30 days
        end if
        if (GetClimFile() /= '(None)') then
            if (GetSimulation_ToDayNr() > GetClimRecord_ToDayNr()) then
                call SetSimulation_ToDayNr(GetClimRecord_ToDayNr())
            end if
            if (GetSimulation_ToDayNr() < GetClimRecord_FromDayNr()) then
                call SetSimulation_ToDayNr(GetClimRecord_FromDayNr())
            end if
        end if
    case(.false.)
        if (GetSimulation_FromDayNr() > GetCrop_Day1()) then
            call SetSimulation_FromDayNr(GetCrop_Day1())
        end if
        call SetSimulation_ToDayNr(GetCrop_DayN())
        if ((GetClimFile() /= '(None)') .and. &
            ((GetSimulation_FromDayNr() <= GetClimRecord_FromDayNr()) &
             .or. (GetSimulation_FromDayNr() >= GetClimRecord_ToDayNr()))) then
            call SetSimulation_FromDayNr(GetClimRecord_FromDayNr())
            call SetSimulation_ToDayNr(GetSimulation_FromDayNr() + 30) ! 30 days
       end if
    end select

    ! adjust initial depth and quality of the groundwater when required
    if ((.not. GetSimulParam_ConstGwt()) .and. &
        (IniSimFromDayNr /= GetSimulation_FromDayNr())) then
        if (GetGroundWaterFile() == '(None)') then
            FullFileName = GetPathNameProg()// 'GroundWater.AqC'
        else
            FullFileName = GetGroundWaterFileFull()
        end if
        ! initialize ZiAqua and ECiAqua
        ZiAqua_tmp = GetZiAqua()
        ECiAqua_tmp = GetECiAqua()
        call LoadGroundWater(FullFileName, GetSimulation_FromDayNr(), &
                 ZiAqua_tmp, ECiAqua_tmp)
        call SetZiAqua(ZiAqua_tmp)
        call SetECiAqua(ECiAqua_tmp)
        Compartment_temp = GetCompartment()
        call CalculateAdjustedFC((GetZiAqua()/100._sp), Compartment_temp)
        call SetCompartment(Compartment_temp)
        if (GetSimulation_IniSWC_AtFC()) then
            call ResetSWCToFC
        end if
    end if
end subroutine AdjustSimPeriod


subroutine ResetSWCToFC()

    integer(int32) :: layeri, Loci, compi, celli

    call SetSimulation_IniSWC_AtDepths(.false.)
    if (GetZiAqua() < 0) then ! no ground water table
        call SetSimulation_IniSWC_NrLoc(GetSoil_NrSoilLayers())
        do layeri = 1, GetSoil_NrSoilLayers()
            call SetSimulation_IniSWC_Loc_i(layeri, &
                                            GetSoilLayer_Thickness(layeri))
            call SetSimulation_IniSWC_VolProc_i(layeri, &
                                                GetSoilLayer_FC(layeri))
            call SetSimulation_IniSWC_SaltECe_i(layeri, 0._sp)
        end do
        do layeri = (GetSoil_NrSoilLayers() + 1), max_No_compartments
            call SetSimulation_IniSWC_Loc_i(layeri, undef_double)
            call SetSimulation_IniSWC_VolProc_i(layeri, undef_double)
            call SetSimulation_IniSWC_SaltECe_i(layeri, undef_double)
        end do
    else
        call SetSimulation_IniSWC_NrLoc(int(GetNrCompartments(), kind=int8))
        do Loci = 1, GetSimulation_IniSWC_NrLoc()
            call SetSimulation_IniSWC_Loc_i(Loci, &
                                            GetCompartment_Thickness(Loci))
            call SetSimulation_IniSWC_VolProc_i(Loci, &
                                                GetCompartment_FCadj(Loci))
            call SetSimulation_IniSWC_SaltECe_i(Loci, 0.0_sp)
        end do
    end if
    do compi = 1, GetNrCompartments()
        call SetCompartment_Theta(compi, GetCompartment_FCadj(compi)/100._sp)
        call SetSimulation_ThetaIni_i(compi, GetCompartment_Theta(compi))
        do celli = 1, GetSoilLayer_SCP1(GetCompartment_Layer(compi))
            ! salinity in cells
            call SetCompartment_Salt(compi, celli, 0.0_sp)
            call SetCompartment_Depo(compi, celli, 0.0_sp)
        end do
    end do
end subroutine ResetSWCToFC


subroutine LoadCrop(FullName)
    character(len=*), intent(in) :: FullName

    integer :: fhandle
    integer(int32) :: XX, YY
    real(sp) :: VersionNr
    integer(int8) :: TempShortInt, perenperiod_onsetOcc_temp
    integer(int8) :: perenperiod_endOcc_temp
    integer(int32) :: TempInt, perenperiod_onsetFD_temp
    integer(int32) :: perenperiod_onsetFM_temp, perenperiod_onsetLSP_temp
    integer(int32) :: perenperiod_onsetPV_temp, perenperiod_endLD_temp
    integer(int32) :: perenperiod_endLM_temp, perenperiod_extrayears_temp
    integer(int32) :: perenperiod_endLSP_temp, perenperiod_endPV_temp
    real(sp) :: TempDouble
    logical :: TempBoolean
    real(sp) :: perenperiod_onsetTV_temp, perenperiod_endTV_temp
    real(sp) :: Crop_SmaxTop_temp, Crop_SmaxBot_temp
    character(len=1024) :: CropDescriptionLocal

    open(newunit=fhandle, file=trim(FullName), status='old', action='read')
    read(fhandle, '(a)') CropDescriptionLocal
    call SetCropDescription(trim(CropDescriptionLocal))
    read(fhandle, *) VersionNr ! AquaCrop version
    read(fhandle, *)  ! Protected or Open file

    ! subkind
    read(fhandle, *) XX
    select case (XX)
        case(1)
            call SetCrop_subkind(subkind_Vegetative)
        case(2)
            call SetCrop_subkind(subkind_Grain)
        case(3)
            call SetCrop_subkind(subkind_Tuber)
        case(4)
            call SetCrop_subkind(subkind_Forage)
    end select

    ! type of planting
    read(fhandle, *) XX
    select case (XX)
        case(1)
            call SetCrop_Planting(plant_Seed)
        case(0)
            call SetCrop_Planting(plant_Transplant)
        case(-9)
            call SetCrop_Planting(plant_Regrowth)
        case default
            call SetCrop_Planting(plant_Seed)
    end select

    ! mode
    read(fhandle, *) XX
    if (XX == 0) then
        call SetCrop_ModeCycle(ModeCycle_GDDays)
    else
        call SetCrop_ModeCycle(ModeCycle_CalendarDays)
    end if

    ! adjustment p to ETo
    read(fhandle, *) YY
    if (YY == 0) then
        call SetCrop_pMethod(pMethod_NoCorrection)
    elseif (YY == 1) then
        call SetCrop_pMethod(pMethod_FAOCorrection)
    end if

    ! temperatures controlling crop development
    read(fhandle, *) TempDouble
    call SetCrop_Tbase(TempDouble)
    read(fhandle, *) TempDouble
    call SetCrop_Tupper(TempDouble)

    ! required growing degree days to complete the crop cycle
    ! (is identical as to maturity)
    read(fhandle, *) TempInt
    call SetCrop_GDDaysToHarvest(TempInt)

    ! water stress
    read(fhandle, *) TempDouble
    call SetCrop_pLeafDefUL(TempDouble)
    read(fhandle, *) TempDouble
    call SetCrop_pLeafDefLL(TempDouble)
    read(fhandle, *) TempDouble
    call SetCrop_KsShapeFactorLeaf(TempDouble)
    read(fhandle, *) TempDouble
    call SetCrop_pdef(TempDouble)
    read(fhandle, *) TempDouble
    call SetCrop_KsShapeFactorStomata(TempDouble)
    read(fhandle, *) TempDouble
    call SetCrop_pSenescence(TempDouble)
    read(fhandle, *) TempDouble
    call SetCrop_KsShapeFactorSenescence(TempDouble)
    read(fhandle, *) TempInt
    call SetCrop_SumEToDelaySenescence(TempInt)
    read(fhandle, *) TempDouble
    call SetCrop_pPollination(TempDouble)
    read(fhandle, *) TempInt
    call SetCrop_AnaeroPoint(TempInt)

    ! soil fertility/salinity stress
    read(fhandle, *) TempShortInt   ! Soil fertility stress at calibration (%)
    call SetCrop_StressResponse_Stress(TempShortInt)
    read(fhandle, *) TempDouble     ! Shape factor for the response of Canopy
                                    ! Growth Coefficient to soil
                                    ! fertility/salinity stress
    call SetCrop_StressResponse_ShapeCGC(TempDouble)
    read(fhandle, *) TempDouble     ! Shape factor for the response of Maximum
                                    ! Canopy Cover to soil
                                    ! fertility/salinity stress
    call SetCrop_StressResponse_ShapeCCX(TempDouble)
    read(fhandle, *) TempDouble     ! Shape factor for the response of Crop
                                    ! Water Producitity to soil
                                    ! fertility stress
    call SetCrop_StressResponse_ShapeWP(TempDouble)
    read(fhandle, *) TempDouble     ! Shape factor for the response of Decline
                                    ! of Canopy Cover to soil
                                    ! fertility/salinity stress
    call SetCrop_StressResponse_ShapeCDecline(TempDouble)

    if (roundc(VersionNr*10, mold=1) >= 40) then
    ! UPDATE required for Version 4.0 and next
        read(fhandle, *)  ! Shape factor for the response of Stomatal Closure
                          ! to soil salinity stress NO LONGER VALID
    end if

    ! continue with soil fertility/salinity stress
    if ((GetCrop_StressResponse_ShapeCGC() > 24.9_sp) &
                .and. (GetCrop_StressResponse_ShapeCCX() > 24.9_sp) &
                .and. (GetCrop_StressResponse_ShapeWP() > 24.9_sp) &
                .and. (GetCrop_StressResponse_ShapeCDecline() > 24.9_sp)) then
        call SetCrop_StressResponse_Calibrated(.false.)
    else
        call SetCrop_StressResponse_Calibrated(.true.)
    end if

    ! temperature stress
    read(fhandle, *) TempShortInt   ! Minimum air temperature below which
                                    ! pollination starts to fail
                                    ! (cold stress) (degC)
    call SetCrop_Tcold(TempShortInt)
    read(fhandle, *) TempShortInt   ! Maximum air temperature above which
                                    ! pollination starts to fail
                                    ! (heat stress) (degC)
    call SetCrop_Theat(TempShortInt)
    read(fhandle, *) TempDouble     ! Minimum growing degrees required for full
                                    ! biomass production (degC - day)
    call SetCrop_GDtranspLow(TempDouble)

    ! salinity stress (Version 3.2 and higher)
    ! -----  UPDATE salinity stress
    if (roundc(VersionNr*10, mold=1) < 32) then
    ! UPDATE required for Version 3.0 and 3.1
        call SetCrop_ECemin(2_int8)     ! upper threshold ECe
        call SetCrop_ECemax(15_int8)    ! lower threhsold ECe
    else
        read(fhandle, *) TempShortInt       ! upper threshold ECe
        call SetCrop_ECemin(TempShortInt)   ! upper threshold ECe
        read(fhandle, *) TempShortInt       ! lower threhsold ECe
        call SetCrop_ECemax(TempShortInt)   ! upper threshold ECe
        read(fhandle, *) ! WAS shape factor of the Ks(salinity) - soil
                         ! saturation extract (ECe) relationship
    end if
    ! -----  UPDATE salinity stress (Version 5.1 and higher)
    if (roundc(VersionNr*10, mold=1) < 51) then
    ! UPDATE required for previous versions
        call SetCrop_CCsaltDistortion(25_int8)   ! distortion canopy cover for
                                            ! simulation of effect of
                                            ! salinity stress (%)
        call SetCrop_ResponseECsw(100) ! Response of Ks stomata to ECsw:
                                       ! From 0 (none) to +200 (very strong)
    else
        read(fhandle, *) TempShortInt
        call SetCrop_CCsaltDistortion(TempShortInt)
        read(fhandle, *) TempInt
        call SetCrop_ResponseECsw(TempInt)
    end if

    ! evapotranspiration
    read(fhandle, *) TempDouble
    call SetCrop_KcTop(TempDouble)
    read(fhandle, *) TempDouble
    call SetCrop_KcDecline(TempDouble)
    read(fhandle, *) TempDouble
    call SetCrop_RootMin(TempDouble)
    read(fhandle, *) TempDouble
    call SetCrop_RootMax(TempDouble)
    if (GetCrop_RootMin() > GetCrop_RootMax()) then
        call SetCrop_RootMin(GetCrop_RootMax()) ! security for sine function
    end if
    read(fhandle, *) TempShortInt
    call SetCrop_RootShape(TempShortInt)
    read(fhandle, *) TempDouble
    call SetCrop_SmaxTopQuarter(TempDouble)
    read(fhandle, *) TempDouble
    call SetCrop_SmaxBotQuarter(TempDouble)
    Crop_SmaxTop_temp = GetCrop_SmaxTop()
    Crop_SmaxBot_temp = GetCrop_SmaxBot()
    call DeriveSmaxTopBottom(GetCrop_SmaxTopQuarter(), &
                             GetCrop_SmaxBotQuarter(), &
                             Crop_SmaxTop_temp, Crop_SmaxBot_temp)
    call SetCrop_SmaxTop(Crop_SmaxTop_temp)
    call SetCrop_SmaxBot(Crop_SmaxBot_temp)
    read(fhandle, *) TempInt
    call SetCrop_CCEffectEvapLate(TempInt)

    ! crop development
    read(fhandle, *) TempDouble
    call SetCrop_SizeSeedling(TempDouble)
    if (roundc(VersionNr*10, mold=1) < 50) then
    ! UPDATE required for Version not yet 5.0
        call SetCrop_SizePlant(GetCrop_SizeSeedling())
    else
        read(fhandle, *) TempDouble ! Canopy size of individual plant
                                    ! (re-growth) at 1st day (cm2)
        call SetCrop_SizePlant(TempDouble)
    end if
    read(fhandle, *) TempInt
    call SetCrop_PlantingDens(TempInt)
    call SetCrop_CCo((GetCrop_PlantingDens()/10000._sp) &
                        * (GetCrop_SizeSeedling()/10000._sp))
    call SetCrop_CCini((GetCrop_PlantingDens()/10000._sp) &
                        * (GetCrop_SizePlant()/10000._sp))
    read(fhandle, *) TempDouble
    call SetCrop_CGC(TempDouble)

    read(fhandle, *) TempShortInt   ! Number of years at which CCx declines
                                    ! to 90 % of its value due to
                                    ! self-thinning - for Perennials
    call SetCrop_YearCCx(TempShortInt)
    read(fhandle, *) TempDouble     ! Shape factor of the decline of CCx over
                                    ! the years due to self-thinning
                                    ! for Perennials
    call SetCrop_CCxRoot(TempDouble)
    read(fhandle, *)

    read(fhandle, *) TempDouble
    call SetCrop_CCx(TempDouble)
    read(fhandle, *) TempDouble
    call SetCrop_CDC(TempDouble)
    read(fhandle, *) TempInt
    call SetCrop_DaysToGermination(TempInt)
    read(fhandle, *) TempInt
    call SetCrop_DaysToMaxRooting(TempInt)
    read(fhandle, *) TempInt
    call SetCrop_DaysToSenescence(TempInt)
    read(fhandle, *) TempInt
    call SetCrop_DaysToHarvest(TempInt)
    read(fhandle, *) TempInt
    call SetCrop_DaysToFlowering(TempInt)
    read(fhandle, *) TempInt
    call SetCrop_LengthFlowering(TempInt)
    ! -----  UPDATE crop development for Version 3.1
    ! leafy vegetable crop has an Harvest Index which builds up starting from sowing
    if ((GetCrop_subkind() == subkind_Vegetative) &
            .or. (GetCrop_subkind() == subkind_Forage)) then
        call SetCrop_DaysToFlowering(0)
        call SetCrop_LengthFlowering(0)
    end if

    ! Crop.DeterminancyLinked
    read(fhandle, *) XX
    select case (XX)
        case(1)
            call SetCrop_DeterminancyLinked(.true.)
        case default
            call SetCrop_DeterminancyLinked(.false.)
    end select

    ! Potential excess of fruits (%) and building up HI
    if ((GetCrop_subkind() == subkind_Vegetative) &
            .or. (GetCrop_subkind() == subkind_Forage)) then
        read(fhandle, *)  ! PercCycle no longer considered
        call SetCrop_fExcess(int(undef_int, kind=int16))
    else
        read(fhandle, *) TempInt
        call SetCrop_fExcess(int(TempInt, kind=int16))
    end if
    read(fhandle, *) TempInt
    call SetCrop_DaysToHIo(TempInt)

    ! yield response to water
    read(fhandle, *) TempDouble
    call SetCrop_WP(TempDouble)
    read(fhandle, *) TempInt
    call SetCrop_WPy(TempInt)
    ! adaptation to elevated CO2 (Version 3.2 and higher)
    ! -----  UPDATE Crop performance under elevated atmospheric CO2 concentration (%)
    if (roundc(VersionNr*10, mold=1) < 32) then ! UPDATE required for Version 3.0 and 3.1
        call SetCrop_AdaptedToCO2(50_int8)
    else
        read(fhandle, *) TempShortInt
        call SetCrop_AdaptedToCO2(TempShortInt)
    end if
    read(fhandle, *) TempInt
    call SetCrop_HI(TempInt)
    read(fhandle, *) TempShortInt
    call SetCrop_HIincrease(TempShortInt)   ! possible increase (%) of HI due
                                            ! to water stress before flowering
    read(fhandle, *) TempDouble
    call SetCrop_aCoeff(TempDouble)     ! coefficient describing impact of
                                        ! restricted vegetative growth at
                                        ! flowering on HI
    read(fhandle, *) TempDouble
    call SetCrop_bCoeff(TempDouble)     ! coefficient describing impact of
                                        ! stomatal closure at flowering on HI
    read(fhandle, *) TempShortInt
    call SetCrop_DHImax(TempShortInt)   ! allowable maximum increase (%) of
                                        ! specified HI
    ! -----  UPDATE yield response to water for Version 3.1
    ! leafy vegetable crop has an Harvest Index (default is 85 %)
    if ((roundc(VersionNr*10, mold=1) == 30) &
                .and. ((GetCrop_subkind() == subkind_Vegetative) &
                        .or. (GetCrop_subkind() == subkind_Forage))) then
        if (GetCrop_HI() == undef_int) then
            call SetCrop_HI(85)
        end if
    end if

    ! growing degree days
    read(fhandle, *) TempInt
    call SetCrop_GDDaysToGermination(TempInt)
    read(fhandle, *) TempInt
    call SetCrop_GDDaysToMaxRooting(TempInt)
    read(fhandle, *) TempInt
    call SetCrop_GDDaysToSenescence(TempInt)
    read(fhandle, *) TempInt
    call SetCrop_GDDaysToHarvest(TempInt)
    read(fhandle, *) TempInt
    call SetCrop_GDDaysToFlowering(TempInt)
    read(fhandle, *) TempInt
    call SetCrop_GDDLengthFlowering(TempInt)
    read(fhandle, *) TempDouble
    call SetCrop_GDDCGC(TempDouble)
    read(fhandle, *) TempDouble
    call SetCrop_GDDCDC(TempDouble)
    read(fhandle, *) TempInt
    call SetCrop_GDDaysToHIo(TempInt)

    ! -----  UPDATE yield response to water for Version 3.1
    ! leafy vegetable crop has an Harvest Index which builds up
    ! starting from sowing
    if ((GetCrop_ModeCycle() == ModeCycle_GDDays) &
            .and. ((GetCrop_subkind() == subkind_Vegetative) &
                .or. (GetCrop_subkind() == subkind_Forage))) then
        call SetCrop_GDDaysToFlowering(0)
        call SetCrop_GDDLengthFlowering(0)
    end if

    ! extra version 6.2
    if (roundc(VersionNr*10, mold=1) < 62) then
    ! UPDATE required for Version 6.2
        call SetCrop_DryMatter(int(undef_int, kind=int8)) ! undefined
    else
        read(fhandle, *) TempShortInt
        call SetCrop_DryMatter(TempShortInt) ! dry matter content (%)
                                             ! of fresh yield
    end if

    ! extra version 7.0
    if (roundc(VersionNr*10, mold=1) < 62) then
    ! UPDATE required for Version 7.0
        call SetCrop_RootMinYear1(GetCrop_RootMin()) ! regrowth not yet possible
        TempBoolean = (GetCrop_Planting() == plant_Seed)
        call SetCrop_SownYear1(TempBoolean)  ! type of planting first year
        ! transfer of assimilates
        call SetCrop_Assimilates_On(.false.) ! Transfer of assimilates between
                                             ! root system and above ground
                                             ! parts is NOT considered
        call SetCrop_Assimilates_Period(0)
        call SetCrop_Assimilates_Stored(0_int8)
        call SetCrop_Assimilates_Mobilized(0_int8)
    else
        read(fhandle, *) TempDouble
        call SetCrop_RootMinYear1(TempDouble) ! Minimum rooting depth in first
                                              ! year in meter (for regrowth)
        read(fhandle, *) XX
        select case (XX)
            case(1)
                call SetCrop_SownYear1(.true.)  ! crop is sown in 1 st year
                                                ! (for perennials)
            case default
                call SetCrop_SownYear1(.false.) ! crop is transplanted in
        end select                                ! 1st year (for regrowth)
        ! transfer of assimilates
        read(fhandle, *) XX
        select case (XX)
            case(1)
                call SetCrop_Assimilates_On(.true.)
                                                ! Transfer of assimilates from
                                                ! above ground parts to root
                                                ! system is considered
            case default
                call SetCrop_Assimilates_On(.false.)
                                                ! Transfer of assimilates from
                                                ! above ground parts to root
                                                ! system is NOT considered
        end select
        read(fhandle, *) TempInt
        call SetCrop_Assimilates_Period(TempInt)
                                                ! Number of days at end of season
                                                ! during which assimilates are
                                                ! stored in root system
        read(fhandle, *) TempShortInt
        call SetCrop_Assimilates_Stored(TempShortInt)
                                                ! Percentage of assimilates,
                                                ! transferred to root system
                                                ! at last day of season
        read(fhandle, *) TempShortInt
        call SetCrop_Assimilates_Mobilized(TempShortInt)
                                        ! Percentage of stored
                                        ! assimilates, transferred to above
                                        ! ground parts in next season
    end if

    if (GetCrop_subkind() == subkind_Forage) then
        ! data for the determination of the growing period
        ! 1. Title
        do XX = 1, 3
            read(fhandle, *)
        end do
        ! 2. ONSET
        read(fhandle, *) XX
        if (XX == 0) then
            call SetPerennialPeriod_GenerateOnset(.false.) ! onset is fixed on
                                                           ! a specific day
        else
            ! onset is generated by an air temperature criterion
            call SetPerennialPeriod_GenerateOnset(.true.)
            select case (XX)
                case(12)
                    call SetPerennialPeriod_OnsetCriterion(AirTCriterion_TMeanPeriod)
                            ! Criterion: mean air temperature
                case(13)
                    call SetPerennialPeriod_OnsetCriterion(AirTCriterion_GDDPeriod)
                            ! Criterion: growing-degree days
                case default
                    call SetPerennialPeriod_GenerateOnset(.false.)
            end select
        end if
        read(fhandle, *) perenperiod_onsetFD_temp
        call SetPerennialPeriod_OnsetFirstDay(perenperiod_onsetFD_temp)
        read(fhandle, *) perenperiod_onsetFM_temp
        call SetPerennialPeriod_OnsetFirstMonth(perenperiod_onsetFM_temp)
        read(fhandle, *) perenperiod_onsetLSP_temp
        call  SetPerennialPeriod_OnsetLengthSearchPeriod(perenperiod_onsetLSP_temp)
        read(fhandle, *) perenperiod_onsetTV_temp ! Mean air temperature
                                                  ! or Growing-degree days
        call SetPerennialPeriod_OnsetThresholdValue(perenperiod_onsetTV_temp)
        read(fhandle, *) perenperiod_onsetPV_temp ! number of succesive days
        call SetPerennialPeriod_OnsetPeriodValue(perenperiod_onsetPV_temp)
        read(fhandle, *) perenperiod_onsetOcc_temp  ! number of occurrence
        call SetPerennialPeriod_OnsetOccurrence(perenperiod_onsetOcc_temp)
        if (GetPerennialPeriod_OnsetOccurrence() > 3_int8) then
            call SetPerennialPeriod_OnsetOccurrence(3_int8)
        end if
        ! 3. END of growing period
        read(fhandle, *) XX
        if (XX == 0) then
            call SetPerennialPeriod_GenerateEnd(.false.)  ! end is fixed on a
                                                          ! specific day
        else
            ! end is generated by an air temperature criterion
            call SetPerennialPeriod_GenerateEnd(.true.)
            select case (XX)
                case(62)
                    call SetPerennialPeriod_EndCriterion(AirTCriterion_TMeanPeriod)
                                                ! Criterion: mean air temperature
                case(63)
                    call SetPerennialPeriod_EndCriterion(AirTCriterion_GDDPeriod)
                                                ! Criterion: growing-degree days
                case default
                    call SetPerennialPeriod_GenerateEnd(.false.)
            end select
        end if
        read(fhandle, *) perenperiod_endLD_temp
        call SetPerennialPeriod_EndLastDay(perenperiod_endLD_temp)
        read(fhandle, *) perenperiod_endLM_temp
        call SetPerennialPeriod_EndLastMonth(perenperiod_endLM_temp)
        read(fhandle, *) perenperiod_extrayears_temp
        call SetPerennialPeriod_ExtraYears(perenperiod_extrayears_temp)
        read(fhandle, *) perenperiod_endLSP_temp
        call SetPerennialPeriod_EndLengthSearchPeriod(perenperiod_endLSP_temp)
        read(fhandle, *) perenperiod_endTV_temp ! Mean air temperature
                                                ! or Growing-degree days
        call SetPerennialPeriod_EndThresholdValue(perenperiod_endTV_temp)
        read(fhandle, *) perenperiod_endPV_temp ! number of succesive days
        call SetPerennialPeriod_EndPeriodValue(perenperiod_endPV_temp)
        read(fhandle, *) perenperiod_endOcc_temp ! number of occurrence
        call SetPerennialPeriod_EndOccurrence(perenperiod_endOcc_temp)
        if (GetPerennialPeriod_EndOccurrence() > 3_int8) then
            call SetPerennialPeriod_EndOccurrence(3_int8)
        end if
    end if
    close(fhandle)
    ! maximum rooting depth in given soil profile
    call SetSoil_RootMax(RootMaxInSoilProfile(GetCrop_RootMax(),&
                                              GetSoil_NrSoilLayers(),&
                                              GetSoilLayer()))

    ! copy to CropFileSet
    call SetCropFileSet_DaysFromSenescenceToEnd(GetCrop_DaysToHarvest() &
                                                - GetCrop_DaysToSenescence())
    call SetCropFileSet_DaysToHarvest(GetCrop_DaysToHarvest())
    if (GetCrop_ModeCycle() == ModeCycle_GDDays) then
        call SetCropFileSet_GDDaysFromSenescenceToEnd(GetCrop_GDDaysToHarvest() &
                                                  - GetCrop_GDDaysToSenescence())
        call SetCropFileSet_GDDaysToHarvest(GetCrop_GDDaysToHarvest())
    else
        call SetCropFileSet_GDDaysFromSenescenceToEnd(undef_int)
        call SetCropFileSet_GDDaysToHarvest(undef_int)
    end if

end subroutine LoadCrop


real(sp) function SeasonalSumOfKcPot(TheDaysToCCini, TheGDDaysToCCini, L0, L12, &
                                     L123, L1234, GDDL0, GDDL12, GDDL123, &
                                     GDDL1234, CCo, CCx, CGC, GDDCGC, CDC, &
                                     GDDCDC, KcTop, KcDeclAgeing, &
                                     CCeffectProcent, Tbase, Tupper, TDayMin, &
                                     TDayMax, GDtranspLow, CO2i, TheModeCycle, ReferenceClimate)
    integer(int32), intent(in) :: TheDaysToCCini
    integer(int32), intent(in) :: TheGDDaysToCCini
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: L12
    integer(int32), intent(in) :: L123
    integer(int32), intent(in) :: L1234
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDL12
    integer(int32), intent(in) :: GDDL123
    integer(int32), intent(in) :: GDDL1234
    real(sp), intent(in) :: CCo
    real(sp), intent(in) :: CCx
    real(sp), intent(in) :: CGC
    real(sp), intent(in) :: GDDCGC
    real(sp), intent(in) :: CDC
    real(sp), intent(in) :: GDDCDC
    real(sp), intent(in) :: KcTop
    real(sp), intent(in) :: KcDeclAgeing
    real(sp), intent(in) :: CCeffectProcent
    real(sp), intent(in) :: Tbase
    real(sp), intent(in) :: Tupper
    real(sp), intent(in) :: TDayMin
    real(sp), intent(in) :: TDayMax
    real(sp), intent(in) :: GDtranspLow
    real(sp), intent(in) :: CO2i
    integer(intEnum), intent(in) :: TheModeCycle
    logical :: ReferenceClimate

    integer(int32), parameter :: EToStandard = 5
    real(sp) :: SumGDD, GDDi, SumKcPot, SumGDDforPlot, SumGDDfromDay1
    real(sp) :: Tndayi, Txdayi, CCi, CCxWitheredForB, TpotForB, EpotTotForB
    real(sp) :: CCinitial, DayFraction, GDDayFraction
    integer(int32) :: DayCC, Tadj, GDDTadj
    integer :: fhandle
    integer(int32) :: Dayi
    integer(int32) :: rc
    logical :: GrowthON
    integer :: i

    ! 1. Open Temperature file
    if ((GetTemperatureFile() /= '(None)') .and. &
        (GetTemperatureFile() /= '(External)')) then
        if (ReferenceClimate .eqv. .true.) then
            open(newunit=fhandle, file=trim(GetPathNameSimul()//'TCropReference.SIM'), &
                 status='old', action='read')
        else
            open(newunit=fhandle, file=trim(GetPathNameSimul()//'TCrop.SIM'), &
                 status='old', action='read')
        end if
    end if

    ! 2. Initialise global settings
    call SetSimulation_DelayedDays(0) ! required for CalculateETpot
    SumKcPot = 0._sp
    SumGDDforPlot = real(undef_int, kind=sp)
    SumGDD = real(undef_int, kind=sp)
    SumGDDfromDay1 = 0._sp
    GrowthON = .false.
    GDDTadj = undef_int
    DayFraction = real(undef_int, kind=sp)
    GDDayFraction = real(undef_int, kind=sp)
    ! 2.bis Initialise 1st day
    if (TheDaysToCCini /= 0) then
        ! regrowth
        if (TheDaysToCCini == undef_int) then
            ! CCx on 1st day
            Tadj = L12 - L0
            if (TheModeCycle == modeCycle_GDDays) then
                GDDTadj = GDDL12 - GDDL0
                SumGDD = GDDL12
            end if
            CCinitial = CCx
        else
            ! CC on 1st day is < CCx
            Tadj = TheDaysToCCini
            DayCC = Tadj + L0
            if (TheModeCycle == modeCycle_GDDays) then
                GDDTadj = TheGDDaysToCCini
                SumGDD = GDDL0 + TheGDDaysToCCini
                SumGDDforPlot = SumGDD
            end if
            CCinitial = CanopyCoverNoStressSF(DayCC, L0, L123, L1234, GDDL0, &
                                              GDDL123, GDDL1234, CCo, CCx, &
                                              CGC, CDC, GDDCGC, GDDCDC, &
                                              SumGDDforPlot, TheModeCycle, &
                                              (0_int8), (0_int8))
        end if
        ! Time reduction for days between L12 and L123
        DayFraction = (L123-L12)/real(Tadj + L0 + (L123-L12), kind=sp)
        if (TheModeCycle == modeCycle_GDDays) then
            GDDayFraction = (GDDL123-GDDL12)/real(GDDTadj + GDDL0 + &
                                                  (GDDL123-GDDL12), kind=sp)
        end if
    else
        ! sowing or transplanting
        Tadj = 0
        if (TheModeCycle == modeCycle_GDDays) then
            GDDTadj = 0
            SumGDD = 0._sp
        end if
        CCinitial = CCo
    end if

    ! 3. Calculate Sum
    i = 0
    do Dayi = 1, L1234
        ! 3.1 calculate growing degrees for the day
        if (GetTemperatureFile() == '(None)') then
            GDDi = DegreesDay(Tbase, Tupper, TDayMin, TDayMax, &
                                    GetSimulParam_GDDMethod())
        elseif (GetTemperatureFile() == '(External)') then
            i = i + 1
            if (i == size(GetTminCropReferenceRun())) then
                i = 1
            endif
            Tndayi = real(GetTminCropReferenceRun_i(i),kind=sp)
            Txdayi = real(GetTmaxCropReferenceRun_i(i),kind=sp)
            GDDi = DegreesDay(Tbase, Tupper, Tndayi, Txdayi, &
                                    GetSimulParam_GDDMethod())
        else
            read(fhandle, *, iostat=rc) Tndayi, Txdayi
            if ((rc == iostat_end) .and. (ReferenceClimate .eqv. .true.)) then
                close(fhandle)
                open(newunit=fhandle, file=trim(GetPathNameSimul()//'TCropReference.SIM'), &
                    status='old', action='read')
                read(fhandle, *, iostat=rc) Tndayi, Txdayi
            end if
            GDDi = DegreesDay(Tbase, Tupper, Tndayi, Txdayi, &
                                    GetSimulParam_GDDMethod())
        end if
        if (TheModeCycle == modeCycle_GDDays) then
            SumGDD = SumGDD + GDDi
            SumGDDfromDay1 = SumGDDfromDay1 + GDDi
        end if

        ! 3.2 calculate CCi
        if (GrowthON .eqv. .false.) then
            ! not yet canopy development
            CCi = 0._sp
            DayCC = Dayi
            if (TheDaysToCCini /= 0) then
                ! regrowth on 1st day
                CCi = CCinitial
                GrowthON = .true.
            else
                ! wait for day of germination or recover of transplant
                if (TheModeCycle == modeCycle_CalendarDays) then
                    if (Dayi == (L0+1)) then
                        CCi = CCinitial
                        GrowthON = .true.
                    end if
                else
                    if (SumGDD > GDDL0) then
                        CCi = CCinitial
                        GrowthON = .true.
                    end if
                end if
            end if
        else
            if (TheDaysToCCini == 0) then
                DayCC = Dayi
            else
                DayCC = Dayi + Tadj + L0 ! adjusted time scale
                if (DayCC > L1234) then
                    DayCC = L1234 ! special case where L123 > L1234
                end if
                if (DayCC > L12) then
                    if (Dayi <= L123) then
                        DayCC = L12 + roundc(DayFraction & ! slow down
                                            * (Dayi+Tadj+L0 - L12), mold=1)
                    else
                        DayCC = Dayi ! switch time scale
                    end if
                end if
            end if
            if (TheModeCycle == modeCycle_GDDays) then
                if (TheGDDaysToCCini == 0) then
                    SumGDDforPlot = SumGDDfromDay1
                else
                    SumGDDforPlot = SumGDD
                    if (SumGDDforPlot > GDDL1234) then
                        SumGDDforPlot = GDDL1234 ! special case
                                                 ! where L123 > L1234
                    end if
                    if (SumGDDforPlot > GDDL12) then
                        if (SumGDDfromDay1 <= GDDL123) then
                            SumGDDforPlot = GDDL12 + roundc(GDDayFraction &
                                   * (SumGDDfromDay1+GDDTadj+GDDL0 - GDDL12), &
                                      mold=1) ! slow down
                        else
                            SumGDDforPlot = SumGDDfromDay1 ! switch time scale
                        end if
                    end if
                end if
            end if
            CCi = CanopyCoverNoStressSF(DayCC, L0, L123, L1234, GDDL0, &
                                        GDDL123, GDDL1234, CCo, CCx, &
                                        CGC, CDC, GDDCGC, GDDCDC, &
                                        SumGDDforPlot, TheModeCycle, &
                                        (0_int8), (0_int8))
        end if

        ! 3.3 calculate CCxWithered
        CCxWitheredForB = CCi
        if (Dayi >= L12) then
            CCxWitheredForB = CCx
        end if

        ! 3.4 Calculate Tpot + Adjust for Low temperature
        ! (no transpiration)
        if (CCi > 0.0001_sp) then
            call CalculateETpot(DayCC, L0, L12, L123, L1234, (0), CCi, &
                           real(EToStandard, kind=sp), KcTop, &
                           KcDeclAgeing, CCx, CCxWitheredForB, &
                           CCeffectProcent, CO2i, &
                           GDDi, GDtranspLow, TpotForB, EpotTotForB)
        else
            TpotForB = 0._sp
        end if

        ! 3.5 Sum of Sum Of KcPot
        SumKcPot = SumKcPot + (TpotForB/EToStandard)
    end do

    ! 5. Close Temperature file
    if ((GetTemperatureFile() /= '(None)') .and. &
        (GetTemperatureFile() /= '(External)')) then
        close(fhandle)
    end if

    ! 6. final sum
    SeasonalSumOfKcPot = SumKcPot
end function SeasonalSumOfKcPot


real(sp) function HarvestIndexDay(DAP, DaysToFlower, HImax, dHIdt, CCi, &
                                  CCxadjusted, TheCCxWithered, &
                                  PercCCxHIfinal, TempPlanting, &
                                  PercentLagPhase, HIfinal)
    integer(int32), intent(in) :: DAP
    integer(int32), intent(in) :: DaysToFlower
    integer(int32), intent(in) :: HImax
    real(sp), intent(in) :: dHIdt
    real(sp), intent(in) :: CCi
    real(sp), intent(in) :: CCxadjusted
    real(sp), intent(in) :: TheCCxWithered
    integer(int8), intent(in) :: PercCCxHIfinal
    integer(intEnum), intent(in) :: TempPlanting
    integer(int8), intent(inout) :: PercentLagPhase
    integer(int32), intent(inout) :: HIfinal


    integer(int32), parameter :: HIo = 1
    real(sp) :: HIGC, HIday, HIGClinear, dHIdt_local
    integer(int32) :: t, tMax, tSwitch

    dHIdt_local = dHIdt
    t = DAP - GetSimulation_DelayedDays() - DaysToFlower
    ! Simulation.WPyON := false;
    PercentLagPhase = 0_int8
    if (t <= 0) then
        HIday = 0._sp
    else
        if ((GetCrop_Subkind() == subkind_Vegetative) &
                            .and. (TempPlanting == plant_Regrowth)) then
            dHIdt_local = 100._sp
        end if
        if ((GetCrop_Subkind() == subkind_Forage) &
                            .and. (TempPlanting == plant_Regrowth)) then
            dHIdt_local = 100._sp
        end if
        if (dHIdt_local > 99._sp) then
            HIday = HImax
            PercentLagPhase = 100_int8
        else
            HIGC = HarvestIndexGrowthCoefficient(real(HImax, kind=sp), &
                                                 dHIdt_local)
            call GetDaySwitchToLinear(HImax, dHIdt_local, HIGC, &
                                      tSwitch, HIGClinear)
            if (t < tSwitch) then
                PercentLagPhase = roundc(100._sp &
                                     * (t/real(tSwitch, kind=sp)), mold=1_int8)
                HIday = (HIo*HImax)/ (HIo+(HImax-HIo)*exp(-HIGC*t))
            else
                PercentLagPhase = 100_int8
                if ((GetCrop_subkind() == subkind_Tuber) &
                            .or. (GetCrop_subkind() == subkind_Vegetative) &
                            .or. (GetCrop_subkind() == subkind_Forage)) then
                    ! continue with logistic equation
                    HIday = (HIo*HImax)/ (HIo+(HImax-HIo)*exp(-HIGC*t))
                    if (HIday >= 0.9799_sp*HImax) then
                        HIday = HImax
                    end if
                else
                    ! switch to linear increase
                    HIday = (HIo*HImax)/ (HIo+(HImax-HIo)*exp(-HIGC*tSwitch))
                    HIday = Hiday + HIGClinear*(t-tSwitch)
                end if
            end if
            if (HIday > HImax) then
                HIday = HImax
            end if
            if (HIday <= (HIo + 0.4_sp)) then
                HIday = 0._sp
            end if
            if ((HImax - HIday) < 0.4_sp) then
                HIday = HImax
            end if
        end if

        ! adjust HIfinal if required for inadequate photosynthesis (unsufficient green canopy)
        tMax = roundc(HImax/dHIdt_local, mold=1)
        if ((HIfinal == HImax) .and. (t <= tmax) &
                              .and. ((CCi+epsilon(0._sp)) <= (PercCCxHIfinal/100._sp)) &
                              .and. (TheCCxWithered > epsilon(0._sp)) &
                              .and. (CCi < TheCCxWithered) &
                              .and. (GetCrop_subkind() /= subkind_Vegetative) &
                              .and. (GetCrop_subkind() /= subkind_Forage)) then
            HIfinal = roundc(HIday, mold=1)
        end if
        if (HIday > HIfinal) then
            HIday = HIfinal
        end if
    end if
    HarvestIndexDay = HIday

end function HarvestIndexDay


subroutine CompleteCropDescription()

    logical :: CGCisGiven
    integer(int32) :: FertStress
    integer(int8) :: RedCGC_temp, RedCCX_temp
    integer(int32) :: Crop_DaysToSenescence_temp
    integer(int32), dimension(4) :: Crop_Length_temp
    integer(int32) :: Crop_DaysToFullCanopy_temp
    real(sp) :: Crop_CGC_temp
    integer(int32) :: Crop_DaysToFullCanopySF_temp

    if ((GetCrop_subkind() == subkind_Vegetative) &
                    .or. (GetCrop_subkind() == subkind_Forage)) then
        if (GetCrop_DaysToHIo() > 0) then
            if (GetCrop_DaysToHIo() > GetCrop_DaysToHarvest()) then
                call SetCrop_dHIdt(GetCrop_HI()&
                                /real(GetCrop_DaysToHarvest(), kind=sp))
            else
                call SetCrop_dHIdt(GetCrop_HI()&
                                /real(GetCrop_DaysToHIo(), kind=sp))
            end if
            if (GetCrop_dHIdt() > 100._sp) then
                call SetCrop_dHIdt(100._sp)
            end if
        else
            call SetCrop_dHIdt(100._sp)
        end if
    else
        !  grain or tuber crops
        if (GetCrop_DaysToHIo() > 0) then
            call SetCrop_dHIdt(GetCrop_HI()&
                                /real(GetCrop_DaysToHIo(), kind=sp))
        else
            call SetCrop_dHIdt(real(undef_int, kind=sp))
        end if
    end if
    if (GetCrop_ModeCycle() == modeCycle_CalendarDays) then
        call SetCrop_DaysToCCini(TimeToCCini(GetCrop_Planting(), &
                                             GetCrop_PlantingDens(), &
                                             GetCrop_SizeSeedling(), &
                                             GetCrop_SizePlant(), &
                                             GetCrop_CCx(), GetCrop_CGC()))
        call SetCrop_DaysToFullCanopy(DaysToReachCCwithGivenCGC(&
                (0.98_sp * GetCrop_CCx()), GetCrop_CCo(), &
                     GetCrop_CCx(), GetCrop_CGC(), GetCrop_DaysToGermination()))
        if (GetManagement_FertilityStress() /= 0_int32) then
            FertStress = GetManagement_FertilityStress()
            Crop_DaysToFullCanopySF_temp = GetCrop_DaysToFullCanopySF()
            RedCGC_temp = GetSimulation_EffectStress_RedCGC()
            RedCCX_temp = GetSimulation_EffectStress_RedCCX()
            call TimeToMaxCanopySF(GetCrop_CCo(), GetCrop_CGC(), GetCrop_CCx(), &
                                   GetCrop_DaysToGermination(), &
                                   GetCrop_DaysToFullCanopy(), &
                                   GetCrop_DaysToSenescence(), &
                                   GetCrop_DaysToFlowering(), &
                                   GetCrop_LengthFlowering(), &
                                   GetCrop_DeterminancyLinked(), &
                                   Crop_DaysToFullCanopySF_temp, RedCGC_temp, &
                                   RedCCX_temp, FertStress)
            call SetManagement_FertilityStress(FertStress)
            call SetSimulation_EffectStress_RedCGC(RedCGC_temp)
            call SetSimulation_EffectStress_RedCCX(RedCCX_temp)
            call SetCrop_DaysToFullCanopySF(Crop_DaysToFullCanopySF_temp)
        else
            call SetCrop_DaysToFullCanopySF(GetCrop_DaysToFullCanopy())
        end if
        call SetCrop_GDDaysToCCini(undef_int)
        call SetCrop_GDDaysToGermination(undef_int)
        call SetCrop_GDDaysToFullCanopy(undef_int)
        call SetCrop_GDDaysToFullCanopySF(undef_int)
        call SetCrop_GDDaysToFlowering(undef_int)
        call SetCrop_GDDLengthFlowering(undef_int)
        call SetCrop_GDDaysToSenescence(undef_int)
        call SetCrop_GDDaysToHarvest(undef_int)
        call SetCrop_GDDaysToMaxRooting(undef_int)
        call SetCrop_GDDCGC(real(undef_int, kind=sp))
        call SetCrop_GDDCDC(real(undef_int, kind=sp))
    else
        call SetCrop_GDDaysToCCini(TimeToCCini(GetCrop_Planting(), &
                           GetCrop_PlantingDens(), GetCrop_SizeSeedling(), &
                           GetCrop_SizePlant(), &
                           GetCrop_CCx(), GetCrop_GDDCGC()))
        call SetCrop_DaysToCCini(TimeToCCini(GetCrop_Planting(), &
                           GetCrop_PlantingDens(), GetCrop_SizeSeedling(), &
                           GetCrop_SizePlant(), GetCrop_CCx(), GetCrop_CGC()))
        call SetCrop_GDDaysToFullCanopy(DaysToReachCCwithGivenCGC(&
                    (0.98_sp * GetCrop_CCx()), GetCrop_CCo(), GetCrop_CCx(), &
                    GetCrop_GDDCGC(), GetCrop_GDDaysToGermination()))
        ! Crop.GDDaysToFullCanopySF is determined in RUN or ManagementUnit if required
    end if

    CGCisGiven = .true. ! required to adjust Crop.DaysToFullCanopy (does not exist)
    Crop_DaysToSenescence_temp = GetCrop_DaysToSenescence()
    Crop_Length_temp = GetCrop_Length()
    Crop_DaysToFullCanopy_temp = GetCrop_DaysToFullCanopy()
    Crop_CGC_temp = GetCrop_CGC()
    call DetermineLengthGrowthStages(GetCrop_CCo(), GetCrop_CCx(), &
                                     GetCrop_CDC(), &
                                     GetCrop_DaysToGermination(), &
                                     GetCrop_DaysToHarvest(), CGCisGiven, &
                                     GetCrop_DaysToCCini(), GetCrop_Planting(), &
                                     Crop_DaysToSenescence_temp, &
                                     Crop_Length_temp, &
                                     Crop_DaysToFullCanopy_temp, Crop_CGC_temp)
    call SetCrop_DaysToSenescence(Crop_DaysToSenescence_temp)
    call SetCrop_Length(Crop_Length_temp)
    call SetCrop_DaysToFullCanopy(Crop_DaysToFullCanopy_temp)
    call SetCrop_CGC(Crop_CGC_temp)

    call SetCrop_CCoAdjusted(GetCrop_CCo())
    call SetCrop_CCxAdjusted(GetCrop_CCx())
    call SetCrop_CCxWithered(GetCrop_CCx())
    call SetSumWaBal_Biomass(0._sp)
    call SetSumWaBal_BiomassPot(0._sp)
    call SetSumWaBal_BiomassUnlim(0._sp)
    call SetSumWaBal_BiomassTot(0._sp) ! crop and weeds (for soil fertility stress)
    call SetSumWaBal_YieldPart(0._sp)
    call SetSimulation_EvapLimitON(.false.)
end subroutine CompleteCropDescription


subroutine NoManagementOffSeason()

    integer(int32) :: Nri

    call SetOffSeasonDescription('No specific off-season conditions')
    ! mulches
    call SetManagement_SoilCoverBefore(0_int8)
    call SetManagement_SoilCoverAfter(0_int8)
    call SetManagement_EffectMulchOffS(50_int8)
    ! off-season irrigation
    call SetSimulParam_IrriFwOffSeason(100_int8)
    call SetIrriECw_PreSeason(0.0_sp) ! dS/m
    do Nri = 1, 5
        call SetIrriBeforeSeason_DayNr(Nri, 0)
        call SetIrriBeforeSeason_Param(Nri, 0)
    end do
    call SetIrriECw_PostSeason(0.0_sp) ! dS/m
    do Nri = 1, 5
        call SetIrriAfterSeason_DayNr(Nri, 0)
        call SetIrriAfterSeason_Param(Nri, 0)
    end do
end subroutine NoManagementOffSeason


subroutine LoadOffSeason(FullName)
    character(len=*), intent(in) :: FullName

    integer :: fhandle
    integer(int32) :: Nri, NrEvents1, NrEvents2
    character(len=:), allocatable :: ParamString
    real(sp) :: Par1, Par2
    real(sp) :: VersionNr
    real(sp) :: PreSeason_in
    real(sp) :: PostSeason_in
    integer(int8) :: TempShortInt, simul_irri_of
    character(len=1025) :: OffSeasonDescr_temp

    logical :: file_exists

    inquire(file=trim(FullName), exist=file_exists)
    if (file_exists) then
        open(newunit=fhandle, file=trim(FullName), status='old', action='read')
    else
        write(*,*) 'LoadOffSeason file not found'
    end if
    read(fhandle, '(a)') OffSeasonDescr_temp
    call SetOffSeasonDescription(trim(OffSeasonDescr_temp))
    read(fhandle, *) VersionNr ! AquaCrop Version
    ! mulches
    read(fhandle, *) TempShortInt
    call SetManagement_SoilCoverBefore(TempShortInt)
    read(fhandle, *) TempShortInt
    call SetManagement_SoilCoverAfter(TempShortInt)
    read(fhandle, *) TempShortInt
    call SetManagement_EffectMulchOffS(TempShortInt)

    ! irrigation events - initialise
    do Nri = 1, 5
        call SetIrriBeforeSeason_DayNr(Nri, 0)
        call SetIrriBeforeSeason_Param(Nri, 0)
        call SetIrriAfterSeason_DayNr(Nri, 0)
        call SetIrriAfterSeason_Param(Nri, 0)
    end do
    read(fhandle, *) NrEvents1 ! number of irrigation events BEFORE growing period
    if (roundc(10*VersionNr, mold=1) < 32) then ! irrigation water quality BEFORE growing period
        call SetIrriECw_PreSeason(0.0_sp)
    else
        read(fhandle, *) PreSeason_in
        call SetIrriECw_PreSeason(PreSeason_in)
    end if
    read(fhandle, *) NrEvents2 ! number of irrigation events AFTER growing period
    if (roundc(10*VersionNr, mold=1) < 32) then ! irrigation water quality AFTER growing period
        call SetIrriECw_PostSeason(0.0_sp)
    else
        read(fhandle, *) PostSeason_in
        call SetIrriECw_PostSeason(PostSeason_in)
    end if
    read(fhandle, *) simul_irri_of ! percentage of soil surface wetted
    call SetSimulParam_IrriFwOffSeason(simul_irri_of)
    ! irrigation events - get events before and after season
    if ((NrEvents1 > 0) .or. (NrEvents2 > 0)) then
        do Nri = 1, 3
            read(fhandle, *) ! title
        end do
    end if
    if (NrEvents1 > 0) then
        do Nri = 1, NrEvents1
            ! events BEFORE growing period
            read(fhandle, '(a)') ParamString
            call SplitStringInTwoParams(ParamString, Par1, Par2)
            call SetIrriBeforeSeason_DayNr(Nri, roundc(Par1, mold=1))
            call SetIrriBeforeSeason_Param(Nri, roundc(Par2, mold=1))
        end do
    end if
    if (NrEvents2 > 0) then
        do Nri = 1, NrEvents2
            ! events AFTER growing period
            read(fhandle, '(a)') ParamString
            call SplitStringInTwoParams(ParamString, Par1, Par2)
            call SetIrriAfterSeason_DayNr(Nri, roundc(Par1, mold=1))
            call SetIrriAfterSeason_Param(Nri, roundc(Par2, mold=1))
        end do
    end if
    close(fhandle)
end subroutine LoadOffSeason


subroutine AdjustThetaInitial(PrevNrComp, PrevThickComp, PrevVolPrComp, PrevECdSComp)
    integer(int8), intent(in) :: PrevNrComp
    real(sp), dimension(max_No_compartments), intent(in) :: PrevThickComp
    real(sp), dimension(max_No_compartments), intent(in) :: PrevVolPrComp
    real(sp), dimension(max_No_compartments), intent(in) :: PrevECdSComp

    integer(int32) :: layeri, compi
    real(sp) :: TotDepthC, TotDepthL, Total
    type(CompartmentIndividual), &
                        dimension(max_No_compartments) :: Compartment_temp

    ! 1. Actual total depth of compartments
    TotDepthC = 0._sp
    do compi = 1, GetNrCompartments()
        TotDepthC = TotDepthC + GetCompartment_Thickness(compi)
    end do

    ! 2. Stretch thickness of bottom soil layer if required
    TotDepthL = 0._sp
    do layeri = 1, GetSoil_NrSoilLayers()
        TotDepthL = TotDepthL + GetSoilLayer_Thickness(layeri)
    end do
    if (TotDepthC > TotDepthL) then
        call SetSoilLayer_Thickness(int(GetSoil_NrSoilLayers(), kind=int32), &
               GetSoilLayer_Thickness(int(GetSoil_NrSoilLayers(), kind=int32)) &
                                                        + (TotDepthC - TotDepthL))
    end if

    ! 3. Assign a soil layer to each soil compartment
    Compartment_temp = GetCompartment()
    call DesignateSoilLayerToCompartments(GetNrCompartments(), &
                             int(GetSoil_NrSoilLayers(), kind=int32), &
                             Compartment_temp)
    call SetCompartment(Compartment_temp)

    ! 4. Adjust initial Soil Water Content of soil compartments
    if (GetSimulation_ResetIniSWC()) then
        if (GetSimulation_IniSWC_AtDepths()) then
            Compartment_temp = GetCompartment()
            call TranslateIniPointsToSWProfile(GetSimulation_IniSWC_NrLoc(), &
                                               GetSimulation_IniSWC_Loc(), &
                                               GetSimulation_IniSWC_VolProc(), &
                                               GetSimulation_IniSWC_SaltECe(), &
                                               GetNrCompartments(), &
                                               Compartment_temp)
            call SetCompartment(Compartment_temp)
        else
            Compartment_temp = GetCompartment()
            call TranslateIniLayersToSWProfile(GetSimulation_IniSWC_NrLoc(), &
                                               GetSimulation_IniSWC_Loc(), &
                                               GetSimulation_IniSWC_VolProc(), &
                                               GetSimulation_IniSWC_SaltECe(), &
                                               GetNrCompartments(), &
                                               Compartment_temp)
            call SetCompartment(Compartment_temp)
        end if
    else
        Compartment_temp = GetCompartment()
        call TranslateIniLayersToSWProfile(PrevNrComp, PrevThickComp, &
                                           PrevVolPrComp, PrevECdSComp, &
                                           GetNrCompartments(), &
                                           Compartment_temp)
        call SetCompartment(Compartment_temp)
    end if

    ! 5. Adjust watercontent in soil layers and determine ThetaIni
    Total = 0._sp
    do layeri = 1, GetSoil_NrSoilLayers()
        call SetSoilLayer_WaterContent(layeri, 0._sp)
    end do
    do compi = 1, GetNrCompartments()
        call SetSimulation_ThetaIni_i(compi, GetCompartment_Theta(compi))
        call SetSoilLayer_WaterContent(GetCompartment_Layer(compi), &
                         GetSoilLayer_WaterContent(GetCompartment_Layer(compi)) &
                            + GetSimulation_ThetaIni_i(compi)*100._sp*10._sp &
                                * GetCompartment_Thickness(compi))
    end do
    do layeri = 1, GetSoil_NrSoilLayers()
        Total = Total + GetSoilLayer_WaterContent(layeri)
    end do
    call SetTotalWaterContent_BeginDay(Total)
end subroutine AdjustThetaInitial


subroutine LoadClim(FullName, ClimateDescription, ClimateRecord)
    character(len=*), intent(in) :: FullName
    character(len=*), intent(inout) :: ClimateDescription
    type(rep_clim), intent(inout) :: ClimateRecord

    integer :: fhandle
    integer(int32) :: Ni, rc

    logical :: file_exists

    inquire(file=trim(FullName), exist=file_exists)
    if (file_exists) then
        open(newunit=fhandle, file=trim(FullName), status='old', action='read', &
        iostat=rc)
    else
        write(*,*) 'Climate file not found: ' // trim(FullName)
        return
    end if

    read(fhandle, '(a)', iostat=rc) ClimateDescription
    read(fhandle, *, iostat=rc) Ni
    if (Ni == 1) then
        ClimateRecord%DataType = datatype_Daily
    elseif (Ni == 2) then
        ClimateRecord%DataType = datatype_Decadely
    else
        ClimateRecord%DataType = datatype_Monthly
    end if
    read(fhandle, *, iostat=rc) ClimateRecord%FromD
    read(fhandle, *, iostat=rc) ClimateRecord%FromM
    read(fhandle, *, iostat=rc) ClimateRecord%FromY
    read(fhandle, *, iostat=rc)
    read(fhandle, *, iostat=rc)
    read(fhandle, *, iostat=rc)
    ClimateRecord%NrObs = 0
    read(fhandle, *, iostat=rc)
    do while (rc /= iostat_end)
        ClimateRecord%NrObs = ClimateRecord%NrObs + 1
        read(fhandle, *, iostat=rc)
    end do
    close(fhandle)
    call CompleteClimateDescription(ClimateRecord)
end subroutine LoadClim


subroutine LoadGroundWater(FullName, AtDayNr, Zcm, ECdSm)
    character(len=*), intent(in) :: FullName
    integer(int32), intent(in) :: AtDayNr
    integer(int32), intent(inout) :: Zcm
    real(sp), intent(inout) :: ECdSm

    integer :: fhandle
    integer(int32) :: i, dayi, monthi, yeari, Year1Gwt, rc
    integer(int32) :: DayNr1Gwt, DayNr1, DayNr2, DayNrN, AtDayNr_local
    character(len=255) :: StringREAD
    real(sp) :: DayDouble, Z1, EC1, Z2, EC2, ZN, ECN
    logical :: TheEnd

    logical :: file_exists

    AtDayNr_local = AtDayNr
    ! initialize
    TheEnd = .false.
    Year1Gwt = 1901
    DayNr1 = 1
    DayNr2 = 1

    inquire(file=trim(FullName), exist=file_exists)
    if (file_exists) then
        open(newunit=fhandle, file=trim(FullName), status='old', action='read')
    else
        write(*,*) 'Groundwater file not found'
    end if
    read(fhandle, '(a)') GroundWaterDescription
    read(fhandle, *) ! AquaCrop Version

    ! mode groundwater table
    read(fhandle, *) i
    select case (i)
        case(0)
            ! no groundwater table
            Zcm = undef_int
            ECdSm = real(undef_int, kind=sp)
            call SetSimulParam_ConstGwt(.true.)
            TheEnd = .true.
        case(1)
            ! constant groundwater table
            call SetSimulParam_ConstGwt(.true.)
        case default
            call SetSimulParam_ConstGwt(.false.)
    end select

    ! first day of observations (only for variable groundwater table)
    if (.not. GetSimulParam_ConstGwt()) then
        read(fhandle, *) dayi
        read(fhandle, *) monthi
        read(fhandle, *) Year1Gwt
        call DetermineDayNr(dayi, monthi, Year1Gwt, DayNr1Gwt)
    end if

    ! single observation (Constant Gwt) or first observation (Variable Gwt)
    if (i > 0) then
        ! groundwater table is present
        read(fhandle, *)
        read(fhandle, *)
        read(fhandle, *)
        read(fhandle, '(a)', iostat=rc) StringREAD
        call SplitStringInThreeParams(StringREAD, DayDouble, Z2, EC2)
        if ((i == 1) .or. (rc == iostat_end)) then
            ! Constant groundwater table or single observation
            Zcm = roundc(100._sp*Z2, mold=1)
            ECdSm = EC2
            TheEnd = .true.
        else
            DayNr2 = DayNr1Gwt + roundc(DayDouble, mold=1) - 1
        end if
    end if

    ! other observations
    if (.not. TheEnd) then
        ! variable groundwater table with more than 1 observation
        ! adjust AtDayNr
        call DetermineDate(AtDayNr_local, dayi, monthi, yeari)
        if ((yeari == 1901) .and. (Year1Gwt /= 1901)) then
            ! Make AtDayNr defined
            call DetermineDayNr(dayi, monthi, Year1Gwt, AtDayNr_local)
        end if
        if ((yeari /= 1901) .and. (Year1Gwt == 1901)) then
            ! Make AtDayNr undefined
            call DetermineDayNr(dayi, monthi, Year1Gwt, AtDayNr_local)
        end if
        ! get observation at AtDayNr
        if (Year1Gwt /= 1901) then
            ! year is defined
            if (AtDayNr_local <= DayNr2) then
                Zcm = roundc(100._sp*Z2, mold=1)
                ECdSm = EC2
            else
                do while (.not. TheEnd)
                    DayNr1 = DayNr2
                    Z1 = Z2
                    EC1 = EC2
                    read(fhandle, '(a)', iostat=rc) StringREAD
                    call SplitStringInThreeParams(StringREAD, DayDouble, &
                                                                Z2, EC2)
                    DayNr2 = DayNr1Gwt + roundc(DayDouble, mold=1) - 1
                    if (AtDayNr_local <= DayNr2) then
                        call FindValues(AtDayNr_local, DayNr1, DayNr2, Z1, &
                                        EC1, Z2, EC2, Zcm, ECdSm)
                        TheEnd = .true.
                    end if
                    if ((rc == iostat_end) .and. (.not. TheEnd)) then
                        Zcm = roundc(100._sp*Z2, mold=1)
                        ECdSm = EC2
                        TheEnd = .true.
                    end if
                end do
            end if
        else
            ! year is undefined
            if (AtDayNr_local <= DayNr2) then
                DayNr2 = DayNr2 + 365
                AtDayNr_local = AtDayNr_local + 365
                do while (rc /= iostat_end)
                    read(fhandle, '(a)') StringREAD
                    call SplitStringInThreeParams(StringREAD, DayDouble, &
                                                                Z1, EC1)
                    DayNr1 = DayNr1Gwt + roundc(DayDouble, mold=1) - 1
                end do
                call FindValues(AtDayNr_local, DayNr1, DayNr2, Z1, EC1, Z2, EC2, &
                                                               Zcm, ECdSm)
            else
                DayNrN = DayNr2 + 365
                ZN = Z2
                ECN = EC2
                do while (.not. TheEnd)
                    DayNr1 = DayNr2
                    Z1 = Z2
                    EC1 = EC2
                    read(fhandle, '(a)', iostat=rc) StringREAD
                    call SplitStringInThreeParams(StringREAD, DayDouble, &
                                                                Z2, EC2)
                    DayNr2 = DayNr1Gwt + roundc(DayDouble, mold=1) - 1
                    if (AtDayNr_local <= DayNr2) then
                        call FindValues(AtDayNr_local, DayNr1, DayNr2, Z1, EC1, &
                                                     Z2, EC2, Zcm, ECdSm)
                        TheEnd = .true.
                    end if
                    if ((rc == iostat_end) .and. (.not. TheEnd)) then
                        call FindValues(AtDayNr_local, DayNr2, DayNrN, Z2, EC2, &
                                                     ZN, ECN, Zcm, ECdSm)
                        TheEnd = .true.
                    end if
                end do
            end if
        end if
        ! variable groundwater table with more than 1 observation
    end if
    Close(fhandle)


    contains


    subroutine FindValues(AtDayNr, DayNr1, DayNr2, Z1, EC1, Z2, EC2, &
                                                         Zcm, ECdSm)
        integer(int32), intent(in) :: AtDayNr
        integer(int32), intent(in) :: DayNr1
        integer(int32), intent(in) :: DayNr2
        real(sp), intent(in) :: Z1
        real(sp), intent(in) :: EC1
        real(sp), intent(in) :: Z2
        real(sp), intent(in) :: EC2
        integer(int32), intent(inout) :: Zcm
        real(sp), intent(inout) :: ECdSm

        Zcm = roundc(100._sp * (Z1 + (Z2-Z1) &
                    * real(AtDayNr-DayNr1, kind=sp) &
                    / real(DayNr2-Daynr1, kind=sp)), mold=1)
        ECdSm = EC1 + (EC2-EC1) &
                        * real(AtDayNr-DayNr1, kind=sp) &
                        / real(DayNr2-Daynr1, kind=sp)
        end subroutine FindValues

end subroutine LoadGroundWater


subroutine AdjustClimRecordTo(CDayN)
    integer(int32), intent(in) :: CDayN

    integer(int32) :: dayi, monthi, yeari
    integer(int32) :: ToDayNr_tmp

    call DetermineDate(CDayN, dayi, monthi, yeari)
    call SetClimRecord_ToD(31)
    call SetClimRecord_ToM(12)
    call SetClimRecord_ToY(yeari)
    call DetermineDayNr(GetClimRecord_ToD(), GetClimRecord_ToM(), &
                        GetClimRecord_ToY(), ToDayNr_tmp)
    call SetClimRecord_ToDayNr(ToDayNr_tmp)
end subroutine AdjustClimRecordTo


subroutine TranslateIniLayersToSWProfile(NrLay, LayThickness, LayVolPr, LayECdS, NrComp, Comp)
    integer(int8), intent(in) :: NrLay
    real(sp), dimension(max_No_compartments), intent(in) :: LayThickness
    real(sp), dimension(max_No_compartments), intent(in) :: LayVolPr
    real(sp), dimension(max_No_compartments), intent(in) :: LayECdS
    integer(int32), intent(in) :: NrComp
    type(CompartmentIndividual), dimension(max_No_compartments), intent(inout) :: Comp

    integer(int8) :: Compi, Layeri, i
    real(sp) :: SDLay, SDComp, FracC
    logical :: GoOn

    ! from specific layers to Compartments
    do Compi = 1, int(NrComp, kind=int8)
        Comp(Compi)%Theta = 0._sp
        Comp(Compi)%WFactor = 0._sp  ! used for ECe in this procedure
    end do
    Compi = 0_int8
    SDComp = 0._sp
    Layeri = 1_int8
    SDLay = LayThickness(1)
    GoOn = .true.
    do while (Compi < NrComp)
        FracC = 0._sp
        Compi = Compi + 1_int8
        SDComp = SDComp + Comp(compi)%Thickness
        if (SDLay >= SDComp) then
            Comp(Compi)%Theta = Comp(Compi)%Theta + (1._sp-FracC) &
                                                *LayVolPr(Layeri)/100._sp
            Comp(Compi)%WFactor = Comp(Compi)%WFactor + (1._sp-FracC) &
                                                *LayECdS(Layeri)
        else
            ! go to next layer
            do while ((SDLay < SDComp) .and. GoOn)
                ! finish handling previous layer
                FracC = (SDLay - (SDComp-Comp(Compi)%Thickness)) &
                                /(Comp(Compi)%Thickness) - FracC
                Comp(Compi)%Theta = Comp(Compi)%Theta + &
                                    FracC*LayVolPr(Layeri)/100._sp
                Comp(Compi)%WFactor = Comp(Compi)%WFactor + FracC*LayECdS(Layeri)
                FracC = (SDLay - (SDComp-Comp(Compi)%Thickness)) &
                                    /(Comp(Compi)%Thickness)
                ! add next layer
                if (Layeri < NrLay) then
                    Layeri = Layeri + 1_int8
                    SDLay = SDLay + LayThickness(Layeri)
                else
                    GoOn = .false.
                end if
            end do
            Comp(Compi)%Theta = Comp(Compi)%Theta + (1._sp-FracC) &
                                                    *LayVolPr(Layeri)/100._sp
            Comp(Compi)%WFactor = Comp(Compi)%WFactor + (1._sp-FracC) &
                                                    *LayECdS(Layeri)
        end if
        ! next Compartment
    end do
    if (.not. GoOn) then
        do i = (Compi+1_int8), int(NrComp, kind=int8)
            Comp(i)%Theta = LayVolPr(NrLay)/100._sp
            Comp(i)%WFactor = LayECdS(NrLay)
        end do
    end if

    ! final check of SWC
    do Compi = 1_int8, int(NrComp, kind=int8)
        if (Comp(Compi)%Theta > &
                (GetSoilLayer_SAT(Comp(compi)%Layer))/100._sp) then
            Comp(Compi)%Theta = (GetSoilLayer_SAT(Comp(compi)%Layer)) &
                                                            /100._sp
        end if
    end do
    ! salt distribution in cellls
    do Compi = 1_int8, int(NrComp, kind=int8)
        call DetermineSaltContent(Comp(Compi)%WFactor, Comp(Compi))
    end do
end subroutine TranslateIniLayersToSWProfile


subroutine TranslateIniPointsToSWProfile(NrLoc, LocDepth, LocVolPr, LocECdS, &
                                                               NrComp, Comp)
    integer(int8), intent(in) :: NrLoc
    real(sp), dimension(max_No_compartments), intent(in) :: LocDepth
    real(sp), dimension(max_No_compartments), intent(in) :: LocVolPr
    real(sp), dimension(max_No_compartments), intent(in) :: LocECdS
    integer(int32), intent(in) :: NrComp
    type(CompartmentIndividual), dimension(max_No_compartments), &
                                              intent(inout) :: Comp

    integer(int32) :: Compi, Loci
    real(sp) :: TotD, Depthi, D1, D2, Th1, Th2, DTopComp, &
                ThTopComp, ThBotComp
    real(sp) :: EC1, EC2, ECTopComp, ECBotComp
    logical :: AddComp, TheEnd

    TotD = 0
    do Compi = 1, NrComp
        Comp(Compi)%Theta = 0._sp
        Comp(Compi)%WFactor = 0._sp ! used for salt in (10*VolSat*dZ * EC)
        TotD = TotD + Comp(Compi)%Thickness
    end do
    Compi = 0
    Depthi = 0._sp
    AddComp = .true.
    Th2 = LocVolPr(1)
    EC2 = LocECds(1)
    D2 = 0._sp
    Loci = 0
    do while ((Compi < NrComp) .or. &
                ((Compi == NrComp) .and. (AddComp .eqv. .false.)))
        ! upper and lower boundaries location
        D1 = D2
        Th1 = Th2
        EC1 = EC2
        if (Loci < NrLoc) then
            Loci = Loci + 1
            D2 = LocDepth(Loci)
            Th2 = LocVolPr(Loci)
            EC2 = LocECdS(Loci)
        else
            D2 = TotD
        end if
        ! transfer water to compartment (SWC in mm) and salt in (10*VolSat*dZ * EC)
        TheEnd = .false.
        DTopComp = D1  ! Depthi is the bottom depth
        ThBotComp = Th1
        ECBotComp = EC1
        loop: do
            ThTopComp = ThBotComp
            ECTopComp = ECBotComp
            if (AddComp) then
                Compi = Compi + 1
                Depthi = Depthi + Comp(Compi)%Thickness
            end if
            if (Depthi < D2) then
                ThBotComp = Th1 + (Th2-Th1)*(Depthi-D1)/real(D2-D1, kind=sp)
                Comp(Compi)%Theta = Comp(Compi)%Theta &
                                   + 10._sp*(Depthi-DTopComp) &
                                           * ((ThTopComp+ThBotComp)/2._sp)
                ECBotComp = EC1 + (EC2-EC1)*(Depthi-D1)/real(D2-D1, kind=sp)
                Comp(Compi)%WFactor = Comp(Compi)%WFactor &
                                      + (10._sp*(Depthi-DTopComp) &
                                       *GetSoilLayer_SAT(Comp(Compi)%Layer)) &
                                       *((ECTopComp+ECbotComp)/2._sp)
                AddComp = .true.
                DTopComp = Depthi
                if (Compi == NrComp) then
                    TheEnd = .true.
                end if
            else
                ThBotComp = Th2
                ECBotComp = EC2
                Comp(Compi)%Theta = Comp(Compi)%Theta &
                                  + 10._sp*(D2-DTopComp) &
                                          * ((ThTopComp+ThBotComp)/2._sp)
                Comp(Compi)%WFactor = Comp(Compi)%WFactor &
                                      + (10._sp*(D2-DTopComp) &
                                         *GetSoilLayer_SAT(Comp(Compi)%Layer)) &
                                         *((ECTopComp+ECbotComp)/2._sp)
                if (abs(Depthi - D2) < epsilon(0._sp)) then
                    AddComp = .true.
                else
                    AddComp = .false.
                end if
                TheEnd = .true.
            end if
            if (TheEnd) exit loop
        end do loop
    end do

    do Compi = 1, NrComp
        ! from mm(water) to theta and final check
        Comp(Compi)%Theta = Comp(Compi)%Theta/(1000._sp*Comp(Compi)%Thickness)
        if (Comp(Compi)%Theta &
                    > (GetSoilLayer_SAT(Comp(compi)%Layer))/100._sp) then
            Comp(Compi)%Theta = (GetSoilLayer_SAT(Comp(compi)%Layer))/100._sp
        end if
        if (Comp(Compi)%Theta < 0._sp) then
            Comp(Compi)%Theta = 0._sp
        end if
    end do

    do Compi = 1, NrComp
        ! from (10*VolSat*dZ * EC) to ECe and distribution in cellls
        Comp(Compi)%WFactor = Comp(Compi)%WFactor &
                                /(10._sp*Comp(Compi)%Thickness &
                                    *GetSoilLayer_SAT(Comp(Compi)%Layer))
        call DetermineSaltContent(Comp(Compi)%WFactor, Comp(Compi))
    end do
end subroutine TranslateIniPointsToSWProfile


real(sp) function CCiniTotalFromTimeToCCini(TempDaysToCCini, TempGDDaysToCCini, &
                                            L0, L12, L12SF, L123, L1234, GDDL0, &
                                            GDDL12, GDDL12SF, GDDL123, &
                                            GDDL1234, CCo, CCx, CGC, GDDCGC, &
                                            CDC, GDDCDC, RatDGDD, SFRedCGC, &
                                            SFRedCCx, SFCDecline, fWeed, &
                                            TheModeCycle)
    integer(int32), intent(in) :: TempDaysToCCini
    integer(int32), intent(in) :: TempGDDaysToCCini
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: L12
    integer(int32), intent(in) :: L12SF
    integer(int32), intent(in) :: L123
    integer(int32), intent(in) :: L1234
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDL12
    integer(int32), intent(in) :: GDDL12SF
    integer(int32), intent(in) :: GDDL123
    integer(int32), intent(in) :: GDDL1234
    real(sp), intent(in) :: CCo
    real(sp), intent(in) :: CCx
    real(sp), intent(in) :: CGC
    real(sp), intent(in) :: GDDCGC
    real(sp), intent(in) :: CDC
    real(sp), intent(in) :: GDDCDC
    real(sp), intent(in) :: RatDGDD
    integer(int8), intent(in) :: SFRedCGC
    integer(int8), intent(in) :: SFRedCCx
    real(sp), intent(in) :: SFCDecline
    real(sp), intent(in) :: fWeed
    integer(intEnum), intent(in) :: TheModeCycle


    integer(int32) :: DayCC
    real(sp) :: SumGDDforCCini, TempCCini
    integer(int32) :: Tadj, GDDTadj

    if (TempDaysToCCini /= 0) then
        ! regrowth
        SumGDDforCCini = real(undef_int, kind=sp)
        GDDTadj = undef_int
        ! find adjusted calendar and GDD time
        if (TempDaysToCCini == undef_int) then
            ! CCx on 1st day
            Tadj = L12 - L0
            if (TheModeCycle == modeCycle_GDDays) then
                GDDTadj = GDDL12 - GDDL0
            end if
        else
            ! CC on 1st day is < CCx
            Tadj = TempDaysToCCini
            if (TheModeCycle == modeCycle_GDDays) then
                GDDTadj = TempGDDaysToCCini
            end if
        end if
        ! calculate CCini with adjusted time
        DayCC = L0 + Tadj
        if (TheModeCycle == modeCycle_GDDays) then
            SumGDDforCCini = GDDL0 + GDDTadj
        end if
        TempCCini = CCiNoWaterStressSF(DayCC, L0, L12SF, L123, L1234, GDDL0, &
                                       GDDL12SF, GDDL123, GDDL1234, &
                                       (CCo*fWeed), (CCx*fWeed), CGC, GDDCGC, &
                                       (CDC*(fWeed*CCx+2.29_sp)/(CCx+2.29_sp)), &
                                       (GDDCDC*(fWeed*CCx+2.29_sp)/(CCx+2.29_sp)), &
                                       SumGDDforCCini, RatDGDD, SFRedCGC, &
                                       SFRedCCx, SFCDecline, TheModeCycle)
        ! correction for fWeed is already in TempCCini (since DayCC > 0);
    else
        TempCCini = (CCo*fWeed) ! sowing or transplanting
    end if

    CCiniTotalFromTimeToCCini = TempCCini
end function CCiniTotalFromTimeToCCini


subroutine AdjustCropYearToClimFile(CDay1, CDayN)
    integer(int32), intent(inout) :: CDay1
    integer(int32), intent(inout) :: CDayN

    integer(int32) :: dayi, monthi, yeari
    character(len=:), allocatable :: temp_str

    call DetermineDate(CDay1, dayi, monthi, yeari)
    if (GetClimFile() == '(None)') then
        yeari = 1901  ! yeari = 1901 if undefined year
    else
        yeari = GetClimRecord_FromY() ! yeari = 1901 if undefined year
    end if
    call DetermineDayNr(dayi, monthi, yeari, CDay1)
    temp_str = EndGrowingPeriod(CDay1, CDayN)
end subroutine AdjustCropYearToClimFile


function EndGrowingPeriod(Day1, DayN) result(EndGrowingPeriod_out)
    integer(int32), intent(in) :: Day1
    integer(int32), intent(inout) :: DayN

    integer(int32) :: dayi, monthi, yeari
    character(len=2) :: Strday
    character(len=:), allocatable :: StrMonth
    character(len=:), allocatable :: EndGrowingPeriod_out

    ! This function determines Crop.DayN and the string
    DayN = Day1 + GetCrop_DaysToHarvest() - 1
    if (DayN < Day1) then
        DayN = Day1
    end if
    call DetermineDate(DayN, dayi, monthi, yeari)
    write(Strday, '(i2)') dayi
    StrMonth = NameMonth(monthi)
    EndGrowingPeriod_out = Strday // ' ' // StrMonth // '  '
end function EndGrowingPeriod


subroutine LoadInitialConditions(SWCiniFileFull, IniSurfaceStorage)
    character(len=*), intent(in) :: SWCiniFileFull
    real(sp), intent(inout) :: IniSurfaceStorage

    integer :: fhandle
    integer(int32) :: i
    character(len=1025) :: StringParam, swcinidescr_temp
    real(sp) :: VersionNr
    real(sp) :: CCini_temp, Bini_temp, Zrini_temp, ECStorageIni_temp
    integer(int8) :: NrLoc_temp
    real(sp) :: Loc_i_temp, VolProc_i_temp, SaltECe_i_temp

    ! IniSWCRead attribute of the function was removed to fix a
    ! bug occurring when the function was called in TempProcessing.pas
    ! Keep in mind that this could affect the graphical interface
    open(newunit=fhandle, file=trim(SWCiniFileFull), status='old', action='read')
    read(fhandle, '(a)') swcinidescr_temp
    call SetSWCiniDescription(swcinidescr_temp)
    read(fhandle, *) VersionNr ! AquaCrop Version
    if (roundc(10*VersionNr, mold=1) < 41) then ! initial CC at start of simulation period
        call SetSimulation_CCini(real(undef_int, kind=sp))
    else
        read(fhandle, *) CCini_temp
        call SetSimulation_CCini(CCini_temp)
    end if
    if (roundc(10*VersionNr, mold=1) < 41) then ! B produced before start of simulation period
        call SetSimulation_Bini(0.000_sp)
    else
        read(fhandle, *) Bini_temp
        call SetSimulation_Bini(Bini_temp)
    end if
    if (roundc(10*VersionNr, mold=1) < 41) then ! initial rooting depth at start of simulation period
        call SetSimulation_Zrini(real(undef_int, kind=sp))
    else
        read(fhandle, *) Zrini_temp
        call SetSimulation_Zrini(Zrini_temp)
    end if
    read(fhandle, *) IniSurfaceStorage
    if (roundc(10*VersionNr, mold=1) < 32) then ! EC of the ini surface storage
        call SetSimulation_ECStorageIni(0._sp)
    else
        read(fhandle, *) ECStorageIni_temp
        call SetSimulation_ECStorageIni(ECStorageIni_temp)
    end if
    read(fhandle, *) i
    if (i == 1) then
        call SetSimulation_IniSWC_AtDepths(.true.)
    else
        call SetSimulation_IniSWC_AtDepths(.false.)
    end if
    read(fhandle, *) NrLoc_temp
    call SetSimulation_IniSWC_NrLoc(NrLoc_temp)
    read(fhandle, *)
    read(fhandle, *)
    read(fhandle, *)
    do i = 1, GetSimulation_IniSWC_NrLoc()
        read(fhandle, '(a)') StringParam
        Loc_i_temp = GetSimulation_IniSWC_Loc_i(i)
        VolProc_i_temp = GetSimulation_IniSWC_VolProc_i(i)
        if (roundc(10*VersionNr, mold=1) < 32) then ! ECe at the locations
            call SplitStringInTwoParams(StringParam, Loc_i_temp, &
                                                    VolProc_i_temp)
            call SetSimulation_IniSWC_SaltECe_i(i, 0._sp)
        else
            SaltECe_i_temp = GetSimulation_IniSWC_SaltECe_i(i)
            call SplitStringInThreeParams(StringParam, Loc_i_temp, &
                                          VolProc_i_temp, SaltECe_i_temp)
            call SetSimulation_IniSWC_SaltECe_i(i, SaltECe_i_temp)
        end if
        call SetSimulation_IniSWC_Loc_i(i, Loc_i_temp)
        call SetSimulation_IniSWC_VolProc_i(i, VolProc_i_temp)
    end do
    close(fhandle)
    call SetSimulation_IniSWC_AtFC(.false.)
end subroutine LoadInitialConditions


subroutine AdjustSizeCompartments(CropZx)
    real(sp), intent(in) :: CropZx

    integer(int32) :: i , compi
    real(sp) :: TotDepthC, fAdd
    integer(int8) :: PrevNrComp
    real(sp), dimension(max_No_Compartments) :: PrevThickComp, &
                                                PrevVolPrComp, &
                                                PrevECdSComp

    ! 1. Save intial soil water profile (required when initial soil
    ! water profile is NOT reset at start simulation - see 7.)
    PrevNrComp = int(GetNrCompartments(), kind=int8)
    do compi = 1, PrevNrComp
        PrevThickComp(compi) = GetCompartment_Thickness(compi)
        PrevVolPrComp(compi) = 100._sp*GetCompartment_Theta(compi)
    end do

    ! 2. Actual total depth of compartments
    TotDepthC = 0._sp
    do i = 1, GetNrCompartments()
        TotDepthC = TotDepthC + GetCompartment_Thickness(i)
    end do

    ! 3. Increase number of compartments (if less than 12)
    if (GetNrCompartments() < 12) then
        loop: do
            call SetNrCompartments(GetNrCompartments() + 1)
            if ((CropZx - TotDepthC) > GetSimulParam_CompDefThick()) then
                call SetCompartment_Thickness(GetNrCompartments(), &
                                              GetSimulParam_CompDefThick())
            else
                call SetCompartment_Thickness(GetNrCompartments(), &
                                              CropZx - TotDepthC)
            end if
            TotDepthC = TotDepthC &
                        + GetCompartment_Thickness(GetNrCompartments())
            if ((GetNrCompartments() == max_No_compartments) &
                        .or. ((TotDepthC + 0.00001) >= CropZx)) exit loop
        end do loop
    end if

    ! 4. Adjust size of compartments (if total depth of compartments < rooting depth)
    if ((TotDepthC + 0.00001_sp) < CropZx) then
        call SetNrCompartments(12)
        fAdd = (CropZx/0.1_sp - 12._sp)/78._sp
        do i = 1, 12
            call SetCompartment_Thickness(i, 0.1 * (1._sp + i*fAdd))
            call SetCompartment_Thickness(i, 0.05 &
                    * real(roundc(GetCompartment_Thickness(i) &
                                    * 20._sp, mold=1), kind=sp))
        end do
        TotDepthC = 0._sp
        do i = 1, GetNrCompartments()
            TotDepthC = TotDepthC + GetCompartment_Thickness(i)
        end do
        if (CropZx - TotDepthC > ac_zero_threshold) then
            loop2: do
                call SetCompartment_Thickness(12, &
                                        GetCompartment_Thickness(12) &
                                                          + 0.05_sp)
                TotDepthC = TotDepthC + 0.05_sp
                if (TotDepthC >= CropZx) exit loop2
            end do loop2
        else
            do while ((TotDepthC - 0.04999999_sp) >= CropZx)
                call SetCompartment_Thickness(12, &
                                        GetCompartment_Thickness(12) &
                                                          - 0.05_sp)
                TotDepthC = TotDepthC - 0.05_sp
            end do
        end if
    end if

    ! 5. Adjust soil water content and theta initial
    call AdjustThetaInitial(PrevNrComp, PrevThickComp, &
                            PrevVolPrComp, PrevECdSComp)
end subroutine AdjustSizeCompartments


subroutine CheckForKeepSWC(RunWithKeepSWC, ConstZrxForRun)
    !! @NOTE This procedure will try to read from the soil profile file.
    !! If this file does not exist, the necessary information is gathered
    !! from the attributes of the Soil global variable instead.
    logical, intent(out) :: RunWithKeepSWC
    real(sp), intent(out) :: ConstZrxForRun

    integer(int32) :: fhandlex, i, Runi, TotalNrOfRuns
    character(len=1025) :: FileName, FullFileName
    real(sp) :: Zrni, Zrxi, ZrSoili
    real(sp) :: VersionNrCrop
    integer(int8) :: TheNrSoilLayers
    type(SoilLayerIndividual), dimension(max_SoilLayers) :: TheSoilLayer
    character(len=:), allocatable :: PreviousProfFilefull
    logical :: has_external

    ! 1. Initial settings
    RunWithKeepSWC = .false.
    ConstZrxForRun = real(undef_int, kind=sp)

    ! 2. Look for restrictive soil layer
    ! restricted to run 1 since with KeepSWC,
    ! the soil file has to be common between runs
    PreviousProfFilefull = GetProfFilefull() ! keep name soil file
                                             ! (to restore after check)

    Filename = ProjectInput(1)%Soil_Filename
    has_external = Filename == '(External)'

    if (has_external) then
        ! Note: here we use the AquaCrop version number and assume that
        ! the same version can be used in finalizing the soil settings.
        call LoadProfileProcessing(ProjectInput(1)%VersionNr)
    elseif (trim(FileName) == '(None)') then
        FullFileName = GetPathNameSimul() // 'DEFAULT.SOL'
        call LoadProfile(FullFileName)
    else
        FullFileName = ProjectInput(1)%Soil_Directory // FileName
        call LoadProfile(FullFileName)
    end if

    TheNrSoilLayers = GetSoil_NrSoilLayers()
    TheSoilLayer = GetSoilLayer()

    ! 3. Check if runs with KeepSWC exist
    Runi = 1
    TotalNrOfRuns = GetNumberSimulationRuns()
    do while ((RunWithKeepSWC .eqv. .false.) .and. (Runi <= TotalNrOfRuns))
        if (ProjectInput(Runi)%SWCIni_Filename == 'KeepSWC') then
            RunWithKeepSWC = .true.
        end if
        Runi = Runi + 1
    end do

    if (RunWithKeepSWC .eqv. .false.) then
        ConstZrxForRun = real(undef_int, kind=sp) ! reset
    end if

    ! 4. Look for maximum root zone depth IF RunWithKeepSWC
    if (RunWithKeepSWC .eqv. .true.) then
        Runi = 1
        do while (Runi <= TotalNrOfRuns)
            ! Obtain maximum rooting depth from the crop file
            FullFileName = ProjectInput(Runi)%Crop_Directory // &
                           ProjectInput(Runi)%Crop_Filename

            open(newunit=fhandlex, file=trim(FullFileName), &
                                    status='old', action ='read')
            read(fhandlex, *) ! description
            read(fhandlex, *) VersionNrCrop
            if (roundc(VersionNrCrop*10, mold=1) <= 31) then
                do i = 1, 29
                    read(fhandlex, *)  ! no Salinity stress
                        ! (No Reponse Stomata + ECemin + ECemax + ShapeKsSalt)
                end do
            else
                if (roundc(VersionNrCrop*10, mold=1) <= 50) then
                    do i = 1, 32
                        read(fhandlex, *) ! no distortion to salinity and
                                          ! response to ECsw factor
                    end do
                else
                    do i = 1, 34
                        read(fhandlex, *)
                    end do
                end if
            end if
            read(fhandlex, *) Zrni ! minimum rooting depth
            read(fhandlex, *) Zrxi ! maximum rooting depth
            close(fhandlex)

            ZrSoili = RootMaxInSoilProfile(Zrxi, TheNrSoilLayers, &
                                           TheSoilLayer)
            if (real(ZrSoili, kind=sp) > ConstZrxForRun) then
                ConstZrxForRun = ZrSoili
            end if

            Runi = Runi + 1
        end do
    end if

    ! 5. Reload existing soil file
    if (.not. has_external) then
        call SetProfFilefull(PreviousProfFilefull)
        call LoadProfile(GetProfFilefull())
    end if
end subroutine CheckForKeepSWC


subroutine InitializeGlobalStrings()
    !! Initializes all allocatable strings which are global variables
    !! themselves or attributes of derived type global variables.

    call SetCalendarDescription('')
    call SetCalendarFile('')
    call SetCalendarFileFull('')
    call SetClimateDescription('')
    call SetClimateFile('')
    call SetClimateFileFull('')
    call SetClimDescription('')
    call SetClimFile('')
    call SetClimRecord_FromString('')
    call SetClimRecord_ToString('')
    call SetCO2Description('')
    call SetCO2File('')
    call SetCO2FileFull('')
    call SetCropDescription('')
    call SetCropFile('')
    call SetCropFileFull('')
    call SetEToDescription('')
    call SetEToFile('')
    call SetEToFileFull('')
    call SetEToRecord_FromString('')
    call SetEToRecord_ToString('')
    call SetFullFileNameProgramParameters('')
    call SetGroundwaterDescription('')
    call SetGroundWaterFile('')
    call SetGroundWaterFilefull('')
    call SetIrriDescription('')
    call SetIrriFile('')
    call SetIrriFileFull('')
    call SetManDescription('')
    call SetManFile('')
    call SetManFilefull('')
    call SetMultipleProjectDescription('')
    call SetMultipleProjectFile('')
    call SetMultipleProjectFileFull('')
    call SetObservationsDescription('')
    call SetObservationsFile('')
    call SetObservationsFilefull('')
    call SetOffSeasonDescription('')
    call SetOffSeasonFile('')
    call SetOffSeasonFilefull('')
    call SetOutputName('')
    call SetPathNameList('')
    call SetPathNameOutp('')
    call SetPathNameParam('')
    call SetPathNameProg('')
    call SetPathNameSimul('')
    call SetProfDescription('')
    call SetProfFile('')
    call SetProfFilefull('')
    call SetProjectDescription('')
    call SetProjectFile('')
    call SetProjectFileFull('')
    call SetRainDescription('')
    call SetRainFile('')
    call SetRainFileFull('')
    call SetRainRecord_FromString('')
    call SetRainRecord_ToString('')
    call SetSimulation_Storage_CropString('')
    call SetSWCiniDescription('')
    call SetSWCiniFile('')
    call SetSWCiniFileFull('')
    call SetTemperatureDescription('')
    call SetTemperatureFile('')
    call SetTemperatureFileFull('')
    call SetTemperatureRecord_FromString('')
    call SetTemperatureRecord_ToString('')
    call SetTnxReference365DaysFile('')
    call SetTnxReference365DaysFileFull('')
    call SetTnxReferenceFile('')
    call SetTnxReferenceFileFull('')
end subroutine InitializeGlobalStrings


!! Global variables section !!

function GetIrriFile() result(str)
    !! Getter for the "IrriFile" global variable.
    character(len=:), allocatable :: str

    str = IrriFile
end function GetIrriFile


subroutine SetIrriFile(str)
    !! Setter for the "IrriFile" global variable.
    character(len=*), intent(in) :: str

    IrriFile = str
end subroutine SetIrriFile


function GetIrriFileFull() result(str)
    !! Getter for the "IrriFileFull" global variable.
    character(len=:), allocatable :: str

    str = IrriFileFull
end function GetIrriFileFull


subroutine SetIrriFileFull(str)
    !! Setter for the "IrriFileFull" global variable.
    character(len=*), intent(in) :: str

    IrriFileFull = str
end subroutine SetIrriFileFull


function GetClimateFile() result(str)
    !! Getter for the "ClimateFile" global variable.
    character(len=:), allocatable :: str

    str = ClimateFile
end function GetClimateFile


subroutine SetClimateFile(str)
    !! Setter for the "ClimateFile" global variable.
    character(len=*), intent(in) :: str

    ClimateFile = str
end subroutine SetClimateFile


function GetClimateFileFull() result(str)
    !! Getter for the "ClimateFileFull" global variable.
    character(len=:), allocatable :: str

    str = ClimateFileFull
end function GetClimateFileFull


subroutine SetClimateFileFull(str)
    !! Setter for the "ClimateFileFull" global variable.
    character(len=*), intent(in) :: str

    ClimateFileFull = str
end subroutine SetClimateFileFull


function GetClimateDescription() result(str)
    !! Getter for the "ClimateDescription" global variable.
    character(len=:), allocatable :: str

    str = ClimateDescription
end function GetClimateDescription


subroutine SetClimateDescription(str)
    !! Setter for the "ClimateDescription" global variable.
    character(len=*), intent(in) :: str

    ClimateDescription = str
end subroutine SetClimateDescription


function GetClimFile() result(str)
    !! Getter for the "ClimFile" global variable.
    character(len=:), allocatable :: str

    str = ClimFile
end function GetClimFile


subroutine SetClimFile(str)
    !! Setter for the "ClimFile" global variable.
    character(len=*), intent(in) :: str

    ClimFile = str
end subroutine SetClimFile


function GetSWCiniFile() result(str)
    !! Getter for the "SWCiniFile" global variable.
    character(len=:), allocatable :: str

    str = SWCiniFile
end function GetSWCiniFile


subroutine SetSWCiniFile(str)
    !! Setter for the "SWCiniFile" global variable.
    character(len=*), intent(in) :: str

    SWCiniFile = str
end subroutine SetSWCiniFile


function GetSWCiniDescription() result(str)
    !! Getter for the "SWCiniDescription" global variable.
    character(len=:), allocatable :: str

    str = SWCiniDescription
end function GetSWCiniDescription


subroutine SetSWCiniDescription(str)
    !! Setter for the "SWCiniDescription" global variable.
    character(len=*), intent(in) :: str

    SWCiniDescription = str
end subroutine SetSWCiniDescription


function GetSWCiniFileFull() result(str)
    !! Getter for the "SWCiniFileFull" global variable.
    character(len=:), allocatable :: str

    str = SWCiniFileFull
end function GetSWCiniFileFull


subroutine SetSWCiniFileFull(str)
    !! Setter for the "SWCiniFileFull" global variable.
    character(len=*), intent(in) :: str

    SWCiniFileFull = str
end subroutine SetSWCiniFileFull


function GetPathNameProg() result(str)
    !! Getter for the "PathNameProg" global variable.
    character(len=:), allocatable :: str

    str = PathNameProg
end function GetPathNameProg


subroutine SetPathNameProg(str)
    !! Setter for the "PathNameProg" global variable.
    character(len=*), intent(in) :: str

    PathNameProg = str
end subroutine SetPathNameProg


function GetPathNameOutp() result(str)
    !! Getter for the "PathNameOutp" global variable.
    character(len=:), allocatable :: str

    str = PathNameOutp
end function GetPathNameOutp


subroutine SetPathNameOutp(str)
    !! Setter for the "PathNameOutp" global variable.
    character(len=*), intent(in) :: str

    PathNameOutp = str
end subroutine SetPathNameOutp


function GetPathNameSimul() result(str)
    !! Getter for the "PathNameSimul" global variable.
    character(len=:), allocatable :: str

    str = PathNameSimul
end function GetPathNameSimul


subroutine SetPathNameSimul(str)
    !! Setter for the "PathNameSimul" global variable.
    character(len=*), intent(in) :: str

    PathNameSimul = str
end subroutine SetPathNameSimul


function GetProjectFile() result(str)
    !! Getter for the "ProjectFile" global variable.
    character(len=:), allocatable :: str

    str = ProjectFile
end function GetProjectFile


subroutine SetProjectFile(str)
    !! Setter for the "ProjectFile" global variable.
    character(len=*), intent(in) :: str

    ProjectFile = str
end subroutine SetProjectFile


function GetProjectFileFull() result(str)
    !! Getter for the "ProjectFileFull" global variable.
    character(len=:), allocatable :: str

    str = ProjectFileFull
end function GetProjectFileFull


subroutine SetProjectFileFull(str)
    !! Setter for the "ProjectFileFull" global variable.
    character(len=*), intent(in) :: str

    ProjectFileFull = str
end subroutine SetProjectFileFull


function GetMultipleProjectFile() result(str)
    !! Getter for the "MultipleProjectFile" global variable.
    character(len=:), allocatable :: str

    str = MultipleProjectFile
end function GetMultipleProjectFile


subroutine SetMultipleProjectFile(str)
    !! Setter for the "MultipleProjectFile" global variable.
    character(len=*), intent(in) :: str

    MultipleProjectFile = str
end subroutine SetMultipleProjectFile


function GetMultipleProjectFileFull() result(str)
    !! Getter for the "MultipleProjectFileFull" global variable.
    character(len=:), allocatable :: str

    str = MultipleProjectFileFull
end function GetMultipleProjectFileFull


subroutine SetMultipleProjectFileFull(str)
    !! Setter for the "MultipleProjectFileFull" global variable.
    character(len=*), intent(in) :: str

    MultipleProjectFileFull = str
end subroutine SetMultipleProjectFileFull


function GetFullFileNameProgramParameters() result(str)
    !! Getter for the "FullFileNameProgramParameters" global variable.
    character(len=:), allocatable :: str

    str = FullFileNameProgramParameters
end function GetFullFileNameProgramParameters


subroutine SetFullFileNameProgramParameters(str)
    !! Setter for the "FullFileNameProgramParameters" global variable.
    character(len=*), intent(in) :: str

    FullFileNameProgramParameters = str
end subroutine SetFullFileNameProgramParameters


logical function LeapYear(Year)
    integer(int32), intent(in) :: Year

    LeapYear = .false.
    if (frac(Year/4._sp) <= 0.01_sp) then
        LeapYear = .true.
    end if


    contains


    real(sp) function frac(val)
        real(sp), intent(in) :: val

        frac = val - floor(val)
    end function frac
end function LeapYear


subroutine ComposeOutputFileName(TheProjectFileName)
    character(len=*), intent(in) :: TheProjectFileName

    character(len=len(Trim(TheProjectFileName))) :: TempString
    character(len=:), allocatable :: TempString2
    integer(int32) :: i

    TempString = Trim(TheProjectFileName)
    i = len(TempString)
    TempString2 = TempString(1:i-4)
    call SetOutputName(TempString2)
end subroutine ComposeOutputFileName


subroutine GetFileForProgramParameters(TheFullFileNameProgram, FullFileNameProgramParameters)
    character(len=*), intent(in) :: TheFullFileNameProgram
    character(len=*), intent(inout) :: FullFileNameProgramParameters

    integer(int32) :: TheLength
    character(len=:), allocatable :: TheExtension

    FullFileNameProgramParameters = ''
    TheLength = len(TheFullFileNameProgram)
    TheExtension = TheFullFileNameProgram(TheLength-2:TheLength) ! PRO or PRM

    FullFileNameProgramParameters = TheFullFileNameProgram(1:TheLength-3)
    if (TheExtension == 'PRO') then
        FullFileNameProgramParameters = Trim(FullFileNameProgramParameters)//'PP1'
    else
        FullFileNameProgramParameters = Trim(FullFileNameProgramParameters)//'PPn'
    end if
end subroutine GetFileForProgramParameters


subroutine GlobalZero(SumWabal)
    type(rep_sum), intent(inout) :: SumWabal

    integer(int32) :: i

    SumWabal%Epot = 0.0_sp
    SumWabal%Tpot = 0.0_sp
    SumWabal%Rain = 0.0_sp
    SumWabal%Irrigation = 0.0_sp
    SumWabal%Infiltrated = 0.0_sp
    SumWabal%Runoff = 0.0_sp
    SumWabal%Drain = 0.0_sp
    SumWabal%Eact = 0.0_sp
    SumWabal%Tact = 0.0_sp
    SumWabal%TrW = 0.0_sp
    SumWabal%ECropCycle = 0.0_sp
    SumWabal%Biomass = 0._sp
    SumWabal%BiomassPot = 0._sp
    SumWabal%BiomassUnlim = 0._sp
    SumWabal%BiomassTot = 0._sp ! crop and weeds (for soil fertility stress)
    SumWabal%YieldPart = 0._sp
    SumWabal%SaltIn = 0._sp
    SumWabal%SaltOut = 0._sp
    SumWabal%CRwater = 0._sp
    SumWabal%CRsalt = 0._sp
    call SetTotalWaterContent_BeginDay(0._sp)

    do i =1, GetNrCompartments()
        call SetTotalWaterContent_BeginDay(GetTotalWaterContent_BeginDay() + &
        GetCompartment_theta(i)*1000._sp*GetCompartment_Thickness(i))
    end do
end subroutine GlobalZero


subroutine LoadProjectDescription(DescriptionOfProject)
    character(len=*), intent(out) :: DescriptionOfProject

    DescriptionOfProject = ProjectInput(1)%Description
end subroutine LoadProjectDescription


subroutine CheckFilesInProject(Runi, AllOK, FileOK)
    integer(int32), intent(in) :: Runi
    logical, intent(out) :: AllOK
    type(rep_FileOK), intent(out) :: FileOK

    logical :: FileOK_tmp
    
    AllOK = .true.
    FileOK_tmp = .true.

    ! Check the 14 files
    associate(input => ProjectInput(Runi))
    call check_file(input%Climate_Directory, input%Climate_Filename)
    FileOK%Climate_Filename = FileOK_tmp
    call check_file(input%Temperature_Directory, input%Temperature_Filename)
    FileOK%Temperature_Filename = FileOK_tmp
    call check_file(input%ETo_Directory, input%ETo_Filename)
    FileOK%ETo_Filename = FileOK_tmp
    call check_file(input%Rain_Directory, input%Rain_Filename)
    FileOK%Rain_Filename = FileOK_tmp
    call check_file(input%CO2_Directory, input%CO2_Filename)
    FileOK%CO2_Filename = FileOK_tmp
    call check_file(input%Calendar_Directory, input%Calendar_Filename)
    FileOK%Calendar_Filename = FileOK_tmp
    call check_file(input%Crop_Directory, input%Crop_Filename)
    FileOK%Crop_Filename = FileOK_tmp
    call check_file(input%Irrigation_Directory, input%Irrigation_Filename)
    FileOK%Irrigation_Filename = FileOK_tmp
    call check_file(input%Management_Directory, input%Management_Filename)
    FileOK%Management_Filename = FileOK_tmp
    call check_file(input%GroundWater_Directory, input%GroundWater_Filename)
    FileOK%GroundWater_Filename = FileOK_tmp
    call check_file(input%Soil_Directory, input%Soil_Filename)
    FileOK%Soil_Filename = FileOK_tmp

    if (ProjectInput(Runi)%SWCIni_Filename /= 'KeepSWC') then
        call check_file(input%SWCIni_Directory, input%SWCIni_Filename)
        FileOK%SWCIni_Filename = FileOK_tmp
    end if

    call check_file(input%OffSeason_Directory, input%OffSeason_Filename)
    FileOK%OffSeason_Filename = FileOK_tmp
    call check_file(input%Observations_Directory, input%Observations_Filename)
    FileOK%Observations_Filename = FileOK_tmp
    end associate


    contains


    subroutine check_file(directory, filename)
        ! Sets AllOK to false if expected file does not exist.
        character(len=*), intent(in) :: directory
        character(len=*), intent(in) :: filename

        if (filename /= '(None)') then
            if (.not. FileExists(directory // filename)) then
                AllOK = .false.
                FileOK_tmp = .false.
            else
                FileOK_tmp = .true.
            end if
        end if
    end subroutine check_file
end subroutine CheckFilesInProject


real(sp) function ActualRootingDepth(DAP, L0, LZmax, L1234, GDDL0, GDDLZmax, &
                                     SumGDD, Zmin, Zmax, ShapeFactor,&
                                     TypeDays)
    integer(int32), intent(in) :: DAP
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: LZmax
    integer(int32), intent(in) :: L1234
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDLZmax
    real(sp), intent(in) :: SumGDD
    real(sp), intent(in) :: Zmin
    real(sp), intent(in) :: Zmax
    integer(int8), intent(in) :: ShapeFactor
    integer(intEnum), intent(in) :: TypeDays

    real(sp) :: Zini, Zr
    integer(int32) :: VirtualDay, T0, rootmax_rounded, zmax_rounded

    select case (TypeDays)
    case (modeCycle_GDDays)
        Zr = ActualRootingDepthGDDays(DAP, L1234, GDDL0, GDDLZmax, SumGDD, &
                                      Zmin, Zmax)
    case default
        Zr = ActualRootingDepthDays(DAP, L0, LZmax, L1234, Zmin, Zmax)
    end select

    ! restrictive soil layer
    call SetSimulation_SCor(1._sp)

    rootmax_rounded = roundc(real(GetSoil_RootMax()*1000, kind=sp), &
                             mold=rootmax_rounded)
    zmax_rounded = roundc(Zmax*1000, mold=zmax_rounded)

    if (rootmax_rounded < zmax_rounded) then
        call ZrAdjustedToRestrictiveLayers(Zr, GetSoil_NrSoilLayers(), &
                                           GetSoilLayer(), Zr)
    end if

    ActualRootingDepth = Zr


    contains


    real(sp) function ActualRootingDepthDays(DAP, L0, LZmax, L1234, Zmin, Zmax)
        integer(int32), intent(in) :: DAP
        integer(int32), intent(in) :: L0
        integer(int32), intent(in) :: LZmax
        integer(int32), intent(in) :: L1234
        real(sp), intent(in) :: Zmin
        real(sp), intent(in) :: Zmax

        ! Actual rooting depth at the end of Dayi
        VirtualDay = DAP - GetSimulation_DelayedDays()

        if ((VirtualDay < 1) .or. (VirtualDay > L1234)) then
            ActualRootingDepthDays = 0
        elseif (VirtualDay >= LZmax) then
            ActualRootingDepthDays = Zmax
        elseif (Zmin < Zmax) then
            Zini = ZMin * (GetSimulParam_RootPercentZmin()/100._sp)
            T0 = roundc(L0/2._sp, mold=T0)

            if (LZmax <= T0) then
                Zr = Zini + (Zmax-Zini)*VirtualDay*1._sp/LZmax
            elseif (VirtualDay <= T0) then
                Zr = Zini
            else
                Zr = Zini + (Zmax-Zini) &
                     * TimeRootFunction(real(VirtualDay, kind=sp), ShapeFactor,&
                                        real(LZmax, kind=sp), real(T0, kind=sp))
            end if

            if (Zr > ZMin) then
                ActualRootingDepthDays = Zr
            else
                ActualRootingDepthDays = ZMin
            end if
        else
            ActualRootingDepthDays = ZMax
        end if
    end function ActualRootingDepthDays


    real(sp) function ActualRootingDepthGDDays(DAP, L1234, GDDL0, GDDLZmax, &
                                               SumGDD, Zmin, Zmax)
        integer(int32), intent(in) :: DAP
        integer(int32), intent(in) :: L1234
        integer(int32), intent(in) :: GDDL0
        integer(int32), intent(in) :: GDDLZmax
        real(sp), intent(in) :: SumGDD
        real(sp), intent(in) :: Zmin
        real(sp), intent(in) :: Zmax

        real(sp) :: GDDT0

        ! after sowing the crop has roots even when SumGDD = 0
        VirtualDay = DAP - GetSimulation_DelayedDays()

        if ((VirtualDay < 1) .or. (VirtualDay > L1234)) then
            ActualRootingDepthGDDays = 0
        elseif (SumGDD >= GDDLZmax) then
            ActualRootingDepthGDDays = Zmax
        elseif (Zmin < Zmax) then
            Zini = ZMin * (GetSimulParam_RootPercentZmin()/100._sp)
            GDDT0 = GDDL0/2._sp

            if (GDDLZmax <= GDDT0) then
                Zr = Zini + (Zmax-Zini)*SumGDD/GDDLZmax
            else
                if (SumGDD <= GDDT0) then
                    Zr = Zini
                else
                    Zr = Zini + (Zmax-Zini) &
                         * TimeRootFunction(SumGDD, ShapeFactor, &
                                            real(GDDLZmax, kind=sp), GDDT0)
                end if
            end if

            if (Zr > ZMin) then
                ActualRootingDepthGDDays = Zr
            else
                ActualRootingDepthGDDays = ZMin
            end if
        else
            ActualRootingDepthGDDays = ZMax
        end if
    end function ActualRootingDepthGDDays
end function ActualRootingDepth


subroutine DesignateSoilLayerToCompartments(NrCompartments, NrSoilLayers, &
                                            Compartment)
    integer(int32), intent(in) :: NrCompartments
    integer(int32), intent(in) :: NrSoilLayers
    type(CompartmentIndividual), dimension(max_No_compartments), &
                                        intent(inout) :: Compartment

    integer(int32) :: i, layeri, compi
    real(sp) :: depth, depthi
    logical :: finished, NextLayer

    depth = 0._sp
    depthi = 0._sp
    layeri = 1
    compi = 1

    outer_loop: do
        depth = depth + GetSoilLayer_Thickness(layeri)
        inner_loop: do
            depthi = depthi + Compartment(compi)%Thickness/2._sp

            if (depthi <= depth) then
                Compartment(compi)%Layer = layeri
                NextLayer = .false.
                depthi = depthi + Compartment(compi)%Thickness/2._sp
                compi = compi + 1
                finished = (compi > NrCompartments)
            else
                depthi = depthi - Compartment(compi)%Thickness/2._sp
                NextLayer = .true.
                layeri = layeri + 1
                finished = (layeri > NrSoilLayers)
            end if

            if (finished .or. NextLayer) exit inner_loop
            end do inner_loop

        if (finished) exit outer_loop
        end do outer_loop

    do i = compi, NrCompartments
        Compartment(i)%Layer = NrSoilLayers
    end do

    do i = (NrCompartments+1), max_No_compartments
        Compartment(i)%Thickness = undef_double
    end do
end subroutine DesignateSoilLayerToCompartments


subroutine specify_soil_layer(NrCompartments, NrSoilLayers, SoilLayer, &
                              Compartment, TotalWaterContent)
    integer(int32), intent(in) :: NrCompartments
    integer(int32), intent(in) :: NrSoilLayers
    type(SoilLayerIndividual), dimension(max_SoilLayers), &
                                                intent(inout) :: SoilLayer
    type(CompartmentIndividual), dimension(max_No_compartments), &
                                                intent(inout) :: Compartment
    type(rep_Content), intent(inout) :: TotalWaterContent

    integer(int32) :: layeri, compi, celli
    real(sp) :: Total

    call DesignateSoilLayerToCompartments(NrCompartments, NrSoilLayers, &
                                          Compartment)

    ! Set soil layers and compartments at Field Capacity and determine Watercontent (mm)
    ! No salinity in soil layers and compartmens
    ! Absence of ground water table (FCadj = FC)
    Total = 0._sp
    do layeri = 1, NrSoilLayers
        SoilLayer(layeri)%WaterContent = 0._sp
    end do

    do compi = 1, NrCompartments
        Compartment(compi)%Theta = SoilLayer(Compartment(compi)%Layer)%FC/100._sp
        Compartment(compi)%FCadj = SoilLayer(Compartment(compi)%Layer)%FC
        Compartment(compi)%DayAnaero = 0

        do celli = 1, SoilLayer(Compartment(compi)%Layer)%SCP1
            ! salinity in cells
            Compartment(compi)%Salt(celli) = 0.0_sp
            Compartment(compi)%Depo(celli) = 0.0_sp
        end do

        call SetSimulation_ThetaIni_i(compi, Compartment(compi)%Theta)
        call SetSimulation_ECeIni_i(compi, 0._sp) ! initial soil salinity in dS/m

        SoilLayer(Compartment(compi)%Layer)%WaterContent = &
                SoilLayer(Compartment(compi)%Layer)%WaterContent &
                + GetSimulation_ThetaIni_i(compi)*100._sp &
                             *10._sp*Compartment(compi)%Thickness
    end do

    do layeri = 1, NrSoilLayers
        Total = Total + SoilLayer(layeri)%WaterContent
    end do
    call SetTotalWaterContent_BeginDay(Total)

    ! initial soil water content and no salts
    call DeclareInitialCondAtFCandNoSalt()

    ! Number of days with RootZone Anaerobic Conditions
    call SetSimulation_DayAnaero(0_int8)
end subroutine specify_soil_layer


subroutine Calculate_Saltmobility(layer, SaltDiffusion, Macro, Mobil)
    integer(int32), intent(in) :: layer
    integer(int8), intent(in) :: SaltDiffusion
    integer(int8), intent(in) :: Macro
    real(sp), dimension(11), intent(inout) :: Mobil

    integer(int32) :: i, CelMax
    real(sp) :: Mix, a, b, xi, yi, UL

    Mix = SaltDiffusion/100._sp ! global salt mobility expressed as a fraction
    UL = GetSoilLayer_UL(layer) * 100._sp ! upper limit in VOL% of SC cell

    ! 1. convert Macro (vol%) in SaltCelNumber
    if (Macro > UL) then
        CelMax = GetSoilLayer_SCP1(layer)
    else
        CelMax = roundc((Macro/UL)*GetSoilLayer_SC(layer), mold=1)
    end if

    if (CelMax <= 0) then
        CelMax = 1
    end if

    ! 2. find a and b
    if (Mix < 0.5_sp) then
        a = Mix * 2._sp
        b = exp(10._sp*(0.5_sp-Mix)*log(10._sp))
    else
        a = 2._sp * (1._sp - Mix)
        b = exp(10._sp*(Mix-0.5_sp)*log(10._sp))
    end if

    ! 3. calculate mobility for cells = 1 to Macro
    do i = 1, CelMax-1
        xi = i *1._sp / (CelMax-1)

        if (Mix > 0._sp) then
            if (Mix < 0.5_sp) then
                yi = exp(log(a)+xi*log(b))
                Mobil(i) = (yi-a)/(a*b-a)
            elseif ((Mix >= 0.5_sp - epsilon(0.0_sp)) &
                    .and. (Mix <= 0.5_sp + epsilon(0.0_sp)))  then
                Mobil(i) = xi
            elseif (Mix < 1._sp) then
                yi = exp(log(a)+(1._sp-xi)*log(b))
                Mobil(i) = 1._sp - (yi-a)/(a*b-a)
            else
                Mobil(i) = 1._sp
            end if
        else
            Mobil(i) = 0._sp
        end if
    end do

    ! 4. Saltmobility between Macro and SAT
    do i = CelMax, GetSoilLayer_SCP1(layer)
        Mobil(i) = 1._sp
    end do
end subroutine Calculate_Saltmobility


subroutine CompleteProfileDescription()

    integer(int32) :: i
    type(rep_Content) :: TotalWaterContent_temp
    type(CompartmentIndividual), &
                dimension(max_No_compartments) :: Compartment_temp
    type(SoilLayerIndividual) :: soillayer_i_temp
    type(SoilLayerIndividual), dimension(max_SoilLayers) :: soillayer_temp

    do i = GetSoil_NrSoilLayers()+1, max_SoilLayers
        soillayer_i_temp = GetSoilLayer_i(i)
        call set_layer_undef(soillayer_i_temp)
        call SetSoilLayer_i(i, soillayer_i_temp)
    end do

    call SetSimulation_ResetIniSWC(.true.) ! soil water content and soil salinity
    TotalWaterContent_temp = GetTotalWaterContent()
    Compartment_temp = GetCompartment()
    soillayer_temp = GetSoilLayer()
    call specify_soil_layer(GetNrCompartments(), &
                            int(GetSoil_NrSoilLayers(), kind=int32), &
                            soillayer_temp, Compartment_temp, &
                            TotalWaterContent_temp)
    call SetSoilLayer(soillayer_temp)
    call SetTotalWaterContent(TotalWaterContent_temp)
    call SetCompartment(Compartment_temp)
end subroutine CompleteProfileDescription


subroutine LoadProfile(FullName)
    !! Reads in soil data from the given file.
    !! Further initializations happen via a call to LoadProfileProcessing().
    character(len=*), intent(in) :: FullName

    integer :: fhandle
    integer(int32) :: i
    character(len=3) :: blank
    real(sp) :: VersionNr
    integer(int8) :: TempShortInt
    character(len=1024) :: ProfDescriptionLocal
    real(sp) :: thickness_temp, SAT_temp, FC_temp, WP_temp, infrate_temp
    real(sp) :: cra_temp, crb_temp
    character(len=25) :: description_temp
    integer(int8) :: penetrability_temp, gravelm_temp

    open(newunit=fhandle, file=trim(FullName), status='old', action='read')
    read(fhandle, '(a)') ProfDescriptionLocal
    call SetProfDescription(trim(ProfDescriptionLocal))
    read(fhandle, *) VersionNr  ! AquaCrop version
    read(fhandle, *) TempShortInt
    call SetSoil_CNvalue(TempShortInt)
    read(fhandle, *) TempShortInt
    call SetSoil_REW(TempShortInt)
    read(fhandle, *) TempShortInt
    call SetSoil_NrSoilLayers(TempShortInt)
    read(fhandle, *) ! depth of restrictive soil layer which is no longer applicable
    read(fhandle, *)
    read(fhandle, *)

    ! Load characteristics of each soil layer
    do i = 1, GetSoil_NrSoilLayers()
        ! Parameters for capillary rise missing in Versions 3.0 and 3.1
        if (roundc(VersionNr*10, mold=1) < 40) then
            read(fhandle, *) thickness_temp, SAT_temp, FC_temp, &
                             WP_temp, infrate_temp, blank, description_temp
            call SetSoilLayer_Thickness(i, thickness_temp)
            call SetSoilLayer_SAT(i, SAT_temp)
            call SetSoilLayer_FC(i, FC_temp)
            call SetSoilLayer_WP(i, WP_temp)
            call SetSoilLayer_InfRate(i, infrate_temp)
            call SetSoilLayer_Description(i, description_temp)
            ! Default values for Penetrability and Gravel
            call SetSoilLayer_Penetrability(i, 100_int8)
            call SetSoilLayer_GravelMass(i, 0_int8)
            ! determine volume gravel
            call SetSoilLayer_GravelVol(i, 0._sp)
        else
            if (roundc(VersionNr*10, mold=1) < 60) then
                            ! UPDATE required for Version 6.0
                read(fhandle, *) thickness_temp, SAT_temp, FC_temp, &
                                 WP_temp, infrate_temp, cra_temp, &
                                 crb_temp, blank, description_temp
                call SetSoilLayer_Thickness(i, thickness_temp)
                call SetSoilLayer_SAT(i, SAT_temp)
                call SetSoilLayer_FC(i, FC_temp)
                call SetSoilLayer_WP(i, WP_temp)
                call SetSoilLayer_InfRate(i, infrate_temp)
                call SetSoilLayer_CRa(i, cra_temp)
                call SetSoilLayer_CRb(i, crb_temp)
                call SetSoilLayer_Description(i, description_temp)
                ! Default values for Penetrability and Gravel
                call SetSoilLayer_Penetrability(i, 100_int8)
                call SetSoilLayer_GravelMass(i, 0_int8)
                ! determine volume gravel
                call SetSoilLayer_GravelVol(i, 0._sp)
            else
                read(fhandle, *) thickness_temp, SAT_temp, FC_temp, WP_temp, &
                                 infrate_temp, penetrability_temp, &
                                 gravelm_temp, cra_temp, crb_temp, &
                                 description_temp
                call SetSoilLayer_Thickness(i, thickness_temp)
                call SetSoilLayer_SAT(i, SAT_temp)
                call SetSoilLayer_FC(i, FC_temp)
                call SetSoilLayer_WP(i, WP_temp)
                call SetSoilLayer_InfRate(i, infrate_temp)
                call SetSoilLayer_Penetrability(i, penetrability_temp)
                call SetSoilLayer_GravelMass(i, gravelm_temp)
                call SetSoilLayer_CRa(i, cra_temp)
                call SetSoilLayer_CRb(i, crb_temp)
                call SetSoilLayer_Description(i, description_temp)
                ! determine volume gravel
                call SetSoilLayer_GravelVol(i, &
                            FromGravelMassToGravelVolume(GetSoilLayer_SAT(i), &
                                                    GetSoilLayer_GravelMass(i)))
            end if
        end if
    end do

    close(fhandle)
    call LoadProfileProcessing(VersionNr)
end subroutine LoadProfile


subroutine LoadProfileProcessing(VersionNr)
    !! Further initializations after soil profile attributes have been set
    !! (e.g. via a call to LoadProfile()).
    real(sp), intent(in) :: VersionNr
        !! AquaCrop Version (e.g. 7.0)

    integer(int32) :: i
    real(sp) :: cra_temp, crb_temp, dx_temp
    real(sp), dimension(11) :: saltmob_temp

    call SetSimulation_SurfaceStorageIni(0.0_sp)
    call SetSimulation_ECStorageIni(0.0_sp)

    do i = 1, GetSoil_NrSoilLayers()
        ! determine drainage coefficient
        call SetSoilLayer_tau(i, TauFromKsat(GetSoilLayer_InfRate(i)))

        ! determine number of salt cells based on infiltration rate
        if (GetSoilLayer_InfRate(i) <= 112._sp) then
            call SetSoilLayer_SCP1(i, 11_int8)
        else
            call SetSoilLayer_SCP1(i, &
                            roundc(1.6_sp + 1000._sp/GetSoilLayer_InfRate(i), &
                                   mold=1_int8))
            if (GetSoilLayer_SCP1(i) < 2_int8) then
                call SetSoilLayer_SCP1(i, 2_int8)
            end if
        end if

        ! determine parameters for soil salinity
        call SetSoilLayer_SC(i, GetSoilLayer_SCP1(i) - 1_int8)
        call SetSoilLayer_Macro(i, roundc(GetSoilLayer_FC(i), mold=1_int8))
        call SetSoilLayer_UL(i, ((GetSoilLayer_SAT(i))/100._sp) &
                    * (GetSoilLayer_SC(i)/(GetSoilLayer_SC(i)+2._sp))) ! m3/m3
        dx_temp = (GetSoilLayer_UL(i))/GetSoilLayer_SC(i)
        call SetSoilLayer_Dx(i, dx_temp)  ! m3/m3

        saltmob_temp = GetSoilLayer_SaltMobility(i)
        call Calculate_SaltMobility(i, GetSimulParam_SaltDiff(), &
                                    GetSoilLayer_Macro(i), saltmob_temp)
        call SetSoilLayer_SaltMobility(i, saltmob_temp)


        ! determine default parameters for capillary rise if missing
        call SetSoilLayer_SoilClass(i, NumberSoilClass(GetSoilLayer_SAT(i), &
                                    GetSoilLayer_FC(i), GetSoilLayer_WP(i), &
                                    GetSoilLayer_InfRate(i)))

        if (roundc(VersionNr*10, mold=1) < 40) then
            cra_temp = GetSoilLayer_CRa(i)
            crb_temp = GetSoilLayer_CRb(i)
            call DetermineParametersCR(GetSoilLayer_SoilClass(i), &
                                       GetSoilLayer_InfRate(i), &
                                       cra_temp, crb_temp)
            call SetSoilLayer_CRa(i, cra_temp)
            call SetSoilLayer_CRb(i, crb_temp)
        end if
    end do

    call DetermineNrandThicknessCompartments()
    call SetSoil_RootMax(RootMaxInSoilProfile(GetCrop_RootMax(), &
                                              GetSoil_NrSoilLayers(), &
                                              GetSoilLayer()))
end subroutine LoadProfileProcessing


subroutine DetermineRootZoneWC(RootingDepth, ZtopSWCconsidered)
    real(sp), intent(in) :: RootingDepth
    logical, intent(inout) :: ZtopSWCconsidered


    real(sp) :: CumDepth, Factor, frac_value, DrRel, DZtopRel, TopSoilInMeter
    integer(int32) :: compi

    ! calculate SWC in root zone
    CumDepth = 0._sp
    compi = 0
    call SetRootZoneWC_Actual(0._sp)
    call SetRootZoneWC_FC(0._sp)
    call SetRootZoneWC_WP(0._sp)
    call SetRootZoneWC_SAT(0._sp)
    call SetRootZoneWC_Leaf(0._sp)
    call SetRootZoneWC_Thresh(0._sp)
    call SetRootZoneWC_Sen(0._sp)
    loop: do
        compi = compi + 1
        CumDepth = CumDepth + GetCompartment_Thickness(compi)
        if (CumDepth <= RootingDepth) then
            Factor = 1._sp
        else
            frac_value = RootingDepth - (CumDepth - GetCompartment_Thickness(compi))
            if (frac_value > 0._sp) then
                Factor = frac_value/GetCompartment_Thickness(compi)
            else
                Factor = 0._sp
            end if
        end if
        call SetRootZoneWC_Actual(GetRootZoneWC_Actual() + Factor * 1000._sp * &
             GetCompartment_Theta(compi) * GetCompartment_Thickness(compi)* &
             (1._sp - GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._sp))
        call SetRootZoneWC_FC(GetRootZoneWC_FC() + Factor * 10._sp * &
             GetSoilLayer_FC(GetCompartment_Layer(compi)) * &
             GetCompartment_Thickness(compi) * (1._sp - &
             GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._sp))
        call SetRootZoneWC_Leaf(GetRootZoneWC_Leaf() + Factor * 10._sp * &
             GetCompartment_Thickness(compi) * &
             (GetSoilLayer_FC(GetCompartment_Layer(compi)) - GetCrop_pLeafAct() &
             * (GetSoilLayer_FC(GetCompartment_Layer(compi))- &
             GetSoilLayer_WP(GetCompartment_Layer(compi)))) * (1._sp &
             - GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._sp))
        call SetRootZoneWC_Thresh(GetRootZoneWC_Thresh() + Factor * 10._sp * &
             GetCompartment_Thickness(compi) * &
             (GetSoilLayer_FC(GetCompartment_Layer(compi)) - GetCrop_pActStom() &
             * (GetSoilLayer_FC(GetCompartment_Layer(compi))- &
             GetSoilLayer_WP(GetCompartment_Layer(compi)))) * (1._sp - &
             GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._sp))
        call SetRootZoneWC_Sen(GetRootZoneWC_Sen() + Factor * 10._sp * &
             GetCompartment_Thickness(compi) * &
             (GetSoilLayer_FC(GetCompartment_Layer(compi)) - GetCrop_pSenAct() &
             * (GetSoilLayer_FC(GetCompartment_Layer(compi))- &
             GetSoilLayer_WP(GetCompartment_Layer(compi)))) * (1._sp - &
             GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._sp))
        call SetRootZoneWC_WP(GetRootZoneWC_WP() + Factor * 10._sp * &
             GetSoilLayer_WP(GetCompartment_Layer(compi))* &
             GetCompartment_Thickness(compi)* (1._sp - &
             GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._sp))
        call SetRootZoneWC_SAT(GetRootZoneWC_SAT()+ Factor * 10._sp * &
             GetSoilLayer_SAT(GetCompartment_Layer(compi)) * &
             GetCompartment_Thickness(compi) * (1._sp - &
             GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._sp))
        if ((CumDepth >= RootingDepth) .or. (compi == NrCompartments)) exit loop
    end do loop
    ! calculate SWC in top soil (top soil in meter = SimulParam.ThicknessTopSWC/100)
    if ((RootingDepth*100._sp) <= GetSimulParam_ThicknessTopSWC()) then
        call SetRootZoneWC_ZtopAct(GetRootZoneWC_Actual())
        call SetRootZoneWC_ZtopFC(GetRootZoneWC_FC())
        call SetRootZoneWC_ZtopWP(GetRootZoneWC_WP())
        call SetRootZoneWC_ZtopThresh(GetRootZoneWC_Thresh())
    else
        CumDepth = 0._sp
        compi = 0
        call SetRootZoneWC_ZtopAct(0._sp)
        call SetRootZoneWC_ZtopFC(0._sp)
        call SetRootZoneWC_ZtopWP(0._sp)
        call SetRootZoneWC_ZtopThresh(0._sp)
        TopSoilInMeter = GetSimulParam_ThicknessTopSWC()/100._sp
        loop_2: do
            compi = compi + 1
            CumDepth = CumDepth + GetCompartment_Thickness(compi)
            if ((CumDepth*100._sp) <= GetSimulParam_ThicknessTopSWC()) then
                Factor = 1._sp
            else
                frac_value = TopSoilInMeter - (CumDepth - GetCompartment_Thickness(compi))
                if (frac_value > 0._sp) then
                    Factor = frac_value/GetCompartment_Thickness(compi)
                else
                    Factor = 0._sp
                end if
            end if
            call SetRootZoneWC_ZtopAct(GetRootZoneWC_ZtopAct() + Factor * &
                 1000._sp * GetCompartment_Theta(compi) * &
                 GetCompartment_Thickness(compi) * (1._sp - &
                 GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._sp))
            call SetRootZoneWC_ZtopFC(GetRootZoneWC_ZtopFC() + Factor * 10._sp &
                 * GetSoilLayer_FC(GetCompartment_Layer(compi)) * &
                 GetCompartment_Thickness(compi) * (1._sp - &
                 GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._sp))
            call SetRootZoneWC_ZtopWP(GetRootZoneWC_ZtopWP() + Factor * 10._sp &
                 * GetSoilLayer_WP(GetCompartment_Layer(compi)) * &
                 GetCompartment_Thickness(compi) * (1._sp - &
                 GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._sp))
            call SetRootZoneWC_ZtopThresh(GetRootZoneWC_ZtopThresh() + Factor * &
                 10._sp * GetCompartment_Thickness(compi) * &
                 (GetSoilLayer_FC(GetCompartment_Layer(compi))  - &
                 GetCrop_pActStom() * (GetSoilLayer_FC(GetCompartment_Layer(compi)) &
                 -GetSoilLayer_WP(GetCompartment_Layer(compi)))) * (1._sp - &
                 GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._sp))
            if ((CumDepth >= TopSoilInMeter) .or. (compi == NrCompartments)) exit loop_2
        end do loop_2
    end if

    ! Relative depletion in rootzone and in top soil
    if (roundc(1000._sp*(GetRootZoneWc_FC() - GetRootZoneWc_WP()), mold=1) > 0) then
        DrRel = (GetRootZoneWc_FC() - GetRootZoneWC_Actual())/(GetRootZoneWc_FC() &
                                                          - GetRootZoneWc_WP())
    else
        DrRel = 0._sp
    end if
    if (roundc(1000._sp*(GetRootZoneWC_ZtopFC() - GetRootZoneWc_ZtopWP()),mold=1) > 0) then
        DZtopRel = (GetRootZoneWC_ZtopFC() - GetRootZoneWc_ZtopAct())/ &
                              (GetRootZoneWC_ZtopFC() - GetRootZoneWc_ZtopWP())
    else
        DZtopRel = 0._sp
    end if

    ! Zone in soil profile considered for determining stress response
    if (DZtopRel < DrRel) then
        ZtopSWCconsidered = .true.  ! top soil is relative wetter than root zone
    else
        ZtopSWCconsidered = .false.
    end if
end subroutine DetermineRootZoneWC


subroutine CalculateETpot(DAP, L0, L12, L123, LHarvest, DayLastCut, CCi, &
                          EToVal, KcVal, KcDeclineVal, CCx, CCxWithered, &
                          CCeffectProcent, CO2i, GDDayi, TempGDtranspLow, &
                          TpotVal, EpotVal)
    integer(int32), intent(in) :: DAP
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: L12
    integer(int32), intent(in) :: L123
    integer(int32), intent(in) :: LHarvest
    integer(int32), intent(in) :: DayLastCut
    real(sp), intent(in) :: CCi
    real(sp), intent(in) :: EToVal
    real(sp), intent(in) :: KcVal
    real(sp), intent(in) :: KcDeclineVal
    real(sp), intent(in) :: CCx
    real(sp), intent(in) :: CCxWithered
    real(sp), intent(in) :: CCeffectProcent
    real(sp), intent(in) :: CO2i
    real(sp), intent(in) :: GDDayi
    real(sp), intent(in) :: TempGDtranspLow
    real(sp), intent(inout) :: TpotVal
    real(sp), intent(inout) :: EpotVal

    real(sp) :: KcVal_local
    real(sp) :: EpotMin, EpotMax, CCiAdjusted, Multiplier, KsTrCold
    integer(int32) :: VirtualDay

    ! CalculateETpot
    VirtualDay = DAP - GetSimulation_DelayedDays()
    if (((VirtualDay < L0) .and. (roundc(100._sp*CCi, mold=1) == 0)) &
                          .or. (VirtualDay > LHarvest)) then
        ! To handlle Forage crops: Round(100*CCi) = 0
        TpotVal = 0._sp
        EpotVal = GetSimulParam_KcWetBare()*EToVal
    else
        ! Correction for micro-advection
        CCiAdjusted = 1.72_sp*CCi - 1._sp*(CCi*CCi) + 0.30_sp*(CCi*CCi*CCi)
        if (CCiAdjusted < epsilon(1._sp)) then
            CCiAdjusted = 0._sp
        end if
        if (CCiAdjusted > 1._sp) then
            CCiAdjusted = 1._sp
        end if

        ! Correction for ageing effects - is a function of calendar days
        if ((VirtualDay-DayLastCut) > (L12+5)) then
            KcVal_local = KcVal - (VirtualDay-DayLastCut-(L12+5._sp)) &
                        * (KcDeclineVal/100._sp)*CCxWithered
        else
            KcVal_local = KcVal
        end if

        ! Correction for elevated atmospheric CO2 concentration
        if (CO2i > 369.41_sp) then
            KcVal_local = KcVal_local * &
                    (1._sp - 0.05_sp * (CO2i-369.41_sp)/(550._sp-369.41_sp))
        end if

        ! Correction for Air temperature stress
        if ((CCiAdjusted <= ac_zero_threshold) &
            .or. (roundc(GDDayi, mold=1) < 0)) then
            KsTrCold = 1._sp
        else
            KsTrCold = KsTemperature(0._sp, TempGDtranspLow, GDDayi)
        end if

        ! First estimate of Epot and Tpot
        TpotVal = CCiAdjusted * KsTrCold * KcVal_local * EToVal
        EpotVal = GetSimulParam_KcWetBare() * (1._sp - CCiAdjusted) * EToVal

        ! Maximum Epot with withered canopy as a result of (early) senescence
        EpotMax = GetSimulParam_KcWetBare() * EToVal * &
                        (1._sp - CCxWithered * CCEffectProcent/100._sp)

        ! Correction Epot for dying crop in late-season stage
        if ((VirtualDay > L123) .and. (CCx > epsilon(1._sp))) then
            if (CCi > (CCx/2._sp)) then
                ! not yet full effect
                if (CCi > CCx) then
                    Multiplier = 0._sp  ! no effect
                else
                    Multiplier = (CCx-CCi)/(CCx/2._sp)
                end if
            else
                Multiplier = 1._sp ! full effect
            end if
            EpotVal = EpotVal * (1._sp - CCx * (CCEffectProcent/100._sp) * Multiplier)
            EpotMin = GetSimulParam_KcWetBare() &
                      * (1._sp - 1.72_sp*CCx + 1._sp*(CCx*CCx) &
                            - 0.30_sp*(CCx*CCx*CCx)) * EToVal
            if (EpotMin < epsilon(1._sp)) then
                EpotMin = 0._sp
            end if
            if (EpotVal < EpotMin) then
                EpotVal = EpotMin
            end if
            if (EpotVal > EpotMax) then
                EpotVal = EpotMax
            end if
        end if

        ! Correction for canopy senescence before late-season stage
        if (GetSimulation_EvapLimitON()) then
            if (EpotVal > EpotMax) then
                EpotVal = EpotMax
            end if
        end if

        ! Correction for drop in photosynthetic capacity of a dying green canopy
        if (CCi < CCxWithered) then
            if ((CCxWithered > 0.01_sp) .and. (CCi > 0.001_sp)) then
                TpotVal = TpotVal &
                           * exp(GetSimulParam_ExpFsen() &
                                * log(CCi/CCxWithered))
            end if
        end if
    end if
end subroutine CalculateETpot



subroutine LoadProgramParametersProject(FullFileNameProgramParameters)
    character(len=*), intent(in) :: FullFileNameProgramParameters

    integer :: fhandle
    integer(int32) :: i, simul_RpZmi, simul_lowox
    integer(int8) :: simul_ed, effrainperc, effrainshow, effrainrootE, &
                     simul_saltdiff, simul_saltsolub, simul_root, simul_pCCHIf, &
                     simul_SFR, simul_TAWg, simul_beta, simul_Tswc, simul_GDD, &
                     simul_EZma
    real(sp) :: simul_rod, simul_kcWB, simul_RZEma, simul_pfao, simul_expFsen, &
                simul_Tmi, simul_Tma
    logical :: file_exists

    inquire(file=trim(FullFileNameProgramParameters), exist=file_exists)
    if (file_exists) then
        ! load set of program parameters
        open(newunit=fhandle, file=trim(FullFileNameProgramParameters), &
             status='old', action='read')
        ! crop
        read(fhandle, *) simul_ed ! evaporation decline factor in stage 2
        call SetSimulParam_EvapDeclineFactor(simul_ed)
        read(fhandle, *) simul_kcWB ! Kc wet bare soil [-]
        call SetSimulParam_KcWetBare(simul_kcWB)
        read(fhandle, *) simul_pCCHIf ! CC threshold below which HI no longer
                                      ! increase(% of 100)
        call SetSimulParam_PercCCxHIfinal(simul_pCCHIf)
        read(fhandle, *) simul_RpZmi ! Starting depth of root sine function
                                     ! (% of Zmin)
        call SetSimulParam_RootPercentZmin(simul_RpZmi)
        read(fhandle, *) simul_RZEma ! cm/day
        call SetSimulParam_MaxRootZoneExpansion(simul_RZEma)
        call SetSimulParam_MaxRootZoneExpansion(5.00_sp) ! fixed at 5 cm/day
        read(fhandle, *) simul_SFR ! Shape factor for effect water stress
                                   ! on rootzone expansion
        call SetSimulParam_KsShapeFactorRoot(simul_SFR)
        read(fhandle, *) simul_TAWg  ! Soil water content (% TAW) required
                                     ! at sowing depth for germination
        call SetSimulParam_TAWGermination(simul_TAWg)
        read(fhandle, *) simul_pfao ! Adjustment factor for FAO-adjustment
                                    ! soil water depletion (p) for various ET
        call SetSimulParam_pAdjFAO(simul_pfao)
        read(fhandle, *) simul_lowox ! number of days for full effect of
                                     ! deficient aeration
        call SetSimulParam_DelayLowOxygen(simul_lowox)
        read(fhandle, *) simul_expFsen ! exponent of senescence factor
                                       ! adjusting drop in photosynthetic
                                       ! activity of dying crop
        call SetSimulParam_ExpFsen(simul_expFsen)
        read(fhandle, *) simul_beta ! Decrease (percentage) of p(senescence)
                                    ! once early canopy senescence is triggered
        call SetSimulParam_Beta(simul_beta)
        read(fhandle, *) simul_Tswc  ! Thickness top soil (cm) in which soil
                                     ! water depletion has to be determined
        call SetSimulParam_ThicknessTopSWC(simul_Tswc)
        ! field
        read(fhandle, *) simul_EZma ! maximum water extraction depth by soil
                                    ! evaporation [cm]
        call SetSimulParam_EvapZmax(simul_EZma)
        ! soil
        read(fhandle, *) simul_rod ! considered depth (m) of soil profile for
                                   ! calculation of mean soil water content
        call SetSimulParam_RunoffDepth(simul_rod)
        read(fhandle, *) i   ! correction CN for Antecedent Moisture Class
        if (i == 1) then
            call SetSimulParam_CNcorrection(.true.)
        else
            call SetSimulParam_CNcorrection(.false.)
        end if
        read(fhandle, *) simul_saltdiff ! salt diffusion factor (%)
        read(fhandle, *) simul_saltsolub ! salt solubility (g/liter)
        read(fhandle, *) simul_root ! shape factor capillary rise factor
        call SetSimulParam_SaltDiff(simul_saltdiff)
        call SetSimulParam_SaltSolub(simul_saltsolub)
        call SetSimulParam_RootNrDF(simul_root)
        call SetSimulParam_IniAbstract(5_int8) ! fixed in Version 5.0 cannot be &
                                     ! changed since linked with equations for
                                     ! CN AMCII and CN converions
        ! Temperature
        read(fhandle, *) simul_Tmi   ! Default minimum temperature (degC) if no &
                                     ! temperature file is specified
        call SetSimulParam_Tmin(simul_Tmi)
        read(fhandle, *) simul_Tma   ! Default maximum temperature (degC) if no &
                                     ! temperature file is specified
        call SetSimulParam_Tmax(simul_Tma)
        read(fhandle, *) simul_GDD ! Default method for GDD calculations
        call SetSimulParam_GDDMethod(simul_GDD)
        if (GetSimulParam_GDDMethod() > 3_int8) then
            call SetSimulParam_GDDMethod(3_int8)
        end if
        if (GetSimulParam_GDDMethod()< 1_int8) then
            call SetSimulParam_GDDMethod(3_int8)
        end if
        ! Rainfall
        read(fhandle, *) i
        select case (i)
            case (0)
                call SetSimulParam_EffectiveRain_Method(EffectiveRainMethod_Full)
            case (1)
                call SetSimulParam_EffectiveRain_Method(EffectiveRainMethod_USDA)
            case (2)
                call SetSimulParam_EffectiveRain_Method(EffectiveRainMethod_Percentage)
        end select
        read(fhandle, *) effrainperc ! IF Method is Percentage
        call SetSimulParam_EffectiveRain_PercentEffRain(effrainperc)
        read(fhandle, *) effrainshow  ! For estimation of surface run-off
        call SetSimulParam_EffectiveRain_ShowersInDecade(effrainshow)
        read(fhandle, *) effrainrootE ! For reduction of soil evaporation
        call SetSimulParam_EffectiveRain_RootNrEvap(effrainrootE)
        ! close
        Close(fhandle)
    else
        ! take the default set of program parameters
        call ReadSoilSettings
        call ReadRainfallSettings
        call ReadCropSettingsParameters
        call ReadFieldSettingsParameters
        call ReadTemperatureSettingsParameters
    end if
end subroutine LoadProgramParametersProject


subroutine ReadCropSettingsParameters()

    integer :: fhandle
    character(len=:), allocatable :: FullName
    integer(int8) :: simul_ed, simul_pCCHIf, simul_SFR, simul_TAWg, &
                     simul_beta, simul_Tswc
    real(sp) :: simul_kcWB, simul_RZEma, simul_pfao, simul_expFsen
    integer(int32) :: simul_RpZmi, simul_lowox

    FullName = GetPathNameSimul() // 'Crop.PAR'
    open(newunit=fhandle, file=trim(FullName), status='old', action='read')
    read(fhandle, *) simul_ed ! evaporation decline factor in stage 2
    call SetSimulParam_EvapDeclineFactor(simul_ed)
    read(fhandle, *) simul_kcWB ! Kc wet bare soil [-]
    call SetSimulParam_KcWetBare(simul_kcWB)
    read(fhandle, *) simul_pCCHIf ! CC threshold below which HI no
                                  ! longer increase(% of 100)
    call SetSimulParam_PercCCxHIfinal(simul_pCCHIf)
    read(fhandle, *) simul_RpZmi ! Starting depth of root sine function
                                 ! (% of Zmin)
    call SetSimulParam_RootPercentZmin(simul_RpZmi)
    read(fhandle, *) simul_RZEma ! cm/day
    call SetSimulParam_MaxRootZoneExpansion(simul_RZEma)
    call SetSimulParam_MaxRootZoneExpansion(5.00_sp) ! fixed at 5 cm/day
    read(fhandle, *) simul_SFR ! Shape factor for effect water stress on
                               ! rootzone expansion
    call SetSimulParam_KsShapeFactorRoot(simul_SFR)
    read(fhandle, *) simul_TAWg  ! Soil water content (% TAW) required
                                 ! at sowing depth for germination
    call SetSimulParam_TAWGermination(simul_TAWg)
    read(fhandle, *) simul_pfao ! Adjustment factor for FAO-adjustment soil
                                ! water depletion (p) for various ET
    call SetSimulParam_pAdjFAO(simul_pfao)
    read(fhandle, *) simul_lowox ! number of days for full effect of
                                 ! deficient aeration
    call SetSimulParam_DelayLowOxygen(simul_lowox)
    read(fhandle, *) simul_expFsen ! exponent of senescence factor adjusting
                               ! drop in photosynthetic activity of dying crop
    call SetSimulParam_ExpFsen(simul_expFsen)
    read(fhandle, *) simul_beta ! Decrease (percentage) of p(senescence) once
                                ! early canopy senescence is triggered
    call SetSimulParam_Beta(simul_beta)
    read(fhandle, *) simul_Tswc ! Thickness top soil (cm) in which soil water
                                ! depletion has to be determined
    call SetSimulParam_ThicknessTopSWC(simul_Tswc)
    close(fhandle)
end subroutine ReadCropSettingsParameters


subroutine ReadFieldSettingsParameters()

    integer :: fhandle
    character(len=:), allocatable :: FullName
    integer(int8) :: simul_evmax

    FullName = GetPathNameSimul() // 'Field.PAR'
    open(newunit=fhandle, file=trim(FullName), status='old', action='read')
    read(fhandle, *) simul_evmax ! maximum water extraction depth by
                                 ! soil evaporation [cm]
    call SetSimulParam_EvapZmax(simul_evmax)
    close(fhandle)
end subroutine ReadFieldSettingsParameters


subroutine ReadTemperatureSettingsParameters()

    integer :: fhandle
    character(len=:), allocatable :: FullName
    integer(int8) :: simul_GDD
    real(sp) :: simul_Tmi, simul_Tma

    FullName = GetPathNameSimul() // 'Temperature.PAR'
    open(newunit=fhandle, file=trim(FullName), status='old', action='read')
    read(fhandle, *)
    read(fhandle, *) simul_Tmi ! Default minimum temperature (degC) if no
                               ! temperature file is specified
    call SetSimulParam_Tmin(simul_Tmi)
    read(fhandle, *) simul_Tma ! Default maximum temperature (degC) if no
                               ! temperature file is specified
    call SetSimulParam_Tmax(simul_Tma)
    read(fhandle, *) simul_GDD ! Default method for GDD calculations
    call SetSimulParam_GDDMethod(simul_GDD)
    if (GetSimulParam_GDDMethod() > 3_int8) then
        call SetSimulParam_GDDMethod(3_int8)
    end if
    if (GetSimulParam_GDDMethod() < 1_int8) then
        call SetSimulParam_GDDMethod(1_int8)
    end if
    close(fhandle)
end subroutine ReadTemperatureSettingsParameters



subroutine CompleteClimateDescription(ClimateRecord)
    type(rep_clim), intent(inout) :: ClimateRecord

    character(len=2) :: dayStr
    character(len=4) ::  yearStr
    integer(int32) :: Deci

    call DetermineDayNr(ClimateRecord%FromD, ClimateRecord%FromM, &
                        ClimateRecord%FromY, ClimateRecord%FromDayNr)
    select case (ClimateRecord%DataType)
        case(datatype_daily)
            ClimateRecord%ToDayNr = ClimateRecord%FromDayNr &
                                    + ClimateRecord%NrObs - 1
            call DetermineDate(ClimateRecord%ToDayNr, ClimateRecord%ToD, &
                               ClimateRecord%ToM, ClimateRecord%ToY)
        case(datatype_decadely)
            Deci = roundc((ClimateRecord%FromD+9)/10._sp, mold=1) &
                            + ClimateRecord%NrObs - 1
            ClimateRecord%ToM = ClimateRecord%FromM
            ClimateRecord%ToY = ClimateRecord%FromY
            do while (Deci > 3)
                Deci = Deci - 3
                ClimateRecord%ToM = ClimateRecord%ToM + 1
                if (ClimateRecord%ToM > 12) then
                    ClimateRecord%ToM = 1
                    ClimateRecord%ToY = ClimateRecord%ToY  + 1
                end if
            end do
            ClimateRecord%ToD = 10
            if (Deci == 2) then
                ClimateRecord%ToD = 20
            end if
            if (Deci == 3) then
                ClimateRecord%ToD = DaysInMonth(ClimateRecord%ToM)
                if ((ClimateRecord%ToM == 2) &
                            .and. LeapYear(ClimateRecord%ToY)) then
                    ClimateRecord%ToD = ClimateRecord%ToD + 1
                end if
            end if
            call DetermineDayNr(ClimateRecord%ToD, ClimateRecord%ToM, &
                                ClimateRecord%ToY, ClimateRecord%ToDayNr)
        case(datatype_monthly)
            ClimateRecord%ToY = ClimateRecord%FromY
            ClimateRecord%ToM = ClimateRecord%FromM + ClimateRecord%NrObs - 1
            do while (ClimateRecord%ToM > 12)
                ClimateRecord%ToY = ClimateRecord%ToY + 1
                ClimateRecord%ToM = ClimateRecord%ToM - 12
            end do
            ClimateRecord%ToD = DaysInMonth(ClimateRecord%ToM)
            if ((ClimateRecord%ToM == 2) &
                        .and. LeapYear(ClimateRecord%ToY)) then
                ClimateRecord%ToD = ClimateRecord%ToD + 1
            end if
            call DetermineDayNr(ClimateRecord%ToD, ClimateRecord%ToM, &
                                ClimateRecord%ToY, ClimateRecord%ToDayNr)
    end select
    write(dayStr, '(i2)') ClimateRecord%FromD
    if (ClimateRecord%FromY == 1901) then
        yearStr = ''
    else
        write(yearStr, '(i4)') ClimateRecord%FromY
    end if
    ClimateRecord%FromString = dayStr // ' ' // &
                               NameMonth(ClimateRecord%FromM) // &
                               ' ' // yearStr
    write(dayStr, '(i2)') ClimateRecord%ToD
    if (ClimateRecord%FromY == 1901) then
        yearStr = ''
    else
        write(yearStr, '(i4)') ClimateRecord%ToY
    end if
    ClimateRecord%ToString = dayStr // ' ' // NameMonth(ClimateRecord%ToM) // &
                             ' ' // yearStr
end subroutine CompleteClimateDescription

integer(int32) function SumCalendarDaysReferenceTnx(ValGDDays, RefCropDay1,&
                                        StartDayNr, Tbase, Tupper,&
                                        TDayMin, TDayMax)
    integer(int32), intent(in) :: ValGDDays
    integer(int32), intent(in) :: RefCropDay1
    integer(int32), intent(in) :: StartDayNr
    real(sp), intent(in) :: Tbase
    real(sp), intent(in) :: Tupper
    real(sp), intent(in) :: TDayMin
    real(sp), intent(in) :: TDayMax

    integer(int32) :: i
    integer(int32) :: NrCDays
    real(sp) :: RemainingGDDays, DayGDD
    real(sp) :: TDayMin_loc, TDayMax_loc

    TDayMin_loc = TDayMin
    TDayMax_loc = TDayMax

    NrCdays = 0
    if (ValGDDays > 0) then
        if (GetTnxReferenceFile() == '(None)') then
            ! given average Tmin and Tmax
            DayGDD = DegreesDay(Tbase, Tupper, &
                       TDayMin_loc, TDayMax_loc, GetSimulParam_GDDMethod())
            if (abs(DayGDD) < epsilon(1._sp)) then
                NrCDays = 0
            else
                NrCDays = roundc(ValGDDays/DayGDD, mold=1_int32)
            end if
        else
            ! Get TCropReference: mean daily Tnx (365 days) from RefCropDay1 onwards
            ! determine corresponding calendar days
            RemainingGDDays = ValGDDays

            ! TminCropReference and TmaxCropReference arrays contain the TemperatureFilefull data
            i = StartDayNr - RefCropDay1

            ! For crops with very long cycles (unrealistic but avoid crashing)
            if (i >= 365) then
                i = i - 365
            end if

            do while (RemainingGDDays > 0.1_sp)
                i = i + 1
                if (i == size(GetTminCropReferenceRun())) then
                    i = 1
                end if
                TDayMin_loc = GetTminCropReferenceRun_i(i)
                TDayMax_loc = GetTmaxCropReferenceRun_i(i)

                DayGDD = DegreesDay(Tbase, Tupper, TDayMin_loc, &
                                    TDayMax_loc, &
                                    GetSimulParam_GDDMethod())
                if (DayGDD > RemainingGDDays) then
                    if (roundc((DayGDD-RemainingGDDays)/RemainingGDDays,mold=1) >= 1) then
                        NrCDays = NrCDays + 1
                    end if
                else
                    NrCDays = NrCDays + 1
                end if
                RemainingGDDays = RemainingGDDays - DayGDD
            end do
        end if
    end if
    SumCalendarDaysReferenceTnx = NrCDays
end function SumCalendarDaysReferenceTnx


!! Global variables section !!


function GetOutputName() result(str)
    !! Getter for the "OutputName" global variable.
    character(len=:), allocatable :: str

    str = OutputName
end function GetOutputName


subroutine SetOutputName(str)
    !! Setter for the "OutputName" global variable.
    character(len=*), intent(in) :: str

    OutputName = str
end subroutine SetOutputName


function GetCO2File() result(str)
    !! Getter for the "CO2File" global variable.
    character(len=:), allocatable :: str

    str = CO2File
end function GetCO2File


subroutine SetCO2File(str)
    !! Setter for the "CO2File" global variable.
    character(len=*), intent(in) :: str

    CO2File = str
end subroutine SetCO2File


function GetCO2FileFull() result(str)
    !! Getter for the "CO2FileFull" global variable.
    character(len=:), allocatable :: str

    str = CO2FileFull
end function GetCO2FileFull


subroutine SetCO2FileFull(str)
    !! Setter for the "CO2FileFull" global variable.
    character(len=*), intent(in) :: str

    CO2FileFull = str
end subroutine SetCO2FileFull


function GetCO2Description() result(str)
    !! Getter for the "CO2Description" global variable.
    character(len=:), allocatable :: str

    str = CO2Description
end function GetCO2Description


subroutine SetCO2Description(str)
    !! Setter for the "CO2Description" global variable.
    character(len=*), intent(in) :: str

    CO2Description = str
end subroutine SetCO2Description


type(rep_RootZoneWC) function GetRootZoneWC()
    !! Getter for the "RootZoneWC" global variable.

    GetRootZoneWC = RootZoneWC
end function GetRootZoneWC


real(sp) function GetRootZoneWC_Actual()
    !! Getter for the "Rootzonewc" global variable.

    GetRootZoneWC_Actual = RootZoneWC%Actual
end function GetRootZoneWC_Actual


real(sp) function GetRootZoneWC_FC()
    !! Getter for the "Rootzonewc" global variable.

    GetRootZoneWC_FC = RootZoneWC%FC
end function GetRootZoneWC_FC


real(sp) function GetRootZoneWC_WP()
    !! Getter for the "Rootzonewc" global variable.

    GetRootZoneWC_WP = RootZoneWC%WP
end function GetRootZoneWC_WP


real(sp) function GetRootZoneWC_SAT()
    !! Getter for the "Rootzonewc" global variable.

    GetRootZoneWC_SAT = RootZoneWC%SAT
end function GetRootZoneWC_SAT


real(sp) function GetRootZoneWC_Leaf()
    !! Getter for the "Rootzonewc" global variable.

    GetRootZoneWC_Leaf = RootZoneWC%Leaf
end function GetRootZoneWC_Leaf


real(sp) function GetRootZoneWC_Thresh()
    !! Getter for the "Rootzonewc" global variable.

    GetRootZoneWC_Thresh = RootZoneWC%Thresh
end function GetRootZoneWC_Thresh


real(sp) function GetRootZoneWC_Sen()
    !! Getter for the "Rootzonewc" global variable.

    GetRootZoneWC_Sen = RootZoneWC%Sen
end function GetRootZoneWC_Sen


real(sp) function GetRootZoneWC_ZtopAct()
    !! Getter for the "Rootzonewc" global variable.

    GetRootZoneWC_ZtopAct = RootZoneWC%ZtopAct
end function GetRootZoneWC_ZtopAct


real(sp) function GetRootZoneWC_ZtopFC()
    !! Getter for the "Rootzonewc" global variable.

    GetRootZoneWC_ZtopFC = RootZoneWC%ZtopFC
end function GetRootZoneWC_ZtopFC


real(sp) function GetRootZoneWC_ZtopWP()
    !! Getter for the "Rootzonewc" global variable.

    GetRootZoneWC_ZtopWP = RootZoneWC%ZtopWP
end function GetRootZoneWC_ZtopWP


real(sp) function GetRootZoneWC_ZtopThresh()
    !! Getter for the "Rootzonewc" global variable.

    GetRootZoneWC_ZtopThresh = RootZoneWC%ZtopThresh
end function GetRootZoneWC_ZtopThresh


subroutine SetRootZoneWC_Actual(Actual)
    !! Setter for the "RootZoneWC" global variable.
    real(sp), intent(in) :: Actual

    RootZoneWC%Actual = Actual
end subroutine SetRootZoneWC_Actual


subroutine SetRootZoneWC_FC(FC)
    !! Setter for the "RootZoneWC" global variable.
    real(sp), intent(in) :: FC

    RootZoneWC%FC = FC
end subroutine SetRootZoneWC_FC


subroutine SetRootZoneWC_WP(WP)
    !! Setter for the "RootZoneWC" global variable.
    real(sp), intent(in) :: WP

    RootZoneWC%WP = WP
end subroutine SetRootZoneWC_WP


subroutine SetRootZoneWC_SAT(SAT)
    !! Setter for the "RootZoneWC" global variable.
    real(sp), intent(in) :: SAT

    RootZoneWC%SAT = SAT
end subroutine SetRootZoneWC_SAT


subroutine SetRootZoneWC_Leaf(Leaf)
    !! Setter for the "RootZoneWC" global variable.
    real(sp), intent(in) :: Leaf

    RootZoneWC%Leaf = Leaf
end subroutine SetRootZoneWC_Leaf


subroutine SetRootZoneWC_Thresh(Thresh)
    !! Setter for the "RootZoneWC" global variable.
    real(sp), intent(in) :: Thresh

    RootZoneWC%Thresh = Thresh
end subroutine SetRootZoneWC_Thresh


subroutine SetRootZoneWC_Sen(Sen)
    !! Setter for the "RootZoneWC" global variable.
    real(sp), intent(in) :: Sen

    RootZoneWC%Sen = Sen
end subroutine SetRootZoneWC_Sen


subroutine SetRootZoneWC_ZtopAct(ZtopAct)
    !! Setter for the "RootZoneWC" global variable.
    real(sp), intent(in) :: ZtopAct

    RootZoneWC%ZtopAct = ZtopAct
end subroutine SetRootZoneWC_ZtopAct


subroutine SetRootZoneWC_ZtopFC(ZtopFC)
    !! Setter for the "RootZoneWC" global variable.
    real(sp), intent(in) :: ZtopFC

    RootZoneWC%ZtopFC = ZtopFC
end subroutine SetRootZoneWC_ZtopFC


subroutine SetRootZoneWC_ZtopWP(ZtopWP)
    !! Setter for the "RootZoneWC" global variable.
    real(sp), intent(in) :: ZtopWP

    RootZoneWC%ZtopWP = ZtopWP
end subroutine SetRootZoneWC_ZtopWP


subroutine SetRootZoneWC_ZtopThresh(ZtopThresh)
    !! Setter for the "RootZoneWC" global variable.
    real(sp), intent(in) :: ZtopThresh

    RootZoneWC%ZtopThresh = ZtopThresh
end subroutine SetRootZoneWC_ZtopThresh


function GetCalendarFile() result(str)
    !! Getter for the "CalendarFile" global variable.
    character(len=:), allocatable :: str

    str = CalendarFile
end function GetCalendarFile


subroutine SetCalendarFile(str)
    !! Setter for the "CalendarFile" global variable.
    character(len=*), intent(in) :: str

    CalendarFile = str
end subroutine SetCalendarFile


function GetCalendarFileFull() result(str)
    !! Getter for the "CalendarFileFull" global variable.
    character(len=:), allocatable :: str

    str = CalendarFileFull
end function GetCalendarFileFull


subroutine SetCalendarFileFull(str)
    !! Setter for the "CalendarFileFull" global variable.
    character(len=*), intent(in) :: str

    CalendarFileFull = str
end subroutine SetCalendarFileFull


function GetCalendarDescription() result(str)
    !! Getter for the "CalendarDescription" global variable.
    character(len=:), allocatable :: str

    str = CalendarDescription
end function GetCalendarDescription


subroutine SetCalendarDescription(str)
    !! Setter for the "CalendarDescription" global variable.
    character(len=*), intent(in) :: str

    CalendarDescription = str
end subroutine SetCalendarDescription


function GetCropFile() result(str)
    !! Getter for the "CropFile" global variable.
    character(len=:), allocatable :: str

    str = CropFile
end function GetCropFile


subroutine SetCropFile(str)
    !! Setter for the "CropFile" global variable.
    character(len=*), intent(in) :: str

    CropFile = str
end subroutine SetCropFile


function GetCropFileFull() result(str)
    !! Getter for the "CropFile" global variable.
    character(len=:), allocatable :: str

    str = CropFileFull
end function GetCropFileFull


subroutine SetCropFileFull(str)
    !! Setter for the "CropFile" global variable.
    character(len=*), intent(in) :: str

    CropFileFull = str
end subroutine SetCropFileFull


function GetCropDescription() result(str)
    !! Getter for the "CropDescription" global variable.
    character(len=:), allocatable :: str

    str = CropDescription
end function GetCropDescription


subroutine SetCropDescription(str)
    !! Setter for the "CropDescription" global variable.
    character(len=*), intent(in) :: str

    CropDescription = str
end subroutine SetCropDescription


type(rep_IrriECw) function GetIrriECw()
    !! Getter for the "IrriECw" global variable.

    GetIrriECw = IrriECw
end function GetIrriECw


real(sp) function GetIrriECw_PreSeason()
    !! Getter for the "IrriECw" global variable.

    GetIrriECw_PreSeason = IrriECw%PreSeason
end function GetIrriECw_PreSeason


subroutine SetIrriECw(IrriECw_in)
    !! Setter for the "IrriECw" global variable.
    type(rep_IrriECw), intent(in) :: IrriECw_in

    IrriECw = IrriECw_in
end subroutine SetIrriECw


subroutine SetIrriECw_PreSeason(PreSeason)
    !! Setter for the "IrriECw" global variable.
    real(sp), intent(in) :: PreSeason

    IrriECw%PreSeason = PreSeason
end subroutine SetIrriECw_PreSeason


real(sp) function GetIrriECw_PostSeason()
    !! Getter for the "IrriECw" global variable.

    GetIrriECw_PostSeason = IrriECw%PostSeason
end function GetIrriECw_PostSeason


subroutine SetIrriECw_PostSeason(PostSeason)
    !! Setter for the "IrriECw" global variable.
    real(sp), intent(in) :: PostSeason

    IrriECw%PostSeason = PostSeason
end subroutine SetIrriECw_PostSeason


function GetProfFile() result(str)
    !! Getter for the "ProfFile" global variable.
    character(len=:), allocatable :: str

    str = ProfFile
end function GetProfFile


subroutine SetProfFile(str)
    !! Setter for the "ProfFile" global variable.
    character(len=*), intent(in) :: str

    ProfFile = str
end subroutine SetProfFile


function GetProfFilefull() result(str)
    !! Getter for the "ProfFilefull" global variable.
    character(len=:), allocatable :: str

    str = ProfFilefull
end function GetProfFilefull


subroutine SetProfFilefull(str)
    !! Setter for the "ProfFilefull" global variable.
    character(len=*), intent(in) :: str

    ProfFilefull = str
end subroutine SetProfFilefull


function GetProfDescription() result(str)
    !! Getter for the "ProfDescription" global variable.
    character(len=:), allocatable :: str

    str = ProfDescription
end function GetProfDescription


subroutine SetProfDescription(str)
    !! Setter for the "ProfDescription" global variable.
    character(len=*), intent(in) :: str

    ProfDescription = str
end subroutine SetProfDescription


function GetManFile() result(str)
    !! Getter for the "ManFile" global variable.
    character(len=:), allocatable :: str

    str = ManFile
end function GetManFile


subroutine SetManFile(str)
    !! Setter for the "ManFile" global variable.
    character(len=*), intent(in) :: str

    ManFile = str
end subroutine SetManFile


function GetManFilefull() result(str)
    !! Getter for the "ManFilefull" global variable.
    character(len=:), allocatable :: str

    str = ManFilefull
end function GetManFilefull


subroutine SetManFilefull(str)
    !! Setter for the "ManFilefull" global variable.
    character(len=*), intent(in) :: str

    ManFilefull = str
end subroutine SetManFilefull


function GetOffSeasonFile() result(str)
    !! Getter for the "OffSeasonFile" global variable.
    character(len=:), allocatable :: str

    str = OffSeasonFile
end function GetOffSeasonFile


subroutine SetOffSeasonFile(str)
    !! Setter for the "OffSeasonFile" global variable.
    character(len=*), intent(in) :: str

    OffSeasonFile = str
end subroutine SetOffSeasonFile


function GetOffSeasonFilefull() result(str)
    !! Getter for the "OffSeasonFilefull" global variable.
    character(len=:), allocatable :: str

    str = OffSeasonFilefull
end function GetOffSeasonFilefull


subroutine SetOffSeasonFilefull(str)
    !! Setter for the "OffSeasonFilefull" global variable.
    character(len=*), intent(in) :: str

    OffSeasonFilefull = str
end subroutine SetOffSeasonFilefull


function GetObservationsFile() result(str)
    !! Getter for the "ObservationsFile" global variable.
    character(len=:), allocatable :: str

    str = ObservationsFile
end function GetObservationsFile


subroutine SetObservationsFile(str)
    !! Setter for the "ObservationsFile" global variable.
    character(len=*), intent(in) :: str

    ObservationsFile = str
end subroutine SetObservationsFile


function GetObservationsFilefull() result(str)
    !! Getter for the "ObservationsFilefull" global variable.
    character(len=:), allocatable :: str

    str = ObservationsFilefull
end function GetObservationsFilefull


subroutine SetObservationsFilefull(str)
    !! Setter for the "ObservationsFilefull" global variable.
    character(len=*), intent(in) :: str

    ObservationsFilefull = str
end subroutine SetObservationsFilefull


function GetObservationsDescription() result(str)
    !! Getter for the "ObservationsDescription" global variable.
    character(len=:), allocatable :: str

    str = ObservationsDescription
end function GetObservationsDescription


subroutine SetObservationsDescription(str)
    !! Setter for the "ObservationsDescription" global variable.
    character(len=*), intent(in) :: str

    ObservationsDescription = str
end subroutine SetObservationsDescription


function GetGroundWaterFile() result(str)
    !! Getter for the "GroundWaterFile" global variable.
    character(len=:), allocatable :: str

    str = GroundWaterFile
end function GetGroundWaterFile


subroutine SetGroundWaterFile(str)
    !! Setter for the "GroundWaterFile" global variable.
    character(len=*), intent(in) :: str

    GroundWaterFile = str
end subroutine SetGroundWaterFile


function GetGroundWaterFilefull() result(str)
    !! Getter for the "GroundWaterFilefull" global variable.
    character(len=:), allocatable :: str

    str = GroundWaterFilefull
end function GetGroundWaterFilefull


subroutine SetGroundWaterFilefull(str)
    !! Setter for the "GroundWaterFilefull" global variable.
    character(len=*), intent(in) :: str

    GroundWaterFilefull = str
end subroutine SetGroundWaterFilefull


type(rep_CropFileSet) function GetCropFileSet()
    !! Getter for the "CropFileSet" global variable.

    GetCropFileSet = CropFileSet
end function GetCropFileSet


subroutine SetCropFileSet_DaysFromSenescenceToEnd(DaysFromSenescenceToEnd)
    !! Setter for the "CropFileSet" global variable.
    integer(int32), intent(in) :: DaysFromSenescenceToEnd

    CropFileSet%DaysFromSenescenceToEnd = DaysFromSenescenceToEnd
end subroutine SetCropFileSet_DaysFromSenescenceToEnd


subroutine SetCropFileSet_DaysToHarvest(DaysToHarvest)
    !! Setter for the "CropFileSet" global variable.
    integer(int32), intent(in) :: DaysToHarvest

    CropFileSet%DaysToHarvest = DaysToHarvest
end subroutine SetCropFileSet_DaysToHarvest


subroutine SetCropFileSet_GDDaysFromSenescenceToEnd(GDDaysFromSenescenceToEnd)
    !! Setter for the "CropFileSet" global variable.
    integer(int32), intent(in) :: GDDaysFromSenescenceToEnd

    CropFileSet%GDDaysFromSenescenceToEnd = GDDaysFromSenescenceToEnd
end subroutine SetCropFileSet_GDDaysFromSenescenceToEnd


subroutine SetCropFileSet_GDDaysToHarvest(GDDaysToHarvest)
    !! Setter for the "CropFileSet" global variable.
    integer(int32), intent(in) :: GDDaysToHarvest

    CropFileSet%GDDaysToHarvest = GDDaysToHarvest
end subroutine SetCropFileSet_GDDaysToHarvest


function GetEToFile() result(str)
    !! Getter for the "EToFile" global variable.
    character(len=:), allocatable :: str

    str = EToFile
end function GetEToFile


subroutine SetEToFile(str)
    !! Setter for the "EToFile" global variable.
    character(len=*), intent(in) :: str

    EToFile = str
end subroutine SetEToFile


function GetEToFileFull() result(str)
    !! Getter for the "EToFileFull" global variable.
    character(len=:), allocatable :: str

    str = EToFileFull
end function GetEToFileFull


subroutine SetEToFileFull(str)
    !! Setter for the "EToFileFull" global variable.
    character(len=*), intent(in) :: str

    EToFileFull = str
end subroutine SetEToFileFull


function GetEToDescription() result(str)
    !! Getter for the "EToDescription" global variable.
    character(len=:), allocatable :: str

    str = EToDescription
end function GetEToDescription


subroutine SetEToDescription(str)
    !! Setter for the "EToDescription" global variable.
    character(len=*), intent(in) :: str

    EToDescription = str
end subroutine SetEToDescription


function GetRainFile() result(str)
    !! Getter for the "RainFile" global variable.
    character(len=:), allocatable :: str

    str = RainFile
end function GetRainFile


subroutine SetRainFile(str)
    !! Setter for the "RainFile" global variable.
    character(len=*), intent(in) :: str

    RainFile = str
end subroutine SetRainFile


function GetRainFileFull() result(str)
    !! Getter for the "RainFileFull" global variable.
    character(len=:), allocatable :: str

    str = RainFileFull
end function GetRainFileFull


subroutine SetRainFileFull(str)
    !! Setter for the "RainFileFull" global variable.
    character(len=*), intent(in) :: str

    RainFileFull = str
end subroutine SetRainFileFull


function GetRainDescription() result(str)
    !! Getter for the "RainDescription" global variable.
    character(len=:), allocatable :: str

    str = RainDescription
end function GetRainDescription


subroutine SetRainDescription(str)
    !! Setter for the "RainDescription" global variable.
    character(len=*), intent(in) :: str

    RainDescription = str
end subroutine SetRainDescription


type(rep_Manag) function GetManagement()
    !! Getter for the "Management" global variable.

    GetManagement = Management
end function GetManagement


integer(int8) function GetManagement_Mulch()
    !! Getter for the "Management" global variable.

    GetManagement_Mulch = Management%Mulch
end function GetManagement_Mulch


integer(int8) function GetManagement_SoilCoverBefore()
    !! Getter for the "Management" global variable.

    GetManagement_SoilCoverBefore = Management%SoilCoverBefore
end function GetManagement_SoilCoverBefore


integer(int8) function GetManagement_SoilCoverAfter()
    !! Getter for the "Management" global variable.

    GetManagement_SoilCoverAfter = Management%SoilCoverAfter
end function GetManagement_SoilCoverAfter


integer(int8) function GetManagement_EffectMulchOffS()
    !! Getter for the "Management" global variable.

    GetManagement_EffectMulchOffS = Management%EffectMulchOffS
end function GetManagement_EffectMulchOffS


integer(int8) function GetManagement_EffectMulchInS()
    !! Getter for the "Management" global variable.

    GetManagement_EffectMulchInS = Management%EffectMulchInS
end function GetManagement_EffectMulchInS


integer(int32) function GetManagement_FertilityStress()
    !! Getter for the "Management" global variable.

    GetManagement_FertilityStress = Management%FertilityStress
end function GetManagement_FertilityStress


real(sp) function GetManagement_BundHeight()
    !! Getter for the "Management" global variable.

    GetManagement_BundHeight = Management%BundHeight
end function GetManagement_BundHeight


logical function GetManagement_RunoffOn()
    !! Getter for the "Management" global variable.

    GetManagement_RunoffOn = Management%RunoffOn
end function GetManagement_RunoffOn


integer(int32) function GetManagement_CNcorrection()
    !! Getter for the "Management" global variable.

    GetManagement_CNcorrection = Management%CNcorrection
end function GetManagement_CNcorrection


integer(int8) function GetManagement_WeedRC()
    !! Getter for the "Management" global variable.

    GetManagement_WeedRC = Management%WeedRC
end function GetManagement_WeedRC


integer(int32) function GetManagement_WeedDeltaRC()
    !! Getter for the "Management" global variable.

    GetManagement_WeedDeltaRC = Management%WeedDeltaRC
end function GetManagement_WeedDeltaRC


real(sp) function GetManagement_WeedShape()
    !! Getter for the "Management" global variable.

    GetManagement_WeedShape = Management%WeedShape
end function GetManagement_WeedShape


integer(int8) function GetManagement_WeedAdj()
    !! Getter for the "Management" global variable.

    GetManagement_WeedAdj = Management%WeedAdj
end function GetManagement_WeedAdj


type(rep_Cuttings) function GetManagement_Cuttings()
    !! Setter for the "Management" global variable.

    GetManagement_Cuttings = Management%Cuttings
end function GetManagement_Cuttings


subroutine SetManagement_Mulch(Mulch)
    !! Setter for the "Management" global variable.
    integer(int8), intent(in) :: Mulch

    Management%Mulch = Mulch
end subroutine SetManagement_Mulch


subroutine SetManagement_SoilCoverBefore(SoilCoverBefore)
    !! Setter for the "Management" global variable.
    integer(int8), intent(in) :: SoilCoverBefore

    Management%SoilCoverBefore = SoilCoverBefore
end subroutine SetManagement_SoilCoverBefore


subroutine SetManagement_SoilCoverAfter(SoilCoverAfter)
    !! Setter for the "Management" global variable.
    integer(int8), intent(in) :: SoilCoverAfter

    Management%SoilCoverAfter = SoilCoverAfter
end subroutine SetManagement_SoilCoverAfter


subroutine SetManagement_EffectMulchOffS(EffectMulchOffS)
    !! Setter for the "Management" global variable.
    integer(int8), intent(in) :: EffectMulchOffS

    Management%EffectMulchOffS = EffectMulchOffS
end subroutine SetManagement_EffectMulchOffS


subroutine SetManagement_EffectMulchInS(EffectMulchInS)
    !! Setter for the "Management" global variable.
    integer(int8), intent(in) :: EffectMulchInS

    Management%EffectMulchInS = EffectMulchInS
end subroutine SetManagement_EffectMulchInS


subroutine SetManagement_FertilityStress(FertilityStress)
    !! Setter for the "Management" global variable.
    integer(int32), intent(in) :: FertilityStress

    Management%FertilityStress = FertilityStress
end subroutine SetManagement_FertilityStress


subroutine SetManagement_BundHeight(BundHeight)
    !! Setter for the "Management" global variable.
    real(sp), intent(in) :: BundHeight

    Management%BundHeight = BundHeight
end subroutine SetManagement_BundHeight


subroutine SetManagement_RunoffOn(RunoffOn)
    !! Setter for the "Management" global variable.
    logical, intent(in) :: RunoffOn

    Management%RunoffOn = RunoffOn
end subroutine SetManagement_RunoffOn


subroutine SetManagement_CNcorrection(CNcorrection)
    !! Setter for the "Management" global variable.
    integer(int32), intent(in) :: CNcorrection

    Management%CNcorrection = CNcorrection
end subroutine SetManagement_CNcorrection


subroutine SetManagement_WeedRC(WeedRC)
    !! Setter for the "Management" global variable.
    integer(int8), intent(in) :: WeedRC

    Management%WeedRC = WeedRC
end subroutine SetManagement_WeedRC


subroutine SetManagement_WeedDeltaRC(WeedDeltaRC)
    !! Setter for the "Management" global variable.
    integer(int32), intent(in) :: WeedDeltaRC

    Management%WeedDeltaRC = WeedDeltaRC
end subroutine SetManagement_WeedDeltaRC


subroutine SetManagement_WeedShape(WeedShape)
    !! Setter for the "Management" global variable.
    real(sp), intent(in) :: WeedShape

    Management%WeedShape = WeedShape
end subroutine SetManagement_WeedShape


subroutine SetManagement_WeedAdj(WeedAdj)
    !! Setter for the "Management" global variable.
    integer(int8), intent(in) :: WeedAdj

    Management%WeedAdj = WeedAdj
end subroutine SetManagement_WeedAdj


subroutine SetManagement_Cuttings(Cuttings)
    !! Setter for the "Management" global variable.
    type(rep_Cuttings), intent(in) :: Cuttings

    Management%Cuttings = Cuttings
end subroutine SetManagement_Cuttings


logical function GetManagement_Cuttings_Considered()
    !! Getter for the "Cuttings" global variable.

    GetManagement_Cuttings_Considered = Cuttings%Considered
end function GetManagement_Cuttings_Considered


integer(int32) function GetManagement_Cuttings_CCcut()
    !! Getter for the "Cuttings" global variable.

    GetManagement_Cuttings_CCcut = Cuttings%CCcut
end function GetManagement_Cuttings_CCcut


integer(int32) function GetManagement_Cuttings_Day1()
    !! Getter for the "Cuttings" global variable.

    GetManagement_Cuttings_Day1 = Cuttings%Day1
end function GetManagement_Cuttings_Day1


integer(int32) function GetManagement_Cuttings_NrDays()
    !! Setter for the "Cuttings" global variable.

    GetManagement_Cuttings_NrDays = Cuttings%NrDays
end function GetManagement_Cuttings_NrDays


logical function GetManagement_Cuttings_Generate()
    !! Getter for the "Cuttings" global variable.

    GetManagement_Cuttings_Generate = Cuttings%Generate
end function GetManagement_Cuttings_Generate


integer(intEnum) function GetManagement_Cuttings_Criterion()
    !! Getter for the "Cuttings" global variable.

    GetManagement_Cuttings_Criterion = Cuttings%Criterion
end function GetManagement_Cuttings_Criterion


logical function GetManagement_Cuttings_HarvestEnd()
    !! Getter for the "Cuttings" global variable.

    GetManagement_Cuttings_HarvestEnd = Cuttings%HarvestEnd
end function GetManagement_Cuttings_HarvestEnd


integer(int32) function GetManagement_Cuttings_FirstDayNr()
    !! Getter for the "Cuttings" global variable.

    GetManagement_Cuttings_FirstDayNr = Cuttings%FirstDayNr
end function GetManagement_Cuttings_FirstDayNr


subroutine SetManagement(Management_in)
    !! Setter for the "Management" global variable.
    type(rep_Manag), intent(in) :: Management_in

    Management = Management_in
end subroutine SetManagement


subroutine SetManagement_Cuttings_Considered(Considered)
    !! Setter for the "Cuttings" global variable.
    logical, intent(in) :: Considered

    Cuttings%Considered = Considered
end subroutine SetManagement_Cuttings_Considered


subroutine SetManagement_Cuttings_CCcut(CCcut)
    !! Setter for the "Cuttings" global variable.
    integer(int32), intent(in) :: CCcut

    Cuttings%CCcut = CCcut
end subroutine SetManagement_Cuttings_CCcut


subroutine SetManagement_Cuttings_Day1(Day1)
    !! Setter for the "Cuttings" global variable.
    integer(int32), intent(in) :: Day1

    Cuttings%Day1 = Day1
end subroutine SetManagement_Cuttings_Day1


subroutine SetManagement_Cuttings_NrDays(NrDays)
    !! Setter for the "Cuttings" global variable.
    integer(int32), intent(in) :: NrDays

    Cuttings%NrDays = NrDays
end subroutine SetManagement_Cuttings_NrDays


subroutine SetManagement_Cuttings_Generate(Generate)
    !! Setter for the "Cuttings" global variable.
    logical, intent(in) :: Generate

    Cuttings%Generate = Generate
end subroutine SetManagement_Cuttings_Generate


subroutine SetManagement_Cuttings_Criterion(Criterion)
    !! Setter for the "Cuttings" global variable.
    integer(intEnum), intent(in) :: Criterion

    Cuttings%Criterion = Criterion
end subroutine SetManagement_Cuttings_Criterion


subroutine SetManagement_Cuttings_HarvestEnd(HarvestEnd)
    !! Setter for the "Cuttings" global variable.
    logical, intent(in) :: HarvestEnd

    Cuttings%HarvestEnd = HarvestEnd
end subroutine SetManagement_Cuttings_HarvestEnd


subroutine SetManagement_Cuttings_FirstDayNr(FirstDayNr)
    !! Setter for the "Cuttings" global variable.
    integer(int32), intent(in) :: FirstDayNr

    Cuttings%FirstDayNr = FirstDayNr
end subroutine SetManagement_Cuttings_FirstDayNr



function Geteffectiverain() result(effectiverain_out)
    !! Getter for the "effectiverain" global variable.
    type(rep_EffectiveRain) :: effectiverain_out

    effectiverain_out = effectiverain
end function Geteffectiverain


function Geteffectiverain_Method() result(Method)
    !! Getter for the "Method" attribute of the "effectiverain" global variable.
    integer(intEnum) :: Method

    Method = effectiverain%Method
end function Geteffectiverain_Method


function Geteffectiverain_PercentEffRain() result(PercentEffRain)
    !! Getter for the "PercentEffRain" attribute of the "effectiverain" global variable.
    integer(int8) :: PercentEffRain

    PercentEffRain = effectiverain%PercentEffRain
end function Geteffectiverain_PercentEffRain


function Geteffectiverain_ShowersInDecade() result(ShowersInDecade)
    !! Getter for the "ShowersInDecade" attribute of the "effectiverain" global variable.
    integer(int8) :: ShowersInDecade

    ShowersInDecade = effectiverain%ShowersInDecade
end function Geteffectiverain_ShowersInDecade


function Geteffectiverain_RootNrEvap() result(RootNrEvap)
    !! Getter for the "RootNrEvap" attribute of the "effectiverain" global variable.
    integer(int8) :: RootNrEvap

    RootNrEvap = effectiverain%RootNrEvap
end function Geteffectiverain_RootNrEvap


subroutine Seteffectiverain(effectiverain_in)
    !! Setter for the "effectiverain" global variable.
    type(rep_EffectiveRain), intent(in) :: effectiverain_in

    effectiverain = effectiverain_in
end subroutine Seteffectiverain


subroutine Seteffectiverain_Method(Method)
    !! Setter for the "Method" attribute of the "effectiverain" global variable.
    integer(intEnum), intent(in) :: Method

    effectiverain%Method = Method
end subroutine Seteffectiverain_Method


subroutine Seteffectiverain_PercentEffRain(PercentEffRain)
    !! Setter for the "PercentEffRain" attribute of the "effectiverain" global variable.
    integer(int8), intent(in) :: PercentEffRain

    effectiverain%PercentEffRain = PercentEffRain
end subroutine Seteffectiverain_PercentEffRain


subroutine Seteffectiverain_ShowersInDecade(ShowersInDecade)
    !! Setter for the "ShowersInDecade" attribute of the "effectiverain" global variable.
    integer(int8), intent(in) :: ShowersInDecade

    effectiverain%ShowersInDecade = ShowersInDecade
end subroutine Seteffectiverain_ShowersInDecade


subroutine Seteffectiverain_RootNrEvap(RootNrEvap)
    !! Setter for the "RootNrEvap" attribute of the "effectiverain" global variable.
    integer(int8), intent(in) :: RootNrEvap

    effectiverain%RootNrEvap = RootNrEvap
end subroutine Seteffectiverain_RootNrEvap


function GetSimulParam() result(SimulParam_out)
    !! Getter for the "simulparam" global variable.
    type(rep_param) :: SimulParam_out

    SimulParam_out = simulparam
end function GetSimulParam


function GetSimulParam_EvapDeclineFactor() result(EvapDeclineFactor)
    !! Getter for the "EvapDeclineFactor" attribute of the "simulparam" global variable.
    integer(int8) :: EvapDeclineFactor

    EvapDeclineFactor = simulparam%EvapDeclineFactor
end function GetSimulParam_EvapDeclineFactor


type(rep_EffectiveRain) function GetSimulParam_EffectiveRain()
    !! Setter for the "EffectiveRain" global variable.

    GetSimulParam_EffectiveRain = SimulParam%EffectiveRain
end function GetSimulParam_EffectiveRain


function GetSimulParam_KcWetBare() result(KcWetBare)
    !! Getter for the "KcWetBare" attribute of the "simulparam" global variable.
    real(sp) :: KcWetBare

    KcWetBare = simulparam%KcWetBare
end function GetSimulParam_KcWetBare


function GetSimulParam_PercCCxHIfinal() result(PercCCxHIfinal)
    !! Getter for the "PercCCxHIfinal" attribute of the "simulparam" global variable.
    integer(int8) :: PercCCxHIfinal

    PercCCxHIfinal = simulparam%PercCCxHIfinal
end function GetSimulParam_PercCCxHIfinal


function GetSimulParam_RootPercentZmin() result(RootPercentZmin)
    !! Getter for the "RootPercentZmin" attribute of the "simulparam" global variable.
    integer(int32) :: RootPercentZmin

    RootPercentZmin = simulparam%RootPercentZmin
end function GetSimulParam_RootPercentZmin


function GetSimulParam_MaxRootZoneExpansion() result(MaxRootZoneExpansion)
    !! Getter for the "MaxRootZoneExpansion" attribute of the "simulparam" global variable.
    real(sp) :: MaxRootZoneExpansion

    MaxRootZoneExpansion = simulparam%MaxRootZoneExpansion
end function GetSimulParam_MaxRootZoneExpansion


function GetSimulParam_KsShapeFactorRoot() result(KsShapeFactorRoot)
    !! Getter for the "KsShapeFactorRoot" attribute of the "simulparam" global variable.
    integer(int8) :: KsShapeFactorRoot

    KsShapeFactorRoot = simulparam%KsShapeFactorRoot
end function GetSimulParam_KsShapeFactorRoot


function GetSimulParam_TAWGermination() result(TAWGermination)
    !! Getter for the "TAWGermination" attribute of the "simulparam" global variable.
    integer(int8) :: TAWGermination

    TAWGermination = simulparam%TAWGermination
end function GetSimulParam_TAWGermination


function GetSimulParam_pAdjFAO() result(pAdjFAO)
    !! Getter for the "pAdjFAO" attribute of the "simulparam" global variable.
    real(sp) :: pAdjFAO

    pAdjFAO = simulparam%pAdjFAO
end function GetSimulParam_pAdjFAO


function GetSimulParam_DelayLowOxygen() result(DelayLowOxygen)
    !! Getter for the "DelayLowOxygen" attribute of the "simulparam" global variable.
    integer(int32) :: DelayLowOxygen

    DelayLowOxygen = simulparam%DelayLowOxygen
end function GetSimulParam_DelayLowOxygen


function GetSimulParam_ExpFsen() result(ExpFsen)
    !! Getter for the "ExpFsen" attribute of the "simulparam" global variable.
    real(sp) :: ExpFsen

    ExpFsen = simulparam%ExpFsen
end function GetSimulParam_ExpFsen


function GetSimulParam_Beta() result(Beta)
    !! Getter for the "Beta" attribute of the "simulparam" global variable.
    integer(int8) :: Beta

    Beta = simulparam%Beta
end function GetSimulParam_Beta


function GetSimulParam_ThicknessTopSWC() result(ThicknessTopSWC)
    !! Getter for the "ThicknessTopSWC" attribute of the "simulparam" global variable.
    integer(int8) :: ThicknessTopSWC

    ThicknessTopSWC = simulparam%ThicknessTopSWC
end function GetSimulParam_ThicknessTopSWC


function GetSimulParam_EvapZmax() result(EvapZmax)
    !! Getter for the "EvapZmax" attribute of the "simulparam" global variable.
    integer(int8) :: EvapZmax

    EvapZmax = simulparam%EvapZmax
end function GetSimulParam_EvapZmax


function GetSimulParam_RunoffDepth() result(RunoffDepth)
    !! Getter for the "RunoffDepth" attribute of the "simulparam" global variable.
    real(sp) :: RunoffDepth

    RunoffDepth = simulparam%RunoffDepth
end function GetSimulParam_RunoffDepth


function GetSimulParam_CNcorrection() result(CNcorrection)
    !! Getter for the "CNcorrection" attribute of the "simulparam" global variable.
    logical :: CNcorrection

    CNcorrection = simulparam%CNcorrection
end function GetSimulParam_CNcorrection


function GetSimulParam_Tmin() result(Tmin)
    !! Getter for the "Tmin" attribute of the "simulparam" global variable.
    real(sp) :: Tmin

    Tmin = simulparam%Tmin
end function GetSimulParam_Tmin


function GetSimulParam_Tmax() result(Tmax)
    !! Getter for the "Tmax" attribute of the "simulparam" global variable.
    real(sp) :: Tmax

    Tmax = simulparam%Tmax
end function GetSimulParam_Tmax


function GetSimulParam_GDDMethod() result(GDDMethod)
    !! Getter for the "GDDMethod" attribute of the "simulparam" global variable.
    integer(int8) :: GDDMethod

    GDDMethod = simulparam%GDDMethod
end function GetSimulParam_GDDMethod


function GetSimulParam_PercRAW() result(PercRAW)
    !! Getter for the "PercRAW" attribute of the "simulparam" global variable.
    integer(int32) :: PercRAW

    PercRAW = simulparam%PercRAW
end function GetSimulParam_PercRAW


function GetSimulParam_CompDefThick() result(CompDefThick)
    !! Getter for the "CompDefThick" attribute of the "simulparam" global variable.
    real(sp) :: CompDefThick

    CompDefThick = simulparam%CompDefThick
end function GetSimulParam_CompDefThick


function GetSimulParam_CropDay1() result(CropDay1)
    !! Getter for the "CropDay1" attribute of the "simulparam" global variable.
    integer(int32) :: CropDay1

    CropDay1 = simulparam%CropDay1
end function GetSimulParam_CropDay1


function GetSimulParam_Tbase() result(Tbase)
    !! Getter for the "Tbase" attribute of the "simulparam" global variable.
    real(sp) :: Tbase

    Tbase = simulparam%Tbase
end function GetSimulParam_Tbase


function GetSimulParam_Tupper() result(Tupper)
    !! Getter for the "Tupper" attribute of the "simulparam" global variable.
    real(sp) :: Tupper

    Tupper = simulparam%Tupper
end function GetSimulParam_Tupper


function GetSimulParam_IrriFwInSeason() result(IrriFwInSeason)
    !! Getter for the "IrriFwInSeason" attribute of the "simulparam" global variable.
    integer(int8) :: IrriFwInSeason

    IrriFwInSeason = simulparam%IrriFwInSeason
end function GetSimulParam_IrriFwInSeason


function GetSimulParam_IrriFwOffSeason() result(IrriFwOffSeason)
    !! Getter for the "IrriFwOffSeason" attribute of the "simulparam" global variable.
    integer(int8) :: IrriFwOffSeason

    IrriFwOffSeason = simulparam%IrriFwOffSeason
end function GetSimulParam_IrriFwOffSeason


function GetSimulParam_SaltDiff() result(SaltDiff)
    !! Getter for the "SaltDiff" attribute of the "simulparam" global variable.
    integer(int8) :: SaltDiff

    SaltDiff = simulparam%SaltDiff
end function GetSimulParam_SaltDiff


function GetSimulParam_SaltSolub() result(SaltSolub)
    !! Getter for the "SaltSolub" attribute of the "simulparam" global variable.
    integer(int8) :: SaltSolub

    SaltSolub = simulparam%SaltSolub
end function GetSimulParam_SaltSolub


function GetSimulParam_ConstGwt() result(ConstGwt)
    !! Getter for the "ConstGwt" attribute of the "simulparam" global variable.
    logical :: ConstGwt

    ConstGwt = simulparam%ConstGwt
end function GetSimulParam_ConstGwt


function GetSimulParam_RootNrDF() result(RootNrDF)
    !! Getter for the "RootNrDF" attribute of the "simulparam" global variable.
    integer(int8) :: RootNrDF

    RootNrDF = simulparam%RootNrDF
end function GetSimulParam_RootNrDF


function GetSimulParam_IniAbstract() result(IniAbstract)
    !! Getter for the "IniAbstract" attribute of the "simulparam" global variable.
    integer(int8) :: IniAbstract

    IniAbstract = simulparam%IniAbstract
end function GetSimulParam_IniAbstract


subroutine SetSimulParam(SimulParam_in)
    !! Setter for the "simulparam" global variable.
    type(rep_param), intent(in) :: SimulParam_in

    simulparam = SimulParam_in
end subroutine SetSimulParam


subroutine SetSimulParam_EvapDeclineFactor(EvapDeclineFactor)
    !! Setter for the "EvapDeclineFactor" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: EvapDeclineFactor

    simulparam%EvapDeclineFactor = EvapDeclineFactor
end subroutine SetSimulParam_EvapDeclineFactor


subroutine SetSimulParam_EffectiveRain(EffectiveRain)
    !! Setter for the "EffectiveRain" global variable.
    type(rep_EffectiveRain), intent(in) :: EffectiveRain

    SimulParam%EffectiveRain = EffectiveRain
end subroutine SetSimulParam_EffectiveRain


subroutine SetSimulParam_KcWetBare(KcWetBare)
    !! Setter for the "KcWetBare" attribute of the "simulparam" global variable.
    real(sp), intent(in) :: KcWetBare

    simulparam%KcWetBare = KcWetBare
end subroutine SetSimulParam_KcWetBare


subroutine SetSimulParam_PercCCxHIfinal(PercCCxHIfinal)
    !! Setter for the "PercCCxHIfinal" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: PercCCxHIfinal

    simulparam%PercCCxHIfinal = PercCCxHIfinal
end subroutine SetSimulParam_PercCCxHIfinal


subroutine SetSimulParam_RootPercentZmin(RootPercentZmin)
    !! Setter for the "RootPercentZmin" attribute of the "simulparam" global variable.
    integer(int32), intent(in) :: RootPercentZmin

    simulparam%RootPercentZmin = RootPercentZmin
end subroutine SetSimulParam_RootPercentZmin


subroutine SetSimulParam_MaxRootZoneExpansion(MaxRootZoneExpansion)
    !! Setter for the "MaxRootZoneExpansion" attribute of the "simulparam" global variable.
    real(sp), intent(in) :: MaxRootZoneExpansion

    simulparam%MaxRootZoneExpansion = MaxRootZoneExpansion
end subroutine SetSimulParam_MaxRootZoneExpansion


subroutine SetSimulParam_KsShapeFactorRoot(KsShapeFactorRoot)
    !! Setter for the "KsShapeFactorRoot" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: KsShapeFactorRoot

    simulparam%KsShapeFactorRoot = KsShapeFactorRoot
end subroutine SetSimulParam_KsShapeFactorRoot


subroutine SetSimulParam_TAWGermination(TAWGermination)
    !! Setter for the "TAWGermination" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: TAWGermination

    simulparam%TAWGermination = TAWGermination
end subroutine SetSimulParam_TAWGermination


subroutine SetSimulParam_pAdjFAO(pAdjFAO)
    !! Setter for the "pAdjFAO" attribute of the "simulparam" global variable.
    real(sp), intent(in) :: pAdjFAO

    simulparam%pAdjFAO = pAdjFAO
end subroutine SetSimulParam_pAdjFAO


subroutine SetSimulParam_DelayLowOxygen(DelayLowOxygen)
    !! Setter for the "DelayLowOxygen" attribute of the "simulparam" global variable.
    integer(int32), intent(in) :: DelayLowOxygen

    simulparam%DelayLowOxygen = DelayLowOxygen
end subroutine SetSimulParam_DelayLowOxygen


subroutine SetSimulParam_ExpFsen(ExpFsen)
    !! Setter for the "ExpFsen" attribute of the "simulparam" global variable.
    real(sp), intent(in) :: ExpFsen

    simulparam%ExpFsen = ExpFsen
end subroutine SetSimulParam_ExpFsen


subroutine SetSimulParam_Beta(Beta)
    !! Setter for the "Beta" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: Beta

    simulparam%Beta = Beta
end subroutine SetSimulParam_Beta


subroutine SetSimulParam_ThicknessTopSWC(ThicknessTopSWC)
    !! Setter for the "ThicknessTopSWC" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: ThicknessTopSWC

    simulparam%ThicknessTopSWC = ThicknessTopSWC
end subroutine SetSimulParam_ThicknessTopSWC


subroutine SetSimulParam_EvapZmax(EvapZmax)
    !! Setter for the "EvapZmax" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: EvapZmax

    simulparam%EvapZmax = EvapZmax
end subroutine SetSimulParam_EvapZmax


subroutine SetSimulParam_RunoffDepth(RunoffDepth)
    !! Setter for the "RunoffDepth" attribute of the "simulparam" global variable.
    real(sp), intent(in) :: RunoffDepth

    simulparam%RunoffDepth = RunoffDepth
end subroutine SetSimulParam_RunoffDepth


subroutine SetSimulParam_CNcorrection(CNcorrection)
    !! Setter for the "CNcorrection" attribute of the "simulparam" global variable.
    logical, intent(in) :: CNcorrection

    simulparam%CNcorrection = CNcorrection
end subroutine SetSimulParam_CNcorrection


subroutine SetSimulParam_Tmin(Tmin)
    !! Setter for the "Tmin" attribute of the "simulparam" global variable.
    real(sp), intent(in) :: Tmin

    simulparam%Tmin = Tmin
end subroutine SetSimulParam_Tmin


subroutine SetSimulParam_Tmax(Tmax)
    !! Setter for the "Tmax" attribute of the "simulparam" global variable.
    real(sp), intent(in) :: Tmax

    simulparam%Tmax = Tmax
end subroutine SetSimulParam_Tmax


subroutine SetSimulParam_GDDMethod(GDDMethod)
    !! Setter for the "GDDMethod" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: GDDMethod

    simulparam%GDDMethod = GDDMethod
end subroutine SetSimulParam_GDDMethod


subroutine SetSimulParam_PercRAW(PercRAW)
    !! Setter for the "PercRAW" attribute of the "simulparam" global variable.
    integer(int32), intent(in) :: PercRAW

    simulparam%PercRAW = PercRAW
end subroutine SetSimulParam_PercRAW


subroutine SetSimulParam_CompDefThick(CompDefThick)
    !! Setter for the "CompDefThick" attribute of the "simulparam" global variable.
    real(sp), intent(in) :: CompDefThick

    simulparam%CompDefThick = CompDefThick
end subroutine SetSimulParam_CompDefThick


subroutine SetSimulParam_CropDay1(CropDay1)
    !! Setter for the "CropDay1" attribute of the "simulparam" global variable.
    integer(int32), intent(in) :: CropDay1

    simulparam%CropDay1 = CropDay1
end subroutine SetSimulParam_CropDay1


subroutine SetSimulParam_Tbase(Tbase)
    !! Setter for the "Tbase" attribute of the "simulparam" global variable.
    real(sp), intent(in) :: Tbase

    simulparam%Tbase = Tbase
end subroutine SetSimulParam_Tbase


subroutine SetSimulParam_Tupper(Tupper)
    !! Setter for the "Tupper" attribute of the "simulparam" global variable.
    real(sp), intent(in) :: Tupper

    simulparam%Tupper = Tupper
end subroutine SetSimulParam_Tupper


subroutine SetSimulParam_IrriFwInSeason(IrriFwInSeason)
    !! Setter for the "IrriFwInSeason" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: IrriFwInSeason

    simulparam%IrriFwInSeason = IrriFwInSeason
end subroutine SetSimulParam_IrriFwInSeason


subroutine SetSimulParam_IrriFwOffSeason(IrriFwOffSeason)
    !! Setter for the "IrriFwOffSeason" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: IrriFwOffSeason

    simulparam%IrriFwOffSeason = IrriFwOffSeason
end subroutine SetSimulParam_IrriFwOffSeason


subroutine SetSimulParam_SaltDiff(SaltDiff)
    !! Setter for the "SaltDiff" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: SaltDiff

    simulparam%SaltDiff = SaltDiff
end subroutine SetSimulParam_SaltDiff


subroutine SetSimulParam_SaltSolub(SaltSolub)
    !! Setter for the "SaltSolub" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: SaltSolub

    simulparam%SaltSolub = SaltSolub
end subroutine SetSimulParam_SaltSolub


subroutine SetSimulParam_ConstGwt(ConstGwt)
    !! Setter for the "ConstGwt" attribute of the "simulparam" global variable.
    logical, intent(in) :: ConstGwt

    simulparam%ConstGwt = ConstGwt
end subroutine SetSimulParam_ConstGwt


subroutine SetSimulParam_RootNrDF(RootNrDF)
    !! Setter for the "RootNrDF" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: RootNrDF

    simulparam%RootNrDF = RootNrDF
end subroutine SetSimulParam_RootNrDF


subroutine SetSimulParam_IniAbstract(IniAbstract)
    !! Setter for the "IniAbstract" attribute of the "simulparam" global variable.
    integer(int8), intent(in) :: IniAbstract

    simulparam%IniAbstract = IniAbstract
end subroutine SetSimulParam_IniAbstract


integer(intEnum) function GetSimulParam_EffectiveRain_Method()
    !! Getter for the "EffectiveRain" global variable.

    GetSimulParam_EffectiveRain_Method = EffectiveRain%Method
end function GetSimulParam_EffectiveRain_Method


integer(int8) function GetSimulParam_EffectiveRain_PercentEffRain()
    !! Getter for the "EffectiveRain" global variable.

    GetSimulParam_EffectiveRain_PercentEffRain = EffectiveRain%PercentEffRain
end function GetSimulParam_EffectiveRain_PercentEffRain


integer(int8) function GetSimulParam_EffectiveRain_ShowersInDecade()
    !! Getter for the "EffectiveRain" global variable.

    GetSimulParam_EffectiveRain_ShowersInDecade = EffectiveRain%ShowersInDecade
end function GetSimulParam_EffectiveRain_ShowersInDecade


integer(int8) function GetSimulParam_EffectiveRain_RootNrEvap()
    !! Getter for the "EffectiveRain" global variable.

    GetSimulParam_EffectiveRain_RootNrEvap = EffectiveRain%RootNrEvap
end function GetSimulParam_EffectiveRain_RootNrEvap


subroutine SetSimulParam_EffectiveRain_Method(Method)
    !! Setter for the "EffectiveRain" global variable.
    integer(intEnum), intent(in) :: Method

    EffectiveRain%Method = Method
end subroutine SetSimulParam_EffectiveRain_Method


subroutine SetSimulParam_EffectiveRain_PercentEffRain(PercentEffRain)
    !! Setter for the "EffectiveRain" global variable.
    integer(int8), intent(in) :: PercentEffRain

    EffectiveRain%PercentEffRain = PercentEffRain
end subroutine SetSimulParam_EffectiveRain_PercentEffRain


subroutine SetSimulParam_EffectiveRain_ShowersInDecade(ShowersInDecade)
    !! Setter for the "EffectiveRain" global variable.
    integer(int8), intent(in) :: ShowersInDecade

    EffectiveRain%ShowersInDecade= ShowersInDecade
end subroutine SetSimulParam_EffectiveRain_ShowersInDecade


subroutine SetSimulParam_EffectiveRain_RootNrEvap(RootNrEvap)
    !! Setter for the "EffectiveRain" global variable.
    integer(int8), intent(in) :: RootNrEvap

    EffectiveRain%RootNrEvap= RootNrEvap
end subroutine SetSimulParam_EffectiveRain_RootNrEvap


type(rep_sum) function GetSumWaBal()
    !! Getter for the "SymWaBal" global variable.

    GetSumWaBal = SumWaBal
end function GetSumWaBal


real(sp) function GetSumWaBal_Epot()
    !! Getter for the "SumWaBal" global variable.

     GetSumWaBal_Epot = SumWaBal%Epot
end function GetSumWaBal_Epot


real(sp) function GetSumWaBal_Tpot()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_Tpot = SumWaBal%Tpot
end function GetSumWaBal_Tpot


real(sp) function GetSumWaBal_Rain()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_Rain = SumWaBal%Rain
end function GetSumWaBal_Rain


real(sp) function GetSumWaBal_Irrigation()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_Irrigation = SumWaBal%Irrigation
end function GetSumWaBal_Irrigation


real(sp) function GetSumWaBal_Infiltrated()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_Infiltrated = SumWaBal%Infiltrated
end function GetSumWaBal_Infiltrated


real(sp) function GetSumWaBal_Runoff()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_Runoff = SumWaBal%Runoff
end function GetSumWaBal_Runoff


real(sp) function GetSumWaBal_Drain()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_Drain = SumWaBal%Drain
end function GetSumWaBal_Drain


real(sp) function GetSumWaBal_Eact()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_Eact = SumWaBal%Eact
end function GetSumWaBal_Eact


real(sp) function GetSumWaBal_Tact()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_Tact = SumWaBal%Tact
end function GetSumWaBal_Tact


real(sp) function GetSumWaBal_TrW()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_TrW = SumWaBal%TrW
end function GetSumWaBal_TrW


real(sp) function GetSumWaBal_ECropCycle()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_ECropCycle = SumWaBal%ECropCycle
end function GetSumWaBal_ECropCycle


real(sp) function GetSumWaBal_CRwater()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_CRwater = SumWaBal%CRwater
end function GetSumWaBal_CRwater


real(sp) function GetSumWaBal_Biomass()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_Biomass = SumWaBal%Biomass
end function GetSumWaBal_Biomass


real(sp) function GetSumWaBal_YieldPart()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_YieldPart = SumWaBal%YieldPart
end function GetSumWaBal_YieldPart


real(sp) function GetSumWaBal_BiomassPot()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_BiomassPot = SumWaBal%BiomassPot
end function GetSumWaBal_BiomassPot


real(sp) function GetSumWaBal_BiomassUnlim()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_BiomassUnlim = SumWaBal%BiomassUnlim
end function GetSumWaBal_BiomassUnlim


real(sp) function GetSumWaBal_BiomassTot()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_BiomassTot = SumWaBal%BiomassTot
end function GetSumWaBal_BiomassTot


real(sp) function GetSumWaBal_SaltIn()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_SaltIn = SumWaBal%SaltIn
end function GetSumWaBal_SaltIn


real(sp) function GetSumWaBal_SaltOut()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_SaltOut = SumWaBal%SaltOut
end function GetSumWaBal_SaltOut


real(sp) function GetSumWaBal_CRSalt()
    !! Getter for the "SumWaBal" global variable.

    GetSumWaBal_CRSalt = SumWaBal%CRSalt
end function GetSumWaBal_CRSalt


subroutine SetSumWaBal(SumWaBal_in)
    !! Setter for the "SumWaBal" global variable.
    type(rep_sum), intent(in) :: SumWaBal_in

    SumWaBal = SumWaBal_in
end subroutine SetSumWaBal


subroutine SetSumWaBal_Epot(Epot)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: Epot

    SumWaBal%Epot = Epot
end subroutine SetSumWaBal_Epot


subroutine SetSumWaBal_Tpot(Tpot)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: Tpot

    SumWaBal%Tpot = Tpot
end subroutine SetSumWaBal_Tpot


subroutine SetSumWaBal_Rain(Rain)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: Rain

    SumWaBal%Rain = Rain
end subroutine SetSumWaBal_Rain


subroutine SetSumWaBal_Irrigation(Irrigation)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: Irrigation

    SumWaBal%Irrigation = Irrigation
end subroutine SetSumWaBal_Irrigation


subroutine SetSumWaBal_Infiltrated(Infiltrated)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: Infiltrated

    SumWaBal%Infiltrated = Infiltrated
end subroutine SetSumWaBal_Infiltrated


subroutine SetSumWaBal_Runoff(Runoff)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: Runoff

    SumWaBal%Runoff = Runoff
end subroutine SetSumWaBal_Runoff


subroutine SetSumWaBal_Drain(Drain)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: Drain

    SumWaBal%Drain = Drain
end subroutine SetSumWaBal_Drain


subroutine SetSumWaBal_Eact(Eact)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: Eact

    SumWaBal%Eact = Eact
end subroutine SetSumWaBal_Eact


subroutine SetSumWaBal_Tact(Tact)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: Tact

    SumWaBal%Tact = Tact
end subroutine SetSumWaBal_Tact


subroutine SetSumWaBal_TrW(TrW)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: TrW

    SumWaBal%TrW = TrW
end subroutine SetSumWaBal_TrW


subroutine SetSumWaBal_ECropCycle(ECropCycle)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: ECropCycle

    SumWaBal%ECropCycle = ECropCycle
end subroutine SetSumWaBal_ECropCycle


subroutine SetSumWaBal_CRwater(CRwater)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: CRwater

    SumWaBal%CRwater = CRwater
end subroutine SetSumWaBal_CRwater


subroutine SetSumWaBal_Biomass(Biomass)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: Biomass

    SumWaBal%Biomass = Biomass
end subroutine SetSumWaBal_Biomass


subroutine SetSumWaBal_YieldPart(YieldPart)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: YieldPart

    SumWaBal%YieldPart = YieldPart
end subroutine SetSumWaBal_YieldPart


subroutine SetSumWaBal_BiomassPot(BiomassPot)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: BiomassPot

    SumWaBal%BiomassPot = BiomassPot
end subroutine SetSumWaBal_BiomassPot


subroutine SetSumWaBal_BiomassUnlim(BiomassUnlim)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: BiomassUnlim

    SumWaBal%BiomassUnlim = BiomassUnlim
end subroutine SetSumWaBal_BiomassUnlim


subroutine SetSumWaBal_BiomassTot(BiomassTot)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: BiomassTot

    SumWaBal%BiomassTot = BiomassTot
end subroutine SetSumWaBal_BiomassTot


subroutine SetSumWaBal_SaltIn(SaltIn)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: SaltIn

    SumWaBal%SaltIn = SaltIn
end subroutine SetSumWaBal_SaltIn


subroutine SetSumWaBal_SaltOut(SaltOut)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: SaltOut

    SumWaBal%SaltOut = SaltOut
end subroutine SetSumWaBal_SaltOut


subroutine SetSumWaBal_CRSalt(CRSalt)
    !! Setter for the "SumWaBal" global variable.
    real(sp), intent(in) :: CRSalt

    SumWaBal%CRSalt = CRSalt
end subroutine SetSumWaBal_CRSalt


type(rep_soil) function GetSoil()
    !! Getter for the "Soil" global variable.

    GetSoil = Soil
end function GetSoil


integer(int8) function GetSoil_REW()
    !! Getter for "REW" attribute of the "soil" global variable.

    GetSoil_REW = soil%REW
end function GetSoil_REW


real(sp) function GetSoil_RootMax()
    !! Getter for "RootMax" attribute of the "soil" global variable.

    GetSoil_Rootmax = soil%RootMax
end function GetSoil_RootMax


integer(int8) function GetSoil_NrSoilLayers()
    !! Getter for "NrSoilLayers" attribute of the "soil" global variable.

    GetSoil_NrSoilLayers = soil%NrSoilLayers
end function GetSoil_NrSoilLayers


integer(int8) function GetSoil_CNvalue()
    !! Getter for "CNvalue" attribute of the "soil" global variable.

    GetSoil_CNvalue = soil%CNvalue
end function GetSoil_CNvalue


subroutine SetSoil(Soil_in)
    !! Setter for the "Soil" global variable.
    type(rep_soil), intent(in) :: Soil_in

    Soil = Soil_in
end subroutine SetSoil


subroutine SetSoil_REW(REW)
    !! Setter for the "Soil" global variable.
    integer(int8), intent(in) :: REW

    Soil%REW = REW
end subroutine SetSoil_REW


subroutine SetSoil_NrSoilLayers(NrSoilLayers)
    !! Setter for the "Soil" global variable.
    integer(int8), intent(in) :: NrSoilLayers

    Soil%NrSoilLayers = NrSoilLayers
end subroutine SetSoil_NrSoilLayers


subroutine SetSoil_CNvalue(CNvalue)
    !! Setter for the "Soil" global variable.
    integer(int8), intent(in) :: CNvalue

    Soil%CNvalue = CNvalue
end subroutine SetSoil_CNvalue


subroutine SetSoil_RootMax(RootMax)
    !! Setter for the "Soil" global variable.
    real(sp), intent(in) :: RootMax

    Soil%RootMax = RootMax
end subroutine SetSoil_RootMax


function GetCrop_StressResponse() result(StressResponse)
    !! Getter for the "StressResponse" attribute of the "crop" global variable.
    type(rep_Shapes) :: StressResponse

    StressResponse = crop%StressResponse
end function GetCrop_StressResponse


function GetCrop_StressResponse_Stress() result(Stress)
    !! Getter for the "Stress" attribute of the "StressResponse" attribute of the "crop" global variable.
    integer(int8) :: Stress

    Stress = crop%StressResponse%Stress
end function GetCrop_StressResponse_Stress


function GetCrop_StressResponse_ShapeCGC() result(ShapeCGC)
    !! Getter for the "ShapeCGC" attribute of the "StressResponse" attribute of the "crop" global variable.
    real(sp) :: ShapeCGC

    ShapeCGC = crop%StressResponse%ShapeCGC
end function GetCrop_StressResponse_ShapeCGC


function GetCrop_StressResponse_ShapeCCX() result(ShapeCCX)
    !! Getter for the "ShapeCCX" attribute of the "StressResponse" attribute of the "crop" global variable.
    real(sp) :: ShapeCCX

    ShapeCCX = crop%StressResponse%ShapeCCX
end function GetCrop_StressResponse_ShapeCCX


function GetCrop_StressResponse_ShapeWP() result(ShapeWP)
    !! Getter for the "ShapeWP" attribute of the "StressResponse" attribute of the "crop" global variable.
    real(sp) :: ShapeWP

    ShapeWP = crop%StressResponse%ShapeWP
end function GetCrop_StressResponse_ShapeWP


function GetCrop_StressResponse_ShapeCDecline() result(ShapeCDecline)
    !! Getter for the "ShapeCDecline" attribute of the "StressResponse" attribute of the "crop" global variable.
    real(sp) :: ShapeCDecline

    ShapeCDecline = crop%StressResponse%ShapeCDecline
end function GetCrop_StressResponse_ShapeCDecline


function GetCrop_StressResponse_Calibrated() result(Calibrated)
    !! Getter for the "Calibrated" attribute of the "StressResponse" attribute of the "crop" global variable.
    logical :: Calibrated

    Calibrated = crop%StressResponse%Calibrated
end function GetCrop_StressResponse_Calibrated


subroutine SetCrop_StressResponse(StressResponse)
    !! Setter for the "StressResponse" attribute of the "crop" global variable.
    type(rep_Shapes), intent(in) :: StressResponse

    crop%StressResponse = StressResponse
end subroutine SetCrop_StressResponse


subroutine SetCrop_StressResponse_Stress(Stress)
    !! Setter for the "Stress" attribute of the "StressResponse" attribute of the "crop" global variable.
    integer(int8), intent(in) :: Stress

    crop%StressResponse%Stress = Stress
end subroutine SetCrop_StressResponse_Stress


subroutine SetCrop_StressResponse_ShapeCGC(ShapeCGC)
    !! Setter for the "ShapeCGC" attribute of the "StressResponse" attribute of the "crop" global variable.
    real(sp), intent(in) :: ShapeCGC

    crop%StressResponse%ShapeCGC = ShapeCGC
end subroutine SetCrop_StressResponse_ShapeCGC


subroutine SetCrop_StressResponse_ShapeCCX(ShapeCCX)
    !! Setter for the "ShapeCCX" attribute of the "StressResponse" attribute of the "crop" global variable.
    real(sp), intent(in) :: ShapeCCX

    crop%StressResponse%ShapeCCX = ShapeCCX
end subroutine SetCrop_StressResponse_ShapeCCX


subroutine SetCrop_StressResponse_ShapeWP(ShapeWP)
    !! Setter for the "ShapeWP" attribute of the "StressResponse" attribute of the "crop" global variable.
    real(sp), intent(in) :: ShapeWP

    crop%StressResponse%ShapeWP = ShapeWP
end subroutine SetCrop_StressResponse_ShapeWP


subroutine SetCrop_StressResponse_ShapeCDecline(ShapeCDecline)
    !! Setter for the "ShapeCDecline" attribute of the "StressResponse" attribute of the "crop" global variable.
    real(sp), intent(in) :: ShapeCDecline

    crop%StressResponse%ShapeCDecline = ShapeCDecline
end subroutine SetCrop_StressResponse_ShapeCDecline


subroutine SetCrop_StressResponse_Calibrated(Calibrated)
    !! Setter for the "Calibrated" attribute of the "StressResponse" attribute of the "crop" global variable.
    logical, intent(in) :: Calibrated

    crop%StressResponse%Calibrated = Calibrated
end subroutine SetCrop_StressResponse_Calibrated


function GetCrop() result(Crop_out)
    !! Getter for the "crop" global variable.
    type(rep_Crop) :: Crop_out

    Crop_out = crop
end function GetCrop


function GetCrop_subkind() result(subkind)
    !! Getter for the "subkind" attribute of the "crop" global variable.
    integer(intEnum) :: subkind

    subkind = crop%subkind
end function GetCrop_subkind


function GetCrop_ModeCycle() result(ModeCycle)
    !! Getter for the "ModeCycle" attribute of the "crop" global variable.
    integer(intEnum) :: ModeCycle

    ModeCycle = crop%ModeCycle
end function GetCrop_ModeCycle


function GetCrop_Planting() result(Planting)
    !! Getter for the "Planting" attribute of the "crop" global variable.
    integer(intEnum) :: Planting

    Planting = crop%Planting
end function GetCrop_Planting


function GetCrop_pMethod() result(pMethod)
    !! Getter for the "pMethod" attribute of the "crop" global variable.
    integer(intEnum) :: pMethod

    pMethod = crop%pMethod
end function GetCrop_pMethod


function GetCrop_pdef() result(pdef)
    !! Getter for the "pdef" attribute of the "crop" global variable.
    real(sp) :: pdef

    pdef = crop%pdef
end function GetCrop_pdef


function GetCrop_pActStom() result(pActStom)
    !! Getter for the "pActStom" attribute of the "crop" global variable.
    real(sp) :: pActStom

    pActStom = crop%pActStom
end function GetCrop_pActStom


function GetCrop_KsShapeFactorLeaf() result(KsShapeFactorLeaf)
    !! Getter for the "KsShapeFactorLeaf" attribute of the "crop" global variable.
    real(sp) :: KsShapeFactorLeaf

    KsShapeFactorLeaf = crop%KsShapeFactorLeaf
end function GetCrop_KsShapeFactorLeaf


function GetCrop_KsShapeFactorStomata() result(KsShapeFactorStomata)
    !! Getter for the "KsShapeFactorStomata" attribute of the "crop" global variable.
    real(sp) :: KsShapeFactorStomata

    KsShapeFactorStomata = crop%KsShapeFactorStomata
end function GetCrop_KsShapeFactorStomata


function GetCrop_KsShapeFactorSenescence() result(KsShapeFactorSenescence)
    !! Getter for the "KsShapeFactorSenescence" attribute of the "crop" global variable.
    real(sp) :: KsShapeFactorSenescence

    KsShapeFactorSenescence = crop%KsShapeFactorSenescence
end function GetCrop_KsShapeFactorSenescence


function GetCrop_pLeafDefUL() result(pLeafDefUL)
    !! Getter for the "pLeafDefUL" attribute of the "crop" global variable.
    real(sp) :: pLeafDefUL

    pLeafDefUL = crop%pLeafDefUL
end function GetCrop_pLeafDefUL


function GetCrop_pLeafDefLL() result(pLeafDefLL)
    !! Getter for the "pLeafDefLL" attribute of the "crop" global variable.
    real(sp) :: pLeafDefLL

    pLeafDefLL = crop%pLeafDefLL
end function GetCrop_pLeafDefLL


function GetCrop_pLeafAct() result(pLeafAct)
    !! Getter for the "pLeafAct" attribute of the "crop" global variable.
    real(sp) :: pLeafAct

    pLeafAct = crop%pLeafAct
end function GetCrop_pLeafAct


function GetCrop_pSenescence() result(pSenescence)
    !! Getter for the "pSenescence" attribute of the "crop" global variable.
    real(sp) :: pSenescence

    pSenescence = crop%pSenescence
end function GetCrop_pSenescence


function GetCrop_pSenAct() result(pSenAct)
    !! Getter for the "pSenAct" attribute of the "crop" global variable.
    real(sp) :: pSenAct

    pSenAct = crop%pSenAct
end function GetCrop_pSenAct


function GetCrop_pPollination() result(pPollination)
    !! Getter for the "pPollination" attribute of the "crop" global variable.
    real(sp) :: pPollination

    pPollination = crop%pPollination
end function GetCrop_pPollination


function GetCrop_SumEToDelaySenescence() result(SumEToDelaySenescence)
    !! Getter for the "SumEToDelaySenescence" attribute of the "crop" global variable.
    integer(int32) :: SumEToDelaySenescence

    SumEToDelaySenescence = crop%SumEToDelaySenescence
end function GetCrop_SumEToDelaySenescence


function GetCrop_AnaeroPoint() result(AnaeroPoint)
    !! Getter for the "AnaeroPoint" attribute of the "crop" global variable.
    integer(int32) :: AnaeroPoint

    AnaeroPoint = crop%AnaeroPoint
end function GetCrop_AnaeroPoint


function GetCrop_ECemin() result(ECemin)
    !! Getter for the "ECemin" attribute of the "crop" global variable.
    integer(int8) :: ECemin

    ECemin = crop%ECemin
end function GetCrop_ECemin


function GetCrop_ECemax() result(ECemax)
    !! Getter for the "ECemax" attribute of the "crop" global variable.
    integer(int8) :: ECemax

    ECemax = crop%ECemax
end function GetCrop_ECemax


function GetCrop_CCsaltDistortion() result(CCsaltDistortion)
    !! Getter for the "CCsaltDistortion" attribute of the "crop" global variable.
    integer(int8) :: CCsaltDistortion

    CCsaltDistortion = crop%CCsaltDistortion
end function GetCrop_CCsaltDistortion


function GetCrop_ResponseECsw() result(ResponseECsw)
    !! Getter for the "ResponseECsw" attribute of the "crop" global variable.
    integer(int32) :: ResponseECsw

    ResponseECsw = crop%ResponseECsw
end function GetCrop_ResponseECsw


function GetCrop_SmaxTopQuarter() result(SmaxTopQuarter)
    !! Getter for the "SmaxTopQuarter" attribute of the "crop" global variable.
    real(sp) :: SmaxTopQuarter

    SmaxTopQuarter = crop%SmaxTopQuarter
end function GetCrop_SmaxTopQuarter


function GetCrop_SmaxBotQuarter() result(SmaxBotQuarter)
    !! Getter for the "SmaxBotQuarter" attribute of the "crop" global variable.
    real(sp) :: SmaxBotQuarter

    SmaxBotQuarter = crop%SmaxBotQuarter
end function GetCrop_SmaxBotQuarter


function GetCrop_SmaxTop() result(SmaxTop)
    !! Getter for the "SmaxTop" attribute of the "crop" global variable.
    real(sp) :: SmaxTop

    SmaxTop = crop%SmaxTop
end function GetCrop_SmaxTop


function GetCrop_SmaxBot() result(SmaxBot)
    !! Getter for the "SmaxBot" attribute of the "crop" global variable.
    real(sp) :: SmaxBot

    SmaxBot = crop%SmaxBot
end function GetCrop_SmaxBot


function GetCrop_KcTop() result(KcTop)
    !! Getter for the "KcTop" attribute of the "crop" global variable.
    real(sp) :: KcTop

    KcTop = crop%KcTop
end function GetCrop_KcTop


function GetCrop_KcDecline() result(KcDecline)
    !! Getter for the "KcDecline" attribute of the "crop" global variable.
    real(sp) :: KcDecline

    KcDecline = crop%KcDecline
end function GetCrop_KcDecline


function GetCrop_CCEffectEvapLate() result(CCEffectEvapLate)
    !! Getter for the "CCEffectEvapLate" attribute of the "crop" global variable.
    integer(int32) :: CCEffectEvapLate

    CCEffectEvapLate = crop%CCEffectEvapLate
end function GetCrop_CCEffectEvapLate


function GetCrop_Day1() result(Day1)
    !! Getter for the "Day1" attribute of the "crop" global variable.
    integer(int32) :: Day1

    Day1 = crop%Day1
end function GetCrop_Day1


function GetCrop_DayN() result(DayN)
    !! Getter for the "DayN" attribute of the "crop" global variable.
    integer(int32) :: DayN

    DayN = crop%DayN
end function GetCrop_DayN


function GetCrop_Length() result(Length)
    !! Getter for the "Length" attribute of the "crop" global variable.
    integer(int32),dimension(4) :: Length

    Length = crop%Length
end function GetCrop_Length


function GetCrop_RootMin() result(RootMin)
    !! Getter for the "RootMin" attribute of the "crop" global variable.
    real(sp) :: RootMin

    RootMin = crop%RootMin
end function GetCrop_RootMin


function GetCrop_RootMax() result(RootMax)
    !! Getter for the "RootMax" attribute of the "crop" global variable.
    real(sp) :: RootMax

    RootMax = crop%RootMax
end function GetCrop_RootMax


function GetCrop_RootShape() result(RootShape)
    !! Getter for the "RootShape" attribute of the "crop" global variable.
    integer(int8) :: RootShape

    RootShape = crop%RootShape
end function GetCrop_RootShape


function GetCrop_Tbase() result(Tbase)
    !! Getter for the "Tbase" attribute of the "crop" global variable.
    real(sp) :: Tbase

    Tbase = crop%Tbase
end function GetCrop_Tbase


function GetCrop_Tupper() result(Tupper)
    !! Getter for the "Tupper" attribute of the "crop" global variable.
    real(sp) :: Tupper

    Tupper = crop%Tupper
end function GetCrop_Tupper


function GetCrop_Tcold() result(Tcold)
    !! Getter for the "Tcold" attribute of the "crop" global variable.
    integer(int8) :: Tcold

    Tcold = crop%Tcold
end function GetCrop_Tcold


function GetCrop_Theat() result(Theat)
    !! Getter for the "Theat" attribute of the "crop" global variable.
    integer(int8) :: Theat

    Theat = crop%Theat
end function GetCrop_Theat


function GetCrop_GDtranspLow() result(GDtranspLow)
    !! Getter for the "GDtranspLow" attribute of the "crop" global variable.
    real(sp) :: GDtranspLow

    GDtranspLow = crop%GDtranspLow
end function GetCrop_GDtranspLow


function GetCrop_SizeSeedling() result(SizeSeedling)
    !! Getter for the "SizeSeedling" attribute of the "crop" global variable.
    real(sp) :: SizeSeedling

    SizeSeedling = crop%SizeSeedling
end function GetCrop_SizeSeedling


function GetCrop_SizePlant() result(SizePlant)
    !! Getter for the "SizePlant" attribute of the "crop" global variable.
    real(sp) :: SizePlant

    SizePlant = crop%SizePlant
end function GetCrop_SizePlant


function GetCrop_PlantingDens() result(PlantingDens)
    !! Getter for the "PlantingDens" attribute of the "crop" global variable.
    integer(int32) :: PlantingDens

    PlantingDens = crop%PlantingDens
end function GetCrop_PlantingDens


function GetCrop_CCo() result(CCo)
    !! Getter for the "CCo" attribute of the "crop" global variable.
    real(sp) :: CCo

    CCo = crop%CCo
end function GetCrop_CCo


function GetCrop_CCini() result(CCini)
    !! Getter for the "CCini" attribute of the "crop" global variable.
    real(sp) :: CCini

    CCini = crop%CCini
end function GetCrop_CCini


function GetCrop_CGC() result(CGC)
    !! Getter for the "CGC" attribute of the "crop" global variable.
    real(sp) :: CGC

    CGC = crop%CGC
end function GetCrop_CGC


function GetCrop_GDDCGC() result(GDDCGC)
    !! Getter for the "GDDCGC" attribute of the "crop" global variable.
    real(sp) :: GDDCGC

    GDDCGC = crop%GDDCGC
end function GetCrop_GDDCGC


function GetCrop_CCx() result(CCx)
    !! Getter for the "CCx" attribute of the "crop" global variable.
    real(sp) :: CCx

    CCx = crop%CCx
end function GetCrop_CCx


function GetCrop_CDC() result(CDC)
    !! Getter for the "CDC" attribute of the "crop" global variable.
    real(sp) :: CDC

    CDC = crop%CDC
end function GetCrop_CDC


function GetCrop_GDDCDC() result(GDDCDC)
    !! Getter for the "GDDCDC" attribute of the "crop" global variable.
    real(sp) :: GDDCDC

    GDDCDC = crop%GDDCDC
end function GetCrop_GDDCDC


function GetCrop_CCxAdjusted() result(CCxAdjusted)
    !! Getter for the "CCxAdjusted" attribute of the "crop" global variable.
    real(sp) :: CCxAdjusted

    CCxAdjusted = crop%CCxAdjusted
end function GetCrop_CCxAdjusted


function GetCrop_CCxWithered() result(CCxWithered)
    !! Getter for the "CCxWithered" attribute of the "crop" global variable.
    real(sp) :: CCxWithered

    CCxWithered = crop%CCxWithered
end function GetCrop_CCxWithered


function GetCrop_CCoAdjusted() result(CCoAdjusted)
    !! Getter for the "CCoAdjusted" attribute of the "crop" global variable.
    real(sp) :: CCoAdjusted

    CCoAdjusted = crop%CCoAdjusted
end function GetCrop_CCoAdjusted


function GetCrop_DaysToCCini() result(DaysToCCini)
    !! Getter for the "DaysToCCini" attribute of the "crop" global variable.
    integer(int32) :: DaysToCCini

    DaysToCCini = crop%DaysToCCini
end function GetCrop_DaysToCCini


function GetCrop_DaysToGermination() result(DaysToGermination)
    !! Getter for the "DaysToGermination" attribute of the "crop" global variable.
    integer(int32) :: DaysToGermination

    DaysToGermination = crop%DaysToGermination
end function GetCrop_DaysToGermination


function GetCrop_DaysToFullCanopy() result(DaysToFullCanopy)
    !! Getter for the "DaysToFullCanopy" attribute of the "crop" global variable.
    integer(int32) :: DaysToFullCanopy

    DaysToFullCanopy = crop%DaysToFullCanopy
end function GetCrop_DaysToFullCanopy


function GetCrop_DaysToFullCanopySF() result(DaysToFullCanopySF)
    !! Getter for the "DaysToFullCanopySF" attribute of the "crop" global variable.
    integer(int32) :: DaysToFullCanopySF

    DaysToFullCanopySF = crop%DaysToFullCanopySF
end function GetCrop_DaysToFullCanopySF


function GetCrop_DaysToFlowering() result(DaysToFlowering)
    !! Getter for the "DaysToFlowering" attribute of the "crop" global variable.
    integer(int32) :: DaysToFlowering

    DaysToFlowering = crop%DaysToFlowering
end function GetCrop_DaysToFlowering


function GetCrop_LengthFlowering() result(LengthFlowering)
    !! Getter for the "LengthFlowering" attribute of the "crop" global variable.
    integer(int32) :: LengthFlowering

    LengthFlowering = crop%LengthFlowering
end function GetCrop_LengthFlowering


function GetCrop_DaysToSenescence() result(DaysToSenescence)
    !! Getter for the "DaysToSenescence" attribute of the "crop" global variable.
    integer(int32) :: DaysToSenescence

    DaysToSenescence = crop%DaysToSenescence
end function GetCrop_DaysToSenescence


function GetCrop_DaysToHarvest() result(DaysToHarvest)
    !! Getter for the "DaysToHarvest" attribute of the "crop" global variable.
    integer(int32) :: DaysToHarvest

    DaysToHarvest = crop%DaysToHarvest
end function GetCrop_DaysToHarvest


function GetCrop_DaysToMaxRooting() result(DaysToMaxRooting)
    !! Getter for the "DaysToMaxRooting" attribute of the "crop" global variable.
    integer(int32) :: DaysToMaxRooting

    DaysToMaxRooting = crop%DaysToMaxRooting
end function GetCrop_DaysToMaxRooting


function GetCrop_DaysToHIo() result(DaysToHIo)
    !! Getter for the "DaysToHIo" attribute of the "crop" global variable.
    integer(int32) :: DaysToHIo

    DaysToHIo = crop%DaysToHIo
end function GetCrop_DaysToHIo


function GetCrop_GDDaysToCCini() result(GDDaysToCCini)
    !! Getter for the "GDDaysToCCini" attribute of the "crop" global variable.
    integer(int32) :: GDDaysToCCini

    GDDaysToCCini = crop%GDDaysToCCini
end function GetCrop_GDDaysToCCini


function GetCrop_GDDaysToGermination() result(GDDaysToGermination)
    !! Getter for the "GDDaysToGermination" attribute of the "crop" global variable.
    integer(int32) :: GDDaysToGermination

    GDDaysToGermination = crop%GDDaysToGermination
end function GetCrop_GDDaysToGermination


function GetCrop_GDDaysToFullCanopy() result(GDDaysToFullCanopy)
    !! Getter for the "GDDaysToFullCanopy" attribute of the "crop" global variable.
    integer(int32) :: GDDaysToFullCanopy

    GDDaysToFullCanopy = crop%GDDaysToFullCanopy
end function GetCrop_GDDaysToFullCanopy


function GetCrop_GDDaysToFullCanopySF() result(GDDaysToFullCanopySF)
    !! Getter for the "GDDaysToFullCanopySF" attribute of the "crop" global variable.
    integer(int32) :: GDDaysToFullCanopySF

    GDDaysToFullCanopySF = crop%GDDaysToFullCanopySF
end function GetCrop_GDDaysToFullCanopySF


function GetCrop_GDDaysToFlowering() result(GDDaysToFlowering)
    !! Getter for the "GDDaysToFlowering" attribute of the "crop" global variable.
    integer(int32) :: GDDaysToFlowering

    GDDaysToFlowering = crop%GDDaysToFlowering
end function GetCrop_GDDaysToFlowering


function GetCrop_GDDLengthFlowering() result(GDDLengthFlowering)
    !! Getter for the "GDDLengthFlowering" attribute of the "crop" global variable.
    integer(int32) :: GDDLengthFlowering

    GDDLengthFlowering = crop%GDDLengthFlowering
end function GetCrop_GDDLengthFlowering


function GetCrop_GDDaysToSenescence() result(GDDaysToSenescence)
    !! Getter for the "GDDaysToSenescence" attribute of the "crop" global variable.
    integer(int32) :: GDDaysToSenescence

    GDDaysToSenescence = crop%GDDaysToSenescence
end function GetCrop_GDDaysToSenescence


function GetCrop_GDDaysToHarvest() result(GDDaysToHarvest)
    !! Getter for the "GDDaysToHarvest" attribute of the "crop" global variable.
    integer(int32) :: GDDaysToHarvest

    GDDaysToHarvest = crop%GDDaysToHarvest
end function GetCrop_GDDaysToHarvest


function GetCrop_GDDaysToMaxRooting() result(GDDaysToMaxRooting)
    !! Getter for the "GDDaysToMaxRooting" attribute of the "crop" global variable.
    integer(int32) :: GDDaysToMaxRooting

    GDDaysToMaxRooting = crop%GDDaysToMaxRooting
end function GetCrop_GDDaysToMaxRooting


function GetCrop_GDDaysToHIo() result(GDDaysToHIo)
    !! Getter for the "GDDaysToHIo" attribute of the "crop" global variable.
    integer(int32) :: GDDaysToHIo

    GDDaysToHIo = crop%GDDaysToHIo
end function GetCrop_GDDaysToHIo


function GetCrop_WP() result(WP)
    !! Getter for the "WP" attribute of the "crop" global variable.
    real(sp) :: WP

    WP = crop%WP
end function GetCrop_WP


function GetCrop_WPy() result(WPy)
    !! Getter for the "WPy" attribute of the "crop" global variable.
    integer(int32) :: WPy

    WPy = crop%WPy
end function GetCrop_WPy


function GetCrop_AdaptedToCO2() result(AdaptedToCO2)
    !! Getter for the "AdaptedToCO2" attribute of the "crop" global variable.
    integer(int8) :: AdaptedToCO2

    AdaptedToCO2 = crop%AdaptedToCO2
end function GetCrop_AdaptedToCO2


function GetCrop_HI() result(HI)
    !! Getter for the "HI" attribute of the "crop" global variable.
    integer(int32) :: HI

    HI = crop%HI
end function GetCrop_HI


function GetCrop_dHIdt() result(dHIdt)
    !! Getter for the "dHIdt" attribute of the "crop" global variable.
    real(sp) :: dHIdt

    dHIdt = crop%dHIdt
end function GetCrop_dHIdt


function GetCrop_HIincrease() result(HIincrease)
    !! Getter for the "HIincrease" attribute of the "crop" global variable.
    integer(int8) :: HIincrease

    HIincrease = crop%HIincrease
end function GetCrop_HIincrease


function GetCrop_aCoeff() result(aCoeff)
    !! Getter for the "aCoeff" attribute of the "crop" global variable.
    real(sp) :: aCoeff

    aCoeff = crop%aCoeff
end function GetCrop_aCoeff


function GetCrop_bCoeff() result(bCoeff)
    !! Getter for the "bCoeff" attribute of the "crop" global variable.
    real(sp) :: bCoeff

    bCoeff = crop%bCoeff
end function GetCrop_bCoeff


function GetCrop_DHImax() result(DHImax)
    !! Getter for the "DHImax" attribute of the "crop" global variable.
    integer(int8) :: DHImax

    DHImax = crop%DHImax
end function GetCrop_DHImax


function GetCrop_DeterminancyLinked() result(DeterminancyLinked)
    !! Getter for the "DeterminancyLinked" attribute of the "crop" global variable.
    logical :: DeterminancyLinked

    DeterminancyLinked = crop%DeterminancyLinked
end function GetCrop_DeterminancyLinked


function GetCrop_fExcess() result(fExcess)
    !! Getter for the "fExcess" attribute of the "crop" global variable.
    integer(int16) :: fExcess

    fExcess = crop%fExcess
end function GetCrop_fExcess


function GetCrop_DryMatter() result(DryMatter)
    !! Getter for the "DryMatter" attribute of the "crop" global variable.
    integer(int8) :: DryMatter

    DryMatter = crop%DryMatter
end function GetCrop_DryMatter


function GetCrop_RootMinYear1() result(RootMinYear1)
    !! Getter for the "RootMinYear1" attribute of the "crop" global variable.
    real(sp) :: RootMinYear1

    RootMinYear1 = crop%RootMinYear1
end function GetCrop_RootMinYear1


function GetCrop_SownYear1() result(SownYear1)
    !! Getter for the "SownYear1" attribute of the "crop" global variable.
    logical :: SownYear1

    SownYear1 = crop%SownYear1
end function GetCrop_SownYear1


function GetCrop_YearCCx() result(YearCCx)
    !! Getter for the "YearCCx" attribute of the "crop" global variable.
    integer(int8) :: YearCCx

    YearCCx = crop%YearCCx
end function GetCrop_YearCCx


function GetCrop_CCxRoot() result(CCxRoot)
    !! Getter for the "CCxRoot" attribute of the "crop" global variable.
    real(sp) :: CCxRoot

    CCxRoot = crop%CCxRoot
end function GetCrop_CCxRoot


subroutine SetCrop(Crop_in)
    !! Setter for the "crop" global variable.
    type(rep_Crop), intent(in) :: Crop_in

    crop = Crop_in
end subroutine SetCrop


subroutine SetCrop_subkind(subkind)
    !! Setter for the "subkind" attribute of the "crop" global variable.
    integer(intEnum), intent(in) :: subkind

    crop%subkind = subkind
end subroutine SetCrop_subkind


subroutine SetCrop_ModeCycle(ModeCycle)
    !! Setter for the "ModeCycle" attribute of the "crop" global variable.
    integer(intEnum), intent(in) :: ModeCycle

    crop%ModeCycle = ModeCycle
end subroutine SetCrop_ModeCycle


subroutine SetCrop_Planting(Planting)
    !! Setter for the "Planting" attribute of the "crop" global variable.
    integer(intEnum), intent(in) :: Planting

    crop%Planting = Planting
end subroutine SetCrop_Planting


subroutine SetCrop_pMethod(pMethod)
    !! Setter for the "pMethod" attribute of the "crop" global variable.
    integer(intEnum), intent(in) :: pMethod

    crop%pMethod = pMethod
end subroutine SetCrop_pMethod


subroutine SetCrop_pdef(pdef)
    !! Setter for the "pdef" attribute of the "crop" global variable.
    real(sp), intent(in) :: pdef

    crop%pdef = pdef
end subroutine SetCrop_pdef


subroutine SetCrop_pActStom(pActStom)
    !! Setter for the "pActStom" attribute of the "crop" global variable.
    real(sp), intent(in) :: pActStom

    crop%pActStom = pActStom
end subroutine SetCrop_pActStom


subroutine SetCrop_KsShapeFactorLeaf(KsShapeFactorLeaf)
    !! Setter for the "KsShapeFactorLeaf" attribute of the "crop" global variable.
    real(sp), intent(in) :: KsShapeFactorLeaf

    crop%KsShapeFactorLeaf = KsShapeFactorLeaf
end subroutine SetCrop_KsShapeFactorLeaf


subroutine SetCrop_KsShapeFactorStomata(KsShapeFactorStomata)
    !! Setter for the "KsShapeFactorStomata" attribute of the "crop" global variable.
    real(sp), intent(in) :: KsShapeFactorStomata

    crop%KsShapeFactorStomata = KsShapeFactorStomata
end subroutine SetCrop_KsShapeFactorStomata


subroutine SetCrop_KsShapeFactorSenescence(KsShapeFactorSenescence)
    !! Setter for the "KsShapeFactorSenescence" attribute of the "crop" global variable.
    real(sp), intent(in) :: KsShapeFactorSenescence

    crop%KsShapeFactorSenescence = KsShapeFactorSenescence
end subroutine SetCrop_KsShapeFactorSenescence


subroutine SetCrop_pLeafDefUL(pLeafDefUL)
    !! Setter for the "pLeafDefUL" attribute of the "crop" global variable.
    real(sp), intent(in) :: pLeafDefUL

    crop%pLeafDefUL = pLeafDefUL
end subroutine SetCrop_pLeafDefUL


subroutine SetCrop_pLeafDefLL(pLeafDefLL)
    !! Setter for the "pLeafDefLL" attribute of the "crop" global variable.
    real(sp), intent(in) :: pLeafDefLL

    crop%pLeafDefLL = pLeafDefLL
end subroutine SetCrop_pLeafDefLL


subroutine SetCrop_pLeafAct(pLeafAct)
    !! Setter for the "pLeafAct" attribute of the "crop" global variable.
    real(sp), intent(in) :: pLeafAct

    crop%pLeafAct = pLeafAct
end subroutine SetCrop_pLeafAct


subroutine SetCrop_pSenescence(pSenescence)
    !! Setter for the "pSenescence" attribute of the "crop" global variable.
    real(sp), intent(in) :: pSenescence

    crop%pSenescence = pSenescence
end subroutine SetCrop_pSenescence


subroutine SetCrop_pSenAct(pSenAct)
    !! Setter for the "pSenAct" attribute of the "crop" global variable.
    real(sp), intent(in) :: pSenAct

    crop%pSenAct = pSenAct
end subroutine SetCrop_pSenAct


subroutine SetCrop_pPollination(pPollination)
    !! Setter for the "pPollination" attribute of the "crop" global variable.
    real(sp), intent(in) :: pPollination

    crop%pPollination = pPollination
end subroutine SetCrop_pPollination


subroutine SetCrop_SumEToDelaySenescence(SumEToDelaySenescence)
    !! Setter for the "SumEToDelaySenescence" attribute of the "crop" global variable.
    integer(int32), intent(in) :: SumEToDelaySenescence

    crop%SumEToDelaySenescence = SumEToDelaySenescence
end subroutine SetCrop_SumEToDelaySenescence


subroutine SetCrop_AnaeroPoint(AnaeroPoint)
    !! Setter for the "AnaeroPoint" attribute of the "crop" global variable.
    integer(int32), intent(in) :: AnaeroPoint

    crop%AnaeroPoint = AnaeroPoint
end subroutine SetCrop_AnaeroPoint


subroutine SetCrop_ECemin(ECemin)
    !! Setter for the "ECemin" attribute of the "crop" global variable.
    integer(int8), intent(in) :: ECemin

    crop%ECemin = ECemin
end subroutine SetCrop_ECemin


subroutine SetCrop_ECemax(ECemax)
    !! Setter for the "ECemax" attribute of the "crop" global variable.
    integer(int8), intent(in) :: ECemax

    crop%ECemax = ECemax
end subroutine SetCrop_ECemax


subroutine SetCrop_CCsaltDistortion(CCsaltDistortion)
    !! Setter for the "CCsaltDistortion" attribute of the "crop" global variable.
    integer(int8), intent(in) :: CCsaltDistortion

    crop%CCsaltDistortion = CCsaltDistortion
end subroutine SetCrop_CCsaltDistortion


subroutine SetCrop_ResponseECsw(ResponseECsw)
    !! Setter for the "ResponseECsw" attribute of the "crop" global variable.
    integer(int32), intent(in) :: ResponseECsw

    crop%ResponseECsw = ResponseECsw
end subroutine SetCrop_ResponseECsw


subroutine SetCrop_SmaxTopQuarter(SmaxTopQuarter)
    !! Setter for the "SmaxTopQuarter" attribute of the "crop" global variable.
    real(sp), intent(in) :: SmaxTopQuarter

    crop%SmaxTopQuarter = SmaxTopQuarter
end subroutine SetCrop_SmaxTopQuarter


subroutine SetCrop_SmaxBotQuarter(SmaxBotQuarter)
    !! Setter for the "SmaxBotQuarter" attribute of the "crop" global variable.
    real(sp), intent(in) :: SmaxBotQuarter

    crop%SmaxBotQuarter = SmaxBotQuarter
end subroutine SetCrop_SmaxBotQuarter


subroutine SetCrop_SmaxTop(SmaxTop)
    !! Setter for the "SmaxTop" attribute of the "crop" global variable.
    real(sp), intent(in) :: SmaxTop

    crop%SmaxTop = SmaxTop
end subroutine SetCrop_SmaxTop


subroutine SetCrop_SmaxBot(SmaxBot)
    !! Setter for the "SmaxBot" attribute of the "crop" global variable.
    real(sp), intent(in) :: SmaxBot

    crop%SmaxBot = SmaxBot
end subroutine SetCrop_SmaxBot


subroutine SetCrop_KcTop(KcTop)
    !! Setter for the "KcTop" attribute of the "crop" global variable.
    real(sp), intent(in) :: KcTop

    crop%KcTop = KcTop
end subroutine SetCrop_KcTop


subroutine SetCrop_KcDecline(KcDecline)
    !! Setter for the "KcDecline" attribute of the "crop" global variable.
    real(sp), intent(in) :: KcDecline

    crop%KcDecline = KcDecline
end subroutine SetCrop_KcDecline


subroutine SetCrop_CCEffectEvapLate(CCEffectEvapLate)
    !! Setter for the "CCEffectEvapLate" attribute of the "crop" global variable.
    integer(int32), intent(in) :: CCEffectEvapLate

    crop%CCEffectEvapLate = CCEffectEvapLate
end subroutine SetCrop_CCEffectEvapLate


subroutine SetCrop_Day1(Day1)
    !! Setter for the "Day1" attribute of the "crop" global variable.
    integer(int32), intent(in) :: Day1

    crop%Day1 = Day1
end subroutine SetCrop_Day1


subroutine SetCrop_DayN(DayN)
    !! Setter for the "DayN" attribute of the "crop" global variable.
    integer(int32), intent(in) :: DayN

    crop%DayN = DayN
end subroutine SetCrop_DayN


subroutine SetCrop_Length(Length)
    !! Setter for the "Length" attribute of the "crop" global variable.
    integer(int32), dimension(4), intent(in) :: Length

    crop%Length = Length
end subroutine SetCrop_Length


subroutine SetCrop_RootMin(RootMin)
    !! Setter for the "RootMin" attribute of the "crop" global variable.
    real(sp), intent(in) :: RootMin

    crop%RootMin = RootMin
end subroutine SetCrop_RootMin


subroutine SetCrop_RootMax(RootMax)
    !! Setter for the "RootMax" attribute of the "crop" global variable.
    real(sp), intent(in) :: RootMax

    crop%RootMax = RootMax
end subroutine SetCrop_RootMax


subroutine SetCrop_RootShape(RootShape)
    !! Setter for the "RootShape" attribute of the "crop" global variable.
    integer(int8), intent(in) :: RootShape

    crop%RootShape = RootShape
end subroutine SetCrop_RootShape


subroutine SetCrop_Tbase(Tbase)
    !! Setter for the "Tbase" attribute of the "crop" global variable.
    real(sp), intent(in) :: Tbase

    crop%Tbase = Tbase
end subroutine SetCrop_Tbase


subroutine SetCrop_Tupper(Tupper)
    !! Setter for the "Tupper" attribute of the "crop" global variable.
    real(sp), intent(in) :: Tupper

    crop%Tupper = Tupper
end subroutine SetCrop_Tupper


subroutine SetCrop_Tcold(Tcold)
    !! Setter for the "Tcold" attribute of the "crop" global variable.
    integer(int8), intent(in) :: Tcold

    crop%Tcold = Tcold
end subroutine SetCrop_Tcold


subroutine SetCrop_Theat(Theat)
    !! Setter for the "Theat" attribute of the "crop" global variable.
    integer(int8), intent(in) :: Theat

    crop%Theat = Theat
end subroutine SetCrop_Theat


subroutine SetCrop_GDtranspLow(GDtranspLow)
    !! Setter for the "GDtranspLow" attribute of the "crop" global variable.
    real(sp), intent(in) :: GDtranspLow

    crop%GDtranspLow = GDtranspLow
end subroutine SetCrop_GDtranspLow


subroutine SetCrop_SizeSeedling(SizeSeedling)
    !! Setter for the "SizeSeedling" attribute of the "crop" global variable.
    real(sp), intent(in) :: SizeSeedling

    crop%SizeSeedling = SizeSeedling
end subroutine SetCrop_SizeSeedling


subroutine SetCrop_SizePlant(SizePlant)
    !! Setter for the "SizePlant" attribute of the "crop" global variable.
    real(sp), intent(in) :: SizePlant

    crop%SizePlant = SizePlant
end subroutine SetCrop_SizePlant


subroutine SetCrop_PlantingDens(PlantingDens)
    !! Setter for the "PlantingDens" attribute of the "crop" global variable.
    integer(int32), intent(in) :: PlantingDens

    crop%PlantingDens = PlantingDens
end subroutine SetCrop_PlantingDens


subroutine SetCrop_CCo(CCo)
    !! Setter for the "CCo" attribute of the "crop" global variable.
    real(sp), intent(in) :: CCo

    crop%CCo = CCo
end subroutine SetCrop_CCo


subroutine SetCrop_CCini(CCini)
    !! Setter for the "CCini" attribute of the "crop" global variable.
    real(sp), intent(in) :: CCini

    crop%CCini = CCini
end subroutine SetCrop_CCini


subroutine SetCrop_CGC(CGC)
    !! Setter for the "CGC" attribute of the "crop" global variable.
    real(sp), intent(in) :: CGC

    crop%CGC = CGC
end subroutine SetCrop_CGC


subroutine SetCrop_GDDCGC(GDDCGC)
    !! Setter for the "GDDCGC" attribute of the "crop" global variable.
    real(sp), intent(in) :: GDDCGC

    crop%GDDCGC = GDDCGC
end subroutine SetCrop_GDDCGC


subroutine SetCrop_CCx(CCx)
    !! Setter for the "CCx" attribute of the "crop" global variable.
    real(sp), intent(in) :: CCx

    crop%CCx = CCx
end subroutine SetCrop_CCx


subroutine SetCrop_CDC(CDC)
    !! Setter for the "CDC" attribute of the "crop" global variable.
    real(sp), intent(in) :: CDC

    crop%CDC = CDC
end subroutine SetCrop_CDC


subroutine SetCrop_GDDCDC(GDDCDC)
    !! Setter for the "GDDCDC" attribute of the "crop" global variable.
    real(sp), intent(in) :: GDDCDC

    crop%GDDCDC = GDDCDC
end subroutine SetCrop_GDDCDC


subroutine SetCrop_CCxAdjusted(CCxAdjusted)
    !! Setter for the "CCxAdjusted" attribute of the "crop" global variable.
    real(sp), intent(in) :: CCxAdjusted

    crop%CCxAdjusted = CCxAdjusted
end subroutine SetCrop_CCxAdjusted


subroutine SetCrop_CCxWithered(CCxWithered)
    !! Setter for the "CCxWithered" attribute of the "crop" global variable.
    real(sp), intent(in) :: CCxWithered

    crop%CCxWithered = CCxWithered
end subroutine SetCrop_CCxWithered


subroutine SetCrop_CCoAdjusted(CCoAdjusted)
    !! Setter for the "CCoAdjusted" attribute of the "crop" global variable.
    real(sp), intent(in) :: CCoAdjusted

    crop%CCoAdjusted = CCoAdjusted
end subroutine SetCrop_CCoAdjusted


subroutine SetCrop_DaysToCCini(DaysToCCini)
    !! Setter for the "DaysToCCini" attribute of the "crop" global variable.
    integer(int32), intent(in) :: DaysToCCini

    crop%DaysToCCini = DaysToCCini
end subroutine SetCrop_DaysToCCini


subroutine SetCrop_DaysToGermination(DaysToGermination)
    !! Setter for the "DaysToGermination" attribute of the "crop" global variable.
    integer(int32), intent(in) :: DaysToGermination

    crop%DaysToGermination = DaysToGermination
end subroutine SetCrop_DaysToGermination


subroutine SetCrop_DaysToFullCanopy(DaysToFullCanopy)
    !! Setter for the "DaysToFullCanopy" attribute of the "crop" global variable.
    integer(int32), intent(in) :: DaysToFullCanopy

    crop%DaysToFullCanopy = DaysToFullCanopy
end subroutine SetCrop_DaysToFullCanopy


subroutine SetCrop_DaysToFullCanopySF(DaysToFullCanopySF)
    !! Setter for the "DaysToFullCanopySF" attribute of the "crop" global variable.
    integer(int32), intent(in) :: DaysToFullCanopySF

    crop%DaysToFullCanopySF = DaysToFullCanopySF
end subroutine SetCrop_DaysToFullCanopySF


subroutine SetCrop_DaysToFlowering(DaysToFlowering)
    !! Setter for the "DaysToFlowering" attribute of the "crop" global variable.
    integer(int32), intent(in) :: DaysToFlowering

    crop%DaysToFlowering = DaysToFlowering
end subroutine SetCrop_DaysToFlowering


subroutine SetCrop_LengthFlowering(LengthFlowering)
    !! Setter for the "LengthFlowering" attribute of the "crop" global variable.
    integer(int32), intent(in) :: LengthFlowering

    crop%LengthFlowering = LengthFlowering
end subroutine SetCrop_LengthFlowering


subroutine SetCrop_DaysToSenescence(DaysToSenescence)
    !! Setter for the "DaysToSenescence" attribute of the "crop" global variable.
    integer(int32), intent(in) :: DaysToSenescence

    crop%DaysToSenescence = DaysToSenescence
end subroutine SetCrop_DaysToSenescence


subroutine SetCrop_DaysToHarvest(DaysToHarvest)
    !! Setter for the "DaysToHarvest" attribute of the "crop" global variable.
    integer(int32), intent(in) :: DaysToHarvest

    crop%DaysToHarvest = DaysToHarvest
end subroutine SetCrop_DaysToHarvest


subroutine SetCrop_DaysToMaxRooting(DaysToMaxRooting)
    !! Setter for the "DaysToMaxRooting" attribute of the "crop" global variable.
    integer(int32), intent(in) :: DaysToMaxRooting

    crop%DaysToMaxRooting = DaysToMaxRooting
end subroutine SetCrop_DaysToMaxRooting


subroutine SetCrop_DaysToHIo(DaysToHIo)
    !! Setter for the "DaysToHIo" attribute of the "crop" global variable.
    integer(int32), intent(in) :: DaysToHIo

    crop%DaysToHIo = DaysToHIo
end subroutine SetCrop_DaysToHIo


subroutine SetCrop_GDDaysToCCini(GDDaysToCCini)
    !! Setter for the "GDDaysToCCini" attribute of the "crop" global variable.
    integer(int32), intent(in) :: GDDaysToCCini

    crop%GDDaysToCCini = GDDaysToCCini
end subroutine SetCrop_GDDaysToCCini


subroutine SetCrop_GDDaysToGermination(GDDaysToGermination)
    !! Setter for the "GDDaysToGermination" attribute of the "crop" global variable.
    integer(int32), intent(in) :: GDDaysToGermination

    crop%GDDaysToGermination = GDDaysToGermination
end subroutine SetCrop_GDDaysToGermination


subroutine SetCrop_GDDaysToFullCanopy(GDDaysToFullCanopy)
    !! Setter for the "GDDaysToFullCanopy" attribute of the "crop" global variable.
    integer(int32), intent(in) :: GDDaysToFullCanopy

    crop%GDDaysToFullCanopy = GDDaysToFullCanopy
end subroutine SetCrop_GDDaysToFullCanopy


subroutine SetCrop_GDDaysToFullCanopySF(GDDaysToFullCanopySF)
    !! Setter for the "GDDaysToFullCanopySF" attribute of the "crop" global variable.
    integer(int32), intent(in) :: GDDaysToFullCanopySF

    crop%GDDaysToFullCanopySF = GDDaysToFullCanopySF
end subroutine SetCrop_GDDaysToFullCanopySF


subroutine SetCrop_GDDaysToFlowering(GDDaysToFlowering)
    !! Setter for the "GDDaysToFlowering" attribute of the "crop" global variable.
    integer(int32), intent(in) :: GDDaysToFlowering

    crop%GDDaysToFlowering = GDDaysToFlowering
end subroutine SetCrop_GDDaysToFlowering


subroutine SetCrop_GDDLengthFlowering(GDDLengthFlowering)
    !! Setter for the "GDDLengthFlowering" attribute of the "crop" global variable.
    integer(int32), intent(in) :: GDDLengthFlowering

    crop%GDDLengthFlowering = GDDLengthFlowering
end subroutine SetCrop_GDDLengthFlowering


subroutine SetCrop_GDDaysToSenescence(GDDaysToSenescence)
    !! Setter for the "GDDaysToSenescence" attribute of the "crop" global variable.
    integer(int32), intent(in) :: GDDaysToSenescence

    crop%GDDaysToSenescence = GDDaysToSenescence
end subroutine SetCrop_GDDaysToSenescence


subroutine SetCrop_GDDaysToHarvest(GDDaysToHarvest)
    !! Setter for the "GDDaysToHarvest" attribute of the "crop" global variable.
    integer(int32), intent(in) :: GDDaysToHarvest

    crop%GDDaysToHarvest = GDDaysToHarvest
end subroutine SetCrop_GDDaysToHarvest


subroutine SetCrop_GDDaysToMaxRooting(GDDaysToMaxRooting)
    !! Setter for the "GDDaysToMaxRooting" attribute of the "crop" global variable.
    integer(int32), intent(in) :: GDDaysToMaxRooting

    crop%GDDaysToMaxRooting = GDDaysToMaxRooting
end subroutine SetCrop_GDDaysToMaxRooting


subroutine SetCrop_GDDaysToHIo(GDDaysToHIo)
    !! Setter for the "GDDaysToHIo" attribute of the "crop" global variable.
    integer(int32), intent(in) :: GDDaysToHIo

    crop%GDDaysToHIo = GDDaysToHIo
end subroutine SetCrop_GDDaysToHIo


subroutine SetCrop_WP(WP)
    !! Setter for the "WP" attribute of the "crop" global variable.
    real(sp), intent(in) :: WP

    crop%WP = WP
end subroutine SetCrop_WP


subroutine SetCrop_WPy(WPy)
    !! Setter for the "WPy" attribute of the "crop" global variable.
    integer(int32), intent(in) :: WPy

    crop%WPy = WPy
end subroutine SetCrop_WPy


subroutine SetCrop_AdaptedToCO2(AdaptedToCO2)
    !! Setter for the "AdaptedToCO2" attribute of the "crop" global variable.
    integer(int8), intent(in) :: AdaptedToCO2

    crop%AdaptedToCO2 = AdaptedToCO2
end subroutine SetCrop_AdaptedToCO2


subroutine SetCrop_HI(HI)
    !! Setter for the "HI" attribute of the "crop" global variable.
    integer(int32), intent(in) :: HI

    crop%HI = HI
end subroutine SetCrop_HI


subroutine SetCrop_dHIdt(dHIdt)
    !! Setter for the "dHIdt" attribute of the "crop" global variable.
    real(sp), intent(in) :: dHIdt

    crop%dHIdt = dHIdt
end subroutine SetCrop_dHIdt


subroutine SetCrop_HIincrease(HIincrease)
    !! Setter for the "HIincrease" attribute of the "crop" global variable.
    integer(int8), intent(in) :: HIincrease

    crop%HIincrease = HIincrease
end subroutine SetCrop_HIincrease


subroutine SetCrop_aCoeff(aCoeff)
    !! Setter for the "aCoeff" attribute of the "crop" global variable.
    real(sp), intent(in) :: aCoeff

    crop%aCoeff = aCoeff
end subroutine SetCrop_aCoeff


subroutine SetCrop_bCoeff(bCoeff)
    !! Setter for the "bCoeff" attribute of the "crop" global variable.
    real(sp), intent(in) :: bCoeff

    crop%bCoeff = bCoeff
end subroutine SetCrop_bCoeff


subroutine SetCrop_DHImax(DHImax)
    !! Setter for the "DHImax" attribute of the "crop" global variable.
    integer(int8), intent(in) :: DHImax

    crop%DHImax = DHImax
end subroutine SetCrop_DHImax


subroutine SetCrop_DeterminancyLinked(DeterminancyLinked)
    !! Setter for the "DeterminancyLinked" attribute of the "crop" global variable.
    logical, intent(in) :: DeterminancyLinked

    crop%DeterminancyLinked = DeterminancyLinked
end subroutine SetCrop_DeterminancyLinked


subroutine SetCrop_fExcess(fExcess)
    !! Setter for the "fExcess" attribute of the "crop" global variable.
    integer(int16), intent(in) :: fExcess

    crop%fExcess = fExcess
end subroutine SetCrop_fExcess


subroutine SetCrop_DryMatter(DryMatter)
    !! Setter for the "DryMatter" attribute of the "crop" global variable.
    integer(int8), intent(in) :: DryMatter

    crop%DryMatter = DryMatter
end subroutine SetCrop_DryMatter


subroutine SetCrop_RootMinYear1(RootMinYear1)
    !! Setter for the "RootMinYear1" attribute of the "crop" global variable.
    real(sp), intent(in) :: RootMinYear1

    crop%RootMinYear1 = RootMinYear1
end subroutine SetCrop_RootMinYear1


subroutine SetCrop_SownYear1(SownYear1)
    !! Setter for the "SownYear1" attribute of the "crop" global variable.
    logical, intent(in) :: SownYear1

    crop%SownYear1 = SownYear1
end subroutine SetCrop_SownYear1


subroutine SetCrop_YearCCx(YearCCx)
    !! Setter for the "YearCCx" attribute of the "crop" global variable.
    integer(int8), intent(in) :: YearCCx

    crop%YearCCx = YearCCx
end subroutine SetCrop_YearCCx


subroutine SetCrop_CCxRoot(CCxRoot)
    !! Setter for the "CCxRoot" attribute of the "crop" global variable.
    real(sp), intent(in) :: CCxRoot

    crop%CCxRoot = CCxRoot
end subroutine SetCrop_CCxRoot


function GetCrop_Assimilates() result(Assimilates)
    !! Getter for the "Assimilates" attribute of the "crop" global variable.
    type(rep_Assimilates) :: Assimilates

    Assimilates = crop%Assimilates
end function GetCrop_Assimilates


function GetCrop_Assimilates_On() result(On)
    !! Getter for the "On" attribute of the "Assimilates" attribute of the "crop" global variable.
    logical :: On

    On = crop%Assimilates%On
end function GetCrop_Assimilates_On


function GetCrop_Assimilates_Period() result(Period)
    !! Getter for the "Period" attribute of the "Assimilates" attribute of the "crop" global variable.
    integer(int32) :: Period

    Period = crop%Assimilates%Period
end function GetCrop_Assimilates_Period


function GetCrop_Assimilates_Stored() result(Stored)
    !! Getter for the "Stored" attribute of the "Assimilates" attribute of the "crop" global variable.
    integer(int8) :: Stored

    Stored = crop%Assimilates%Stored
end function GetCrop_Assimilates_Stored


function GetCrop_Assimilates_Mobilized() result(Mobilized)
    !! Getter for the "Mobilized" attribute of the "Assimilates" attribute of the "crop" global variable.
    integer(int8) :: Mobilized

    Mobilized = crop%Assimilates%Mobilized
end function GetCrop_Assimilates_Mobilized


subroutine SetCrop_Assimilates(Assimilates)
    !! Setter for the "Assimilates" attribute of the "crop" global variable.
    type(rep_Assimilates), intent(in) :: Assimilates

    crop%Assimilates = Assimilates
end subroutine SetCrop_Assimilates


subroutine SetCrop_Assimilates_On(On)
    !! Setter for the "On" attribute of the "Assimilates" attribute of the "crop" global variable.
    logical, intent(in) :: On

    crop%Assimilates%On = On
end subroutine SetCrop_Assimilates_On


subroutine SetCrop_Assimilates_Period(Period)
    !! Setter for the "Period" attribute of the "Assimilates" attribute of the "crop" global variable.
    integer(int32), intent(in) :: Period

    crop%Assimilates%Period = Period
end subroutine SetCrop_Assimilates_Period


subroutine SetCrop_Assimilates_Stored(Stored)
    !! Setter for the "Stored" attribute of the "Assimilates" attribute of the "crop" global variable.
    integer(int8), intent(in) :: Stored

    crop%Assimilates%Stored = Stored
end subroutine SetCrop_Assimilates_Stored


subroutine SetCrop_Assimilates_Mobilized(Mobilized)
    !! Setter for the "Mobilized" attribute of the "Assimilates" attribute of the "crop" global variable.
    integer(int8), intent(in) :: Mobilized

    crop%Assimilates%Mobilized = Mobilized
end subroutine SetCrop_Assimilates_Mobilized


function GetOnset() result(Onset_out)
    !! Getter for the "onset" global variable.
    type(rep_Onset) :: Onset_out

    Onset_out = onset
end function GetOnset


function GetOnset_GenerateOn() result(GenerateOn)
    !! Getter for the "GenerateOn" attribute of the "onset" global variable.
    logical :: GenerateOn

    GenerateOn = onset%GenerateOn
end function GetOnset_GenerateOn


function GetOnset_GenerateTempOn() result(GenerateTempOn)
    !! Getter for the "GenerateTempOn" attribute of the "onset" global variable.
    logical :: GenerateTempOn

    GenerateTempOn = onset%GenerateTempOn
end function GetOnset_GenerateTempOn


function GetOnset_Criterion() result(Criterion)
    !! Getter for the "Criterion" attribute of the "onset" global variable.
    integer(intEnum) :: Criterion

    Criterion = onset%Criterion
end function GetOnset_Criterion


function GetOnset_AirTCriterion() result(AirTCriterion)
    !! Getter for the "AirTCriterion" attribute of the "onset" global variable.
    integer(intEnum) :: AirTCriterion

    AirTCriterion = onset%AirTCriterion
end function GetOnset_AirTCriterion


function GetOnset_StartSearchDayNr() result(StartSearchDayNr)
    !! Getter for the "StartSearchDayNr" attribute of the "onset" global variable.
    integer(int32) :: StartSearchDayNr

    StartSearchDayNr = onset%StartSearchDayNr
end function GetOnset_StartSearchDayNr


function GetOnset_StopSearchDayNr() result(StopSearchDayNr)
    !! Getter for the "StopSearchDayNr" attribute of the "onset" global variable.
    integer(int32) :: StopSearchDayNr

    StopSearchDayNr = onset%StopSearchDayNr
end function GetOnset_StopSearchDayNr


function GetOnset_LengthSearchPeriod() result(LengthSearchPeriod)
    !! Getter for the "LengthSearchPeriod" attribute of the "onset" global variable.
    integer(int32) :: LengthSearchPeriod

    LengthSearchPeriod = onset%LengthSearchPeriod
end function GetOnset_LengthSearchPeriod


subroutine SetOnset(Onset_in)
    !! Setter for the "onset" global variable.
    type(rep_Onset), intent(in) :: Onset_in

    onset = Onset_in
end subroutine SetOnset


subroutine SetOnset_GenerateOn(GenerateOn)
    !! Setter for the "GenerateOn" attribute of the "onset" global variable.
    logical, intent(in) :: GenerateOn

    onset%GenerateOn = GenerateOn
end subroutine SetOnset_GenerateOn


subroutine SetOnset_GenerateTempOn(GenerateTempOn)
    !! Setter for the "GenerateTempOn" attribute of the "onset" global variable.
    logical, intent(in) :: GenerateTempOn

    onset%GenerateTempOn = GenerateTempOn
end subroutine SetOnset_GenerateTempOn


subroutine SetOnset_Criterion(Criterion)
    !! Setter for the "Criterion" attribute of the "onset" global variable.
    integer(intEnum), intent(in) :: Criterion

    onset%Criterion = Criterion
end subroutine SetOnset_Criterion


subroutine SetOnset_AirTCriterion(AirTCriterion)
    !! Setter for the "AirTCriterion" attribute of the "onset" global variable.
    integer(intEnum), intent(in) :: AirTCriterion

    onset%AirTCriterion = AirTCriterion
end subroutine SetOnset_AirTCriterion


subroutine SetOnset_StartSearchDayNr(StartSearchDayNr)
    !! Setter for the "StartSearchDayNr" attribute of the "onset" global variable.
    integer(int32), intent(in) :: StartSearchDayNr

    onset%StartSearchDayNr = StartSearchDayNr
end subroutine SetOnset_StartSearchDayNr


subroutine SetOnset_StopSearchDayNr(StopSearchDayNr)
    !! Setter for the "StopSearchDayNr" attribute of the "onset" global variable.
    integer(int32), intent(in) :: StopSearchDayNr

    onset%StopSearchDayNr = StopSearchDayNr
end subroutine SetOnset_StopSearchDayNr


subroutine SetOnset_LengthSearchPeriod(LengthSearchPeriod)
    !! Setter for the "LengthSearchPeriod" attribute of the "onset" global variable.
    integer(int32), intent(in) :: LengthSearchPeriod

    onset%LengthSearchPeriod = LengthSearchPeriod
end subroutine SetOnset_LengthSearchPeriod


function GetEndSeason() result(EndSeason_out)
    !! Getter for the "endseason" global variable.
    type(rep_EndSeason) :: EndSeason_out

    EndSeason_out = endseason
end function GetEndSeason


function GetEndSeason_ExtraYears() result(ExtraYears)
    !! Getter for the "ExtraYears" attribute of the "endseason" global variable.
    integer(int32) :: ExtraYears

    ExtraYears = endseason%ExtraYears
end function GetEndSeason_ExtraYears


function GetEndSeason_GenerateTempOn() result(GenerateTempOn)
    !! Getter for the "GenerateTempOn" attribute of the "endseason" global variable.
    logical :: GenerateTempOn

    GenerateTempOn = endseason%GenerateTempOn
end function GetEndSeason_GenerateTempOn


function GetEndSeason_AirTCriterion() result(AirTCriterion)
    !! Getter for the "AirTCriterion" attribute of the "endseason" global variable.
    integer(intEnum) :: AirTCriterion

    AirTCriterion = endseason%AirTCriterion
end function GetEndSeason_AirTCriterion


function GetEndSeason_StartSearchDayNr() result(StartSearchDayNr)
    !! Getter for the "StartSearchDayNr" attribute of the "endseason" global variable.
    integer(int32) :: StartSearchDayNr

    StartSearchDayNr = endseason%StartSearchDayNr
end function GetEndSeason_StartSearchDayNr


function GetEndSeason_StopSearchDayNr() result(StopSearchDayNr)
    !! Getter for the "StopSearchDayNr" attribute of the "endseason" global variable.
    integer(int32) :: StopSearchDayNr

    StopSearchDayNr = endseason%StopSearchDayNr
end function GetEndSeason_StopSearchDayNr


function GetEndSeason_LengthSearchPeriod() result(LengthSearchPeriod)
    !! Getter for the "LengthSearchPeriod" attribute of the "endseason" global variable.
    integer(int32) :: LengthSearchPeriod

    LengthSearchPeriod = endseason%LengthSearchPeriod
end function GetEndSeason_LengthSearchPeriod


subroutine SetEndSeason(EndSeason_in)
    !! Setter for the "endseason" global variable.
    type(rep_EndSeason), intent(in) :: EndSeason_in

    endseason = EndSeason_in
end subroutine SetEndSeason


subroutine SetEndSeason_ExtraYears(ExtraYears)
    !! Setter for the "ExtraYears" attribute of the "endseason" global variable.
    integer(int32), intent(in) :: ExtraYears

    endseason%ExtraYears = ExtraYears
end subroutine SetEndSeason_ExtraYears


subroutine SetEndSeason_GenerateTempOn(GenerateTempOn)
    !! Setter for the "GenerateTempOn" attribute of the "endseason" global variable.
    logical, intent(in) :: GenerateTempOn

    endseason%GenerateTempOn = GenerateTempOn
end subroutine SetEndSeason_GenerateTempOn


subroutine SetEndSeason_AirTCriterion(AirTCriterion)
    !! Setter for the "AirTCriterion" attribute of the "endseason" global variable.
    integer(intEnum), intent(in) :: AirTCriterion

    endseason%AirTCriterion = AirTCriterion
end subroutine SetEndSeason_AirTCriterion


subroutine SetEndSeason_StartSearchDayNr(StartSearchDayNr)
    !! Setter for the "StartSearchDayNr" attribute of the "endseason" global variable.
    integer(int32), intent(in) :: StartSearchDayNr

    endseason%StartSearchDayNr = StartSearchDayNr
end subroutine SetEndSeason_StartSearchDayNr


subroutine SetEndSeason_StopSearchDayNr(StopSearchDayNr)
    !! Setter for the "StopSearchDayNr" attribute of the "endseason" global variable.
    integer(int32), intent(in) :: StopSearchDayNr

    endseason%StopSearchDayNr = StopSearchDayNr
end subroutine SetEndSeason_StopSearchDayNr


subroutine SetEndSeason_LengthSearchPeriod(LengthSearchPeriod)
    !! Setter for the "LengthSearchPeriod" attribute of the "endseason" global variable.
    integer(int32), intent(in) :: LengthSearchPeriod

    endseason%LengthSearchPeriod = LengthSearchPeriod
end subroutine SetEndSeason_LengthSearchPeriod


function GetPerennialPeriod() result(PerennialPeriod_out)
    !! Getter for the "perennialperiod" global variable.
    type(rep_PerennialPeriod) :: PerennialPeriod_out

    PerennialPeriod_out = perennialperiod
end function GetPerennialPeriod


function GetPerennialPeriod_GenerateOnset() result(GenerateOnset)
    !! Getter for the "GenerateOnset" attribute of the "perennialperiod" global variable.
    logical :: GenerateOnset

    GenerateOnset = perennialperiod%GenerateOnset
end function GetPerennialPeriod_GenerateOnset


function GetPerennialPeriod_OnsetCriterion() result(OnsetCriterion)
    !! Getter for the "OnsetCriterion" attribute of the "perennialperiod" global variable.
    integer(intEnum) :: OnsetCriterion

    OnsetCriterion = perennialperiod%OnsetCriterion
end function GetPerennialPeriod_OnsetCriterion


function GetPerennialPeriod_OnsetFirstDay() result(OnsetFirstDay)
    !! Getter for the "OnsetFirstDay" attribute of the "perennialperiod" global variable.
    integer(int32) :: OnsetFirstDay

    OnsetFirstDay = perennialperiod%OnsetFirstDay
end function GetPerennialPeriod_OnsetFirstDay


function GetPerennialPeriod_OnsetFirstMonth() result(OnsetFirstMonth)
    !! Getter for the "OnsetFirstMonth" attribute of the "perennialperiod" global variable.
    integer(int32) :: OnsetFirstMonth

    OnsetFirstMonth = perennialperiod%OnsetFirstMonth
end function GetPerennialPeriod_OnsetFirstMonth


function GetPerennialPeriod_OnsetStartSearchDayNr() result(OnsetStartSearchDayNr)
    !! Getter for the "OnsetStartSearchDayNr" attribute of the "perennialperiod" global variable.
    integer(int32) :: OnsetStartSearchDayNr

    OnsetStartSearchDayNr = perennialperiod%OnsetStartSearchDayNr
end function GetPerennialPeriod_OnsetStartSearchDayNr


function GetPerennialPeriod_OnsetStopSearchDayNr() result(OnsetStopSearchDayNr)
    !! Getter for the "OnsetStopSearchDayNr" attribute of the "perennialperiod" global variable.
    integer(int32) :: OnsetStopSearchDayNr

    OnsetStopSearchDayNr = perennialperiod%OnsetStopSearchDayNr
end function GetPerennialPeriod_OnsetStopSearchDayNr


function GetPerennialPeriod_OnsetLengthSearchPeriod() result(OnsetLengthSearchPeriod)
    !! Getter for the "OnsetLengthSearchPeriod" attribute of the "perennialperiod" global variable.
    integer(int32) :: OnsetLengthSearchPeriod

    OnsetLengthSearchPeriod = perennialperiod%OnsetLengthSearchPeriod
end function GetPerennialPeriod_OnsetLengthSearchPeriod


function GetPerennialPeriod_OnsetThresholdValue() result(OnsetThresholdValue)
    !! Getter for the "OnsetThresholdValue" attribute of the "perennialperiod" global variable.
    real(sp) :: OnsetThresholdValue

    OnsetThresholdValue = perennialperiod%OnsetThresholdValue
end function GetPerennialPeriod_OnsetThresholdValue


function GetPerennialPeriod_OnsetPeriodValue() result(OnsetPeriodValue)
    !! Getter for the "OnsetPeriodValue" attribute of the "perennialperiod" global variable.
    integer(int32) :: OnsetPeriodValue

    OnsetPeriodValue = perennialperiod%OnsetPeriodValue
end function GetPerennialPeriod_OnsetPeriodValue


function GetPerennialPeriod_OnsetOccurrence() result(OnsetOccurrence)
    !! Getter for the "OnsetOccurrence" attribute of the "perennialperiod" global variable.
    integer(int8) :: OnsetOccurrence

    OnsetOccurrence = perennialperiod%OnsetOccurrence
end function GetPerennialPeriod_OnsetOccurrence


function GetPerennialPeriod_GenerateEnd() result(GenerateEnd)
    !! Getter for the "GenerateEnd" attribute of the "perennialperiod" global variable.
    logical :: GenerateEnd

    GenerateEnd = perennialperiod%GenerateEnd
end function GetPerennialPeriod_GenerateEnd


function GetPerennialPeriod_EndCriterion() result(EndCriterion)
    !! Getter for the "EndCriterion" attribute of the "perennialperiod" global variable.
    integer(intEnum) :: EndCriterion

    EndCriterion = perennialperiod%EndCriterion
end function GetPerennialPeriod_EndCriterion


function GetPerennialPeriod_EndLastDay() result(EndLastDay)
    !! Getter for the "EndLastDay" attribute of the "perennialperiod" global variable.
    integer(int32) :: EndLastDay

    EndLastDay = perennialperiod%EndLastDay
end function GetPerennialPeriod_EndLastDay


function GetPerennialPeriod_EndLastMonth() result(EndLastMonth)
    !! Getter for the "EndLastMonth" attribute of the "perennialperiod" global variable.
    integer(int32) :: EndLastMonth

    EndLastMonth = perennialperiod%EndLastMonth
end function GetPerennialPeriod_EndLastMonth


function GetPerennialPeriod_ExtraYears() result(ExtraYears)
    !! Getter for the "ExtraYears" attribute of the "perennialperiod" global variable.
    integer(int32) :: ExtraYears

    ExtraYears = perennialperiod%ExtraYears
end function GetPerennialPeriod_ExtraYears


function GetPerennialPeriod_EndStartSearchDayNr() result(EndStartSearchDayNr)
    !! Getter for the "EndStartSearchDayNr" attribute of the "perennialperiod" global variable.
    integer(int32) :: EndStartSearchDayNr

    EndStartSearchDayNr = perennialperiod%EndStartSearchDayNr
end function GetPerennialPeriod_EndStartSearchDayNr


function GetPerennialPeriod_EndStopSearchDayNr() result(EndStopSearchDayNr)
    !! Getter for the "EndStopSearchDayNr" attribute of the "perennialperiod" global variable.
    integer(int32) :: EndStopSearchDayNr

    EndStopSearchDayNr = perennialperiod%EndStopSearchDayNr
end function GetPerennialPeriod_EndStopSearchDayNr


function GetPerennialPeriod_EndLengthSearchPeriod() result(EndLengthSearchPeriod)
    !! Getter for the "EndLengthSearchPeriod" attribute of the "perennialperiod" global variable.
    integer(int32) :: EndLengthSearchPeriod

    EndLengthSearchPeriod = perennialperiod%EndLengthSearchPeriod
end function GetPerennialPeriod_EndLengthSearchPeriod


function GetPerennialPeriod_EndThresholdValue() result(EndThresholdValue)
    !! Getter for the "EndThresholdValue" attribute of the "perennialperiod" global variable.
    real(sp) :: EndThresholdValue

    EndThresholdValue = perennialperiod%EndThresholdValue
end function GetPerennialPeriod_EndThresholdValue


function GetPerennialPeriod_EndPeriodValue() result(EndPeriodValue)
    !! Getter for the "EndPeriodValue" attribute of the "perennialperiod" global variable.
    integer(int32) :: EndPeriodValue

    EndPeriodValue = perennialperiod%EndPeriodValue
end function GetPerennialPeriod_EndPeriodValue


function GetPerennialPeriod_EndOccurrence() result(EndOccurrence)
    !! Getter for the "EndOccurrence" attribute of the "perennialperiod" global variable.
    integer(int8) :: EndOccurrence

    EndOccurrence = perennialperiod%EndOccurrence
end function GetPerennialPeriod_EndOccurrence


function GetPerennialPeriod_GeneratedDayNrOnset() result(GeneratedDayNrOnset)
    !! Getter for the "GeneratedDayNrOnset" attribute of the "perennialperiod" global variable.
    integer(int32) :: GeneratedDayNrOnset

    GeneratedDayNrOnset = perennialperiod%GeneratedDayNrOnset
end function GetPerennialPeriod_GeneratedDayNrOnset


function GetPerennialPeriod_GeneratedDayNrEnd() result(GeneratedDayNrEnd)
    !! Getter for the "GeneratedDayNrEnd" attribute of the "perennialperiod" global variable.
    integer(int32) :: GeneratedDayNrEnd

    GeneratedDayNrEnd = perennialperiod%GeneratedDayNrEnd
end function GetPerennialPeriod_GeneratedDayNrEnd


subroutine SetPerennialPeriod(PerennialPeriod_in)
    !! Setter for the "perennialperiod" global variable.
    type(rep_PerennialPeriod), intent(in) :: PerennialPeriod_in

    perennialperiod = PerennialPeriod_in
end subroutine SetPerennialPeriod


subroutine SetPerennialPeriod_GenerateOnset(GenerateOnset)
    !! Setter for the "GenerateOnset" attribute of the "perennialperiod" global variable.
    logical, intent(in) :: GenerateOnset

    perennialperiod%GenerateOnset = GenerateOnset
end subroutine SetPerennialPeriod_GenerateOnset


subroutine SetPerennialPeriod_OnsetCriterion(OnsetCriterion)
    !! Setter for the "OnsetCriterion" attribute of the "perennialperiod" global variable.
    integer(intEnum), intent(in) :: OnsetCriterion

    perennialperiod%OnsetCriterion = OnsetCriterion
end subroutine SetPerennialPeriod_OnsetCriterion


subroutine SetPerennialPeriod_OnsetFirstDay(OnsetFirstDay)
    !! Setter for the "OnsetFirstDay" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: OnsetFirstDay

    perennialperiod%OnsetFirstDay = OnsetFirstDay
end subroutine SetPerennialPeriod_OnsetFirstDay


subroutine SetPerennialPeriod_OnsetFirstMonth(OnsetFirstMonth)
    !! Setter for the "OnsetFirstMonth" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: OnsetFirstMonth

    perennialperiod%OnsetFirstMonth = OnsetFirstMonth
end subroutine SetPerennialPeriod_OnsetFirstMonth


subroutine SetPerennialPeriod_OnsetStartSearchDayNr(OnsetStartSearchDayNr)
    !! Setter for the "OnsetStartSearchDayNr" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: OnsetStartSearchDayNr

    perennialperiod%OnsetStartSearchDayNr = OnsetStartSearchDayNr
end subroutine SetPerennialPeriod_OnsetStartSearchDayNr


subroutine SetPerennialPeriod_OnsetStopSearchDayNr(OnsetStopSearchDayNr)
    !! Setter for the "OnsetStopSearchDayNr" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: OnsetStopSearchDayNr

    perennialperiod%OnsetStopSearchDayNr = OnsetStopSearchDayNr
end subroutine SetPerennialPeriod_OnsetStopSearchDayNr


subroutine SetPerennialPeriod_OnsetLengthSearchPeriod(OnsetLengthSearchPeriod)
    !! Setter for the "OnsetLengthSearchPeriod" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: OnsetLengthSearchPeriod

    perennialperiod%OnsetLengthSearchPeriod = OnsetLengthSearchPeriod
end subroutine SetPerennialPeriod_OnsetLengthSearchPeriod


subroutine SetPerennialPeriod_OnsetThresholdValue(OnsetThresholdValue)
    !! Setter for the "OnsetThresholdValue" attribute of the "perennialperiod" global variable.
    real(sp), intent(in) :: OnsetThresholdValue

    perennialperiod%OnsetThresholdValue = OnsetThresholdValue
end subroutine SetPerennialPeriod_OnsetThresholdValue


subroutine SetPerennialPeriod_OnsetPeriodValue(OnsetPeriodValue)
    !! Setter for the "OnsetPeriodValue" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: OnsetPeriodValue

    perennialperiod%OnsetPeriodValue = OnsetPeriodValue
end subroutine SetPerennialPeriod_OnsetPeriodValue


subroutine SetPerennialPeriod_OnsetOccurrence(OnsetOccurrence)
    !! Setter for the "OnsetOccurrence" attribute of the "perennialperiod" global variable.
    integer(int8), intent(in) :: OnsetOccurrence

    perennialperiod%OnsetOccurrence = OnsetOccurrence
end subroutine SetPerennialPeriod_OnsetOccurrence


subroutine SetPerennialPeriod_GenerateEnd(GenerateEnd)
    !! Setter for the "GenerateEnd" attribute of the "perennialperiod" global variable.
    logical, intent(in) :: GenerateEnd

    perennialperiod%GenerateEnd = GenerateEnd
end subroutine SetPerennialPeriod_GenerateEnd


subroutine SetPerennialPeriod_EndCriterion(EndCriterion)
    !! Setter for the "EndCriterion" attribute of the "perennialperiod" global variable.
    integer(intEnum), intent(in) :: EndCriterion

    perennialperiod%EndCriterion = EndCriterion
end subroutine SetPerennialPeriod_EndCriterion


subroutine SetPerennialPeriod_EndLastDay(EndLastDay)
    !! Setter for the "EndLastDay" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: EndLastDay

    perennialperiod%EndLastDay = EndLastDay
end subroutine SetPerennialPeriod_EndLastDay


subroutine SetPerennialPeriod_EndLastMonth(EndLastMonth)
    !! Setter for the "EndLastMonth" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: EndLastMonth

    perennialperiod%EndLastMonth = EndLastMonth
end subroutine SetPerennialPeriod_EndLastMonth


subroutine SetPerennialPeriod_ExtraYears(ExtraYears)
    !! Setter for the "ExtraYears" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: ExtraYears

    perennialperiod%ExtraYears = ExtraYears
end subroutine SetPerennialPeriod_ExtraYears


subroutine SetPerennialPeriod_EndStartSearchDayNr(EndStartSearchDayNr)
    !! Setter for the "EndStartSearchDayNr" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: EndStartSearchDayNr

    perennialperiod%EndStartSearchDayNr = EndStartSearchDayNr
end subroutine SetPerennialPeriod_EndStartSearchDayNr


subroutine SetPerennialPeriod_EndStopSearchDayNr(EndStopSearchDayNr)
    !! Setter for the "EndStopSearchDayNr" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: EndStopSearchDayNr

    perennialperiod%EndStopSearchDayNr = EndStopSearchDayNr
end subroutine SetPerennialPeriod_EndStopSearchDayNr


subroutine SetPerennialPeriod_EndLengthSearchPeriod(EndLengthSearchPeriod)
    !! Setter for the "EndLengthSearchPeriod" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: EndLengthSearchPeriod

    perennialperiod%EndLengthSearchPeriod = EndLengthSearchPeriod
end subroutine SetPerennialPeriod_EndLengthSearchPeriod


subroutine SetPerennialPeriod_EndThresholdValue(EndThresholdValue)
    !! Setter for the "EndThresholdValue" attribute of the "perennialperiod" global variable.
    real(sp), intent(in) :: EndThresholdValue

    perennialperiod%EndThresholdValue = EndThresholdValue
end subroutine SetPerennialPeriod_EndThresholdValue


subroutine SetPerennialPeriod_EndPeriodValue(EndPeriodValue)
    !! Setter for the "EndPeriodValue" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: EndPeriodValue

    perennialperiod%EndPeriodValue = EndPeriodValue
end subroutine SetPerennialPeriod_EndPeriodValue


subroutine SetPerennialPeriod_EndOccurrence(EndOccurrence)
    !! Setter for the "EndOccurrence" attribute of the "perennialperiod" global variable.
    integer(int8), intent(in) :: EndOccurrence

    perennialperiod%EndOccurrence = EndOccurrence
end subroutine SetPerennialPeriod_EndOccurrence


subroutine SetPerennialPeriod_GeneratedDayNrOnset(GeneratedDayNrOnset)
    !! Setter for the "GeneratedDayNrOnset" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: GeneratedDayNrOnset

    perennialperiod%GeneratedDayNrOnset = GeneratedDayNrOnset
end subroutine SetPerennialPeriod_GeneratedDayNrOnset


subroutine SetPerennialPeriod_GeneratedDayNrEnd(GeneratedDayNrEnd)
    !! Setter for the "GeneratedDayNrEnd" attribute of the "perennialperiod" global variable.
    integer(int32), intent(in) :: GeneratedDayNrEnd

    perennialperiod%GeneratedDayNrEnd = GeneratedDayNrEnd
end subroutine SetPerennialPeriod_GeneratedDayNrEnd


type(rep_Content) function GetTotalSaltContent()
    !! Getter for the "TotalSaltContent" global variable.

    GetTotalSaltContent = TotalSaltContent
end function GetTotalSaltContent


real(sp) function GetTotalSaltContent_BeginDay()
    !! Getter for the "TotalSaltContent_BeginDay" global variable.

    GetTotalSaltContent_BeginDay = TotalSaltContent%BeginDay
end function GetTotalSaltContent_BeginDay


real(sp) function GetTotalSaltContent_EndDay()
    !! Getter for the "TotalSaltContent_EndDay" global variable.

    GetTotalSaltContent_EndDay = TotalSaltContent%EndDay
end function GetTotalSaltContent_EndDay


real(sp) function GetTotalSaltContent_ErrorDay()
    !! Getter for the "TotalSaltContent_ErrorDay" global variable.

    GetTotalSaltContent_ErrorDay = TotalSaltContent%ErrorDay
end function GetTotalSaltContent_ErrorDay


type(rep_Content) function GetTotalWaterContent()
    !! Getter for the "TotalWaterContent" global variable.

    GetTotalWaterContent = TotalWaterContent
end function GetTotalWaterContent


real(sp) function GetTotalWaterContent_BeginDay()
    !! Getter for the "TotalWaterContent_BeginDay" global variable.

    GetTotalWaterContent_BeginDay = TotalWaterContent%BeginDay
end function GetTotalWaterContent_BeginDay


real(sp) function GetTotalWaterContent_EndDay()
    !! Getter for the "TotalWaterContent_EndDay" global variable.

    GetTotalWaterContent_EndDay = TotalWaterContent%EndDay
end function GetTotalWaterContent_EndDay


real(sp) function GetTotalWaterContent_ErrorDay()
    !! Getter for the "TotalWaterContent_ErrorDay" global variable.

    GetTotalWaterContent_ErrorDay = TotalWaterContent%ErrorDay
end function GetTotalWaterContent_ErrorDay


subroutine SetTotalSaltContent(TotalSaltContent_in)
    !! Setter for the TotalWaterContent global variable.
    type(rep_content), intent(in) :: TotalSaltContent_in

    TotalSaltContent = TotalSaltContent_in
end subroutine SetTotalSaltContent


subroutine SetTotalSaltContent_BeginDay(BeginDay)
    !! Setter for the "TotalSaltContent" global variable.
    real(sp), intent(in) :: BeginDay

    TotalSaltContent%BeginDay = BeginDay
end subroutine SetTotalSaltContent_BeginDay


subroutine SetTotalSaltContent_EndDay(EndDay)
    !! Setter for the "TotalSaltContent" global variable.
    real(sp), intent(in) :: EndDay

    TotalSaltContent%EndDay = EndDay
end subroutine SetTotalSaltContent_EndDay


subroutine SetTotalSaltContent_ErrorDay(ErrorDay)
    !! Setter for the "TotalSaltContent" global variable.
    real(sp), intent(in) :: ErrorDay

    TotalSaltContent%ErrorDay = ErrorDay
end subroutine SetTotalSaltContent_ErrorDay


subroutine SetTotalWaterContent(TotalWaterContent_in)
    !! Setter for the TotalWaterContent global variable.
    type(rep_content), intent(in) :: TotalWaterContent_in

    TotalWaterContent = TotalWaterContent_in
end subroutine SetTotalWaterContent


subroutine SetTotalWaterContent_BeginDay(BeginDay)
    !! Setter for the "TotalWaterContent" global variable.
    real(sp), intent(in) :: BeginDay

    TotalWaterContent%BeginDay = BeginDay
end subroutine SetTotalWaterContent_BeginDay


subroutine SetTotalWaterContent_EndDay(EndDay)
    !! Setter for the "TotalWaterContent" global variable.
    real(sp), intent(in) :: EndDay

    TotalWaterContent%EndDay = EndDay
end subroutine SetTotalWaterContent_EndDay


subroutine SetTotalWaterContent_ErrorDay(ErrorDay)
    !! Setter for the "TotalWaterContent" global variable.
    real(sp), intent(in) :: ErrorDay

    TotalWaterContent%ErrorDay = ErrorDay
end subroutine SetTotalWaterContent_ErrorDay


type(rep_RootZoneSalt) function GetRootZoneSalt()
    !! Getter for the "RootZoneSalt" global variable.

    GetRootZoneSalt = RootZoneSalt
end function GetRootZoneSalt


real(sp) function GetRootZoneSalt_ECe()
    !! Getter for the "RootZoneSalt" global variable.

    GetRootZoneSalt_ECe = RootZoneSalt%ECe
end function GetRootZoneSalt_ECe


real(sp) function GetRootZoneSalt_ECsw()
    !! Getter for the "RootZoneSalt" global variable.

    GetRootZoneSalt_ECsw = RootZoneSalt%ECsw
end function GetRootZoneSalt_ECsw


real(sp) function GetRootZoneSalt_ECswFC()
    !! Getter for the "RootZoneSalt" global variable.

    GetRootZoneSalt_ECswFC = RootZoneSalt%ECswFC
end function GetRootZoneSalt_ECswFC


real(sp) function GetRootZoneSalt_KsSalt()
    !! Getter for the "RootZoneSalt" global variable.

    GetRootZoneSalt_KsSalt = RootZoneSalt%KsSalt
end function GetRootZoneSalt_KsSalt


subroutine SetRootZoneSalt(RootZoneSalt_in)
    !! Setter for the RootZoneSalt global variable.
    type(rep_rootzonesalt), intent(in) :: RootZoneSalt_in

    RootZoneSalt = RootZoneSalt_in
end subroutine SetRootZoneSalt


subroutine SetRootZoneSalt_ECe(ECe)
    !! Setter for the "RootZoneSalt" global variable.
    real(sp), intent(in) :: ECe

    RootZoneSalt%ECe = ECe
end subroutine SetRootZoneSalt_ECe


subroutine SetRootZoneSalt_ECsw(ECsw)
    !! Setter for the "RootZoneSalt" global variable.
    real(sp), intent(in) :: ECsw

    RootZoneSalt%ECsw = ECsw
end subroutine SetRootZoneSalt_ECsw


subroutine SetRootZoneSalt_ECswFC(ECswFC)
    !! Setter for the "RootZoneSalt" global variable.
    real(sp), intent(in) :: ECswFC

    RootZoneSalt%ECswFC = ECswFC
end subroutine SetRootZoneSalt_ECswFC


subroutine SetRootZoneSalt_KsSalt(KsSalt)
    !! Setter for the "RootZoneSalt" global variable.
    real(sp), intent(in) :: KsSalt

    RootZoneSalt%KsSalt = KsSalt
end subroutine SetRootZoneSalt_KsSalt


integer(intEnum) function GetGenerateTimeMode()
    !! Getter for the "GenerateTimeMode" global variable.

    GetGenerateTimeMode = GenerateTimeMode
end function GetGenerateTimeMode


integer(intEnum) function GetGenerateDepthMode()
    !! Getter for the "GenerateDepthMode" global variable.

    GetGenerateDepthMode = GenerateDepthMode
end function GetGenerateDepthMode


subroutine SetGenerateTimeMode(int_in)
    !! Setter for the "GenerateTimeMode" global variable.
    integer(intEnum), intent(in) :: int_in

    GenerateTimeMode = int_in
end subroutine SetGenerateTimeMode


subroutine SetGenerateDepthMode(int_in)
    !! Setter for the "GenerateDepthMode" global variable.
    integer(intEnum), intent(in) :: int_in

    GenerateDepthMode = int_in
end subroutine SetGenerateDepthMode


integer(intEnum) function GetIrriMode()
    !! Getter for the "IrriMode" global variable.

    GetIrriMode = IrriMode
end function GetIrriMode


integer(intEnum) function GetIrriMethod()
    !! Getter for the "IrriMethod" global variable.

    GetIrriMethod = IrriMethod
end function GetIrriMethod


subroutine SetIrriMode(int_in)
    !! Setter for the "IrriMode" global variable.
    integer(intEnum), intent(in) :: int_in

    IrriMode = int_in
end subroutine SetIrriMode


subroutine SetIrriMethod(int_in)
    !! Setter for the "IrriMethod" global variable.
    integer(intEnum), intent(in) :: int_in

    IrriMethod = int_in
end subroutine SetIrriMethod


function GetTemperatureFile() result(str)
    !! Getter for the "TemperatureFile" global variable.
    character(len=:), allocatable :: str

    str = TemperatureFile
end function GetTemperatureFile


subroutine SetTemperatureFile(str)
    !! Setter for the "TemperatureFile" global variable.
    character(len=*), intent(in) :: str

    TemperatureFile = str
end subroutine SetTemperatureFile


function GetTemperatureFilefull() result(str)
    !! Getter for the "TemperatureFilefull" global variable.
    character(len=:), allocatable :: str

    str = TemperatureFilefull
end function GetTemperatureFilefull


subroutine SetTemperatureFilefull(str)
    !! Setter for the "TemperatureFilefull" global variable.
    character(len=*), intent(in) :: str

    TemperatureFilefull = str
end subroutine SetTemperatureFilefull


function GetTnxReferenceFile() result(str)
    !! Getter for the "TnxReferenceFile" global variable.
    character(len=:), allocatable :: str

    str = TnxReferenceFile
end function GetTnxReferenceFile


subroutine SetTnxReferenceFile(str)
    !! Setter for the "TnxReferenceFile" global variable.
    character(len=*), intent(in) :: str

    TnxReferenceFile = str
end subroutine SetTnxReferenceFile


function GetTnxReferenceFileFull() result(str)
    !! Getter for the "TnxReferenceFileFull" global variable.
    character(len=:), allocatable :: str

    str = TnxReferenceFileFull
end function GetTnxReferenceFileFull


subroutine SetTnxReferenceFileFull(str)
    !! Setter for the "TnxReferenceFileFull" global variable.
    character(len=*), intent(in) :: str

    TnxReferenceFileFull = str
end subroutine SetTnxReferenceFileFull


function GetTnxReference365DaysFile() result(str)
    !! Getter for the "TnxReference365DaysFile" global variable.
    character(len=:), allocatable :: str

    str = TnxReference365DaysFile
end function GetTnxReference365DaysFile


subroutine SetTnxReference365DaysFile(str)
    !! Setter for the "TnxReference365DaysFile" global variable.
    character(len=*), intent(in) :: str

    TnxReference365DaysFile = str
end subroutine SetTnxReference365DaysFile


function GetTnxReference365DaysFileFull() result(str)
    !! Getter for the "TnxReference365DaysFileFull" global variable.
    character(len=:), allocatable :: str

    str = TnxReference365DaysFileFull
end function GetTnxReference365DaysFileFull


subroutine SetTnxReference365DaysFileFull(str)
    !! Setter for the "TnxReference365DaysFileFull" global variable.
    character(len=*), intent(in) :: str

    TnxReference365DaysFileFull = str
end subroutine SetTnxReference365DaysFileFull


function GetTemperatureDescription() result(str)
    !! Getter for the "TemperatureDescription" global variable.
    character(len=:), allocatable :: str

    str = TemperatureDescription
end function GetTemperatureDescription


subroutine SetTemperatureDescription(str)
    !! Setter for the "TemperatureDescription" global variable.
    character(len=*), intent(in) :: str

    TemperatureDescription = str
end subroutine SetTemperatureDescription


function GetClimDescription() result(str)
    !! Getter for the "ClimDescription" global variable.
    character(len=:), allocatable :: str

    str = ClimDescription
end function GetClimDescription


subroutine SetClimDescription(str)
    !! Setter for the "ClimDescription" global variable.
    character(len=*), intent(in) :: str

    ClimDescription = str
end subroutine SetClimDescription


function GetCrop_Length_i(i) result(Length_i)
    !! Getter for the "Length" attribute of "Crop" global variable.
    integer(int32), intent(in) :: i
    integer(int32) :: Length_i

    Length_i = Crop%Length(i)
end function GetCrop_Length_i


subroutine SetCrop_Length_i(i, Length_i)
    !! Setter for the "Length" attribute of "Crop" global variable.
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: Length_i

    Crop%Length(i) = Length_i
end subroutine SetCrop_Length_i


type(rep_clim) function GetTemperatureRecord()
    !! Getter for the "TemperatureRecord" global variable.

    GetTemperatureRecord = TemperatureRecord
end function GetTemperatureRecord


integer(intEnum) function GetTemperatureRecord_DataType()
    !! Getter for the "TemperatureRecord" global variable.

    GetTemperatureRecord_DataType = TemperatureRecord%DataType
end function GetTemperatureRecord_DataType


integer(int32) function GetTemperatureRecord_FromD()
    !! Getter for the "TemperatureRecord" global variable.

    GetTemperatureRecord_FromD = TemperatureRecord%FromD
end function GetTemperatureRecord_FromD


integer(int32) function GetTemperatureRecord_FromM()
    !! Getter for the "TemperatureRecord" global variable.

    GetTemperatureRecord_FromM = TemperatureRecord%FromM
end function GetTemperatureRecord_FromM


integer(int32) function GetTemperatureRecord_FromY()
    !! Getter for the "TemperatureRecord" global variable.

    GetTemperatureRecord_FromY = TemperatureRecord%FromY
end function GetTemperatureRecord_FromY


integer(int32) function GetTemperatureRecord_ToD()
    !! Getter for the "TemperatureRecord" global variable.

    GetTemperatureRecord_ToD = TemperatureRecord%ToD
end function GetTemperatureRecord_ToD


integer(int32) function GetTemperatureRecord_ToM()
    !! Getter for the "TemperatureRecord" global variable.

    GetTemperatureRecord_ToM = TemperatureRecord%ToM
end function GetTemperatureRecord_ToM


integer(int32) function GetTemperatureRecord_ToY()
    !! Getter for the "TemperatureRecord" global variable.

    GetTemperatureRecord_ToY = TemperatureRecord%ToY
end function GetTemperatureRecord_ToY


integer(int32) function GetTemperatureRecord_FromDayNr()
    !! Getter for the "TemperatureRecord" global variable.

    GetTemperatureRecord_FromDayNr = TemperatureRecord%FromDayNr
end function GetTemperatureRecord_FromDayNr


integer(int32) function GetTemperatureRecord_ToDayNr()
    !! Getter for the "TemperatureRecord" global variable.

    GetTemperatureRecord_ToDayNr = TemperatureRecord%ToDayNr
end function GetTemperatureRecord_ToDayNr


function GetTemperatureRecord_FromString() result(str)
    !! Getter for the "TemperatureRecord" global variable.
    character(len=:), allocatable :: str

    str = TemperatureRecord%FromString
end function GetTemperatureRecord_FromString


function GetTemperatureRecord_ToString() result(str)
    !! Getter for the "TemperatureRecord" global variable.
    character(len=:), allocatable :: str

    str = TemperatureRecord%ToString
end function GetTemperatureRecord_ToString


integer(int32) function GetTemperatureRecord_NrObs()
    !! Getter for the "TemperatureRecord" global variable.

    GetTemperatureRecord_NrObs = TemperatureRecord%NrObs
end function GetTemperatureRecord_NrObs


subroutine SetTemperatureRecord(rec_in)
    !! Setter for the "TemperatureRecord" global variable.
    type(rep_clim), intent(in) :: rec_in

    TemperatureRecord = rec_in
end subroutine SetTemperatureRecord


subroutine SetTemperatureRecord_DataType(DataType)
    !! Setter for the "TemperatureRecord" global variable.
    integer(intEnum), intent(in) :: DataType

    TemperatureRecord%DataType = DataType
end subroutine SetTemperatureRecord_DataType


subroutine SetTemperatureRecord_FromD(FromD)
    !! Setter for the "TemperatureRecord" global variable.
    integer(int32), intent(in) :: FromD

    TemperatureRecord%FromD = FromD
end subroutine SetTemperatureRecord_FromD


subroutine SetTemperatureRecord_FromM(FromM)
    !! Setter for the "TemperatureRecord" global variable.
    integer(int32), intent(in) :: FromM

    TemperatureRecord%FromM = FromM
end subroutine SetTemperatureRecord_FromM


subroutine SetTemperatureRecord_FromY(FromY)
    !! Setter for the "TemperatureRecord" global variable.
    integer(int32), intent(in) :: FromY

    TemperatureRecord%FromY = FromY
end subroutine SetTemperatureRecord_FromY


subroutine SetTemperatureRecord_ToD(ToD)
    !! Setter for the "TemperatureRecord" global variable.
    integer(int32), intent(in) :: ToD

    TemperatureRecord%ToD = ToD
end subroutine SetTemperatureRecord_ToD


subroutine SetTemperatureRecord_ToM(ToM)
    !! Setter for the "TemperatureRecord" global variable.
    integer(int32), intent(in) :: ToM

    TemperatureRecord%ToM = ToM
end subroutine SetTemperatureRecord_ToM


subroutine SetTemperatureRecord_TOY(ToY)
    !! Setter for the "TemperatureRecord" global variable.
    integer(int32), intent(in) :: ToY

    TemperatureRecord%ToY = ToY
end subroutine SetTemperatureRecord_ToY


subroutine SetTemperatureRecord_ToDayNr(ToDayNr)
    !! Setter for the "TemperatureRecord" global variable.
    integer(int32), intent(in) :: ToDayNr

    TemperatureRecord%ToDayNr = ToDayNr
end subroutine SetTemperatureRecord_ToDayNr


subroutine SetTemperatureRecord_FromDayNr(FromDayNr)
    !! Setter for the "TemperatureRecord" global variable.
    integer(int32), intent(in) :: FromDayNr

    TemperatureRecord%FromDayNr = FromDayNr
end subroutine SetTemperatureRecord_FromDayNr


subroutine SetTemperatureRecord_NrObs(NrObs)
    !! Setter for the "TemperatureRecord" global variable.
    integer(int32), intent(in) :: NrObs

    TemperatureRecord%NrObs = NrObs
end subroutine SetTemperatureRecord_NrObs


subroutine SetTemperatureRecord_ToString(ToString)
    !! Setter for the "TemperatureRecord" global variable.
    character(len=*), intent(in) :: ToString

    TemperatureRecord%ToString = ToString
end subroutine SetTemperatureRecord_ToString


subroutine SetTemperatureRecord_FromString(FromString)
    !! Setter for the "TemperatureRecord" global variable.
    character(len=*), intent(in) :: FromString

    TemperatureRecord%FromString = FromString
end subroutine SetTemperatureRecord_FromString


function GetIrriAfterSeason_i(i) result(IrriAfterSeason_i)
    !! Getter for individual elements of "IrriAfterSeason" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventInt) :: IrriAfterSeason_i

    IrriAfterSeason_i = IrriAfterSeason(i)
end function GetIrriAfterSeason_i


subroutine SetIrriAfterSeason_i(i, IrriAfterSeason_i)
    !! Setter for individual elements of "IrriAfterSeason" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventInt) :: IrriAfterSeason_i

    IrriAfterSeason(i) = IrriAfterSeason_i
end subroutine SetIrriAfterSeason_i


function GetIrriAfterSeason_DayNr(i) result(DayNr)
    !! Getter for the "DayNr" attribute of the "IrriAfterSeason" global variable.
    integer(int32), intent(in) :: i
    integer(int32) :: DayNr

    DayNr = IrriAfterSeason(i)%DayNr
end function GetIrriAfterSeason_DayNr


function GetIrriAfterSeason_Param(i) result(Param)
    !! Getter for the "Param" attribute of the "IrriAfterSeason" global variable.
    integer(int32), intent(in) :: i
    integer(int32) :: Param

    Param = IrriAfterSeason(i)%Param
end function GetIrriAfterSeason_Param


function GetIrriAfterSeason() result(IrriAfterSeason_out)
    !! Getter for the "IrriAfterSeason" global variable.
    type(rep_DayEventInt), dimension(5) :: IrriAfterSeason_out

    IrriAfterSeason_out = IrriAfterSeason
end function GetIrriAfterSeason


subroutine SetIrriAfterSeason(IrriAfterSeason_in)
    !! Setter for the "IrriAfterSeason" global variable.
    type(rep_DayEventInt), dimension(5), intent(in) :: IrriAfterSeason_in

    IrriAfterSeason = IrriAfterSeason_in
end subroutine SetIrriAfterSeason


subroutine SetIrriAfterSeason_DayNr(i, DayNr)
    !! Setter for the "DayNr" attribute of the "IrriAfterSeason" global variable.
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: DayNr

    IrriAfterSeason(i)%DayNr = DayNr
end subroutine SetIrriAfterSeason_DayNr


subroutine SetIrriAfterSeason_Param(i, Param)
    !! Setter for the "Param" attribute of the "IrriAfterSeason" global variable.
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: Param

    IrriAfterSeason(i)%Param = Param
end subroutine SetIrriAfterSeason_Param


function GetIrriBeforeSeason_i(i) result(IrriBeforeSeason_i)
    !! Getter for individual elements of "IrriBeforeSeason" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventInt) :: IrriBeforeSeason_i

    IrriBeforeSeason_i = IrriBeforeSeason(i)
end function GetIrriBeforeSeason_i


subroutine SetIrriBeforeSeason_i(i, IrriBeforeSeason_i)
    !! Setter for individual elements of "IrriBeforeSeason" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventInt) :: IrriBeforeSeason_i

    IrriBeforeSeason(i) = IrriBeforeSeason_i
end subroutine SetIrriBeforeSeason_i


function GetIrriBeforeSeason_DayNr(i) result(DayNr)
    !! Getter for the "DayNr" attribute of the "IrriBeforeSeason" global variable.
    integer(int32), intent(in) :: i
    integer(int32) :: DayNr

    DayNr = IrriBeforeSeason(i)%DayNr
end function GetIrriBeforeSeason_DayNr


function GetIrriBeforeSeason_Param(i) result(Param)
    !! Getter for the "Param" attribute of the "IrriBeforeSeason" global variable.
    integer(int32), intent(in) :: i
    integer(int32) :: Param

    Param = IrriBeforeSeason(i)%Param
end function GetIrriBeforeSeason_Param


function GetIrriBeforeSeason() result(IrriBeforeSeason_out)
    !! Getter for the "IrriBeforeSeason" global variable.
    type(rep_DayEventInt), dimension(5) :: IrriBeforeSeason_out

    IrriBeforeSeason_out = IrriBeforeSeason
end function GetIrriBeforeSeason


subroutine SetIrriBeforeSeason(IrriBeforeSeason_in)
    !! Setter for the "IrriBeforeSeason" global variable.
    type(rep_DayEventInt), dimension(5), intent(in) :: IrriBeforeSeason_in

    IrriBeforeSeason = IrriBeforeSeason_in
end subroutine SetIrriBeforeSeason


subroutine SetIrriBeforeSeason_DayNr(i, DayNr)
    !! Setter for the "DayNr" attribute of the "IrriBeforeSeason" global variable.
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: DayNr

    IrriBeforeSeason(i)%DayNr = DayNr
end subroutine SetIrriBeforeSeason_DayNr


subroutine SetIrriBeforeSeason_Param(i, Param)
    !! Setter for the "Param" attribute of the "IrriBeforeSeason" global variable.
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: Param

    IrriBeforeSeason(i)%Param = Param
end subroutine SetIrriBeforeSeason_Param


integer(int32) function GetIrriFirstDayNr()
    !! Getter for the "IrriFirstDayNr" global variable.

    GetIrriFirstDayNr = IrriFirstDayNr
end function GetIrriFirstDayNr


subroutine SetIrriFirstDayNr(IrriFirstDayNr_in)
    !! Setter for the "IrriFirstDayNr" global variable.
    integer(int32), intent(in) :: IrriFirstDayNr_in

    IrriFirstDayNr = IrriFirstDayNr_in
end subroutine SetIrriFirstDayNr


type(rep_clim) function GetClimRecord()
    !! Getter for the "ClimRecord" global variable.

    GetClimRecord = ClimRecord
end function GetClimRecord


integer(intEnum) function GetClimRecord_DataType()
    !! Getter for the "ClimRecord" global variable.

    GetClimRecord_DataType = ClimRecord%DataType
end function GetClimRecord_DataType


integer(int32) function GetClimRecord_FromD()
    !! Getter for the "ClimRecord" global variable.

    GetClimRecord_FromD = ClimRecord%FromD
end function GetClimRecord_FromD


integer(int32) function GetClimRecord_FromM()
    !! Getter for the "ClimRecord" global variable.

    GetClimRecord_FromM = ClimRecord%FromM
end function GetClimRecord_FromM


integer(int32) function GetClimRecord_FromY()
    !! Getter for the "ClimRecord" global variable.

    GetClimRecord_FromY = ClimRecord%FromY
end function GetClimRecord_FromY


integer(int32) function GetClimRecord_ToD()
    !! Getter for the "ClimRecord" global variable.

    GetClimRecord_ToD = ClimRecord%ToD
end function GetClimRecord_ToD


integer(int32) function GetClimRecord_ToM()
    !! Getter for the "ClimRecord" global variable.

    GetClimRecord_ToM = ClimRecord%ToM
end function GetClimRecord_ToM


integer(int32) function GetClimRecord_ToY()
    !! Getter for the "ClimRecord" global variable.

    GetClimRecord_ToY = ClimRecord%ToY
end function GetClimRecord_ToY


integer(int32) function GetClimRecord_FromDayNr()
    !! Getter for the "ClimRecord" global variable.

    GetClimRecord_FromDayNr = ClimRecord%FromDayNr
end function GetClimRecord_FromDayNr


integer(int32) function GetClimRecord_ToDayNr()
    !! Getter for the "ClimRecord" global variable.

    GetClimRecord_ToDayNr = ClimRecord%ToDayNr
end function GetClimRecord_ToDayNr


function GetClimRecord_FromString() result(str)
    !! Getter for the "ClimRecord" global variable.
    character(len=:), allocatable :: str

    str = ClimRecord%FromString
end function GetClimRecord_FromString


function GetClimRecord_ToString() result(str)
    !! Getter for the "ClimRecord" global variable.
    character(len=:), allocatable :: str

    str = ClimRecord%ToString
end function GetClimRecord_ToString


integer(int32) function GetClimRecord_NrObs()
    !! Getter for the "ClimRecord" global variable.

    GetClimRecord_NrObs = ClimRecord%NrObs
end function GetClimRecord_NrObs


subroutine SetClimRecord(rec_in)
    !! Setter for the "ClimRecord" global variable.
    type(rep_clim), intent(in) :: rec_in

    ClimRecord = rec_in
end subroutine SetClimRecord


subroutine SetClimRecord_DataType(DataType)
    !! Setter for the "ClimRecord" global variable.
    integer(intEnum), intent(in) :: DataType

    ClimRecord%DataType = DataType
end subroutine SetClimRecord_DataType


subroutine SetClimRecord_FromD(FromD)
    !! Setter for the "ClimRecord" global variable.
    integer(int32), intent(in) :: FromD

    ClimRecord%FromD = FromD
end subroutine SetClimRecord_FromD


subroutine SetClimRecord_FromM(FromM)
    !! Setter for the "ClimRecord" global variable.
    integer(int32), intent(in) :: FromM

    ClimRecord%FromM = FromM
end subroutine SetClimRecord_FromM


subroutine SetClimRecord_FromY(FromY)
    !! Setter for the "ClimRecord" global variable.
    integer(int32), intent(in) :: FromY

    ClimRecord%FromY = FromY
end subroutine SetClimRecord_FromY


subroutine SetClimRecord_ToD(ToD)
    !! Setter for the "ClimRecord" global variable.
    integer(int32), intent(in) :: ToD

    ClimRecord%ToD = ToD
end subroutine SetClimRecord_ToD


subroutine SetClimRecord_ToM(ToM)
    !! Setter for the "ClimRecord" global variable.
    integer(int32), intent(in) :: ToM

    ClimRecord%ToM = ToM
end subroutine SetClimRecord_ToM


subroutine SetClimRecord_TOY(ToY)
    !! Setter for the "ClimRecord" global variable.
    integer(int32), intent(in) :: ToY

    ClimRecord%ToY = ToY
end subroutine SetClimRecord_ToY


subroutine SetClimRecord_ToDayNr(ToDayNr)
    !! Setter for the "ClimRecord" global variable.
    integer(int32), intent(in) :: ToDayNr

    ClimRecord%ToDayNr = ToDayNr
end subroutine SetClimRecord_ToDayNr


subroutine SetClimRecord_FromDayNr(FromDayNr)
    !! Setter for the "ClimRecord" global variable.
    integer(int32), intent(in) :: FromDayNr

    ClimRecord%FromDayNr = FromDayNr
end subroutine SetClimRecord_FromDayNr


subroutine SetClimRecord_NrObs(NrObs)
    !! Setter for the "ClimRecord" global variable.
    integer(int32), intent(in) :: NrObs

    ClimRecord%NrObs = NrObs
end subroutine SetClimRecord_NrObs


subroutine SetClimRecord_ToString(ToString)
    !! Setter for the "ClimRecord" global variable.
    character(len=*), intent(in) :: ToString

    ClimRecord%ToString = ToString
end subroutine SetClimRecord_ToString


subroutine SetClimRecord_FromString(FromString)
    !! Setter for the "ClimRecord" global variable.
    character(len=*), intent(in) :: FromString

    ClimRecord%FromString = FromString
end subroutine SetClimRecord_FromString


type(rep_clim) function GetRainRecord()
    !! Getter for the "RainRecord" global variable.

    GetRainRecord = RainRecord
end function GetRainRecord


integer(intEnum) function GetRainRecord_DataType()
    !! Getter for the "RainRecord" global variable.

    GetRainRecord_DataType = RainRecord%DataType
end function GetRainRecord_DataType


integer(int32) function GetRainRecord_FromD()
    !! Getter for the "RainRecord" global variable.

    GetRainRecord_FromD = RainRecord%FromD
end function GetRainRecord_FromD


integer(int32) function GetRainRecord_FromM()
    !! Getter for the "RainRecord" global variable.

    GetRainRecord_FromM = RainRecord%FromM
end function GetRainRecord_FromM


integer(int32) function GetRainRecord_FromY()
    !! Getter for the "RainRecord" global variable.

    GetRainRecord_FromY = RainRecord%FromY
end function GetRainRecord_FromY


integer(int32) function GetRainRecord_ToD()
    !! Getter for the "RainRecord" global variable.

    GetRainRecord_ToD = RainRecord%ToD
end function GetRainRecord_ToD


integer(int32) function GetRainRecord_ToM()
    !! Getter for the "RainRecord" global variable.

    GetRainRecord_ToM = RainRecord%ToM
end function GetRainRecord_ToM


integer(int32) function GetRainRecord_ToY()
    !! Getter for the "RainRecord" global variable.

    GetRainRecord_ToY = RainRecord%ToY
end function GetRainRecord_ToY


integer(int32) function GetRainRecord_FromDayNr()
    !! Getter for the "RainRecord" global variable.

    GetRainRecord_FromDayNr = RainRecord%FromDayNr
end function GetRainRecord_FromDayNr


integer(int32) function GetRainRecord_ToDayNr()
    !! Getter for the "RainRecord" global variable.

    GetRainRecord_ToDayNr = RainRecord%ToDayNr
end function GetRainRecord_ToDayNr


function GetRainRecord_FromString() result(str)
    !! Getter for the "RainRecord" global variable.
    character(len=:), allocatable :: str

    str = RainRecord%FromString
end function GetRainRecord_FromString


function GetRainRecord_ToString() result(str)
    !! Getter for the "RainRecord" global variable.
    character(len=:), allocatable :: str

    str = RainRecord%ToString
end function GetRainRecord_ToString


integer(int32) function GetRainRecord_NrObs()
    !! Getter for the "RainRecord" global variable.

    GetRainRecord_NrObs = RainRecord%NrObs
end function GetRainRecord_NrObs


subroutine SetRainRecord(rec_in)
    !! Setter for the "RainRecord" global variable.
    type(rep_clim), intent(in) :: rec_in

    RainRecord = rec_in
end subroutine SetRainRecord


subroutine SetRainRecord_DataType(DataType)
    !! Setter for the "RainRecord" global variable.
    integer(intEnum), intent(in) :: DataType

    RainRecord%DataType = DataType
end subroutine SetRainRecord_DataType


subroutine SetRainRecord_FromD(FromD)
    !! Setter for the "RainRecord" global variable.
    integer(int32), intent(in) :: FromD

    RainRecord%FromD = FromD
end subroutine SetRainRecord_FromD


subroutine SetRainRecord_FromM(FromM)
    !! Setter for the "RainRecord" global variable.
    integer(int32), intent(in) :: FromM

    RainRecord%FromM = FromM
end subroutine SetRainRecord_FromM


subroutine SetRainRecord_FromY(FromY)
    !! Setter for the "RainRecord" global variable.
    integer(int32), intent(in) :: FromY

    RainRecord%FromY = FromY
end subroutine SetRainRecord_FromY


subroutine SetRainRecord_ToD(ToD)
    !! Setter for the "RainRecord" global variable.
    integer(int32), intent(in) :: ToD

    RainRecord%ToD = ToD
end subroutine SetRainRecord_ToD


subroutine SetRainRecord_ToM(ToM)
    !! Setter for the "RainRecord" global variable.
    integer(int32), intent(in) :: ToM

    RainRecord%ToM = ToM
end subroutine SetRainRecord_ToM


subroutine SetRainRecord_TOY(ToY)
    !! Setter for the "RainRecord" global variable.
    integer(int32), intent(in) :: ToY

    RainRecord%ToY = ToY
end subroutine SetRainRecord_ToY


subroutine SetRainRecord_ToDayNr(ToDayNr)
    !! Setter for the "RainRecord" global variable.
    integer(int32), intent(in) :: ToDayNr

    RainRecord%ToDayNr = ToDayNr
end subroutine SetRainRecord_ToDayNr


subroutine SetRainRecord_FromDayNr(FromDayNr)
    !! Setter for the "RainRecord" global variable.
    integer(int32), intent(in) :: FromDayNr

    RainRecord%FromDayNr = FromDayNr
end subroutine SetRainRecord_FromDayNr


subroutine SetRainRecord_NrObs(NrObs)
    !! Setter for the "RainRecord" global variable.
    integer(int32), intent(in) :: NrObs

    RainRecord%NrObs = NrObs
end subroutine SetRainRecord_NrObs


subroutine SetRainRecord_ToString(ToString)
    !! Setter for the "RainRecord" global variable.
    character(len=*), intent(in) :: ToString

    RainRecord%ToString = ToString
end subroutine SetRainRecord_ToString


subroutine SetRainRecord_FromString(FromString)
    !! Setter for the "RainRecord" global variable.
    character(len=*), intent(in) :: FromString

    RainRecord%FromString = FromString
end subroutine SetRainRecord_FromString


type(rep_clim) function GetEToRecord()
    !! Getter for the "EToRecord" global variable.

    GetEToRecord = EToRecord
end function GetEToRecord


integer(intEnum) function GetEToRecord_DataType()
    !! Getter for the "EToRecord" global variable.

    GetEToRecord_DataType = EToRecord%DataType
end function GetEToRecord_DataType


integer(int32) function GetEToRecord_FromD()
    !! Getter for the "EToRecord" global variable.

    GetEToRecord_FromD = EToRecord%FromD
end function GetEToRecord_FromD


integer(int32) function GetEToRecord_FromM()
    !! Getter for the "EToRecord" global variable.

    GetEToRecord_FromM = EToRecord%FromM
end function GetEToRecord_FromM


integer(int32) function GetEToRecord_FromY()
    !! Getter for the "EToRecord" global variable.

    GetEToRecord_FromY = EToRecord%FromY
end function GetEToRecord_FromY


integer(int32) function GetEToRecord_ToD()
    !! Getter for the "EToRecord" global variable.

    GetEToRecord_ToD = EToRecord%ToD
end function GetEToRecord_ToD


integer(int32) function GetEToRecord_ToM()
    !! Getter for the "EToRecord" global variable.

    GetEToRecord_ToM = EToRecord%ToM
end function GetEToRecord_ToM


integer(int32) function GetEToRecord_ToY()
    !! Getter for the "EToRecord" global variable.

    GetEToRecord_ToY = EToRecord%ToY
end function GetEToRecord_ToY


integer(int32) function GetEToRecord_FromDayNr()
    !! Getter for the "EToRecord" global variable.

    GetEToRecord_FromDayNr = EToRecord%FromDayNr
end function GetEToRecord_FromDayNr


integer(int32) function GetEToRecord_ToDayNr()
    !! Getter for the "EToRecord" global variable.

    GetEToRecord_ToDayNr = EToRecord%ToDayNr
end function GetEToRecord_ToDayNr


function GetEToRecord_FromString() result(str)
    !! Getter for the "EToRecord" global variable.
    character(len=:), allocatable :: str

    str = EToRecord%FromString
end function GetEToRecord_FromString


function GetEToRecord_ToString() result(str)
    !! Getter for the "EToRecord" global variable.
    character(len=:), allocatable :: str

    str = EToRecord%ToString
end function GetEToRecord_ToString


integer(int32) function GetEToRecord_NrObs()
    !! Getter for the "EToRecord" global variable.

    GetEToRecord_NrObs = EToRecord%NrObs
end function GetEToRecord_NrObs


subroutine SetEToRecord(rec_in)
    !! Setter for the "EToRecord" global variable.
    type(rep_clim), intent(in) :: rec_in

    EToRecord = rec_in
end subroutine SetEToRecord


subroutine SetEToRecord_DataType(DataType)
    !! Setter for the "EToRecord" global variable.
    integer(intEnum), intent(in) :: DataType

    EToRecord%DataType = DataType
end subroutine SetEToRecord_DataType


subroutine SetEToRecord_FromD(FromD)
    !! Setter for the "EToRecord" global variable.
    integer(int32), intent(in) :: FromD

    EToRecord%FromD = FromD
end subroutine SetEToRecord_FromD


subroutine SetEToRecord_FromM(FromM)
    !! Setter for the "EToRecord" global variable.
    integer(int32), intent(in) :: FromM

    EToRecord%FromM = FromM
end subroutine SetEToRecord_FromM


subroutine SetEToRecord_FromY(FromY)
    !! Setter for the "EToRecord" global variable.
    integer(int32), intent(in) :: FromY

    EToRecord%FromY = FromY
end subroutine SetEToRecord_FromY


subroutine SetEToRecord_ToD(ToD)
    !! Setter for the "EToRecord" global variable.
    integer(int32), intent(in) :: ToD

    EToRecord%ToD = ToD
end subroutine SetEToRecord_ToD


subroutine SetEToRecord_ToM(ToM)
    !! Setter for the "EToRecord" global variable.
    integer(int32), intent(in) :: ToM

    EToRecord%ToM = ToM
end subroutine SetEToRecord_ToM


subroutine SetEToRecord_TOY(ToY)
    !! Setter for the "EToRecord" global variable.
    integer(int32), intent(in) :: ToY

    EToRecord%ToY = ToY
end subroutine SetEToRecord_ToY


subroutine SetEToRecord_ToDayNr(ToDayNr)
    !! Setter for the "EToRecord" global variable.
    integer(int32), intent(in) :: ToDayNr

    EToRecord%ToDayNr = ToDayNr
end subroutine SetEToRecord_ToDayNr


subroutine SetEToRecord_FromDayNr(FromDayNr)
    !! Setter for the "EToRecord" global variable.
    integer(int32), intent(in) :: FromDayNr

    EToRecord%FromDayNr = FromDayNr
end subroutine SetEToRecord_FromDayNr


subroutine SetEToRecord_NrObs(NrObs)
    !! Setter for the "EToRecord" global variable.
    integer(int32), intent(in) :: NrObs

    EToRecord%NrObs = NrObs
end subroutine SetEToRecord_NrObs


subroutine SetEToRecord_ToString(ToString)
    !! Setter for the "EToRecord" global variable.
    character(len=*), intent(in) :: ToString

    EToRecord%ToString = ToString
end subroutine SetEToRecord_ToString


subroutine SetEToRecord_FromString(FromString)
    !! Setter for the "EToRecord" global variable.
    character(len=*), intent(in) :: FromString

    EToRecord%FromString = FromString
end subroutine SetEToRecord_FromString


function GetSimulation() result(Simulation_out)
    !! Getter for the "simulation" global variable.
    type(rep_sim) :: Simulation_out

    Simulation_out = simulation
end function GetSimulation


function GetSimulation_FromDayNr() result(FromDayNr)
    !! Getter for the "FromDayNr" attribute of the "simulation" global variable.
    integer(int32) :: FromDayNr

    FromDayNr = simulation%FromDayNr
end function GetSimulation_FromDayNr


function GetSimulation_ToDayNr() result(ToDayNr)
    !! Getter for the "ToDayNr" attribute of the "simulation" global variable.
    integer(int32) :: ToDayNr

    ToDayNr = simulation%ToDayNr
end function GetSimulation_ToDayNr


function GetSimulation_ThetaIni_i(i) result(ThetaIni_i)
    !! Getter for the "ThetaIni" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: ThetaIni_i

    ThetaIni_i = simulation%ThetaIni(i)
end function GetSimulation_ThetaIni_i


function GetSimulation_ECeIni_i(i) result(ECeIni_i)
    !! Getter for the "ECeIni" attribute of the "simulation" global variable.
    integer(int32) :: i
    real(sp) :: ECeIni_i

    ECeIni_i = simulation%ECeIni(i)
end function GetSimulation_ECeIni_i


function GetSimulation_SurfaceStorageIni() result(SurfaceStorageIni)
    !! Getter for the "SurfaceStorageIni" attribute of the "simulation" global variable.
    real(sp) :: SurfaceStorageIni

    SurfaceStorageIni = simulation%SurfaceStorageIni
end function GetSimulation_SurfaceStorageIni


function GetSimulation_ECStorageIni() result(ECStorageIni)
    !! Getter for the "ECStorageIni" attribute of the "simulation" global variable.
    real(sp) :: ECStorageIni

    ECStorageIni = simulation%ECStorageIni
end function GetSimulation_ECStorageIni


function GetSimulation_CCini() result(CCini)
    !! Getter for the "CCini" attribute of the "simulation" global variable.
    real(sp) :: CCini

    CCini = simulation%CCini
end function GetSimulation_CCini


function GetSimulation_Bini() result(Bini)
    !! Getter for the "Bini" attribute of the "simulation" global variable.
    real(sp) :: Bini

    Bini = simulation%Bini
end function GetSimulation_Bini


function GetSimulation_Zrini() result(Zrini)
    !! Getter for the "Zrini" attribute of the "simulation" global variable.
    real(sp) :: Zrini

    Zrini = simulation%Zrini
end function GetSimulation_Zrini


function GetSimulation_LinkCropToSimPeriod() result(LinkCropToSimPeriod)
    !! Getter for the "LinkCropToSimPeriod" attribute of the "simulation" global variable.
    logical :: LinkCropToSimPeriod

    LinkCropToSimPeriod = simulation%LinkCropToSimPeriod
end function GetSimulation_LinkCropToSimPeriod


function GetSimulation_ResetIniSWC() result(ResetIniSWC)
    !! Getter for the "ResetIniSWC" attribute of the "simulation" global variable.
    logical :: ResetIniSWC

    ResetIniSWC = simulation%ResetIniSWC
end function GetSimulation_ResetIniSWC


function GetSimulation_InitialStep() result(InitialStep)
    !! Getter for the "InitialStep" attribute of the "simulation" global variable.
    integer(int32) :: InitialStep

    InitialStep = simulation%InitialStep
end function GetSimulation_InitialStep


function GetSimulation_EvapLimitON() result(EvapLimitON)
    !! Getter for the "EvapLimitON" attribute of the "simulation" global variable.
    logical :: EvapLimitON

    EvapLimitON = simulation%EvapLimitON
end function GetSimulation_EvapLimitON


function GetSimulation_EvapWCsurf() result(EvapWCsurf)
    !! Getter for the "EvapWCsurf" attribute of the "simulation" global variable.
    real(sp) :: EvapWCsurf

    EvapWCsurf = simulation%EvapWCsurf
end function GetSimulation_EvapWCsurf


function GetSimulation_EvapStartStg2() result(EvapStartStg2)
    !! Getter for the "EvapStartStg2" attribute of the "simulation" global variable.
    integer(int8) :: EvapStartStg2

    EvapStartStg2 = simulation%EvapStartStg2
end function GetSimulation_EvapStartStg2


function GetSimulation_EvapZ() result(EvapZ)
    !! Getter for the "EvapZ" attribute of the "simulation" global variable.
    real(sp) :: EvapZ

    EvapZ = simulation%EvapZ
end function GetSimulation_EvapZ


function GetSimulation_HIfinal() result(HIfinal)
    !! Getter for the "HIfinal" attribute of the "simulation" global variable.
    integer(int32) :: HIfinal

    HIfinal = simulation%HIfinal
end function GetSimulation_HIfinal


function GetSimulation_DelayedDays() result(DelayedDays)
    !! Getter for the "DelayedDays" attribute of the "simulation" global variable.
    integer(int32) :: DelayedDays

    DelayedDays = simulation%DelayedDays
end function GetSimulation_DelayedDays


function GetSimulation_Germinate() result(Germinate)
    !! Getter for the "Germinate" attribute of the "simulation" global variable.
    logical :: Germinate

    Germinate = simulation%Germinate
end function GetSimulation_Germinate


function GetSimulation_SumEToStress() result(SumEToStress)
    !! Getter for the "SumEToStress" attribute of the "simulation" global variable.
    real(sp) :: SumEToStress

    SumEToStress = simulation%SumEToStress
end function GetSimulation_SumEToStress


function GetSimulation_SumGDD() result(SumGDD)
    !! Getter for the "SumGDD" attribute of the "simulation" global variable.
    real(sp) :: SumGDD

    SumGDD = simulation%SumGDD
end function GetSimulation_SumGDD


function GetSimulation_SumGDDfromDay1() result(SumGDDfromDay1)
    !! Getter for the "SumGDDfromDay1" attribute of the "simulation" global variable.
    real(sp) :: SumGDDfromDay1

    SumGDDfromDay1 = simulation%SumGDDfromDay1
end function GetSimulation_SumGDDfromDay1


function GetSimulation_SCor() result(SCor)
    !! Getter for the "SCor" attribute of the "simulation" global variable.
    real(sp) :: SCor

    SCor = simulation%SCor
end function GetSimulation_SCor


function GetSimulation_MultipleRun() result(MultipleRun)
    !! Getter for the "MultipleRun" attribute of the "simulation" global variable.
    logical :: MultipleRun

    MultipleRun = simulation%MultipleRun
end function GetSimulation_MultipleRun


function GetSimulation_NrRuns() result(NrRuns)
    !! Getter for the "NrRuns" attribute of the "simulation" global variable.
    integer(int32) :: NrRuns

    NrRuns = simulation%NrRuns
end function GetSimulation_NrRuns


function GetSimulation_MultipleRunWithKeepSWC() result(MultipleRunWithKeepSWC)
    !! Getter for the "MultipleRunWithKeepSWC" attribute of the "simulation" global variable.
    logical :: MultipleRunWithKeepSWC

    MultipleRunWithKeepSWC = simulation%MultipleRunWithKeepSWC
end function GetSimulation_MultipleRunWithKeepSWC


function GetSimulation_MultipleRunConstZrx() result(MultipleRunConstZrx)
    !! Getter for the "MultipleRunConstZrx" attribute of the "simulation" global variable.
    real(sp) :: MultipleRunConstZrx

    MultipleRunConstZrx = simulation%MultipleRunConstZrx
end function GetSimulation_MultipleRunConstZrx


function GetSimulation_IrriECw() result(IrriECw)
    !! Getter for the "IrriECw" attribute of the "simulation" global variable.
    real(sp) :: IrriECw

    IrriECw = simulation%IrriECw
end function GetSimulation_IrriECw


function GetSimulation_DayAnaero() result(DayAnaero)
    !! Getter for the "DayAnaero" attribute of the "simulation" global variable.
    integer(int8) :: DayAnaero

    DayAnaero = simulation%DayAnaero
end function GetSimulation_DayAnaero


function GetSimulation_SalinityConsidered() result(SalinityConsidered)
    !! Getter for the "SalinityConsidered" attribute of the "simulation" global variable.
    logical :: SalinityConsidered

    SalinityConsidered = simulation%SalinityConsidered
end function GetSimulation_SalinityConsidered


function GetSimulation_ProtectedSeedling() result(ProtectedSeedling)
    !! Getter for the "ProtectedSeedling" attribute of the "simulation" global variable.
    logical :: ProtectedSeedling

    ProtectedSeedling = simulation%ProtectedSeedling
end function GetSimulation_ProtectedSeedling


function GetSimulation_SWCtopSoilConsidered() result(SWCtopSoilConsidered)
    !! Getter for the "SWCtopSoilConsidered" attribute of the "simulation" global variable.
    logical :: SWCtopSoilConsidered

    SWCtopSoilConsidered = simulation%SWCtopSoilConsidered
end function GetSimulation_SWCtopSoilConsidered


function GetSimulation_LengthCuttingInterval() result(LengthCuttingInterval)
    !! Getter for the "LengthCuttingInterval" attribute of the "simulation" global variable.
    integer(int32) :: LengthCuttingInterval

    LengthCuttingInterval = simulation%LengthCuttingInterval
end function GetSimulation_LengthCuttingInterval


function GetSimulation_YearSeason() result(YearSeason)
    !! Getter for the "YearSeason" attribute of the "simulation" global variable.
    integer(int8) :: YearSeason

    YearSeason = simulation%YearSeason
end function GetSimulation_YearSeason


function GetSimulation_RCadj() result(RCadj)
    !! Getter for the "RCadj" attribute of the "simulation" global variable.
    integer(int8) :: RCadj

    RCadj = simulation%RCadj
end function GetSimulation_RCadj


function GetSimulation_YearStartCropCycle() result(YearStartCropCycle)
    !! Getter for the "YearStartCropCycle" attribute of the "simulation" global variable.
    integer(int32) :: YearStartCropCycle

    YearStartCropCycle = simulation%YearStartCropCycle
end function GetSimulation_YearStartCropCycle


function GetSimulation_CropDay1Previous() result(CropDay1Previous)
    !! Getter for the "CropDay1Previous" attribute of the "simulation" global variable.
    integer(int32) :: CropDay1Previous

    CropDay1Previous = simulation%CropDay1Previous
end function GetSimulation_CropDay1Previous


subroutine SetSimulation(Simulation_in)
    !! Setter for the "simulation" global variable.
    type(rep_sim), intent(in) :: Simulation_in

    simulation = Simulation_in
end subroutine SetSimulation


subroutine SetSimulation_FromDayNr(FromDayNr)
    !! Setter for the "FromDayNr" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: FromDayNr

    simulation%FromDayNr = FromDayNr
end subroutine SetSimulation_FromDayNr


subroutine SetSimulation_ToDayNr(ToDayNr)
    !! Setter for the "ToDayNr" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: ToDayNr

    simulation%ToDayNr = ToDayNr
end subroutine SetSimulation_ToDayNr


subroutine SetSimulation_ThetaIni_i(i, ThetaIni_i)
    !! Setter for the "ThetaIni" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: ThetaIni_i

    simulation%ThetaIni(i) = ThetaIni_i
end subroutine SetSimulation_ThetaIni_i


subroutine SetSimulation_ECeIni_i(i, ECeIni_i)
    !! Setter for the "ECeIni" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: ECeIni_i

    simulation%ECeIni(i) = ECeIni_i
end subroutine SetSimulation_ECeIni_i


subroutine SetSimulation_SurfaceStorageIni(SurfaceStorageIni)
    !! Setter for the "SurfaceStorageIni" attribute of the "simulation" global variable.
    real(sp), intent(in) :: SurfaceStorageIni

    simulation%SurfaceStorageIni = SurfaceStorageIni
end subroutine SetSimulation_SurfaceStorageIni


subroutine SetSimulation_ECStorageIni(ECStorageIni)
    !! Setter for the "ECStorageIni" attribute of the "simulation" global variable.
    real(sp), intent(in) :: ECStorageIni

    simulation%ECStorageIni = ECStorageIni
end subroutine SetSimulation_ECStorageIni


subroutine SetSimulation_CCini(CCini)
    !! Setter for the "CCini" attribute of the "simulation" global variable.
    real(sp), intent(in) :: CCini

    simulation%CCini = CCini
end subroutine SetSimulation_CCini


subroutine SetSimulation_Bini(Bini)
    !! Setter for the "Bini" attribute of the "simulation" global variable.
    real(sp), intent(in) :: Bini

    simulation%Bini = Bini
end subroutine SetSimulation_Bini


subroutine SetSimulation_Zrini(Zrini)
    !! Setter for the "Zrini" attribute of the "simulation" global variable.
    real(sp), intent(in) :: Zrini

    simulation%Zrini = Zrini
end subroutine SetSimulation_Zrini


subroutine SetSimulation_LinkCropToSimPeriod(LinkCropToSimPeriod)
    !! Setter for the "LinkCropToSimPeriod" attribute of the "simulation" global variable.
    logical, intent(in) :: LinkCropToSimPeriod

    simulation%LinkCropToSimPeriod = LinkCropToSimPeriod
end subroutine SetSimulation_LinkCropToSimPeriod


subroutine SetSimulation_ResetIniSWC(ResetIniSWC)
    !! Setter for the "ResetIniSWC" attribute of the "simulation" global variable.
    logical, intent(in) :: ResetIniSWC

    simulation%ResetIniSWC = ResetIniSWC
end subroutine SetSimulation_ResetIniSWC


subroutine SetSimulation_InitialStep(InitialStep)
    !! Setter for the "InitialStep" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: InitialStep

    simulation%InitialStep = InitialStep
end subroutine SetSimulation_InitialStep


subroutine SetSimulation_EvapLimitON(EvapLimitON)
    !! Setter for the "EvapLimitON" attribute of the "simulation" global variable.
    logical, intent(in) :: EvapLimitON

    simulation%EvapLimitON = EvapLimitON
end subroutine SetSimulation_EvapLimitON


subroutine SetSimulation_EvapWCsurf(EvapWCsurf)
    !! Setter for the "EvapWCsurf" attribute of the "simulation" global variable.
    real(sp), intent(in) :: EvapWCsurf

    simulation%EvapWCsurf = EvapWCsurf
end subroutine SetSimulation_EvapWCsurf


subroutine SetSimulation_EvapStartStg2(EvapStartStg2)
    !! Setter for the "EvapStartStg2" attribute of the "simulation" global variable.
    integer(int8), intent(in) :: EvapStartStg2

    simulation%EvapStartStg2 = EvapStartStg2
end subroutine SetSimulation_EvapStartStg2


subroutine SetSimulation_EvapZ(EvapZ)
    !! Setter for the "EvapZ" attribute of the "simulation" global variable.
    real(sp), intent(in) :: EvapZ

    simulation%EvapZ = EvapZ
end subroutine SetSimulation_EvapZ


subroutine SetSimulation_HIfinal(HIfinal)
    !! Setter for the "HIfinal" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: HIfinal

    simulation%HIfinal = HIfinal
end subroutine SetSimulation_HIfinal


subroutine SetSimulation_DelayedDays(DelayedDays)
    !! Setter for the "DelayedDays" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: DelayedDays

    simulation%DelayedDays = DelayedDays
end subroutine SetSimulation_DelayedDays


subroutine SetSimulation_Germinate(Germinate)
    !! Setter for the "Germinate" attribute of the "simulation" global variable.
    logical, intent(in) :: Germinate

    simulation%Germinate = Germinate
end subroutine SetSimulation_Germinate


subroutine SetSimulation_SumEToStress(SumEToStress)
    !! Setter for the "SumEToStress" attribute of the "simulation" global variable.
    real(sp), intent(in) :: SumEToStress

    simulation%SumEToStress = SumEToStress
end subroutine SetSimulation_SumEToStress


subroutine SetSimulation_SumGDD(SumGDD)
    !! Setter for the "SumGDD" attribute of the "simulation" global variable.
    real(sp), intent(in) :: SumGDD

    simulation%SumGDD = SumGDD
end subroutine SetSimulation_SumGDD


subroutine SetSimulation_SumGDDfromDay1(SumGDDfromDay1)
    !! Setter for the "SumGDDfromDay1" attribute of the "simulation" global variable.
    real(sp), intent(in) :: SumGDDfromDay1

    simulation%SumGDDfromDay1 = SumGDDfromDay1
end subroutine SetSimulation_SumGDDfromDay1


subroutine SetSimulation_SCor(SCor)
    !! Setter for the "SCor" attribute of the "simulation" global variable.
    real(sp), intent(in) :: SCor

    simulation%SCor = SCor
end subroutine SetSimulation_SCor


subroutine SetSimulation_MultipleRun(MultipleRun)
    !! Setter for the "MultipleRun" attribute of the "simulation" global variable.
    logical, intent(in) :: MultipleRun

    simulation%MultipleRun = MultipleRun
end subroutine SetSimulation_MultipleRun


subroutine SetSimulation_NrRuns(NrRuns)
    !! Setter for the "NrRuns" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: NrRuns

    simulation%NrRuns = NrRuns
end subroutine SetSimulation_NrRuns


subroutine SetSimulation_MultipleRunWithKeepSWC(MultipleRunWithKeepSWC)
    !! Setter for the "MultipleRunWithKeepSWC" attribute of the "simulation" global variable.
    logical, intent(in) :: MultipleRunWithKeepSWC

    simulation%MultipleRunWithKeepSWC = MultipleRunWithKeepSWC
end subroutine SetSimulation_MultipleRunWithKeepSWC


subroutine SetSimulation_MultipleRunConstZrx(MultipleRunConstZrx)
    !! Setter for the "MultipleRunConstZrx" attribute of the "simulation" global variable.
    real(sp), intent(in) :: MultipleRunConstZrx

    simulation%MultipleRunConstZrx = MultipleRunConstZrx
end subroutine SetSimulation_MultipleRunConstZrx


subroutine SetSimulation_IrriECw(IrriECw)
    !! Setter for the "IrriECw" attribute of the "simulation" global variable.
    real(sp), intent(in) :: IrriECw

    simulation%IrriECw = IrriECw
end subroutine SetSimulation_IrriECw


subroutine SetSimulation_DayAnaero(DayAnaero)
    !! Setter for the "DayAnaero" attribute of the "simulation" global variable.
    integer(int8), intent(in) :: DayAnaero

    simulation%DayAnaero = DayAnaero
end subroutine SetSimulation_DayAnaero


subroutine SetSimulation_SalinityConsidered(SalinityConsidered)
    !! Setter for the "SalinityConsidered" attribute of the "simulation" global variable.
    logical, intent(in) :: SalinityConsidered

    simulation%SalinityConsidered = SalinityConsidered
end subroutine SetSimulation_SalinityConsidered


subroutine SetSimulation_ProtectedSeedling(ProtectedSeedling)
    !! Setter for the "ProtectedSeedling" attribute of the "simulation" global variable.
    logical, intent(in) :: ProtectedSeedling

    simulation%ProtectedSeedling = ProtectedSeedling
end subroutine SetSimulation_ProtectedSeedling


subroutine SetSimulation_SWCtopSoilConsidered(SWCtopSoilConsidered)
    !! Setter for the "SWCtopSoilConsidered" attribute of the "simulation" global variable.
    logical, intent(in) :: SWCtopSoilConsidered

    simulation%SWCtopSoilConsidered = SWCtopSoilConsidered
end subroutine SetSimulation_SWCtopSoilConsidered


subroutine SetSimulation_LengthCuttingInterval(LengthCuttingInterval)
    !! Setter for the "LengthCuttingInterval" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: LengthCuttingInterval

    simulation%LengthCuttingInterval = LengthCuttingInterval
end subroutine SetSimulation_LengthCuttingInterval


subroutine SetSimulation_YearSeason(YearSeason)
    !! Setter for the "YearSeason" attribute of the "simulation" global variable.
    integer(int8), intent(in) :: YearSeason

    simulation%YearSeason = YearSeason
end subroutine SetSimulation_YearSeason


subroutine SetSimulation_RCadj(RCadj)
    !! Setter for the "RCadj" attribute of the "simulation" global variable.
    integer(int8), intent(in) :: RCadj

    simulation%RCadj = RCadj
end subroutine SetSimulation_RCadj


subroutine SetSimulation_YearStartCropCycle(YearStartCropCycle)
    !! Setter for the "YearStartCropCycle" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: YearStartCropCycle

    simulation%YearStartCropCycle = YearStartCropCycle
end subroutine SetSimulation_YearStartCropCycle


subroutine SetSimulation_CropDay1Previous(CropDay1Previous)
    !! Setter for the "CropDay1Previous" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: CropDay1Previous

    simulation%CropDay1Previous = CropDay1Previous
end subroutine SetSimulation_CropDay1Previous


function GetSimulation_IniSWC() result(IniSWC)
    !! Getter for the "IniSWC" attribute of the "simulation" global variable.
    type(rep_IniSWC) :: IniSWC

    IniSWC = simulation%IniSWC
end function GetSimulation_IniSWC


function GetSimulation_IniSWC_AtDepths() result(AtDepths)
    !! Getter for the "AtDepths" attribute of the "IniSWC" attribute of the "simulation" global variable.
    logical :: AtDepths

    AtDepths = simulation%IniSWC%AtDepths
end function GetSimulation_IniSWC_AtDepths


function GetSimulation_IniSWC_NrLoc() result(NrLoc)
    !! Getter for the "NrLoc" attribute of the "IniSWC" attribute of the "simulation" global variable.
    integer(int8) :: NrLoc

    NrLoc = simulation%IniSWC%NrLoc
end function GetSimulation_IniSWC_NrLoc


function GetSimulation_IniSWC_Loc_i(i) result(Loc_i)
    !! Getter for the "Loc" attribute of the "IniSWC" attribute of the "simulation" global variable
    integer(int32), intent(in) :: i
    real(sp) :: Loc_i

    Loc_i = simulation%IniSWC%Loc(i)
end function GetSimulation_IniSWC_Loc_i


function GetSimulation_IniSWC_Loc() result(Loc)
    !! Getter for the "Loc" attribute of the "IniSWC" attribute of the "simulation" global variable

    real(sp), dimension(max_No_compartments) :: Loc

    Loc = simulation%IniSWC%Loc
end function GetSimulation_IniSWC_Loc


function GetSimulation_IniSWC_VolProc_i(i) result(VolProc_i)
    !! Getter for the "VolProc" attribute of the "IniSWC" attribute of the "simulation" global variable
    integer(int32), intent(in) :: i
    real(sp) :: VolProc_i

    VolProc_i = simulation%IniSWC%VolProc(i)
end function GetSimulation_IniSWC_VolProc_i


function GetSimulation_IniSWC_VolProc() result(VolProc)
    !! Getter for the "VolProc" attribute of the "IniSWC" attribute of the "simulation" global variable

    real(sp), dimension(max_No_compartments) :: VolProc

    VolProc = simulation%IniSWC%VolProc
end function GetSimulation_IniSWC_VolProc


function GetSimulation_IniSWC_SaltECe_i(i) result(SaltECe_i)
    !! Getter for the "SaltECe" attribute of the "IniSWC" attribute of the "simulation" global variable
    integer(int32), intent(in) :: i
    real(sp) :: SaltECe_i

    SaltECe_i = simulation%IniSWC%SaltECe(i)
end function GetSimulation_IniSWC_SaltECe_i


function GetSimulation_IniSWC_SaltECe() result(SaltECe)
    !! Getter for the "SaltECe" attribute of the "IniSWC" attribute of the "simulation" global variable

    real(sp), dimension(max_No_compartments) :: SaltECe

    SaltECe = simulation%IniSWC%SaltECe
end function GetSimulation_IniSWC_SaltECe


function GetSimulation_IniSWC_AtFC() result(AtFC)
    !! Getter for the "AtFC" attribute of the "IniSWC" attribute of the "simulation" global variable.
    logical :: AtFC

    AtFC = simulation%IniSWC%AtFC
end function GetSimulation_IniSWC_AtFC


subroutine SetSimulation_IniSWC(IniSWC)
    !! Setter for the "IniSWC" attribute of the "simulation" global variable.
    type(rep_IniSWC), intent(in) :: IniSWC

    simulation%IniSWC = IniSWC
end subroutine SetSimulation_IniSWC


subroutine SetSimulation_IniSWC_AtDepths(AtDepths)
    !! Setter for the "AtDepths" attribute of the "IniSWC" attribute of the "simulation" global variable.
    logical, intent(in) :: AtDepths

    simulation%IniSWC%AtDepths = AtDepths
end subroutine SetSimulation_IniSWC_AtDepths


subroutine SetSimulation_IniSWC_NrLoc(NrLoc)
    !! Setter for the "NrLoc" attribute of the "IniSWC" attribute of the "simulation" global variable.
    integer(int8), intent(in) :: NrLoc

    simulation%IniSWC%NrLoc = NrLoc
end subroutine SetSimulation_IniSWC_NrLoc


subroutine SetSimulation_IniSWC_Loc_i(i, Loc_i)
    !! Setter for the "Loc" attribute of the "IniSWC" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: Loc_i

    simulation%IniSWC%Loc(i) = Loc_i
end subroutine SetSimulation_IniSWC_Loc_i


subroutine SetSimulation_IniSWC_VolProc_i(i, VolProc_i)
    !! Setter for the "VolProc" attribute of the "IniSWC" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: VolProc_i

    simulation%IniSWC%VolProc(i) = VolProc_i
end subroutine SetSimulation_IniSWC_VolProc_i


subroutine SetSimulation_IniSWC_SaltECe_i(i, SaltECe_i)
    !! Setter for the "SaltECe" attribute of the "IniSWC" attribute of the "simulation" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: SaltECe_i

    simulation%IniSWC%SaltECe(i) = SaltECe_i
end subroutine SetSimulation_IniSWC_SaltECe_i
!!END ATTEMPT

subroutine SetSimulation_IniSWC_AtFC(AtFC)
    !! Setter for the "AtFC" attribute of the "IniSWC" attribute of the "simulation" global variable.
    logical, intent(in) :: AtFC

    simulation%IniSWC%AtFC = AtFC
end subroutine SetSimulation_IniSWC_AtFC


function GetSimulation_EffectStress() result(EffectStress)
    !! Getter for the "EffectStress" attribute of the "simulation" global variable.
    type(rep_EffectStress) :: EffectStress

    EffectStress = simulation%EffectStress
end function GetSimulation_EffectStress


function GetSimulation_EffectStress_RedCGC() result(RedCGC)
    !! Getter for the "RedCGC" attribute of the "EffectStress" attribute of the "simulation" global variable.
    integer(int8) :: RedCGC

    RedCGC = simulation%EffectStress%RedCGC
end function GetSimulation_EffectStress_RedCGC


function GetSimulation_EffectStress_RedCCX() result(RedCCX)
    !! Getter for the "RedCCX" attribute of the "EffectStress" attribute of the "simulation" global variable.
    integer(int8) :: RedCCX

    RedCCX = simulation%EffectStress%RedCCX
end function GetSimulation_EffectStress_RedCCX


function GetSimulation_EffectStress_RedWP() result(RedWP)
    !! Getter for the "RedWP" attribute of the "EffectStress" attribute of the "simulation" global variable.
    integer(int8) :: RedWP

    RedWP = simulation%EffectStress%RedWP
end function GetSimulation_EffectStress_RedWP


function GetSimulation_EffectStress_CDecline() result(CDecline)
    !! Getter for the "CDecline" attribute of the "EffectStress" attribute of the "simulation" global variable.
    real(sp) :: CDecline

    CDecline = simulation%EffectStress%CDecline
end function GetSimulation_EffectStress_CDecline


function GetSimulation_EffectStress_RedKsSto() result(RedKsSto)
    !! Getter for the "RedKsSto" attribute of the "EffectStress" attribute of the "simulation" global variable.
    integer(int8) :: RedKsSto

    RedKsSto = simulation%EffectStress%RedKsSto
end function GetSimulation_EffectStress_RedKsSto


subroutine SetSimulation_EffectStress(EffectStress)
    !! Setter for the "EffectStress" attribute of the "simulation" global variable.
    type(rep_EffectStress), intent(in) :: EffectStress

    simulation%EffectStress = EffectStress
end subroutine SetSimulation_EffectStress


subroutine SetSimulation_EffectStress_RedCGC(RedCGC)
    !! Setter for the "RedCGC" attribute of the "EffectStress" attribute of the "simulation" global variable.
    integer(int8), intent(in) :: RedCGC

    simulation%EffectStress%RedCGC = RedCGC
end subroutine SetSimulation_EffectStress_RedCGC


subroutine SetSimulation_EffectStress_RedCCX(RedCCX)
    !! Setter for the "RedCCX" attribute of the "EffectStress" attribute of the "simulation" global variable.
    integer(int8), intent(in) :: RedCCX

    simulation%EffectStress%RedCCX = RedCCX
end subroutine SetSimulation_EffectStress_RedCCX


subroutine SetSimulation_EffectStress_RedWP(RedWP)
    !! Setter for the "RedWP" attribute of the "EffectStress" attribute of the "simulation" global variable.
    integer(int8), intent(in) :: RedWP

    simulation%EffectStress%RedWP = RedWP
end subroutine SetSimulation_EffectStress_RedWP


subroutine SetSimulation_EffectStress_CDecline(CDecline)
    !! Setter for the "CDecline" attribute of the "EffectStress" attribute of the "simulation" global variable.
    real(sp), intent(in) :: CDecline

    simulation%EffectStress%CDecline = CDecline
end subroutine SetSimulation_EffectStress_CDecline


subroutine SetSimulation_EffectStress_RedKsSto(RedKsSto)
    !! Setter for the "RedKsSto" attribute of the "EffectStress" attribute of the "simulation" global variable.
    integer(int8), intent(in) :: RedKsSto

    simulation%EffectStress%RedKsSto = RedKsSto
end subroutine SetSimulation_EffectStress_RedKsSto


function GetSimulation_Storage() result(Storage)
    !! Getter for the "Storage" attribute of the "simulation" global variable.
    type(rep_storage) :: Storage

    Storage = simulation%Storage
end function GetSimulation_Storage


function GetSimulation_Storage_Btotal() result(Btotal)
    !! Getter for the "Btotal" attribute of the "Storage" attribute of the "simulation" global variable.
    real(sp) :: Btotal

    Btotal = simulation%Storage%Btotal
end function GetSimulation_Storage_Btotal


function GetSimulation_Storage_CropString() result(CropString)
    !! Getter for the "CropString" attribute of the "Storage" attribute of the "simulation" global variable.
    character(len=len(simulation%storage%CropString)) :: CropString

    CropString = simulation%Storage%CropString
end function GetSimulation_Storage_CropString


function GetSimulation_Storage_Season() result(Season)
    !! Getter for the "Season" attribute of the "Storage" attribute of the "simulation" global variable.
    integer(int8) :: Season

    Season = simulation%Storage%Season
end function GetSimulation_Storage_Season


subroutine SetSimulation_Storage(Storage)
    !! Setter for the "Storage" attribute of the "simulation" global variable.
    type(rep_storage), intent(in) :: Storage

    simulation%Storage = Storage
end subroutine SetSimulation_Storage


subroutine SetSimulation_Storage_Btotal(Btotal)
    !! Setter for the "Btotal" attribute of the "Storage" attribute of the "simulation" global variable.
    real(sp), intent(in) :: Btotal

    simulation%Storage%Btotal = Btotal
end subroutine SetSimulation_Storage_Btotal


subroutine SetSimulation_Storage_CropString(CropString)
    !! Setter for the "CropString" attribute of the "Storage" attribute of the "simulation" global variable.
    character(len=*), intent(in) :: CropString

    simulation%Storage%CropString = CropString
end subroutine SetSimulation_Storage_CropString


subroutine SetSimulation_Storage_Season(Season)
    !! Setter for the "Season" attribute of the "Storage" attribute of the "simulation" global variable.
    integer(int8), intent(in) :: Season

    simulation%Storage%Season = Season
end subroutine SetSimulation_Storage_Season


function GetCompartment() result(Compartment_out)
    !! Getter for "Compartment" global variable.
    type(CompartmentIndividual), dimension(max_No_compartments) :: Compartment_out

    Compartment_out = Compartment
end function GetCompartment


function GetCompartment_i(i) result(Compartment_i)
    !! Getter for individual elements of "Compartment" global variable.
    integer(int32), intent(in) :: i
    type(CompartmentIndividual) :: Compartment_i

    Compartment_i = Compartment(i)
end function GetCompartment_i


subroutine SetCompartment_i(i, Compartment_i)
    !! Setter for individual elements of "Compartment" global variable.
    integer(int32), intent(in) :: i
    type(CompartmentIndividual) :: Compartment_i

    Compartment(i) = Compartment_i
end subroutine SetCompartment_i


function GetCompartment_Thickness(i) result(Thickness)
    !! Getter for the "Thickness" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: Thickness

    Thickness = compartment(i)%Thickness
end function GetCompartment_Thickness


function GetCompartment_theta(i) result(theta)
    !! Getter for the "theta" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: theta

    theta = compartment(i)%theta
end function GetCompartment_theta


function GetCompartment_fluxout(i) result(fluxout)
    !! Getter for the "fluxout" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: fluxout

    fluxout = compartment(i)%fluxout
end function GetCompartment_fluxout


function GetCompartment_Layer(i) result(Layer)
    !! Getter for the "Layer" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    integer(int32) :: Layer

    Layer = compartment(i)%Layer
end function GetCompartment_Layer


function GetCompartment_Smax(i) result(Smax)
    !! Getter for the "Smax" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: Smax

    Smax = compartment(i)%Smax
end function GetCompartment_Smax


function GetCompartment_FCadj(i) result(FCadj)
    !! Getter for the "FCadj" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: FCadj

    FCadj = compartment(i)%FCadj
end function GetCompartment_FCadj


function GetCompartment_DayAnaero(i) result(DayAnaero)
    !! Getter for the "DayAnaero" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    integer(int32) :: DayAnaero

    DayAnaero = compartment(i)%DayAnaero
end function GetCompartment_DayAnaero


function GetCompartment_WFactor(i) result(WFactor)
    !! Getter for the "WFactor" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: WFactor

    WFactor = compartment(i)%WFactor
end function GetCompartment_WFactor


function GetCompartment_Salt(i1, i2) result(Salt)
    !! Getter for individual elements of "Salt" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i1
    integer(int32), intent(in) :: i2
    real(sp) :: Salt

    Salt = compartment(i1)%Salt(i2)
end function GetCompartment_Salt


function GetCompartment_Depo(i1, i2) result(Depo)
    !! Getter for individual elements of "Depo" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i1
    integer(int32), intent(in) :: i2
    real(sp) :: Depo

    Depo = compartment(i1)%Depo(i2)
end function GetCompartment_Depo


subroutine SetCompartment(Compartment_in)
    !! Setter for the "compartment" global variable.
    type(CompartmentIndividual), dimension(max_No_compartments), intent(in) :: &
                                                                Compartment_in

    compartment = Compartment_in
end subroutine SetCompartment


subroutine SetCompartment_Thickness(i, Thickness)
    !! Setter for the "Thickness" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: Thickness

    compartment(i)%Thickness = Thickness
end subroutine SetCompartment_Thickness


subroutine SetCompartment_theta(i, theta)
    !! Setter for the "theta" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: theta

    compartment(i)%theta = theta
end subroutine SetCompartment_theta


subroutine SetCompartment_fluxout(i, fluxout)
    !! Setter for the "fluxout" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: fluxout

    compartment(i)%fluxout = fluxout
end subroutine SetCompartment_fluxout


subroutine SetCompartment_Layer(i, Layer)
    !! Setter for the "Layer" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: Layer

    compartment(i)%Layer = Layer
end subroutine SetCompartment_Layer


subroutine SetCompartment_Smax(i, Smax)
    !! Setter for the "Smax" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: Smax

    compartment(i)%Smax = Smax
end subroutine SetCompartment_Smax


subroutine SetCompartment_FCadj(i, FCadj)
    !! Setter for the "FCadj" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: FCadj

    compartment(i)%FCadj = FCadj
end subroutine SetCompartment_FCadj


subroutine SetCompartment_DayAnaero(i, DayAnaero)
    !! Setter for the "DayAnaero" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: DayAnaero

    compartment(i)%DayAnaero = DayAnaero
end subroutine SetCompartment_DayAnaero


subroutine SetCompartment_WFactor(i, WFactor)
    !! Setter for the "WFactor" attribute of the "compartment" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: WFactor

    compartment(i)%WFactor = WFactor
end subroutine SetCompartment_WFactor


subroutine SetCompartment_Salt(i1, i2, Salt)
    !! Setter for individual elements of "Salt" attribute of the "compartment global variable.
    integer(int32), intent(in) :: i1
    integer(int32), intent(in) :: i2
    real(sp), intent(in) :: Salt

    compartment(i1)%Salt(i2) = Salt
end subroutine SetCompartment_Salt


subroutine SetCompartment_Depo(i1, i2, Depo)
    !! Setter for individual elements of "Depo" attribute of the "compartment global variable.
    integer(int32), intent(in) :: i1
    integer(int32), intent(in) :: i2
    real(sp), intent(in) :: Depo

    compartment(i1)%Depo(i2) = Depo
end subroutine SetCompartment_Depo


function GetSoilLayer() result(SoilLayer_out)
    !! Getter for the "soillayer" global variable.
    type(SoilLayerIndividual), dimension(max_SoilLayers) :: SoilLayer_out

    SoilLayer_out = soillayer
end function GetSoilLayer


function GetSoilLayer_i(i) result(SoilLayer_i)
    !! Getter for the "soillayer" global variable (individual element)
    integer(int32), intent(in) :: i
    type(SoilLayerIndividual) :: SoilLayer_i

    SoilLayer_i = soillayer(i)
end function GetSoilLayer_i


function GetSoilLayer_Description(i) result(Description)
    !! Getter for the "Description" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    character(len=25) :: Description

    Description = soillayer(i)%Description
end function GetSoilLayer_Description


function GetSoilLayer_Thickness(i) result(Thickness)
    !! Getter for the "Thickness" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: Thickness

    Thickness = soillayer(i)%Thickness
end function GetSoilLayer_Thickness


function GetSoilLayer_SAT(i) result(SAT)
    !! Getter for the "SAT" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: SAT

    SAT = soillayer(i)%SAT
end function GetSoilLayer_SAT


function GetSoilLayer_FC(i) result(FC)
    !! Getter for the "FC" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: FC

    FC = soillayer(i)%FC
end function GetSoilLayer_FC


function GetSoilLayer_WP(i) result(WP)
    !! Getter for the "WP" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: WP

    WP = soillayer(i)%WP
end function GetSoilLayer_WP


function GetSoilLayer_tau(i) result(tau)
    !! Getter for the "tau" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: tau

    tau = soillayer(i)%tau
end function GetSoilLayer_tau


function GetSoilLayer_InfRate(i) result(InfRate)
    !! Getter for the "InfRate" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: InfRate

    InfRate = soillayer(i)%InfRate
end function GetSoilLayer_InfRate


function GetSoilLayer_Penetrability(i) result(Penetrability)
    !! Getter for the "Penetrability" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    integer(int8) :: Penetrability

    Penetrability = soillayer(i)%Penetrability
end function GetSoilLayer_Penetrability


function GetSoilLayer_GravelMass(i) result(GravelMass)
    !! Getter for the "GravelMass" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    integer(int8) :: GravelMass

    GravelMass = soillayer(i)%GravelMass
end function GetSoilLayer_GravelMass


function GetSoilLayer_GravelVol(i) result(GravelVol)
    !! Getter for the "GravelVol" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: GravelVol

    GravelVol = soillayer(i)%GravelVol
end function GetSoilLayer_GravelVol


function GetSoilLayer_WaterContent(i) result(WaterContent)
    !! Getter for the "WaterContent" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: WaterContent

    WaterContent = soillayer(i)%WaterContent
end function GetSoilLayer_WaterContent


function GetSoilLayer_Macro(i) result(Macro)
    !! Getter for the "Macro" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    integer(int8) :: Macro

    Macro = soillayer(i)%Macro
end function GetSoilLayer_Macro


function GetSoilLayer_SaltMobility(i) result(SaltMobility)
    !! Getter for the "SaltMobility" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), dimension(11) :: SaltMobility

    SaltMobility  = soillayer(i)%SaltMobility
end function GetSoilLayer_SaltMobility


function GetSoilLayer_SaltMobility_i(i1, i2) result(SaltMobility_i)
    !! Getter for the "SaltMobility" attribute of the "soillayer" global variable (individual element)
    integer(int32), intent(in) :: i1, i2
    real(sp) :: SaltMobility_i

    SaltMobility_i = soillayer(i1)%SaltMobility(i2)
end function GetSoilLayer_SaltMobility_i


function GetSoilLayer_SC(i) result(SC)
    !! Getter for the "SC" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    integer(int8) :: SC

    SC = soillayer(i)%SC
end function GetSoilLayer_SC


function GetSoilLayer_SCP1(i) result(SCP1)
    !! Getter for the "SCP1" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    integer(int8) :: SCP1

    SCP1 = soillayer(i)%SCP1
end function GetSoilLayer_SCP1


function GetSoilLayer_UL(i) result(UL)
    !! Getter for the "UL" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: UL

    UL = soillayer(i)%UL
end function GetSoilLayer_UL


function GetSoilLayer_Dx(i) result(Dx)
    !! Getter for the "Dx" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: Dx

    Dx = soillayer(i)%Dx
end function GetSoilLayer_Dx


function GetSoilLayer_SoilClass(i) result(SoilClass)
    !! Getter for the "SoilClass" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    integer(int8) :: SoilClass

    SoilClass = soillayer(i)%SoilClass
end function GetSoilLayer_SoilClass


function GetSoilLayer_CRa(i) result(CRa)
    !! Getter for the "CRa" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: CRa

    CRa = soillayer(i)%CRa
end function GetSoilLayer_CRa


function GetSoilLayer_CRb(i) result(CRb)
    !! Getter for the "CRb" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: CRb

    CRb = soillayer(i)%CRb
end function GetSoilLayer_CRb


subroutine SetSoilLayer(SoilLayer_in)
    !! Setter for the "soillayer" global variable.
    type(SoilLayerIndividual), dimension(max_SoilLayers), intent(in) :: SoilLayer_in

    soillayer = SoilLayer_in
end subroutine SetSoilLayer


subroutine SetSoilLayer_i(i, SoilLayer_i)
    !! Setter for the "soillayer" global variable (individual element)
    integer(int32), intent(in) :: i
    type(SoilLayerIndividual), intent(in) :: SoilLayer_i

    soillayer(i) = SoilLayer_i
end subroutine SetSoilLayer_i


subroutine SetSoilLayer_Description(i, Description)
    !! Setter for the "Description" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    character(len=25), intent(in) :: Description

    soillayer(i)%Description = Description
end subroutine SetSoilLayer_Description


subroutine SetSoilLayer_Thickness(i, Thickness)
    !! Setter for the "Thickness" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: Thickness

    soillayer(i)%Thickness = Thickness
end subroutine SetSoilLayer_Thickness


subroutine SetSoilLayer_SAT(i, SAT)
    !! Setter for the "SAT" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: SAT

    soillayer(i)%SAT = SAT
end subroutine SetSoilLayer_SAT


subroutine SetSoilLayer_FC(i, FC)
    !! Setter for the "FC" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: FC

    soillayer(i)%FC = FC
end subroutine SetSoilLayer_FC


subroutine SetSoilLayer_WP(i, WP)
    !! Setter for the "WP" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: WP

    soillayer(i)%WP = WP
end subroutine SetSoilLayer_WP


subroutine SetSoilLayer_tau(i, tau)
    !! Setter for the "tau" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: tau

    soillayer(i)%tau = tau
end subroutine SetSoilLayer_tau


subroutine SetSoilLayer_InfRate(i, InfRate)
    !! Setter for the "InfRate" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: InfRate

    soillayer(i)%InfRate = InfRate
end subroutine SetSoilLayer_InfRate


subroutine SetSoilLayer_Penetrability(i, Penetrability)
    !! Setter for the "Penetrability" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    integer(int8), intent(in) :: Penetrability

    soillayer(i)%Penetrability = Penetrability
end subroutine SetSoilLayer_Penetrability


subroutine SetSoilLayer_GravelMass(i, GravelMass)
    !! Setter for the "GravelMass" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    integer(int8), intent(in) :: GravelMass

    soillayer(i)%GravelMass = GravelMass
end subroutine SetSoilLayer_GravelMass


subroutine SetSoilLayer_GravelVol(i, GravelVol)
    !! Setter for the "GravelVol" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: GravelVol

    soillayer(i)%GravelVol = GravelVol
end subroutine SetSoilLayer_GravelVol


subroutine SetSoilLayer_WaterContent(i, WaterContent)
    !! Setter for the "WaterContent" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: WaterContent

    soillayer(i)%WaterContent = WaterContent
end subroutine SetSoilLayer_WaterContent


subroutine SetSoilLayer_Macro(i, Macro)
    !! Setter for the "Macro" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    integer(int8), intent(in) :: Macro

    soillayer(i)%Macro = Macro
end subroutine SetSoilLayer_Macro


subroutine SetSoilLayer_SaltMobility(i, SaltMobility)
    !! Setter for the "SaltMobility" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), dimension(11) :: SaltMobility

    soillayer(i)%SaltMobility = SaltMobility
end subroutine SetSoilLayer_SaltMobility


subroutine SetSoilLayer_SaltMobility_i(i1, i2, SaltMobility_i)
    !! Setter for the "SaltMobility" attribute of the "soillayer" global variable (individual element)
    integer(int32), intent(in) :: i1, i2
    real(sp) :: SaltMobility_i

    soillayer(i1)%SaltMobility(i2) = SaltMobility_i
end subroutine SetSoilLayer_SaltMobility_i


subroutine SetSoilLayer_SC(i, SC)
    !! Setter for the "SC" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    integer(int8), intent(in) :: SC

    soillayer(i)%SC = SC
end subroutine SetSoilLayer_SC


subroutine SetSoilLayer_SCP1(i, SCP1)
    !! Setter for the "SCP1" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    integer(int8), intent(in) :: SCP1

    soillayer(i)%SCP1 = SCP1
end subroutine SetSoilLayer_SCP1


subroutine SetSoilLayer_UL(i, UL)
    !! Setter for the "UL" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: UL

    soillayer(i)%UL = UL
end subroutine SetSoilLayer_UL


subroutine SetSoilLayer_Dx(i, Dx)
    !! Setter for the "Dx" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: Dx

    soillayer(i)%Dx = Dx
end subroutine SetSoilLayer_Dx


subroutine SetSoilLayer_SoilClass(i, SoilClass)
    !! Setter for the "SoilClass" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    integer(int8), intent(in) :: SoilClass

    soillayer(i)%SoilClass = SoilClass
end subroutine SetSoilLayer_SoilClass


subroutine SetSoilLayer_CRa(i, CRa)
    !! Setter for the "CRa" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: CRa

    soillayer(i)%CRa = CRa
end subroutine SetSoilLayer_CRa


subroutine SetSoilLayer_CRb(i, CRb)
    !! Setter for the "CRb" attribute of the "soillayer" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: CRb

    soillayer(i)%CRb = CRb
end subroutine SetSoilLayer_CRb


function GetManDescription() result(str)
    !! Getter for the "ManDescription" global variable.
    character(len=:), allocatable :: str

    str = ManDescription
end function GetManDescription


subroutine SetManDescription(str)
    !! Setter for the "ManDescription" global variable.
    character(len=*), intent(in) :: str

    ManDescription = str
end subroutine SetManDescription


integer(int32) function GetNrCompartments()
    !! Getter for the "NrCompartments" global variable.

    GetNrCompartments = NrCompartments
end function GetNrCompartments


subroutine SetNrCompartments(NrCompartments_in)
    !! Setter for the "NrCompartments" global variable.
    integer(int32), intent(in) :: NrCompartments_in

    NrCompartments = NrCompartments_in
end subroutine SetNrCompartments


integer(int32) function GetZiAqua()
    !! Getter for the "ZiAqua" global variable.

    GetZiAqua = ZiAqua
end function GetZiAqua


subroutine SetZiAqua(ZiAqua_in)
    !! Setter for the "ZiAqua" global variable.
    integer(int32), intent(in) :: ZiAqua_in

    ZiAqua = ZiAqua_in
end subroutine SetZiAqua


real(sp) function GetDrain()
    !! Getter for the "Drain" global variable.

    GetDrain = Drain
end function GetDrain


subroutine SetDrain(Drain_in)
    !! Setter for the "Drain" global variable.
    real(sp), intent(in) :: Drain_in

    Drain = Drain_in
end subroutine SetDrain


real(sp) function GetRain()
    !! Getter for the "Rain" global variable.

    GetRain = Rain
end function GetRain


subroutine SetRain(Rain_in)
    !! Setter for the "Rain" global variable.
    real(sp), intent(in) :: Rain_in

    Rain = Rain_in
end subroutine SetRain


real(sp) function GetRunoff()
    !! Getter for the "Runoff" global variable.

    GetRunoff = Runoff
end function GetRunoff


subroutine SetRunoff(Runoff_in)
    !! Setter for the "Runoff" global variable.
    real(sp), intent(in) :: Runoff_in

    Runoff = Runoff_in
end subroutine SetRunoff


real(sp) function GetSurfaceStorage()
    !! Getter for the "SurfaceStorage" global variable.
    !! mm/day
    GetSurfaceStorage = SurfaceStorage
end function GetSurfaceStorage


subroutine SetSurfaceStorage(SurfaceStorage_in)
    !! Setter for the "SurfaceStorage" global variable.
    real(sp), intent(in) :: SurfaceStorage_in

    SurfaceStorage = SurfaceStorage_in
end subroutine SetSurfaceStorage


real(sp) function GetECstorage()
    !! Getter for the "ECstorage" global variable.
    !! * EC surface storage dS/m *
    GetECstorage = ECstorage
end function GetECstorage


subroutine SetECstorage(ECstorage_in)
    !! Setter for the "ECstorage" global variable.
    real(sp), intent(in) :: ECstorage_in

    ECstorage = ECstorage_in
end subroutine SetECstorage


function GetOffSeasonDescription() result(str)
    !! Getter for the "OffSeasonDescription" global variable.
    character(len=:), allocatable :: str

    str = OffSeasonDescription
end function GetOffSeasonDescription


subroutine SetOffSeasonDescription(str)
    !! Setter for the "ManDescription" global variable.
    character(len=*), intent(in) :: str

    OffSeasonDescription = str
end subroutine SetOffSeasonDescription


function GetGroundwaterDescription() result(str)
    !! Getter for the "GroundwaterDescription" global variable.
    character(len=:), allocatable :: str

    str = GroundwaterDescription
end function GetGroundwaterDescription


subroutine SetGroundwaterDescription(str)
    !! Setter for the "ManDescription" global variable.
    character(len=*), intent(in) :: str

    GroundwaterDescription = str
end subroutine SetGroundwaterDescription


integer(int32) function GetDaySubmerged()
    !! Getter for the "DaySubmerged" global variable.

    GetDaySubmerged = DaySubmerged
end function GetDaySubmerged


subroutine SetDaySubmerged(DaySubmerged_in)
    !! Setter for the "DaySubmerged" global variable.
    integer(int32), intent(in) :: DaySubmerged_in

    DaySubmerged = DaySubmerged_in
end subroutine SetDaySubmerged


integer(int32) function GetTnxReferenceYear()
    !! Getter for the "TnxReferenceYear" global variable.

    GetTnxReferenceYear = TnxReferenceYear
end function GetTnxReferenceYear


subroutine SetTnxReferenceYear(TnxReferenceYear_in)
    !! Setter for the "TnxReferenceYear" global variable.
    integer(int32), intent(in) :: TnxReferenceYear_in

    TnxReferenceYear = TnxReferenceYear_in
end subroutine SetTnxReferenceYear


real(sp) function GetCRsalt()
    !! Getter for the "CRsalt" global variable.

    GetCRsalt = CRsalt
end function GetCRsalt


subroutine SetCRsalt(CRsalt_in)
    !! Setter for the "CRsalt" global variable.
    real(sp), intent(in) :: CRsalt_in

    CRsalt = CRsalt_in
end subroutine SetCRsalt


integer(int32) function GetMaxPlotNew()
    !! Getter for the "MaxPlotNew" global variable.

    GetMaxPlotNew = MaxPlotNew
end function GetMaxPlotNew


subroutine SetMaxPlotNew(MaxPlotNew_in)
    !! Setter for the "MaxPlotNew" global variable.
    integer(int32), intent(in) :: MaxPlotNew_in

    MaxPlotNew = MaxPlotNew_in
end subroutine SetMaxPlotNew


integer(int8) function GetMaxPlotTr()
    !! Getter for the "MaxPlotTr" global variable.

    GetMaxPlotTr = MaxPlotTr
end function GetMaxPlotTr


subroutine SetMaxPlotTr(MaxPlotTr_in)
    !! Setter for the "MaxPlotTr" global variable.
    integer(int8), intent(in) :: MaxPlotTr_in

    MaxPlotTr = MaxPlotTr_in
end subroutine SetMaxPlotTr


integer(int8) function GetOutputAggregate()
    !! Getter for the "OutputAggregate" global variable.

    GetOutputAggregate = OutputAggregate
end function GetOutputAggregate


subroutine SetOutputAggregate(OutputAggregate_in)
    !! Setter for the "OutputAggregate" global variable.
    integer(int8), intent(in) :: OutputAggregate_in

    OutputAggregate = OutputAggregate_in
end subroutine SetOutputAggregate


logical function GetEvapoEntireSoilSurface()
    !! Getter for the "EvapoEntireSoilSurface" global variable.

    GetEvapoEntireSoilSurface = EvapoEntireSoilSurface
end function GetEvapoEntireSoilSurface


subroutine SetEvapoEntireSoilSurface(EvapoEntireSoilSurface_in)
    !! Setter for the "EvapoEntireSoilSurface" global variable.
    logical, intent(in) :: EvapoEntireSoilSurface_in

    EvapoEntireSoilSurface = EvapoEntireSoilSurface_in
end subroutine SetEvapoEntireSoilSurface


logical function GetPreDay()
    !! Getter for the "PreDay" global variable.

    GetPreDay = PreDay
end function GetPreDay


subroutine SetPreDay(PreDay_in)
    !! Setter for the "PreDay" global variable.
    logical, intent(in) :: PreDay_in

    PreDay = PreDay_in
end subroutine SetPreDay


logical function GetOut1Wabal()
    !! Getter for the "Out1Wabal" global variable.

    GetOut1Wabal = Out1Wabal
end function GetOut1Wabal


subroutine SetOut1Wabal(Out1Wabal_in)
    !! Setter for the "Out1Wabal" global variable.
    logical, intent(in) :: Out1Wabal_in

    Out1Wabal = Out1Wabal_in
end subroutine SetOut1Wabal


logical function GetOut2Crop()
    !! Getter for the "Out2Crop" global variable.

    GetOut2Crop = Out2Crop
end function GetOut2Crop


subroutine SetOut2Crop(Out2Crop_in)
    !! Setter for the "Out2Crop" global variable.
    logical, intent(in) :: Out2Crop_in

    Out2Crop = Out2Crop_in
end subroutine SetOut2Crop


logical function GetOut3Prof()
    !! Getter for the "Out3Prof" global variable.

    GetOut3Prof= Out3Prof
end function GetOut3Prof


subroutine SetOut3Prof(Out3Prof_in)
    !! Setter for the "Out3Prof" global variable.
    logical, intent(in) :: Out3Prof_in

    Out3Prof = Out3Prof_in
end subroutine SetOut3Prof


logical function GetOut4Salt()
    !! Getter for the "Out4Salt" global variable.

    GetOut4Salt= Out4Salt
end function GetOut4Salt


subroutine SetOut4Salt(Out4Salt_in)
    !! Setter for the "Out4Salt" global variable.
    logical, intent(in) :: Out4Salt_in

    Out4Salt = Out4Salt_in
end subroutine SetOut4Salt


logical function GetOut5CompWC()
    !! Getter for the "Out5CompWC" global variable.

    GetOut5CompWC= Out5CompWC
end function GetOut5CompWC


subroutine SetOut5CompWC(Out5CompWC_in)
    !! Setter for the "Out5CompWC" global variable.
    logical, intent(in) :: Out5CompWC_in

    Out5CompWC = Out5CompWC_in
end subroutine SetOut5CompWC


logical function GetOut6CompEC()
    !! Getter for the "Out6CompEC" global variable.

    GetOut6CompEC= Out6CompEC
end function GetOut6CompEC


subroutine SetOut6CompEC(Out6CompEC_in)
    !! Setter for the "Out6CompEC" global variable.
    logical, intent(in) :: Out6CompEC_in

    Out6CompEC = Out6CompEC_in
end subroutine SetOut6CompEC


logical function GetOut7Clim()
    !! Getter for the "Out7Clim" global variable.

    GetOut7Clim= Out7Clim
end function GetOut7Clim


subroutine SetOut7Clim(Out7Clim_in)
    !! Setter for the "Out7Clim" global variable.
    logical, intent(in) :: Out7Clim_in

    Out7Clim = Out7Clim_in
end subroutine SetOut7Clim


integer(int8) function GetIniPercTAW()
    !! Getter for the "IniPercTAW" global variable.

    GetIniPercTAW = IniPercTAW
end function GetIniPercTAW


subroutine SetIniPercTAW(IniPercTAW_in)
    !! Setter for the "IniPercTAW" global variable.
    integer(int8), intent(in) :: IniPercTAW_in

    IniPercTAW = IniPercTAW_in
end subroutine SetIniPercTAW


real(sp) function GetECiAqua()
    !! Getter for the "ECiAqua" global variable.

    GetECiAqua = ECiAqua
end function GetECiAqua


subroutine SetECiAqua(ECiAqua_in)
    !! Setter for the "ECiAqua" global variable.
    real(sp), intent(in) :: ECiAqua_in

    ECiAqua = ECiAqua_in
end subroutine SetECiAqua


real(sp) function GetEact()
    !! Getter for the "Eact" global variable.

    GetEact = Eact
end function GetEact


subroutine SetEact(Eact_in)
    !! Setter for the "Eact" global variable.
    real(sp), intent(in) :: Eact_in

    Eact = Eact_in
end subroutine SetEact


real(sp) function GetETo()
    !! Getter for the "ETo" global variable.

    GetETo = ETo
end function GetETo


subroutine SetETo(ETo_in)
    !! Setter for the "ETo" global variable.
    real(sp), intent(in) :: ETo_in

    ETo = ETo_in
end subroutine SetETo


real(sp) function GetIrrigation()
    !! Getter for the "Irrigation" global variable.

    GetIrrigation = Irrigation
end function GetIrrigation


subroutine SetIrrigation(Irrigation_in)
    !! Setter for the "Irrigation" global variable.
    real(sp), intent(in) :: Irrigation_in

    Irrigation = Irrigation_in
end subroutine SetIrrigation


real(sp) function GetInfiltrated()
    !! Getter for the "Infiltrated" global variable.

    GetInfiltrated = Infiltrated
end function GetInfiltrated


subroutine SetInfiltrated(Infiltrated_in)
    !! Setter for the "Infiltrated" global variable.
    real(sp), intent(in) :: Infiltrated_in

    Infiltrated = Infiltrated_in
end subroutine SetInfiltrated


real(sp) function GetCRwater()
    !! Getter for the "CRwater" global variable.

    GetCRwater = CRwater
end function GetCRwater


subroutine SetCRwater(CRwater_in)
    !! Setter for the "CRwater" global variable.
    real(sp), intent(in) :: CRwater_in

    CRwater = CRwater_in
end subroutine SetCRwater


real(sp) function GetCCiActual()
    !! Getter for the "CCiActual" global variable.

    GetCCiActual = CCiActual
end function GetCCiActual


subroutine SetCCiActual(CCiActual_in)
    !! Setter for the "CCiActual" global variable.
    real(sp), intent(in) :: CCiActual_in

    CCiActual = CCiActual_in
end subroutine SetCCiActual


real(sp) function GetEpot()
    !! Getter for the "Epot" global variable.

    GetEpot = Epot
end function GetEpot


subroutine SetEpot(Epot_in)
    !! Setter for the "Epot" global variable.
    real(sp), intent(in) :: Epot_in

    Epot = Epot_in
end subroutine SetEpot


real(sp) function GetTpot()
    !! Getter for the "Tpot" global variable.

    GetTpot = Tpot
end function GetTpot


subroutine SetTpot(Tpot_in)
    !! Setter for the "Tpot" global variable.
    real(sp), intent(in) :: Tpot_in

    Tpot = Tpot_in
end subroutine SetTpot


function GetMultipleProjectDescription() result(str)
    !! Getter for the "MultipleProjectDescription" global variable.
    character(len=:), allocatable :: str

    str = MultipleProjectDescription
end function GetMultipleProjectDescription


subroutine SetMultipleProjectDescription(str)
    !! Setter for the "MultipleProjectDescription" global variable.
    character(len=*), intent(in) :: str

    MultipleProjectDescription = str
end subroutine SetMultipleProjectDescription


function GetProjectDescription() result(str)
    !! Getter for the "ProjectDescription" global variable.
    character(len=:), allocatable :: str

    str = ProjectDescription
end function GetProjectDescription


subroutine SetProjectDescription(str)
    !! Setter for the "ProjectDescription" global variable.
    character(len=*), intent(in) :: str

    ProjectDescription = str
end subroutine SetProjectDescription


real(sp) function GetRootingDepth()
    !! Getter for the "RootingDepth" global variable.

    GetRootingDepth = RootingDepth
end function GetRootingDepth


subroutine SetRootingDepth(RootingDepth_in)
    !! Setter for the "RootingDepth" global variable.
    real(sp), intent(in) :: RootingDepth_in

    RootingDepth = RootingDepth_in
end subroutine SetRootingDepth


real(sp) function GetCCiPrev()
    !! Getter for the "CCiPrev" global variable.

    GetCCiPrev = CCiPrev
end function GetCCiPrev


subroutine SetCCiPrev(CCiPrev_in)
    !! Setter for the "CCiPrev" global variable.
    real(sp), intent(in) :: CCiPrev_in

    CCiPrev = CCiPrev_in
end subroutine SetCCiPrev


real(sp) function GetCCiTopEarlySen()
    !! Getter for the "CCiTopEarlySen" global variable.

    GetCCiTopEarlySen = CCiTopEarlySen
end function GetCCiTopEarlySen


subroutine SetCCiTopEarlySen(CCiTopEarlySen_in)
    !! Setter for the "CCiTopEarlySen" global variable.
    real(sp), intent(in) :: CCiTopEarlySen_in

    CCiTopEarlySen = CCiTopEarlySen_in
end subroutine SetCCiTopEarlySen


real(sp) function GetECDrain()
    !! Getter for the "ECDrain" global variable.

    GetECDrain = ECDrain
end function GetECDrain


subroutine SetECDrain(ECDrain_in)
    !! Setter for the "ECDrain" global variable.
    real(sp), intent(in) :: ECDrain_in

    ECDrain = ECDrain_in
end subroutine SetECDrain


real(sp) function GetSaltInfiltr()
    !! Getter for the "SaltInfiltr" global variable.

    GetSaltInfiltr = SaltInfiltr
end function GetSaltInfiltr


subroutine SetSaltInfiltr(SaltInfiltr_in)
    !! Setter for the "SaltInfiltr" global variable.
    real(sp), intent(in) :: SaltInfiltr_in

    SaltInfiltr = SaltInfiltr_in
end subroutine SetSaltInfiltr


real(sp) function GetTact()
    !! Getter for the "Tact" global variable.

    GetTact = Tact
end function GetTact


subroutine SetTact(Tact_in)
    !! Setter for the "Tact" global variable.
    real(sp), intent(in) :: Tact_in

    Tact = Tact_in
end subroutine SetTact


real(sp) function GetTactWeedInfested()
    !! Getter for the "TactWeedInfested" global variable.

    GetTactWeedInfested = TactWeedInfested
end function GetTactWeedInfested


subroutine SetTactWeedInfested(TactWeedInfested_in)
    !! Setter for the "TactWeedInfested" global variable.
    real(sp), intent(in) :: TactWeedInfested_in

    TactWeedInfested = TactWeedInfested_in
end subroutine SetTactWeedInfested

! TminTnxReference365Days

function GetTminTnxReference365DaysRun() result(TminTnxReference365DaysRun_out)
    !! Getter for the "TminTnxReference365DaysRun" global variable.
    real(sp), dimension(1:365) :: TminTnxReference365DaysRun_out

    TminTnxReference365DaysRun_out = TminTnxReference365DaysRun
end function GetTminTnxReference365DaysRun


function GetTminTnxReference365DaysRun_i(i) result(TminTnxReference365DaysRun_i)
    !! Getter for individual elements of the "GetTminTnxReference365DaysRun" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: TminTnxReference365DaysRun_i

    TminTnxReference365DaysRun_i = TminTnxReference365DaysRun(i)
end function GetTminTnxReference365DaysRun_i


subroutine SetTminTnxReference365DaysRun(TminTnxReference365DaysRun_in)
    !! Setter for the "TminTnxReference365DaysRun" global variable.
    real(sp), dimension(1:365), intent(in) :: TminTnxReference365DaysRun_in

    TminTnxReference365DaysRun = TminTnxReference365DaysRun_in
end subroutine SetTminTnxReference365DaysRun


subroutine SetTminTnxReference365DaysRun_i(i, TminTnxReference365DaysRun_i)
    !! Setter for individual element for the "TminTnxReference365DaysRun" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: TminTnxReference365DaysRun_i

    TminTnxReference365DaysRun(i) = TminTnxReference365DaysRun_i
end subroutine SetTminTnxReference365DaysRun_i

! TmaxTnxReference365Days

function GetTmaxTnxReference365DaysRun() result(TmaxTnxReference365DaysRun_out)
    !! Getter for the "TmaxTnxReference365DaysRun" global variable.
    real(sp), dimension(1:365) :: TmaxTnxReference365DaysRun_out

    TmaxTnxReference365DaysRun_out = TmaxTnxReference365DaysRun
end function GetTmaxTnxReference365DaysRun


function GetTmaxTnxReference365DaysRun_i(i) result(TmaxTnxReference365DaysRun_i)
    !! Getter for individual elements of the "GetTmaxTnxReference365DaysRun" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: TmaxTnxReference365DaysRun_i

    TmaxTnxReference365DaysRun_i = TmaxTnxReference365DaysRun(i)
end function GetTmaxTnxReference365DaysRun_i


subroutine SetTmaxTnxReference365DaysRun(TmaxTnxReference365DaysRun_in)
    !! Setter for the "TmaxTnxReference365DaysRun" global variable.
    real(sp), dimension(1:365), intent(in) :: TmaxTnxReference365DaysRun_in

    TmaxTnxReference365DaysRun = TmaxTnxReference365DaysRun_in
end subroutine SetTmaxTnxReference365DaysRun


subroutine SetTmaxTnxReference365DaysRun_i(i, TmaxTnxReference365DaysRun_i)
    !! Setter for individual element for the "TmaxTnxReference365DaysRun" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: TmaxTnxReference365DaysRun_i

    TmaxTnxReference365DaysRun(i) = TmaxTnxReference365DaysRun_i
end subroutine SetTmaxTnxReference365DaysRun_i


real(sp) function GetTminTnxReference365Days()
    !! Getter for the "TminTnxReference365Days" global variable.

    GetTminTnxReference365Days = TminTnxReference365Days
end function GetTminTnxReference365Days


subroutine SetTminTnxReference365Days(TminTnxReference365Days_in)
    !! Setter for the "TminTnxReference365Days" global variable.
    real(sp), intent(in) :: TminTnxReference365Days_in

    TminTnxReference365Days = TminTnxReference365Days_in
end subroutine SetTminTnxReference365Days


real(sp) function GetTmaxTnxReference365Days()
    !! Getter for the "TmaxTnxReference365Days" global variable.

    GetTmaxTnxReference365Days = TmaxTnxReference365Days
end function GetTmaxTnxReference365Days


subroutine SetTmaxTnxReference365Days(TmaxTnxReference365Days_in)
    !! Setter for the "TmaxTnxReference365Days" global variable.
    real(sp), intent(in) :: TmaxTnxReference365Days_in

    TmaxTnxReference365Days = TmaxTnxReference365Days_in
end subroutine SetTmaxTnxReference365Days

! TminCropReferenceRun

function GetTminCropReferenceRun() result(TminCropReferenceRun_out)
    !! Getter for the "TminCropReferenceRun" global variable.
    real(sp), dimension(1:365) :: TminCropReferenceRun_out

    TminCropReferenceRun_out = TminCropReferenceRun
end function GetTminCropReferenceRun


function GetTminCropReferenceRun_i(i) result(TminCropReferenceRun_i)
    !! Getter for individual elements of the "GetTminCropReferenceRun" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: TminCropReferenceRun_i

    TminCropReferenceRun_i = TminCropReferenceRun(i)
end function GetTminCropReferenceRun_i


subroutine SetTminCropReferenceRun(TminCropReferenceRun_in)
    !! Setter for the "TminCropReferenceRun" global variable.
    real(sp), dimension(1:365), intent(in) :: TminCropReferenceRun_in

    TminCropReferenceRun = TminCropReferenceRun_in
end subroutine SetTminCropReferenceRun


subroutine SetTminCropReferenceRun_i(i, TminCropReferenceRun_i)
    !! Setter for individual element for the "TminCropReferenceRun" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: TminCropReferenceRun_i

    TminCropReferenceRun(i) = TminCropReferenceRun_i
end subroutine SetTminCropReferenceRun_i

! TmaxCropReferenceRun

function GetTmaxCropReferenceRun() result(TmaxCropReferenceRun_out)
    !! Getter for the "TmaxCropReferenceRun" global variable.
    real(sp), dimension(1:365) :: TmaxCropReferenceRun_out

    TmaxCropReferenceRun_out = TmaxCropReferenceRun
end function GetTmaxCropReferenceRun


function GetTmaxCropReferenceRun_i(i) result(TmaxCropReferenceRun_i)
    !! Getter for individual elements of the "GetTmaxCropReferenceRun" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: TmaxCropReferenceRun_i

    TmaxCropReferenceRun_i = TmaxCropReferenceRun(i)
end function GetTmaxCropReferenceRun_i


subroutine SetTmaxCropReferenceRun(TmaxCropReferenceRun_in)
    !! Setter for the "TmaxCropReferenceRun" global variable.
    real(sp), dimension(1:365), intent(in) :: TmaxCropReferenceRun_in

    TmaxCropReferenceRun = TmaxCropReferenceRun_in
end subroutine SetTmaxCropReferenceRun


subroutine SetTmaxCropReferenceRun_i(i, TmaxCropReferenceRun_i)
    !! Setter for individual element for the "TmaxCropReferenceRun" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: TmaxCropReferenceRun_i

    TmaxCropReferenceRun(i) = TmaxCropReferenceRun_i
end subroutine SetTmaxCropReferenceRun_i


real(sp) function GetTminCropReference()
    !! Getter for the "TminCropReference" global variable.

    GetTminCropReference = TminCropReference
end function GetTminCropReference


subroutine SetTminCropReference(TminCropReference_in)
    !! Setter for the "TminCropReference" global variable.
    real(sp), intent(in) :: TminCropReference_in

    TminCropReference = TminCropReference_in
end subroutine SetTminCropReference


real(sp) function GetTmaxCropReference()
    !! Getter for the "TmaxCropReference" global variable.

    GetTmaxCropReference = TmaxCropReference
end function GetTmaxCropReference


subroutine SetTmaxCropReference(TmaxCropReference_in)
    !! Setter for the "TmaxCropReference" global variable.
    real(sp), intent(in) :: TmaxCropReference_in

    TmaxCropReference = TmaxCropReference_in
end subroutine SetTmaxCropReference

! TminRun

function GetTminRun() result(TminRun_out)
    !! Getter for the "TminRun" global variable.
    real(sp), dimension(1:366) :: TminRun_out

    TminRun_out = TminRun
end function GetTminRun


function GetTminRun_i(i) result(TminRun_i)
    !! Getter for individual elements of the "GetTminRun" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: TminRun_i

    TminRun_i = TminRun(i)
    !TminRun_i = real(roundc(10000*real(TminRun(i),kind=sp),mold=int32)/10000._sp,kind=sp)
end function GetTminRun_i


subroutine SetTminRun(TminRun_in)
    !! Setter for the "TminRun" global variable.
    real(sp), dimension(1:366), intent(in) :: TminRun_in

    TminRun = TminRun_in
end subroutine SetTminRun


subroutine SetTminRun_i(i, TminRun_i)
    !! Setter for individual element for the "TminRun" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: TminRun_i

    TminRun(i) = TminRun_i
end subroutine SetTminRun_i

! TmaxRun

function GetTmaxRun() result(TmaxRun_out)
    !! Getter for the "TmaxRun" global variable.
    real(sp), dimension(1:366) :: TmaxRun_out

    TmaxRun_out = TmaxRun
end function GetTmaxRun


function GetTmaxRun_i(i) result(TmaxRun_i)
    !! Getter for individual elements of the "GetTmaxRun" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: TmaxRun_i

    TmaxRun_i = TmaxRun(i)
end function GetTmaxRun_i


subroutine SetTmaxRun(TmaxRun_in)
    !! Setter for the "TmaxRun" global variable.
    real(sp), dimension(1:366), intent(in) :: TmaxRun_in

    TmaxRun = TmaxRun_in
end subroutine SetTmaxRun


subroutine SetTmaxRun_i(i, TmaxRun_i)
    !! Setter for individual element for the "TmaxRun" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: TmaxRun_i

    TmaxRun(i) = TmaxRun_i
end subroutine SetTmaxRun_i


real(sp) function GetTmin()
    !! Getter for the "Tmin" global variable.

    GetTmin = Tmin
end function GetTmin


subroutine SetTmin(Tmin_in)
    !! Setter for the "Tmin" global variable.
    real(sp), intent(in) :: Tmin_in

    Tmin = Tmin_in
end subroutine SetTmin


real(sp) function GetTmax()
    !! Getter for the "Tmax" global variable.

    GetTmax = Tmax
end function GetTmax


subroutine SetTmax(Tmax_in)
    !! Setter for the "Tmax" global variable.
    real(sp), intent(in) :: Tmax_in

    Tmax = Tmax_in
end subroutine SetTmax


! TminTnxReference12MonthsRun

function GetTminTnxReference12MonthsRun() result(TminTnxReference12MonthsRun_out)
    !! Getter for the "TminTnxReference12MonthsRun" global variable.
    real(sp), dimension(1:12) :: TminTnxReference12MonthsRun_out

    TminTnxReference12MonthsRun_out = TminTnxReference12MonthsRun
end function GetTminTnxReference12MonthsRun


function GetTminTnxReference12MonthsRun_i(i) result(TminTnxReference12MonthsRun_i)
    !! Getter for individual elements of the "GetTminTnxReference12MonthsRun" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: TminTnxReference12MonthsRun_i

    TminTnxReference12MonthsRun_i = TminTnxReference12MonthsRun(i)
end function GetTminTnxReference12MonthsRun_i


subroutine SetTminTnxReference12MonthsRun(TminTnxReference12MonthsRun_in)
    !! Setter for the "TminTnxReference12MonthsRun" global variable.
    real(sp), dimension(1:12), intent(in) :: TminTnxReference12MonthsRun_in

    TminTnxReference12MonthsRun = TminTnxReference12MonthsRun_in
end subroutine SetTminTnxReference12MonthsRun


subroutine SetTminTnxReference12MonthsRun_i(i, TminTnxReference12MonthsRun_i)
    !! Setter for individual element for the "TminTnxReference12MonthsRun" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: TminTnxReference12MonthsRun_i

    TminTnxReference12MonthsRun(i) = TminTnxReference12MonthsRun_i
end subroutine SetTminTnxReference12MonthsRun_i

! TmaxTnxReference12MonthsRun

function GetTmaxTnxReference12MonthsRun() result(TmaxTnxReference12MonthsRun_out)
    !! Getter for the "TmaxTnxReference12MonthsRun" global variable.
    real(sp), dimension(1:12) :: TmaxTnxReference12MonthsRun_out

    TmaxTnxReference12MonthsRun_out = TmaxTnxReference12MonthsRun
end function GetTmaxTnxReference12MonthsRun


function GetTmaxTnxReference12MonthsRun_i(i) result(TmaxTnxReference12MonthsRun_i)
    !! Getter for individual elements of the "GetTmaxTnxReference12MonthsRun" global variable.
    integer(int32), intent(in) :: i
    real(sp) :: TmaxTnxReference12MonthsRun_i

    TmaxTnxReference12MonthsRun_i = TmaxTnxReference12MonthsRun(i)
end function GetTmaxTnxReference12MonthsRun_i


subroutine SetTmaxTnxReference12MonthsRun(TmaxTnxReference12MonthsRun_in)
    !! Setter for the "TmaxTnxReference12MonthsRun" global variable.
    real(sp), dimension(1:12), intent(in) :: TmaxTnxReference12MonthsRun_in

    TmaxTnxReference12MonthsRun = TmaxTnxReference12MonthsRun_in
end subroutine SetTmaxTnxReference12MonthsRun


subroutine SetTmaxTnxReference12MonthsRun_i(i, TmaxTnxReference12MonthsRun_i)
    !! Setter for individual element for the "TmaxTnxReference12MonthsRun" global variable.
    integer(int32), intent(in) :: i
    real(sp), intent(in) :: TmaxTnxReference12MonthsRun_i

    TmaxTnxReference12MonthsRun(i) = TmaxTnxReference12MonthsRun_i
end subroutine SetTmaxTnxReference12MonthsRun_i


! Surf0

real(sp) function GetSurf0()
    !! Getter for the "Surf0" global variable.

    GetSurf0 = Surf0
end function GetSurf0


subroutine SetSurf0(Surf0_in)
    !! Setter for the "Surf0" global variable.
    real(sp), intent(in) :: Surf0_in

    Surf0 = Surf0_in
end subroutine SetSurf0


logical function GetOutDaily()
    !! Getter for the OutDaily global variable

    GetOutDaily = OutDaily
end function GetOutDaily


subroutine SetOutDaily(OutDaily_in)
    !! Setter for the OutDaily global variable
    logical, intent(in) :: OutDaily_in

    OutDaily = OutDaily_in
end subroutine SetOutDaily


logical function GetOut8Irri()
    !! Getter for the Out8Irri global variable

    GetOut8Irri = Out8Irri
end function GetOut8Irri


subroutine SetOut8Irri(Out8Irri_in)
    !! Setter for the Out8Irri global variable
    logical, intent(in) :: Out8Irri_in

    Out8Irri = Out8Irri_in
end subroutine SetOut8Irri


function GetPathNameParam() result(str)
    !! Getter for the "PathNameParam" global variable.
    character(len=:), allocatable :: str

    str = PathNameParam
end function GetPathNameParam


subroutine SetPathNameParam(str)
    !! Setter for the "PathNameParam" global variable.
    character(len=*), intent(in) :: str

    PathNameParam = str
end subroutine SetPathNameParam


function GetPathNameList() result(str)
    !! Getter for the "PathNameList" global variable.
    character(len=:), allocatable :: str

    str = PathNameList
end function GetPathNameList


subroutine SetPathNameList(str)
    !! Setter for the "PathNameList" global variable.
    character(len=*), intent(in) :: str

    PathNameList = str
end subroutine SetPathNameList


logical function GetPart1Mult()
    !! Getter for the Part1Mult global variable

    GetPart1Mult = Part1Mult
end function GetPart1Mult


subroutine SetPart1Mult(Part1Mult_in)
    !! Setter for the Part1Mult global variable
    logical, intent(in) :: Part1Mult_in

    Part1Mult = Part1Mult_in
end subroutine SetPart1Mult


logical function GetPart2Eval()
    !! Getter for the Part2Eval global variable

    GetPart2Eval = Part2Eval
end function GetPart2Eval


subroutine SetPart2Eval(Part2Eval_in)
    !! Setter for the Part2Eval global variable
    logical, intent(in) :: Part2Eval_in

    Part2Eval = Part2Eval_in
end subroutine SetPart2Eval

end module ac_global
