!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0     
!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
module NoahMP401_module
!BOP
!
! !MODULE: NoahMP401_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the NoahMP401 1-d variables.
!  The variables specified in the data structure include:
!
!  \begin{description}      
!   \item[n]
!     nest id. unit: -
!   \item[latitude]
!     latitude in decimal degree. unit: rad
!   \item[logitude]
!     longitude in decimal year. unit: rad
!   \item[year]
!     year of the current time step. unit: -
!   \item[month]
!     month of the current time step. unit: -
!   \item[day]
!     day of the current time step. unit: -
!   \item[hour]
!     hour of the current time step. unit: -
!   \item[minute]
!     minute of the current time step. unit: -
!   \item[dz8w]
!     thickness of atmospheric layers. unit: m
!   \item[dt]
!     timestep. unit: s
!   \item[sldpth]
!     thickness of soil layers. unit: m
!   \item[nsoil]
!     number of soil layers. unit: -
!   \item[nsnow]
!     maximum number of snow layers (e.g. 3). unit: -
!   \item[vegetype]
!     vegetation type. unit: -
!   \item[soiltype]
!     soil type. unit: -
!   \item[shdfac\_monthly]
!     monthly values for green vegetation fraction. unit: 
!   \item[tbot]
!     deep soil temperature. unit: K
!   \item[urban\_vegetype]
!     urban land cover type index. unit: -
!   \item[cropcat]
!     crop category. unit: -
!   \item[planting]
!     planting date. unit: -
!   \item[harvest]
!     harvest date. unit: -
!   \item[season\_gdd]
!     growing season GDD. unit: -
!   \item[landuse\_tbl\_name]
!     Noah model landuse parameter table. unit: -
!   \item[soil\_tbl\_name]
!     Noah model soil parameter table. unit: -
!   \item[gen\_tbl\_name]
!     Noah model general parameter table. unit: -
!   \item[noahmp\_tbl\_name]
!     NoahMP parameter table. unit: -
!   \item[landuse\_scheme\_name]
!     Landuse classification scheme. unit: -
!   \item[soil\_scheme\_name]
!     Soil classification scheme. unit: -
!   \item[dveg\_opt]
!     dynamic vegetation, (1-$>$off; 2-$>$on); with opt\_crs=1. unit: -
!   \item[crs\_opt]
!     canopt stomatal resistance (1-$>$Ball-Berry; 2-$>$Jarvis). unit: -
!   \item[btr\_opt]
!     soil moisture factor for stomatal resistance (1-$>$Noah;2-$>$CLM;3-$>$SSiB). unit: -
!   \item[run\_opt]
!     runoff and groundwater (1-$>$SIMGM; 2-$>$SIMTOP; 3-$>$Schaake96; 4-$>$BATS). unit: -
!   \item[sfc\_opt]
!     surface layer drag coeff (CH \& CM) (1-$>$M-O; 2-$>$Chen97). unit: -
!   \item[frz\_opt]
!     supercooled liquid water (1-$>$NY06; 2-$>$Koren99). unit: -
!   \item[inf\_opt]
!     frozen soil permeability (1-$>$NY06; 2-$>$Koren99). unit: -
!   \item[rad\_opt]
!     radiation transfer (1-$>$gap=F(3D,cosz); 2-$>$gap=0; 3-$>$gap=1-Fveg). unit: -
!   \item[alb\_opt]
!     snow surface albedo (1-$>$BATS; 2-$>$CLASS). unit: -
!   \item[snf\_opt]
!     rainfall \& snowfall (1-$>$Jordan91; 2-$>$BATS; 3-$>$Noah). unit: -
!   \item[tbot\_opt]
!     lower boundary of soil temperature. unit: -
!   \item[stc\_opt]
!     snow/soil temperature time scheme. unit: -
!   \item[gla\_opt]
!     glacier option (1-$>$phase change; 2-$>$simple). unit: -
!   \item[rsf\_opt]
!     surface resistance (1-$>$Sakaguchi/Zeng;2-$>$Seller;3-$>$mod Sellers;4-$>$1+snow). unit: -
!   \item[soil\_opt]
!     soil configuration option. unit: -
!   \item[pedo\_opt]
!     soil pedotransfer function option. unit: -
!   \item[crop\_opt]
!     crop model option (0-$>$none; 1-$>$Liu et al.; 2-$>$Gecros). unit: -
!   \item[urban\_opt]
!     urban physics option. unit: -
!   \item[soilcomp]
!     soil sand and clay percentage. unit: -
!   \item[soilcL1]
!     soil texture in layer 1. unit: -
!   \item[soilcL2]
!     soil texture in layer 2. unit: -
!   \item[soilcL3]
!     soil texture in layer 3. unit: -
!   \item[soilcL4]
!     soil texture in layer 4. unit: -
!   \item[tair]
!     air temperature. unit: K
!   \item[psurf]
!     air pressure. unit: Pa
!   \item[wind\_e]
!     U wind component. unit: m/s
!   \item[wind\_n]
!     V wind component. unit: m/s
!   \item[qair]
!     specific humidity. unit: kg/kg
!   \item[swdown]
!     downward solar radiation. unit: W m-2
!   \item[lwdown]
!     downward longwave radiation. unit: W m-2
!   \item[prcp]
!     total precipitation (rainfall+snowfall). unit: mm
!   \item[tsk]
!     surface radiative temperature. unit: K
!   \item[hfx]
!     sensible heat flux. unit: W m-2
!   \item[qfx]
!     latent heat flux. unit: kg s-1 m-2
!   \item[lh]
!     latent heat flux. unit: W m-2
!   \item[grdflx]
!     ground/snow heat flux. unit: W m-2
!   \item[sfcrunoff]
!     accumulated surface runoff. unit: m
!   \item[udrrunoff]
!     accumulated sub-surface runoff. unit: m
!   \item[albedo]
!     total grid albedo. unit: -
!   \item[snowc]
!     snow cover fraction. unit: -
!   \item[smc]
!     volumtric soil moisture. unit: m3/m3
!   \item[sh2o]
!     volumtric liquid soil moisture. unit: m3/m3
!   \item[tslb]
!     soil temperature. unit: K
!   \item[sneqv]
!     snow water equivalent. unit: mm
!   \item[snowh]
!     physical snow depth. unit: m
!   \item[canwat]
!     total canopy water + ice. unit: mm
!   \item[acsnom]
!     accumulated snow melt leaving pack. unit: -
!   \item[acsnow]
!     accumulated snow on grid. unit: mm
!   \item[emiss]
!     surface bulk emissivity. unit: -
!   \item[rs]
!     total stomatal resistance. unit: s/m
!   \item[isnow]
!     actual no. of snow layers. unit: -
!   \item[tv]
!     vegetation leaf temperature. unit: K
!   \item[tg]
!     bulk ground surface temperature. unit: K
!   \item[canice]
!     canopy-intercepted ice. unit: mm
!   \item[canliq]
!     canopy-intercepted liquid water. unit: mm
!   \item[eah]
!     canopy air vapor pressure. unit: Pa
!   \item[tah]
!     canopy air temperature. unit: K
!   \item[cm]
!     bulk momentum drag coefficient. unit: -
!   \item[ch]
!     bulk sensible heat exchange coefficient. unit: -
!   \item[fwet]
!     wetted or snowed fraction of canopy. unit: -
!   \item[sneqvo]
!     snow mass at last time step. unit: mm h2o
!   \item[albold]
!     snow albedo at last time step. unit: -
!   \item[qsnow]
!     snowfall on the ground. unit: mm/s
!   \item[wslake]
!     lake water storage. unit: mm
!   \item[zwt]
!     water table depth. unit: m
!   \item[wa]
!     water in the "aquifer". unit: mm
!   \item[wt]
!     water in aquifer and saturated soil. unit: mm
!   \item[tsno]
!     snow layer temperature. unit: K
!   \item[zss]
!     snow/soil layer depth from snow surface. unit: m
!   \item[snowice]
!     snow layer ice. unit: mm
!   \item[snowliq]
!     snow layer liquid water. unit: mm
!   \item[lfmass]
!     leaf mass. unit: g/m2
!   \item[rtmass]
!     mass of fine roots. unit: g/m2
!   \item[stmass]
!     stem mass. unit: g/m2
!   \item[wood]
!     mass of wood (including woody roots). unit: g/m2
!   \item[stblcp]
!     stable carbon in deep soil. unit: g/m2
!   \item[fastcp]
!     short-lived carbon in shallow soil. unit: g/m2
!   \item[lai]
!     leaf area index. unit: -
!   \item[sai]
!     stem area index. unit: -
!   \item[tauss]
!     snow age factor. unit: -
!   \item[smoiseq]
!     equilibrium volumetric soil moisture content. unit: m3/m3
!   \item[smcwtd]
!     soil moisture content in the layer to the water table when deep. unit: -
!   \item[deeprech]
!     recharge to the water table when deep. unit: -
!   \item[rech]
!     recharge to the water table (diagnostic). unit: -
!   \item[grain]
!     mass of grain XING. unit: g/m2
!   \item[gdd]
!     growing degree days XING (based on 10C). unit: -
!   \item[pgs]
!     growing degree days XING. unit: -
!   \item[gecros\_state]
!     optional gecros crop. unit: -
!   \item[t2mv]
!     2m temperature of vegetation part. unit: K
!   \item[t2mb]
!     2m temperature of bare ground part. unit: K
!   \item[q2mv]
!     2m mixing ratio of vegetation part. unit: -
!   \item[q2mb]
!     2m mixing ratio of bare ground part. unit: -
!   \item[trad]
!     surface radiative temperature. unit: K
!   \item[nee]
!     net ecosys exchange of CO2. unit: g/m2/s CO2
!   \item[gpp]
!     gross primary assimilation of carbon. unit: g/m2/s C
!   \item[npp]
!     net primary productivity of carbon. unit: g/m2/s C
!   \item[fveg]
!     Noah-MP green vegetation fraction. unit: -
!   \item[runsf]
!     surface runoff. unit: mm/s
!   \item[runsb]
!     subsurface runoff. unit: mm/s
!   \item[ecan]
!     evaporation of intercepted water. unit: mm/s
!   \item[edir]
!     soil surface evaporation rate. unit: mm/s
!   \item[etran]
!     transpiration rate. unit: mm/s
!   \item[rainf]
!     raifall. unit: km s-1
!   \item[snowf]
!     snow fall. unit: kg s-1
!   \item[fsa]
!     total absorbed solar radiation. unit: W/m2
!   \item[fira]
!     total net longwave radiation [+ to atm]. unit: W/m2
!   \item[apar]
!     photosyn active energy by canopy. unit: W/m2
!   \item[psn]
!     total photosynthesis [+]. unit: umol co2/m2/s
!   \item[sav]
!     solar radiation absorbed by vegetation. unit: W/m2
!   \item[sag]
!     solar radiation absorbed by ground. unit: W/m2
!   \item[rssun]
!     sunlit leaf stomatal resistance. unit: s/m
!   \item[rssha]
!     shaded leaf stomatal resistance. unit: s/m
!   \item[bgap]
!     between gap fraction. unit: -
!   \item[wgap]
!     within gap fraction. unit: -
!   \item[tgb]
!     bare ground temperature. unit: K
!   \item[tgv]
!     under canopy ground temperature. unit: K
!   \item[chv]
!     sensible heat exchange coefficient vegetated. unit: -
!   \item[chb]
!     sensible heat exchange coefficient bare-ground. unit: -
!   \item[shg]
!     veg ground sensible heat [+ to atm]. unit: W/m2
!   \item[shc]
!     canopy sensible heat [+ to atm]. unit: W/m2
!   \item[shb]
!     bare sensible heat [+ to atm]. unit: W/m2
!   \item[evg]
!     veg ground evaporation [+ to atm]. unit: W/m2
!   \item[evb]
!     bare soil evaporation [+ to atm]. unit: W/m2
!   \item[ghv]
!     veg ground heat flux [+ to soil]. unit: W/m2
!   \item[ghb]
!     bare ground heat flux [+ to soil]. unit: W/m2
!   \item[irg]
!     veg ground net LW radiation [+ to atm]. unit: W/m2
!   \item[irc]
!     canopy net LW radiation [+ to atm]. unit: W/m2
!   \item[irb]
!     bare net LW radiation [+ to atm]. unit: W/m2
!   \item[tr]
!     transpiration [ to atm]. unit: W/m2
!   \item[evc]
!     canopy evaporation heat [to atm]. unit: W/m2
!   \item[chleaf]
!     leaf exchange coefficient. unit: -
!   \item[chuc]
!     under canopy exchange coefficient. unit: -
!   \item[chv2]
!     veg 2m exchange coefficient. unit: -
!   \item[chb2]
!     bare 2m exchange coefficient. unit: -
!   \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  10/25/18: Shugong Wang, Zhuo Wang Initial implementation for LIS 7 and NoahMP401
!
!EOP
   USE MODULE_SF_NOAHMPLSM_401
    implicit none

    INTEGER, PRIVATE, PARAMETER :: MBAND = 2
    INTEGER, PRIVATE, PARAMETER :: NSOIL = 4
    INTEGER, PRIVATE, PARAMETER :: NSTAGE = 8

#if 0
    TYPE noahmp_parameters ! define a NoahMP parameters type
       
       !---------------------------------------
       ! From the veg section of MPTABLE.TBL
       !---------------------------------------

       LOGICAL :: URBAN_FLAG
       INTEGER :: ISWATER
       INTEGER :: ISBARREN
       INTEGER :: ISICE
       INTEGER :: ISCROP
       INTEGER :: EBLFOREST
       
       REAL :: CH2OP              !maximum intercepted h2o per unit lai+sai (mm)
       REAL :: DLEAF              !characteristic leaf dimension (m)
       REAL :: Z0MVT              !momentum roughness length (m)
       REAL :: HVT                !top of canopy (m)
       REAL :: HVB                !bottom of canopy (m)
       REAL :: DEN                !tree density (no. of trunks per m2)
       REAL :: RC                 !tree crown radius (m)
       REAL :: MFSNO              !snowmelt m parameter ()
       REAL :: SAIM(12)           !monthly stem area index, one-sided
       REAL :: LAIM(12)           !monthly leaf area index, one-sided
       REAL :: SLA                !single-side leaf area per Kg [m2/kg]
       REAL :: DILEFC             !coeficient for leaf stress death [1/s]
       REAL :: DILEFW             !coeficient for leaf stress death [1/s]
       REAL :: FRAGR              !fraction of growth respiration  !original was 0.3 
       REAL :: LTOVRC             !leaf turnover [1/s]
       
       REAL :: C3PSN              !photosynthetic pathway: 0. = c4, 1. = c3
       REAL :: KC25               !co2 michaelis-menten constant at 25c (pa)
       REAL :: AKC                !q10 for kc25
       REAL :: KO25               !o2 michaelis-menten constant at 25c (pa)
       REAL :: AKO                !q10 for ko25
       REAL :: VCMX25             !maximum rate of carboxylation at 25c (umol co2/m**2/s)
       REAL :: AVCMX              !q10 for vcmx25
       REAL :: BP                 !minimum leaf conductance (umol/m**2/s)
       REAL :: MP                 !slope of conductance-to-photosynthesis relationship
       REAL :: QE25               !quantum efficiency at 25c (umol co2 / umol photon)
       REAL :: AQE                !q10 for qe25
       REAL :: RMF25              !leaf maintenance respiration at 25c (umol co2/m**2/s)
       REAL :: RMS25              !stem maintenance respiration at 25c (umol co2/kg bio/s)
       REAL :: RMR25              !root maintenance respiration at 25c (umol co2/kg bio/s)
       REAL :: ARM                !q10 for maintenance respiration
       REAL :: FOLNMX             !foliage nitrogen concentration when f(n)=1 (%)
       REAL :: TMIN               !minimum temperature for photosynthesis (k)
       
       REAL :: XL                 !leaf/stem orientation index
       REAL :: RHOL(MBAND)        !leaf reflectance: 1=vis, 2=nir
       REAL :: RHOS(MBAND)        !stem reflectance: 1=vis, 2=nir
       REAL :: TAUL(MBAND)        !leaf transmittance: 1=vis, 2=nir
       REAL :: TAUS(MBAND)        !stem transmittance: 1=vis, 2=nir

       REAL :: MRP                !microbial respiration parameter (umol co2 /kg c/ s)
       REAL :: CWPVT              !empirical canopy wind parameter
       
       REAL :: WRRAT              !wood to non-wood ratio
       REAL :: WDPOOL             !wood pool (switch 1 or 0) depending on woody or not [-]
       REAL :: TDLEF              !characteristic T for leaf freezing [K]
       
       INTEGER :: NROOT              !number of soil layers with root present
       REAL :: RGL                !Parameter used in radiation stress function
       REAL :: RSMIN              !Minimum stomatal resistance [s m-1]
       REAL :: HS                 !Parameter used in vapor pressure deficit function
       REAL :: TOPT               !Optimum transpiration air temperature [K]
       REAL :: RSMAX              !Maximal stomatal resistance [s m-1]
       
       REAL :: SLAREA
       REAL :: EPS(5)
       
!------------------------------------------------------------------------------------------!
! From the rad section of MPTABLE.TBL
!------------------------------------------------------------------------------------------!

       REAL :: ALBSAT(MBAND)       !saturated soil albedos: 1=vis, 2=nir
       REAL :: ALBDRY(MBAND)       !dry soil albedos: 1=vis, 2=nir
       REAL :: ALBICE(MBAND)       !albedo land ice: 1=vis, 2=nir
       REAL :: ALBLAK(MBAND)       !albedo frozen lakes: 1=vis, 2=nir
       REAL :: OMEGAS(MBAND)       !two-stream parameter omega for snow
       REAL :: BETADS              !two-stream parameter betad for snow
       REAL :: BETAIS              !two-stream parameter betad for snow
       REAL :: EG(2)               !emissivity
       
!------------------------------------------------------------------------------------------!
! From the globals section of MPTABLE.TBL
!------------------------------------------------------------------------------------------!
 
       REAL :: CO2          !co2 partial pressure
       REAL :: O2           !o2 partial pressure
       REAL :: TIMEAN       !gridcell mean topgraphic index (global mean)
       REAL :: FSATMX       !maximum surface saturated fraction (global mean)
       REAL :: Z0SNO        !snow surface roughness length (m) (0.002)
       REAL :: SSI          !liquid water holding capacity for snowpack (m3/m3)
       REAL :: SWEMX        !new snow mass to fully cover old snow (mm)
       REAL :: RSURF_SNOW   !surface resistance for snow(s/m)
       
!------------------------------------------------------------------------------------------!
! From the crop section of MPTABLE.TBL
!------------------------------------------------------------------------------------------!
       
       INTEGER :: PLTDAY           ! Planting date
       INTEGER :: HSDAY            ! Harvest date
       REAL :: PLANTPOP         ! Plant density [per ha] - used?
       REAL :: IRRI             ! Irrigation strategy 0= non-irrigation 1=irrigation (no water-stress)
       REAL :: GDDTBASE         ! Base temperature for GDD accumulation [C]
       REAL :: GDDTCUT          ! Upper temperature for GDD accumulation [C]
       REAL :: GDDS1            ! GDD from seeding to emergence
       REAL :: GDDS2            ! GDD from seeding to initial vegetative 
       REAL :: GDDS3            ! GDD from seeding to post vegetative 
       REAL :: GDDS4            ! GDD from seeding to intial reproductive
       REAL :: GDDS5            ! GDD from seeding to pysical maturity 
       INTEGER :: C3C4             ! photosynthetic pathway:  1 = c3 2 = c4
       REAL :: AREF             ! reference maximum CO2 assimulation rate 
       REAL :: PSNRF            ! CO2 assimulation reduction factor(0-1) (caused by non-modeling part,e.g.pest,weeds)
       REAL :: I2PAR            ! Fraction of incoming solar radiation to photosynthetically active radiation
       REAL :: TASSIM0          ! Minimum temperature for CO2 assimulation [C]
       REAL :: TASSIM1          ! CO2 assimulation linearly increasing until temperature reaches T1 [C]
       REAL :: TASSIM2          ! CO2 assmilation rate remain at Aref until temperature reaches T2 [C]
       REAL :: K                ! light extinction coefficient
       REAL :: EPSI             ! initial light use efficiency
       REAL :: Q10MR            ! q10 for maintainance respiration
       REAL :: FOLN_MX          ! foliage nitrogen concentration when f(n)=1 (%)
       REAL :: LEFREEZ          ! characteristic T for leaf freezing [K]
       REAL :: DILE_FC(NSTAGE)  ! coeficient for temperature leaf stress death [1/s]
       REAL :: DILE_FW(NSTAGE)  ! coeficient for water leaf stress death [1/s]
       REAL :: FRA_GR           ! fraction of growth respiration 
       REAL :: LF_OVRC(NSTAGE)  ! fraction of leaf turnover  [1/s]
       REAL :: ST_OVRC(NSTAGE)  ! fraction of stem turnover  [1/s]
       REAL :: RT_OVRC(NSTAGE)  ! fraction of root tunrover  [1/s]
       REAL :: LFMR25           ! leaf maintenance respiration at 25C [umol CO2/m**2  /s]
       REAL :: STMR25           ! stem maintenance respiration at 25C [umol CO2/kg bio/s]
       REAL :: RTMR25           ! root maintenance respiration at 25C [umol CO2/kg bio/s]
       REAL :: GRAINMR25        ! grain maintenance respiration at 25C [umol CO2/kg bio/s]
       REAL :: LFPT(NSTAGE)     ! fraction of carbohydrate flux to leaf
       REAL :: STPT(NSTAGE)     ! fraction of carbohydrate flux to stem
       REAL :: RTPT(NSTAGE)     ! fraction of carbohydrate flux to root
       REAL :: GRAINPT(NSTAGE)  ! fraction of carbohydrate flux to grain
       REAL :: BIO2LAI          ! leaf are per living leaf biomass [m^2/kg]

!------------------------------------------------------------------------------------------!
! From the SOILPARM.TBL tables, as functions of soil category.
!------------------------------------------------------------------------------------------!
       REAL :: BEXP(NSOIL)   !B parameter
       REAL :: SMCDRY(NSOIL) !dry soil moisture threshold where direct evap from top
       !layer ends (volumetric) (not used MB: 20140718)
       REAL :: SMCWLT(NSOIL) !wilting point soil moisture (volumetric)
       REAL :: SMCREF(NSOIL) !reference soil moisture (field capacity) (volumetric)
       REAL :: SMCMAX(NSOIL) !porosity, saturated value of soil moisture (volumetric)
       REAL :: PSISAT(NSOIL) !saturated soil matric potential
       REAL :: DKSAT(NSOIL)  !saturated soil hydraulic conductivity
       REAL :: DWSAT(NSOIL)  !saturated soil hydraulic diffusivity
       REAL :: QUARTZ(NSOIL) !soil quartz content
       REAL :: F1            !soil thermal diffusivity/conductivity coef (not used MB: 20140718)
!------------------------------------------------------------------------------------------!
! From the GENPARM.TBL file
!------------------------------------------------------------------------------------------!
       REAL :: SLOPE       !slope index (0 - 1)
       REAL :: CSOIL       !vol. soil heat capacity [j/m3/K]
       REAL :: ZBOT        !Depth (m) of lower boundary soil temperature
       REAL :: CZIL        !Calculate roughness length of heat
       REAL :: REFDK
       REAL :: REFKDT
       
       REAL :: KDT         !used in compute maximum infiltration rate (in INFIL)
       REAL :: FRZX        !used in compute maximum infiltration rate (in INFIL)
     
    END TYPE noahmp_parameters
#endif
    type, public :: noahmp401dec
        !------------------------------------------------------
        ! forcing
        !------------------------------------------------------
        real               :: tair
        real               :: sfctmp    ! Yeosang Yoon for snow DA
        real               :: psurf
        real               :: wind_e
        real               :: wind_n
        real               :: qair
        real               :: swdown
        real               :: lwdown
        real               :: prcp
        !--------------------------------------------------------
        ! spatial parameter
        !--------------------------------------------------------
        integer            :: vegetype
        integer            :: soiltype
        real               :: tbot
        real               :: planting
        real               :: harvest
        real               :: season_gdd
        real               :: soilcL1
        real               :: soilcL2
        real               :: soilcL3
        real               :: soilcL4
        !----------------------------------------------------------
        ! multilevel spatial parameter
        !----------------------------------------------------------
        real, pointer      :: shdfac_monthly(:)
        real, pointer      :: soilcomp(:)
        !----------------------------------------------------------
        ! state
        !----------------------------------------------------------
        real               :: sfcrunoff
        real               :: udrrunoff
        real, pointer      :: smc(:)
        real, pointer      :: sh2o(:)
        real, pointer      :: tslb(:)
        real               :: sneqv
        real               :: snowh
        real               :: canwat
        real               :: acsnom
        real               :: acsnow
        integer            :: isnow
        real               :: tv
        real               :: tg
        real               :: canice
        real               :: canliq
        real               :: eah
        real               :: tah
        real               :: cm
        real               :: ch
        real               :: fwet
        real               :: sneqvo
        real               :: albold
        real               :: qsnow
        real               :: wslake
        real               :: zwt
        real               :: wa
        real               :: wt
        real, pointer      :: tsno(:)
        real, pointer      :: zss(:)
        real, pointer      :: snowice(:)
        real, pointer      :: snowliq(:)
        real               :: lfmass
        real               :: rtmass
        real               :: stmass
        real               :: wood
        real               :: stblcp
        real               :: fastcp
        real               :: lai
        real               :: sai
        real               :: tauss
        real, pointer      :: smoiseq(:)
        real               :: smcwtd
        real               :: deeprech
        real               :: rech
        real               :: grain
        real               :: gdd
        integer            :: pgs
        real, pointer      :: gecros_state(:)
        !-------------------------------------------------------
        ! output
        !-------------------------------------------------------
        real               :: tsk
!       real               :: fsh
        real               :: hfx
        real               :: qfx
        real               :: lh
        real               :: grdflx
        real               :: albedo
        real               :: snowc
        real               :: emiss
        real               :: rs
        real               :: t2mv
        real               :: t2mb
        real               :: q2mv
        real               :: q2mb
        real               :: trad
        real               :: nee
        real               :: gpp
        real               :: npp
        real               :: fveg
        real               :: runsf
        real               :: runsb
        real               :: ecan
        real               :: edir
        real               :: etran
        real               :: rainf
        real               :: snowf
        real               :: fsa
        real               :: fira
        real               :: apar
        real               :: psn
        real               :: sav
        real               :: sag
        real               :: rssun
        real               :: rssha
        real               :: bgap
        real               :: wgap
        real               :: tgb
        real               :: tgv
        real               :: chv
        real               :: chb
        real               :: shg
        real               :: shc
        real               :: shb
        real               :: evg
        real               :: evb
        real               :: ghv
        real               :: ghb
        real               :: irg
        real               :: irc
        real               :: irb
        real               :: tr
        real               :: evc
        real               :: chleaf
        real               :: chuc
        real               :: chv2
        real               :: chb2

        !EMK for 557WW
        real :: tair_agl_min
        real :: rhmin

        type(noahmp_parameters) :: param
    end type noahmp401dec

end module NoahMP401_module
