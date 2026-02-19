!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

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
!   \item[wtrflx]
!    total water flux. unit: kg m-2 s-1
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
        !ag (05Jan2021)
        ! 2-way coupling parameters
        real               :: rivsto
        real               :: fldsto
        real               :: fldfrc

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

        ! For WRF-HYDRO
	real               :: sfcheadrt
	real               :: infxs1rt
	real               :: soldrain1rt
#ifdef PARFLOW
        real, pointer      :: wtrflx(:)
#endif
 
    end type noahmp401dec

end module NoahMP401_module
