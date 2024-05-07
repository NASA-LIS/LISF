!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module RDHM356_module
!BOP
!
! !MODULE: RDHM356_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the RDHM356 1-d variables.
!  The variables specified in the data structure include:
!
!  \begin{description}      
!   \item[Tair]
!     air temperature. unit: K
!   \item[Psurf]
!     surface air pressure. unit: Pa
!   \item[Wind\_E]
!     eastward wind. unit: m s-1
!   \item[Wind\_N]
!     northward wind. unit: m s-1
!   \item[Qair]
!     near surface specific humidity. unit: kg kg-1
!   \item[Rainf]
!     rainfall rate. unit: kg m-2 s-1
!   \item[Snowf]
!     snowfall rate. unit: kg m-2 s-1
!   \item[Swdown]
!     incident shortwave radiation. unit: W m-2
!   \item[Lwdown]
!     incident longwave radiation. unit: W m-2
!   \item[Tair\_min]
!     daily minimum air temperature. unit: K
!   \item[Tair\_max]
!     daily maximum air temperature. unit: K
!   \item[TempHeight]
!     observation height of temperature of humidity. unit: m
!   \item[WindHeight]
!     observation height of wind. unit: m
!   \item[DT\_SAC\_SNOW17]
!     simulation time interval of SAC model and Snow-17. unit: s
!   \item[DT\_FRZ]
!     simulation time interval of frozen soil model. unit: s
!   \item[FRZ\_VER\_OPT]
!     version number of frozen soil model. 1: old version, 2: new version. unit: -
!   \item[SNOW17\_OPT]
!     option for snow-17. If SNOW17\_OPT=1, use SNOW-17, otherwise, don't use. unit: -
!   \item[NSTYP]
!     number of soil types. unit: -
!   \item[NVTYP]
!     number of vegetation types. unit: -
!   \item[NDINTW]
!     number of desired soil layers for total and liquid soil moisture. unit: -
!   \item[NDSINT]
!     number of desired soil layers for soil temperature. unit: -
!   \item[NORMALIZE]
!     normalization flag for total and liquid soil moisture output (1-normalized, 0-not). unit: -
!   \item[DSINTW]
!     thickness of desired soil layers for liquid and total soil moisture. unit: cm
!   \item[DSINT]
!     thickness of desired soil layers for soil temperature. unit: cm
!   \item[PET\_MON]
!     multiband monthly PET climatology, time series of spatial parameter. unit: mm
!   \item[GRN\_MON]
!     multiband monthly greenness climatology, time series of spatial parameter [-]. unit: -
!   \item[PETADJ\_MON]
!     adjustment of PET for 12 months. unit: -
!   \item[SoilAlb]
!     snow free ALBEDO (default value 0.15). unit: -
!   \item[SnowAlb]
!     snow ALBEDO (default value 0.7). unit: -
!   \item[SOILTYP]
!     Soil type. unit: -
!   \item[VEGETYP]
!     Vegetation type. unit: -
!   \item[UZTWM]
!     upper zone tension water maximum storage. unit: mm
!   \item[UZFWM]
!     upper zone free water maximum storage. unit: mm
!   \item[UZK]
!     upper zone free water latent depletion rate. unit: day$^{-1}$
!   \item[PCTIM]
!     impervious fraction of the watershad area. unit: -
!   \item[ADIMP]
!     additional impervious area. unit: -
!   \item[RIVA]
!     riparian vegetation area. unit: -
!   \item[ZPERC]
!     maximum percolation rate. unit: -
!   \item[REXP]
!     exponent of the percolation equation (percolation parameter). unit: -
!   \item[LZTWM]
!     lower zone tension water maximum storage. unit: mm
!   \item[LZFSM]
!     lower zone supplemental free water (fast) maximum storage. unit: mm
!   \item[LZFPM]
!     lower zone primary free water (slow) maximum storage. unit: mm
!   \item[LZSK]
!     lower zone supplemental free water depletion rate. unit: day$^{-1}$
!   \item[LZPK]
!     lower zone primary free water depletion rate. unit: day$^{-1}$
!   \item[PFREE]
!     fraction percolation from upper to lower free water storage. unit: day$^{-1}$
!   \item[SIDE]
!     ratio of deep recharge to channel base flow. unit: -
!   \item[RSERV]
!     fraction of lower zone free water not transferable to tension water. unit: -
!   \item[EFC]
!     fraction of forest cover. unit: -
!   \item[TBOT]
!     bottom boundary soil temperature. unit: $^\circ$C
!   \item[RSMAX]
!     maximum residual porosity (usually = 0.58). unit: -
!   \item[CKSL]
!     ratio of frozen to non-frozen surface (increase in frozen ground contact, usually = 8 s/m). unit: s/m
!   \item[ZBOT]
!     lower boundary depth (negative value, usually = -2.5 m). unit: m
!   \item[vegRCMIN]
!     minimal stomatal resistance table for SACHTET, 14 values. unit: s/m
!   \item[climRCMIN]
!     climate dependent miminal stomatal resistance for SACHTET, 14 values. unit: s/m
!   \item[RGL]
!     solar radiation threshold table for SACHTET, 14 values. unit: W m-2
!   \item[HS]
!     vapor pressure resistance factor table for SACHTET, 14 values. unit: -
!   \item[LAI]
!     leaf area index table for SACHTET, 14 values. unit: -
!   \item[D50]
!     the depth (cm) table at which 50\% roots are allocated for SACHTET, 14 values. unit: cm
!   \item[CROOT]
!     root distribution parameter table for SACHTET, 14 values. unit: -
!   \item[Z0]
!     roughness coefficient of surface. unit: m
!   \item[CLAY]
!     clay content for SACHTET, 12 values. unit: -
!   \item[SAND]
!     sand content for sACHTET, 12 values. unit: -
!   \item[SATDK]
!     saturated hydraulic conductivityfor SACHTET, 12 values. unit: m s-1
!   \item[CZIL]
!     default=0.12 Zilitinkevich. unit: -
!   \item[FXEXP]
!     FXEXP(fxexp),(default=2.0) bare soil. unit: -
!   \item[vegRCMAX]
!     RCMAX,(default=5000s/m) maximum stomatal resistance. unit: s/m
!   \item[TOPT]
!     TOPT,(default=298K)optimum air. unit: K
!   \item[PC]
!     plant coef. default pc = -1, 0.6 - 0.8. unit: -
!   \item[PET\_OPT]
!     if PET\_OPT = 0, use non Penmann-based ETP;if penpt $>$ 0 empirical Penmann equation; if penpt $<$ 0, use energy based Pennman. unit: -
!   \item[RDST]
!     default=1 means noah option,this constant allows selection of tension water redistribution option, 
!    if rdst = 0 (ohd), use OHD version of SRT subroutine this SRT uses reference gradient instead an actual. 
!    if rdst = 1 ( noah), use Noah version of SRT subroutine. unit: -
!   \item[thresholdRCMIN]
!     this constant allows change of RCMIN (0.5). unit: s/m
!   \item[SFCREF]
!     reference wind speed for PET adjustment (4 m s-1). unit: m/s
!   \item[BAREADJ]
!     Ek-Chen evaporation threshold switch. Bare soil evaporation option changes according to greenness.. unit: -
!   \item[ALON]
!     logitude. unit: -
!   \item[ALAT]
!     latitude. unit: -
!   \item[SCF]
!     snow fall correction factor. unit: -
!   \item[MFMAX]
!     maximum melt factor. unit: mm/(6hr$^\circ$C)
!   \item[MFMIN]
!     minimum melt factor. unit: mm/(6hr$^\circ$C)
!   \item[NMF]
!     maximum negative melt factor. unit: mm/(6hr$^\circ$C)
!   \item[UADJ]
!     the average wind function during rain-on-snow periods. unit: mm/mb
!   \item[SI]
!     areal water-equivalent above which 100 percent areal snow cover. unit: mm
!   \item[MBASE]
!     base temperature for non-rain melt factor. unit: $^\circ$C
!   \item[PXTEMP]
!     temperature which spereates rain from snow. unit: $^\circ$C
!   \item[PLWHC]
!     maximum amount of liquid-water held against gravity drainage. unit: -
!   \item[TIPM]
!     antecedent snow temperature index parameter. unit: -
!   \item[GM]
!     daily ground melt. unit: mm/day
!   \item[ELEV]
!     elevation. unit: m
!   \item[LAEC]
!     snow-rain split temperature. unit: $^\circ$C
!   \item[ADC]
!     multiband Snow-17 curve coordinates. unit: -
!   \item[SNOW17\_SWITCH]
!     switch variable change liquid water freezing version, 0: Victor's version, 1: Eric's version. unit: -
!   \item[UZTWC]
!     upper zone tension water storage content. unit: mm
!   \item[UZFWC]
!     upper zone free water storage content. unit: mm
!   \item[LZTWC]
!     lower zone tension water storage content. unit: mm
!   \item[LZFPC]
!     lower zone primary free water storage content. unit: mm
!   \item[LZFSC]
!     lower zone supplemental free water storage content. unit: mm
!   \item[ADIMC]
!     additional impervious area content. unit: mm
!   \item[TS0]
!     first soil layer temperature. unit: $^\circ$C
!   \item[TS1]
!     second soil layer temperature. unit: $^\circ$C
!   \item[TS2]
!     third soil layer temperature. unit: $^\circ$C
!   \item[TS3]
!     fourth soil layer temperature. unit: $^\circ$C
!   \item[TS4]
!     fifth soil layer temperature. unit: $^\circ$C
!   \item[UZTWH]
!     unfrozen upper zone tension water. unit: mm
!   \item[UZFWH]
!     unfrozen uppeer zone free water. unit: mm
!   \item[LZTWH]
!     unfrozen lower zone tension water. unit: mm
!   \item[LZFSH]
!     unfrozen lower zone supplemental free water. unit: mm
!   \item[LZFPH]
!     unfrozen lower zone primary free water. unit: mm
!   \item[SMC]
!     volumetric content of total soil moisture at each layer. unit: m^3 m-3
!   \item[SH2O]
!     volumetric content of liquid soil moisture at each layer. unit: m^3 m-3
!   \item[WE]
!     snow water equivalent without liquid water. unit: mm
!   \item[LIQW]
!     liquid water in snow. unit: mm
!   \item[NEGHS]
!     negative snow heat. unit: mm
!   \item[TINDEX]
!     antecedent temperature index. unit: $^\circ$C
!   \item[ACCMAX]
!     cumulated snow water including liquid. unit: mm
!   \item[SNDPT]
!     snow depth. unit: cm
!   \item[SNTMP]
!     average snow temperature. unit: $^\circ$C
!   \item[SB]
!     the last highest snow water equivalent before any snow fall. unit: $^\circ$C
!   \item[SBAESC]
!     internal snow state during melt \& new snow fall (checked with Victor). unit: -
!   \item[SBWS]
!     internal snow state during melt \& new snow fall (checked with Victor). unit: -
!   \item[STORAGE]
!     snow liquid water attenuation storage. unit: mm
!   \item[AEADJ]
!     adjusted areal snow cover fraction. unit: -
!   \item[EXLAG]
!     array of lagged liquid water values. unit: -
!   \item[NEXLAG]
!     number of ordinates in lagged liquid water array (EXLAG). unit: -
!   \item[TA\_PREV]
!     air temperature of previous time step. unit: -
!   \item[SWINT]
!     total volumetric soil moisture contents at desired soil layers (can be different from soil layers). unit: -
!   \item[SWHINT]
!     liquid volumetric soil moisture contents at desired soil layers (can be different from soil layers). unit: -
!   \item[TSINT]
!     soil temperature at desired soil layers (can be different from soil layers). unit: -
!   \item[FRZDUP]
!     depth of the upper border of frozen ground from surface. unit: m
!   \item[FRZDBT]
!     depth of the bottom border of frozen ground from surface. unit: m
!   \item[FROST]
!     frost index. unit: -
!   \item[ALBEDO]
!     land surface albedo. unit: -
!   \item[SURF]
!     Qs $<$=$>$ SURF simulated fast runoff (surface runoff). unit: mm s-1
!   \item[GRND]
!     Qsb $<$=$>$ GRND simulated slow runoff (baseflow). unit: mm s-1
!   \item[TET]
!     Evap $<$=$>$ TET simulated actual evapotranspiration. unit: mm s-1
!   \item[EDMND]
!     PotEvap $<$=$>$ EDMND potential evapotranspiration. unit: mm s-1
!   \item[CH]
!     surface layer exchage coefficient for heat and moisture. unit: s/m
!   \item[CM]
!     surface layer exchange coefficient for momentum (drag coefficient). unit: s/m
!   \item[SnowFrac]
!     snow cover fraction, Snow-17. unit: -
!   \item[SWE]
!     snow water equivalent, Snow-17. unit: kg m-2
!   \item[SnowDepth]
!     snow depth, Snow-17. unit: m
!   \item[RM]
!     rain + melt. unit: mm
!   \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  11/5/13: Shugong Wang Initial implementation for LIS 7 and RDHM356
!
!EOP
    implicit none
    private
    type, public :: rdhm356dec
        !-------------------------------------------------------------------------
        ! forcing
        !-------------------------------------------------------------------------
        real               :: Tair
        real               :: Psurf
        real               :: Wind_E
        real               :: Wind_N
        real               :: Qair
        real               :: Rainf
        real               :: Snowf
        real               :: Swdown
        real               :: Lwdown
        !-------------------------------------------------------------------------
        ! spatial parameter
        !-------------------------------------------------------------------------
        real               :: SoilAlb
        real               :: SnowAlb
        integer            :: SOILTYP
        integer            :: VEGETYP
        real               :: UZTWM
        real               :: UZFWM
        real               :: UZK
        real               :: PCTIM
        real               :: ADIMP
        real               :: RIVA
        real               :: ZPERC
        real               :: REXP
        real               :: LZTWM
        real               :: LZFSM
        real               :: LZFPM
        real               :: LZSK
        real               :: LZPK
        real               :: PFREE
        real               :: SIDE
        real               :: RSERV
        real               :: EFC
        real               :: TBOT
        real               :: RSMAX
        real               :: CKSL
        real               :: ZBOT
        real               :: vegRCMIN
        real               :: climRCMIN
        real               :: RGL
        real               :: HS
        real               :: LAI
        real               :: D50
        real               :: CROOT
        real               :: Z0
        real               :: CLAY
        real               :: SAND
        real               :: SATDK
        real               :: ALON
        real               :: ALAT
        real               :: SCF
        real               :: MFMAX
        real               :: MFMIN
        real               :: NMF
        real               :: UADJ
        real               :: SI
        real               :: MBASE
        real               :: PXTEMP
        real               :: PLWHC
        real               :: TIPM
        real               :: GM
        real               :: ELEV
        real               :: LAEC
        !-------------------------------------------------------------------------
        ! dynamic spatial parameter
        !-------------------------------------------------------------------------
        real               :: Tair_min
        real               :: Tair_max
        !-------------------------------------------------------------------------
        ! multiband spatial parameter
        !-------------------------------------------------------------------------
        real               :: PET_MON(12)
        real               :: GRN_MON(12)
        real               :: ADC(11)
        !-------------------------------------------------------------------------
        ! state
        !-------------------------------------------------------------------------
        real               :: UZTWC
        real               :: UZFWC
        real               :: LZTWC
        real               :: LZFPC
        real               :: LZFSC
        real               :: ADIMC
        real               :: TS0
        real               :: TS1
        real               :: TS2
        real               :: TS3
        real               :: TS4
        real               :: UZTWH
        real               :: UZFWH
        real               :: LZTWH
        real               :: LZFSH
        real               :: LZFPH
        real               :: SMC(6)
        real               :: SH2O(6)
        real               :: WE
        real               :: LIQW
        real               :: NEGHS
        real               :: TINDEX
        real               :: ACCMAX
        real               :: SNDPT
        real               :: SNTMP
        real               :: SB
        real               :: SBAESC
        real               :: SBWS
        real               :: STORAGE
        real               :: AEADJ
        real               :: EXLAG(7)
        integer            :: NEXLAG
        real               :: TA_PREV
        real               :: CH
        real               :: CM
        !-------------------------------------------------------------------------
        ! output
        !-------------------------------------------------------------------------
        real, allocatable      :: SWINT(:)
        real, allocatable      :: SWHINT(:)
        real, allocatable      :: TSINT(:)
        real               :: FRZDUP
        real               :: FRZDBT
        real               :: FROST
        real               :: ALBEDO
        real               :: SURF
        real               :: GRND
        real               :: TET
        real               :: EDMND
        real               :: SnowFrac
        real               :: SWE
        real               :: SnowDepth
        real               :: RM
    end type rdhm356dec
end module RDHM356_module
