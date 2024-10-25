!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module RUC37_module
!BOP
!
! !MODULE: RUC37_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the RUC37 1-d variables.
!  The variables specified in the data structure include:
!
!  \begin{description}      
!   \item[year]
!     year. unit: -
!   \item[month]
!     month. unit: -
!   \item[day]
!     day. unit: -
!   \item[hour]
!     hour. unit: -
!   \item[minute]
!     minute. unit: -
!   \item[lwdown]
!     downward longwave radiation flux at surface (w m-2) [forcing]. unit: W m-2
!   \item[swdown]
!     downward shortwave radiation flux at surface (w m-2) [forcing]. unit: W m-2
!   \item[psurf]
!     surface atmospheric pressure (pa) [forcing]. unit: Pa
!   \item[rainf]
!     rainfall rate (kg m-2 s-1) [forcing]. unit: kg m-2 s-1
!   \item[snowf]
!     snowfall rate (kg m-2 s-1) [forcing]. unit: Kg/m2s
!   \item[tair]
!     air temperature (k) [forcing]. unit: K
!   \item[qair]
!     surface specific humidity (kg kg-1) [forcing]. unit: kg kg-1
!   \item[wind\_e]
!     eastward wind speed (m s-1) [forcing]. unit: m s-1
!   \item[wind\_n]
!     northward wind speed (m s-1) [forcing]. unit: m s-1
!   \item[vegetype]
!     vegetation category. unit: -
!   \item[soiltype]
!     soil category. unit: -
!   \item[dt]
!     time step (seconds).. unit: s
!   \item[soil\_layer\_thickness]
!     thicknesses of each soil level (m). unit: m
!   \item[use\_local\_param]
!     .true. to use table values for albbck, shdfac, and z0brd; .false. to use values for albbck, shdfac, and z0brd as set in this driver routine. unit: -
!   \item[use\_2d\_lai\_map]
!     if rdlai2d == .true., then the xlai value that we pass to lsmruc will be used. if rdlai2d == .false., then xlai will be computed within lsmruc, from table minimum and maximum values in vegparm.tbl, and the current green vegetation fraction.. unit: -
!   \item[use\_monthly\_albedo\_map]
!     if usemonalb == .true., then the alb value passed to lsmruc will be used as the background snow-free albedo term.  if usemonalb == .false., then alb will be computed within lsmruc from minimum and maximum values in vegparm.tbl, and the current green vegetation fraction.. unit: -
!   \item[option\_iz0tlnd]
!     option to turn on (iz0tlnd=1) or off (iz0tlnd=0) the vegetation-category-dependent calculation of the zilitinkivich coefficient czil in the sfcdif subroutines.. unit: -
!   \item[option\_sfcdif]
!     option to use previous (sfcdif\_option=0) or updated (sfcdif\_option=1) version of sfcdif subroutine.. unit: -
!   \item[landuse\_tbl\_name]
!     noah model landuse parameter table. unit: -
!   \item[soil\_tbl\_name]
!     noah model soil parameter table. unit: -
!   \item[landuse\_scheme\_name]
!     landuse classification scheme. unit: -
!   \item[soil\_scheme\_name]
!     soil classification scheme. unit: -
!   \item[nsoil]
!     number of soil levels.. unit: -
!   \item[water\_class\_num]
!     number of water category in llanduse classification. unit: -
!   \item[ice\_class\_num]
!     number of ice category in llanduse classification. unit: -
!   \item[albedo\_monthly]
!     monthly values of background (i.e., snow-free) albedo ( fraction [0.0-1.0] ). unit: -
!   \item[shdfac\_monthly]
!     monthly values for green vegetation fraction ( fraction [0.0-1.0] ). unit: -
!   \item[z0brd\_monthly]
!     monthly values for background (i.e., snow-free) roughness length ( m ). unit: m
!   \item[lai\_monthly]
!     monthly values for leaf area index ( dimensionless ). unit: -
!   \item[albbck]
!     background snow-free albedo (0.0-1.0).. unit: -
!   \item[tbot]
!     deep-soil time-invariant temperature (k).  representing sort of a mean annual air temperature.. unit: K
!   \item[snoalb]
!     maximum snow albedo over deep snow (0.0-1.0). unit: -
!   \item[emiss]
!     surface emissivity (0.0 - 1.0).. unit: -
!   \item[ch]
!     exchange coefficient for head and moisture (m s-1).. unit: s/m
!   \item[cm]
!     exchange coefficient for momentum (m s-1).. unit: s/m
!   \item[sneqv]
!     water equivalent of accumulated snow depth (m).. unit: m
!   \item[snowh]
!     physical snow depth (m).. unit: m
!   \item[snowc]
!     fractional snow cover ( fraction [0.0-1.0] ). unit: -
!   \item[canwat]
!     canopy moisture content (kg m-2). unit: kg m-2
!   \item[alb]
!     surface albedo including possible snow-cover effect.  this is set in lsmruc,. unit: -
!   \item[smc]
!     total soil moisture content (m3 m-3). unit: m^3 m-3
!   \item[sho]
!     liquid soil moisture content (m3 m-3). unit: m^3 m-3
!   \item[stc]
!     soil temperature (k). unit: K
!   \item[smfr]
!     soil ice (m3 m-3). unit: m^3 m-3
!   \item[keepfr]
!     frozen soil flag: 0. or 1.
!   \item[tskin]
!     skin temperature (k). unit: K
!   \item[qvg]
!     effective mixing ratio at the surface ( kg kg{-1} ). unit: kg kg-1
!   \item[qsfc]
!     specific humidity at the surface ( kg kg{-1} ). unit: kg kg-1
!   \item[qcg]
!     effective cloud water mixing ratio at the surface ( kg kg{-1} ). unit: kg kg-1
!   \item[qsg]
!     surface water vapor mixing ratio at satration (kg kg-1). unit: kg/kg
!   \item[snt75cm]
!     snow temperature at 7.5 cm depth (k). unit: K
!   \item[tsnav]
!     average snow temperature in k. unit: K
!   \item[soilm]
!     total soil column moisture content, frozen and unfrozen ( m ). unit: m
!   \item[smroot]
!     available soil moisture in the root zone ( fraction [smcwlt-smcmax]. unit: m^3 m-3
!   \item[qh]
!     sensible heat flux ( w m{-2} ). unit: W m-2
!   \item[qle]
!     latent heat flux (evapotranspiration) ( w m{-2} ). unit: W m-2
!   \item[eta]
!     latent heat flux (evapotranspiration) ( kg m{-2} s{-1} ). unit: kg m-2 s-1
!   \item[qg]
!     soil heat flux ( w m{-2} ). unit: W m-2
!   \item[rhosnf]
!     density of frozen precipitation (kg m{-3}). unit: kg/m3
!   \item[precipfr]
!     time-step frozen precipitation (kg m{-2}). unit: kg m-2
!   \item[snowfallac]
!     run total snowfall accumulation (kg m{-2}). unit: kg m-2
!   \item[acsnow]
!     run total frozen precipitation (kg m{-2}). unit: kg m-2
!   \item[sfcevp]
!     run total evaporation flux  (kg m{-2}). unit: kg m-2
!   \item[snomlt]
!     snow melt water ( m ). unit: m
!   \item[dew]
!     dewfall (or frostfall for t$<$273.15) ( m ). unit: m
!   \item[drip]
!     throughfall of precipitation from canopy (kg m{-2} s{-1}). unit: kg m-2 s-1
!   \item[qs]
!     surface runoff, not infiltrating the soil ( m s{-1} ). unit: kg m-2 s-1
!   \item[qsb]
!     subsurface runoff, drainage out the bottom of the last soil layer ( m s{-1} ). unit: kg m-2 s-1
!   \item[snflx]
!     snow heat flux (w/m^2: negative, if downward from surface). unit: W m-2
!   \item[edir]
!     latent heat flux component: direct soil evaporation ( w m{-2} ). unit: W m-2
!   \item[ec]
!     latent heat flux component: canopy water evaporation ( w m{-2} ). unit: W m-2
!   \item[ett]
!     latent heat flux component: total plant transpiration ( w m{-2} ). unit: W m-2
!   \item[esnow]
!     sublimation from snowpack (w m{-2}). unit: W m-2
!   \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  1/15/15: Shugong Wang Initial implementation for LIS 7 and RUC37
!
!EOP
    implicit none
    private
    type, public :: ruc37dec
        !-------------------------------------------------------------------------
        ! forcing
        !-------------------------------------------------------------------------
        real               :: lwdown
        real               :: swdown
        real               :: psurf
        real               :: rainf
        real               :: snowf
        real               :: tair
        real               :: qair
        real               :: wind_e
        real               :: wind_n
        !-------------------------------------------------------------------------
        ! spatial parameter
        !-------------------------------------------------------------------------
        integer            :: vegetype
        integer            :: soiltype
        real               :: albbck
        real               :: tbot
        real               :: snoalb
        !-------------------------------------------------------------------------
        ! soil parameters
        !-------------------------------------------------------------------------
        real :: qmax
        real :: qmin
        real :: psis
        real :: ksat
        real :: bclh
        real :: qwrtz
        real :: wilt
        real :: ref
        !-------------------------------------------------------------------------
        ! multilevel spatial parameter
        !-------------------------------------------------------------------------
        real, pointer      :: albedo_monthly(:)
        real, pointer      :: shdfac_monthly(:)
        real, pointer      :: z0brd_monthly(:)
        real, pointer      :: lai_monthly(:)
        !-------------------------------------------------------------------------
        ! state
        !-------------------------------------------------------------------------
        real               :: emiss
        real               :: ch
        real               :: cm
        real               :: shdfac 
        real               :: snthresh 
        real               :: sneqv
        real               :: snowh
        real               :: snowc
        real               :: canwat
        real               :: alb
        real, pointer      :: smc(:)
        real, pointer      :: sho(:)
        real, pointer      :: stc(:)
        real, pointer      :: smfr(:)
        real, pointer      :: keepfr(:)
        real               :: tskin
        real               :: qvg
        real               :: qsfc
        real               :: qcg
        real               :: qsg
        real               :: snt75cm
        real               :: tsnav
        real               :: soilm
        real               :: smroot
        !-------------------------------------------------------------------------
        ! output
        !-------------------------------------------------------------------------
        real               :: qh
        real               :: qle
        real               :: eta
        real               :: qg
        real               :: rhosnf
        real               :: precipfr
        real               :: snowfallac
        real               :: acsnow
        real               :: sfcevp
        real               :: snomlt
        real               :: dew
        real               :: drip
        real               :: qs
        real               :: qsb
        real               :: snflx
        real               :: edir
        real               :: ec
        real               :: ett
        real               :: esnow
        real               :: shdmin 
        real               :: shdmax 
    end type ruc37dec
end module RUC37_module
