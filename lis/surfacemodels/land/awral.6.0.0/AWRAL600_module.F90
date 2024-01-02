!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module AWRAL600_module
!BOP
!
! !MODULE: AWRAL600_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the AWRAL600 1-d variables.
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
!     year of the currrent time step. unit: -
!   \item[month]
!     month of the current time step. unit: -
!   \item[day]
!     day of the current time step. unit: -
!   \item[hour]
!     hour of the current time step. unit: -
!   \item[minute]
!     minute of the current time step. unit: -
!   \item[Tair]
!     average air temperature. unit: K
!   \item[Swdown]
!     downward shortwave radiation. unit: W/m2
!   \item[Rainf]
!     daily gross precipitation. unit: kg/m2s
!   \item[Qair]
!     actual vapour pressure. unit: kg/kg
!   \item[Wind\_E]
!     2m wind magnitude. unit: m/s
!   \item[Swdirect]
!     expected downwelling shortwave radiation on a cloudless day. unit: W/m2
!   \item[e0]
!     potential evaporation. unit: mm/d
!   \item[etot]
!     actual evapotranspiration. unit: mm/d
!   \item[dd]
!     vertical drainage from the bottom of the deep soil layer. unit: mm
!   \item[s0\_avg]
!     water storage in the surface soil layer. unit: mm
!   \item[ss\_avg]
!     water content of the shallow soil store. unit: mm
!   \item[sd\_avg]
!     water content of the deep soil store. unit: mm
!   \item[qtot]
!     total discharge to stream. unit: mm
!   \item[sr]
!     volume of water in the surface water store. unit: mm
!   \item[sg]
!     groundwater storage in the unconfined aquifer. unit: mm
!   \item[s0]
!     water storage in the surface soil layer for each hru. unit: mm
!   \item[ss]
!     water content of the shallow soil store for each hru. unit: mm
!   \item[sd]
!     water content of the deep soil store for each hru. unit: mm
!   \item[mleaf]
!     leaf biomass. unit: kg/m2
!   \item[slope\_coeff]
!     scaling factor for slope. unit: -
!   \item[pair]
!     air pressure. unit: Pa
!   \item[kr\_coeff]
!     scaling factor for ratio of saturated hydraulic conductivity. unit: -
!   \item[k\_rout]
!     rate coefficient controlling discharge to stream. unit: -
!   \item[kssat]
!     saturated hydraulic conductivity of shallow soil layer. unit: mm/d
!   \item[prefr]
!     reference value for precipitation. unit: mm
!   \item[s0max]
!     maximum storage of the surface soil layer. unit: mm
!   \item[slope]
!     slope of the land surface. unit: %
!   \item[ssmax]
!     maximum storage of the shallow soil layer. unit: mm
!   \item[k\_gw]
!     groundwater drainage coefficient. unit: 1/d
!   \item[kr\_sd]
!     routing delay factor for the deep layer. unit: -
!   \item[kr\_0s]
!     routing delay factor for the surface layer. unit: -
!   \item[k0sat]
!     saturated hydraulic conductivity of surface soil layer. unit: mm/d
!   \item[sdmax]
!     maximum storage of the deep soil layer. unit: mm
!   \item[kdsat]
!     saturated hydraulic conductivity of shallow soil layer. unit: mm/d
!   \item[ne]
!     effective porosity. unit: -
!   \item[height]
!     elevation of a point on the hypsometric curve. unit: m
!   \item[hypsperc]
!     hypsometric curve distribution percentile bins. unit: %
!   \item[alb\_dry]
!     dry soil albedo for each hru. unit: -
!   \item[alb\_wet]
!     wet soil albedo for each hru. unit: -
!   \item[cgsmax]
!     coefficient relating vegetation photosynthetic capacity to maximum stomatal conductance for each hru. unit: m/s
!   \item[er\_frac\_ref]
!     specific ratio of the mean evaporation rate and the mean rainfall intensity during storms for each hru. unit: -
!   \item[fsoilemax]
!     soil evaporation scaling factor corresponding to unlimited soil water supply for each hru. unit: -
!   \item[lairef]
!     reference leaf area index (at which fv = 0.63) for each hru. unit: -
!   \item[rd]
!     rooting depth for each hru. unit: m
!   \item[s\_sls]
!     specific canopy rainfall storage per unit leaf area for each hru. unit: mm
!   \item[sla]
!     specific leaf area for each hru. unit: m2/kg
!   \item[tgrow]
!     characteristic time scale for vegetation growth towards equilibrium for each hru. unit: d
!   \item[tsenc]
!     characteristic time scale for vegetation senescence towards equilibrium for each hru. unit: d
!   \item[ud0]
!     maximum possible root water uptake from the deep soil store for each hru. unit: mm/d
!   \item[us0]
!     maximum possible root water uptake from the shallow soil store for each hru. unit: mm/d
!   \item[vc]
!     vegetation photosynthetic capacity index per unit canopy cover for each hru. unit: -
!   \item[w0lime]
!     limiting the value of the relative soil moisture content of the top soil layer at which evaporation is reduced for each hru. unit: -
!   \item[w0ref\_alb]
!     Reference value of w0 that determines the rate of albedo decrease with wetness for each hru. unit: -
!   \item[wdlimu]
!     water-limiting relative water content of the deep soil store for each hru. unit: -
!   \item[wslimu]
!     water-limiting relative water content of the shallow soil store for each hru. unit: -
!   \item[fhru]
!     fraction of the cell which contains shallow and deep rooted vegetation. unit: -
!   \item[hveg]
!     vegetation height for each hru. unit: -
!   \item[laimax]
!     leaf area index max for each hru. unit: -
!   \item[timesteps]
!     number of daily timesteps. unit: -
!   \item[cells]
!     number of grid cells. unit: -
!   \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  12/18/18: Wendy Sharples, Shugong Wang Initial implementation for LIS 7 and AWRAL600
!
!EOP
    implicit none
    private
    type, public :: awral600dec
        !-------------------------------------------------------------------------
        ! forcing
        !-------------------------------------------------------------------------
        real               :: Tair
        real               :: Swdown
        real               :: Rainf
        real               :: Qair
        real               :: Wind_E
        real               :: Swdirect
        !-------------------------------------------------------------------------
        ! spatial parameter
        !-------------------------------------------------------------------------
        real               :: k_rout
        real               :: kssat
        real               :: prefr
        real               :: s0max
        real               :: slope
        real               :: ssmax
        real               :: k_gw
        real               :: kr_sd
        real               :: kr_0s
        real               :: k0sat
        real               :: sdmax
        real               :: kdsat
        real               :: ne
        real               :: hypsperc
        !-------------------------------------------------------------------------
        ! multilevel spatial parameter
        !-------------------------------------------------------------------------
        real, allocatable      :: height(:)
        real, allocatable      :: fhru(:)
        real, allocatable      :: hveg(:)
        real, allocatable      :: laimax(:)
        !-------------------------------------------------------------------------
        ! state
        !-------------------------------------------------------------------------
        real                   :: sr
        real                   :: sg
        real, allocatable      :: s0(:)
        real, allocatable      :: ss(:)
        real, allocatable      :: sd(:)
        real, allocatable      :: mleaf(:)
        !-------------------------------------------------------------------------
        ! output
        !-------------------------------------------------------------------------
        real               :: e0
        real               :: etot
        real               :: dd
        real               :: s0_avg
        real               :: ss_avg
        real               :: sd_avg
        real               :: qtot
    end type awral600dec
end module AWRAL600_module
