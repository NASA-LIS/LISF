!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module FLake1_module
!BOP
!
! !MODULE: FLake1_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the FLake1 1-d variables.
!  The variables specified in the data structure include:
!
!  \begin{description}      
!   \item[swdown]
!     incoming solar radiation at the surface. unit: W/m2
!   \item[lwdown]
!     incoming longwave radiation at the surface. unit: W/m2
!   \item[wind\_e]
!     eastward wind speed at height\_wind. unit: m/s
!   \item[wind\_n]
!     northward wind speed at height\_wind. unit: m/s
!   \item[tair]
!     air temperature at height\_tq. unit: K
!   \item[qair]
!     specific air humidity at height\_tq. unit: kg/kg
!   \item[psurf]
!     surface air pressure. unit: Pa
!   \item[height\_wind]
!     height where wind speed is measured. unit: m
!   \item[height\_tq]
!     height where temperature and humidity are measured. unit: -
!   \item[flake\_dt]
!     model time step. unit: s
!   \item[lon]
!     longitude of lake center. unit: -
!   \item[lat]
!     latitude of lake center. unit: -
!   \item[depth\_w]
!     lake depth. unit: m
!   \item[fetch]
!     typical wind fetch. unit: m
!   \item[depth\_bs]
!     depth of the thermally active layer of the bottom sediments. unit: m
!   \item[T\_bs]
!     temperature at the outer edge of the thermally active layer of the bottom sediments. unit: K
!   \item[T\_snow]
!     temperature at the air-snow interface. unit: K
!   \item[T\_ice]
!     temperature at the snow-ice interface. unit: K
!   \item[T\_mnw]
!     mean temperature of the water column. unit: K
!   \item[T\_wML]
!     temperature of mixed layer. unit: K
!   \item[T\_bot]
!     temperature at the water-bottom sediment interface. unit: K
!   \item[T\_b1]
!     temperature at the bottom of the upper layer of the sediments. unit: K
!   \item[C\_T]
!     thermocline shape factor. unit: -
!   \item[H\_snow]
!     snow thickness. unit: m
!   \item[H\_ice]
!     ice thickness. unit: m
!   \item[H\_ML]
!     thickness of mixed layer. unit: m
!   \item[H\_B1]
!     thickness of the upper layer of bottom sediments. unit: m
!   \item[T\_sfc]
!     surface temperature. unit: K
!   \item[albedo\_water]
!     water surface albedo with resect to solar radiation. unit: -
!   \item[albedo\_ice]
!     ice surface albedo with respect to the solar radiation. unit: -
!   \item[albedo\_snow]
!     snow surface albedo with respect to the solar radiation. unit: -
!   \item[ufr\_a]
!     friction velocity in air. unit: m/s
!   \item[ufr\_w]
!     friction velocity in surface water. unit: m/s
!   \item[Wconv]
!     convective velocity scale. unit: m/s
!   \item[Q\_se]
!     sensible surface heat flux. unit: W/m2
!   \item[Q\_la]
!     latent surface heat flux. unit: W/m2
!   \item[I\_w]
!     radiation flux through the ice-water or air-water interface. unit: W/m2
!   \item[Q\_lwa]
!     longwave radiation flux from atmosphere. unit: W/m2
!   \item[Q\_lww]
!     longwave radiation flux from water. unit: W/m2
!   \item[Q\_bot]
!     heat flux across water-sediments boundary. unit: W/m2
!   \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  6/4/13: Shugong Wang Initial implementation for LIS 7 and FLake1
!
!EOP
    implicit none
    private
    type, public :: flake1dec
        !-------------------------------------------------------------------------
        ! forcing
        !-------------------------------------------------------------------------
        real :: swdown
        real :: lwdown
        real :: wind_e
        real :: wind_n
        real :: tair
        real :: qair
        real :: psurf
        !-------------------------------------------------------------------------
        ! spatial parameter
        !-------------------------------------------------------------------------
        real :: lon
        real :: lat
        real :: depth_w
        real :: fetch
        real :: depth_bs
        real :: T_bs
        !-------------------------------------------------------------------------
        ! lookup parameter
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        ! state
        !-------------------------------------------------------------------------
        real :: T_snow
        real :: T_ice
        real :: T_mnw
        real :: T_wML
        real :: T_bot
        real :: T_b1
        real :: C_T
        real :: H_snow
        real :: H_ice
        real :: H_ML
        real :: H_B1
        real :: T_sfc
        real :: albedo_water
        real :: albedo_ice
        real :: albedo_snow
        !-------------------------------------------------------------------------
        ! output
        !-------------------------------------------------------------------------
        real :: ufr_a
        real :: ufr_w
        real :: Wconv
        real :: Q_se
        real :: Q_la
        real :: I_w
        real :: Q_lwa
        real :: Q_lww
        real :: Q_bot
    end type flake1dec
end module FLake1_module
