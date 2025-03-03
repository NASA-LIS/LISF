!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: mos_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the 
!  data structure containing the mosaic 1-d variables. The 
!  variables specified in the data structure include: 
!
!  \begin{description}
!   \item[vegt]
!     vegetation type of tile
!   \item[soiltype]
!     soil type of tile
!   \item{vegp}
!     static vegetation parameter values
!   \item{vegip}
!     interpolated monthly parameter values
!   \item{soilp}
!     static soil parameter values
!   \item{lai}
!     leaf area index of the tile
!   \item{green}
!      greeness fraction value of the tile
!   \item{dsai}
!     stem area index of the tile
!   \item{ct}
!     canopy/soil temperature
!   \item{qa}
!     canopy humidity
!   \item{ics}
!     interception canopy storage
!   \item{snow}
!     snow depth
!   \item{sot}
!     deep soil temperature
!   \item{sowet}
!     soil wetness
!   \item{tair}
!     2m air temperature forcing
!   \item{qair}
!     2m specific humidity forcing
!   \item{swdown}
!     downward shortwave forcing
!   \item{lwdown}
!     downward longwave forcing
!   \item{uwind}
!     u-wind component forcing
!   \item{vwind}
!     v-wind component forcing
!   \item{psurf}
!     surface pressure forcing
!   \item{rainf}
!     total rainfall forcing
!   \item{rainf\_c}
!     convective rainfall forcing
!   \item{snowf}
!     total snowfall forcing
!   \item[swnet]
!     net shortwave radiation (W m-2)
!   \item[lwnet]
!     net longwave radiation (W m-2)
!   \item[qle]
!     latent heat flux (W m-2)
!   \item[qh]
!     sensible heat flux (W m-2)
!   \item[qg]
!     ground heat flux (W m-2)
!   \item[snowfall]
!     snowfall (kg m-2 s-1) 
!   \item[pcp]
!     rainfall (kg m-2 s-1) 
!   \item[evap]
!     evapotranspiration (kg m-2 s-1) 
!   \item[qs]
!     surface runoff (kg m-2 s-1) 
!   \item[qsb]
!     subsurface runoff (kg m-2 s-1) 
!   \item[qsm]
!     snow melt (kg m-2 s-1) 
!   \item[avgsurft]
!     average surface temperature (K)
!   \item[albedo]
!     surface albedo (-)
!   \item[swe]
!     snow water equivalent (kg m-2)
!   \item[soilmoist1]
!     soil moisture for layer1 (kg m-2)
!   \item[soilmoist2]
!     soil moisture for layer2 (kg m-2)
!   \item[soilmoist3]
!     soil moisture for layer3 (kg m-2)
!   \item[soilmoist1m]
!     soil moisture for top 1m of soil (kg m-2)
!   \item[soilwet]
!     total column soil wetness (-)
!   \item[soilwetrz]
!     total column soil wetness (-)
!   \item[ecanop]
!     interception evaporation (kg m-2 s-1)
!   \item[tveg]
!     vegetation transpiration (kg m-2 s-1)
!   \item[canopint]
!     total canopy water storage (kg m-2 s-1) 
!   \item[rootmoist]
!     root zone soil moisture (kg m-2)
!   \item[soilm\_prev]
!     soil moisture from the previous model output 
!   \item[swe\_prev]
!     snow water equivalent from the previous model output
!   \item[water1]
!     ??
!   \item[water2]
!     ??
!   \item[water3]
!     ??
!   \item[dtcanal]
!     change in temperature based on analysis
!   \end{description}
!
! !REVISION HISTORY:
!  11 Feb 2002: Jon Gottschalck; Added AVHRR derived variables
!  25 Sep 2007: Sujay Kumar; Upgraded for LIS5.0
!  09 Nov 2007: Chuck Alonge; Removed nslay - unused
!
! !INTERFACE:
module mos_module
!EOP
  implicit none

  PRIVATE
 
  type, public :: mosdec

     integer :: vegt
     integer :: soiltype
     real :: vegp(24)       
     real :: vegip(6)      
     real :: soilp(10)   
     real :: lai       
     real :: green     
     real :: dsai      
     real :: dtcanal
!-------------------------------------------------------------------------
! Mosaic-State Variables
!-------------------------------------------------------------------------
     real :: ct              
     real :: qa              
     real :: ics             
     real :: snow            
     real :: sot             
     real :: sowet(3)      
     
!-------------------------------------------------------------------------
! Mosaic-Forcing Variables
!-------------------------------------------------------------------------
     real :: tair
     real :: qair
     real :: swdown
     real :: lwdown
     real :: uwind
     real :: vwind
     real :: psurf
     real :: rainf
     real :: rainf_c
     real :: snowf

!-----------------------------------------------------------------------
!  Mosaic-Output variables
!-----------------------------------------------------------------------
     real :: swnet
     real :: lwnet
     real :: qle
     real :: qh
     real :: qg
     real :: snowfall
     real :: pcp
     real :: evap
     real :: qs
     real :: qsb
     real :: qsm
     real :: avgsurft
     real :: albedo
     real :: swe
     real :: snowdepth
     real :: snowcover
     real :: snowhf
     real :: subsnow
     real :: soilmoist1
     real :: soilmoist2
     real :: soilmoist3
     real :: soilmoist1m
     real :: soilwet
     real :: soilwetrz
     real :: snod
     real :: snwfrc
     real :: ecanop
     real :: tveg
     real :: esoil 
     real :: rootmoist
     real :: canopint
     real :: acond
     real :: ccond
     real :: obsz
     real :: soilm_prev
     real :: swe_prev
     real :: water1
     real :: water2
     real :: water3

  end type mosdec
end module mos_module












