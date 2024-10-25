!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module templateOpenWater_module
!BOP
!
! !MODULE: templateOpenWater_module.F90
!
! !DESCRIPTION:
!  Declare forcing-only option (templateOpenWater) variables
!
!  \begin{description}
!   \item[forcing]
!     Array of meteorological forcing
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
!   \end{description}
!
! !REVISION HISTORY:
!
! 21 Jul 2004: Sujay Kumar; Initial Specification
! 23 Oct 2007: Kristi Arsenault;  Added code for LISv5.0
! 
!EOP
  implicit none

  type templateOpenWaterdec

!-------------------------------------------------------------------------
! TemplateOpenWater-Forcing Variables
!-------------------------------------------------------------------------
     integer :: dummy
  end type templateOpenWaterdec

end module templateOpenWater_module
