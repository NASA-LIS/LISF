!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module summa1_module
!BOP
!
! !MODULE: summa1_module.F90
!
! !DESCRIPTION:
!  Declare forcing-only option (summa1) variables
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
  use nrtype
  use data_types

  implicit none
  

  type summa1dec

     ! define indices
     integer(i4b)                     :: iVar                       ! index of a model variable 
     integer(i4b)                     :: iStruct                    ! loop through data structures
     integer(i4b)                     :: iGRU
     integer(i4b)                     :: iHRU,jHRU,kHRU             ! index of the hydrologic response unit

!-------------------------------------------------------------------------
! Summa1-Forcing Variables
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

  end type summa1dec

end module summa1_module
