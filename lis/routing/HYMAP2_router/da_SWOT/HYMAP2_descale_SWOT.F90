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
! !ROUTINE: HYMAP2_descale_SWOT
! \label{HYMAP2_descale_SWOT}
!
! !REVISION HISTORY:
! 15 Apr 24: Yeosang Yoon; Initial specification;
!                          copied from HYMAP2_descale_WL
!
! !INTERFACE:
subroutine HYMAP2_descale_SWOT(n, Routing_State, Routing_Incr_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use HYMAP2_routingMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: Routing_State
  type(ESMF_State)       :: Routing_Incr_State
!
! !DESCRIPTION:
!
!  Descales the SWOT state prognostic variables for
!  data assimilation (currently empty)
! 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[Routing\_State] ESMF State container for Routing state variables \newline
!  \end{description}
!EOP
 

end subroutine HYMAP2_descale_SWOT
