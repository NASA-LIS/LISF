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
! !ROUTINE: HYMAP3_descale_WL
! \label{HYMAP3_descale_WL}
!
! !REVISION HISTORY:
! 07 Nov 2019: Sujay Kumar, Initial specification
!
! !INTERFACE:
subroutine HYMAP3_descale_WL(n, Routing_State, Routing_Incr_State)

! !USES:
  use ESMF
  use HYMAP3_routingMod
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify

  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State), intent(in) :: Routing_State
  type(ESMF_State), intent(in) :: Routing_Incr_State
!
! !DESCRIPTION:
!
!  Descales the water level state prognostic variables for
!  data assimilation (currently empty)
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[Routing\_State] ESMF State container for Routing state variables \newline
!  \end{description}
!EOP


end subroutine HYMAP3_descale_WL

