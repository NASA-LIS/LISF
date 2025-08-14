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
! !ROUTINE: noah39_scale_snip
! \label{noah39_scale_snip}
!
! !REVISION HISTORY:
!  17 Jul 2025: Eric Kemp; Initial specification (copied from USAFSI version)
!
! !INTERFACE:
subroutine noah39_scale_snip(n, LSM_State)

! !USES:
  use ESMF
  use noah39_lsmMod

  ! Defaults
  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the snow related state prognostic variables for
!  data assimilation
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP

end subroutine noah39_scale_snip

