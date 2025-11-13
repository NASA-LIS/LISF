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
! !ROUTINE: noah39_transform_snip
! \label{noah39_transform_snip}
!
! !REVISION HISTORY:
!  17 Jul 2025: Eric Kemp; Initial specification (copied from USAFSI version)
!
! !INTERFACE:
subroutine noah39_transform_snip(n, OBS_State)

! !USES:
  use ESMF
  use LIS_logMod,  only : LIS_verify
  use noah39_lsmMod

  ! Defaults
  implicit none

! Arguments
  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine transforms the SNIP state
!  (meters) to the lsm state
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!EOP

  ! Since SNIP is already in meters, no work is needed here.

end subroutine noah39_transform_snip
