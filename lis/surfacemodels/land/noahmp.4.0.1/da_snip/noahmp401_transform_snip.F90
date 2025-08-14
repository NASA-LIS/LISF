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
! !ROUTINE: noahmp401_transform_snip
! \label{noahmp401_transform_snip}
!
! !REVISION HISTORY:
! 18 Jul 2025: Eric Kemp; Initial specification (copied from USAFSI version)
!
! !INTERFACE:
subroutine noahmp401_transform_snip(n,OBS_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  use noahmp401_lsmMod

  ! Defaults
  implicit none

  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine transforms the SNIP state
!  (mm) to the lsm state
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!EOP

  ! Since SNIP is already in meters, no work is needed here.

end subroutine noahmp401_transform_snip
