!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_transform_burnarea
! \label{noahmp401_transform_burnarea}
!
! !REVISION HISTORY:
! 24 Jul 2022: Sujay Kumar: Initial Specification
!
! !INTERFACE:
subroutine noahmp401_transform_burnarea(n,OBS_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  use noahmp401_lsmMod

!EOP
  implicit none

  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine transforms the BURNAREA state
!  (mm) to the lsm state
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!EOP

end subroutine noahmp401_transform_burnarea
