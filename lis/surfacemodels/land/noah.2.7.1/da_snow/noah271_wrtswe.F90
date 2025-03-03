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
! !ROUTINE: noah271_wrtswe
! \label{noah271_wrtswe}
!
! !REVISION HISTORY:
! 13Nov2007: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah271_wrtswe(ftn, n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use noah271_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: ftn
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the tskin related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP


end subroutine noah271_wrtswe

