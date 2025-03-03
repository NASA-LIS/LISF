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
! !ROUTINE: clsmf25_wrtswe
! \label{clsmf25_wrtswe}
!
! !REVISION HISTORY:
! 13Nov2007: Sujay Kumar; Initial Specification
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!
! !INTERFACE:
subroutine clsmf25_wrtswe(ftn, n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use clsmf25_lsmMod

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


end subroutine clsmf25_wrtswe

