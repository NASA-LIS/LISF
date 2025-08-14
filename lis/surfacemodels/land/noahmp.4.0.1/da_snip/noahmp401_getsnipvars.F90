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
! !ROUTINE: noahmp401_getsnipvars
! \label{noahmp401_getsnipvars}
!
! !REVISION HISTORY:
! 18 Jul 2025: Eric Kemp; Initial specification (copied from USAFSI version)
!
! !INTERFACE:
!
subroutine noahmp401_getsnipvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  use noahmp401_lsmMod

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
!
!EOP
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField

  integer                :: t
  integer                :: status
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)

  call ESMF_StateGet(LSM_State, "SWE", sweField, rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State, "Snowdepth", snodField, rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField, localDE=0, farrayPtr=swe, rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField, localDE=0, farrayPtr=snod, rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     swe(t) = noahmp401_struc(n)%noahmp401(t)%sneqv
     snod(t) = noahmp401_struc(n)%noahmp401(t)%snowh
  enddo
end subroutine noahmp401_getsnipvars

