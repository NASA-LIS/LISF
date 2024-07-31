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
! !ROUTINE: noah271_updateswevars
! \label{noah271_updateswevars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine noah271_updateswevars(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use noah271_lsmMod
  use LIS_logMod,   only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
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
!  \item[LSM\_Incr\_State] ESMF State container for LSM state increments \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)       :: sweField, sweIncrField
  type(ESMF_Field)       :: snodField, snodIncrField

  integer                :: t
  integer                :: status
  real, allocatable          :: swe(:), sweincr(:)
  real, allocatable          :: snod(:), snodincr(:)
 
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State,"SWE",sweIncrField,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweIncrField,localDE=0,farrayPtr=sweincr,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     swe(t) = swe(t) + sweincr(t)
  enddo
end subroutine noah271_updateswevars

