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
! !ROUTINE: noahmp36_getalbedovars
! \label{noahmp36_getalbedovars}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp36_getalbedovars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  use noahmp36_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the related state prognostic variables for
!  veg data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)       :: albdField,albiField
  integer                :: t
  integer                :: status
  real, pointer          :: albd(:), albi(:)
 
  call ESMF_StateGet(LSM_State,"ALBD",albdField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(albdField,localDE=0,farrayPtr=albd,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_State,"ALBI",albiField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(albiField,localDE=0,farrayPtr=albi,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     albd(t) = NOAHMP36_struc(n)%noahmp36(t)%albd(1)
!     albi(t) = NOAHMP36_struc(n)%noahmp36(t)%albd(2)
     albi(t) = NOAHMP36_struc(n)%noahmp36(t)%albi(1)
  enddo
  
end subroutine noahmp36_getalbedovars

