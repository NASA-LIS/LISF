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
! !ROUTINE: noahmp36_setalbedovars
! \label{noahmp36_setalbedovars}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp36_setalbedovars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use noahmp36_lsmMod
  use NOAHMP_VEG_PARAMETERS_36
  use LIS_logMod, only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
! 
! !DESCRIPTION:
! 
!  This routine assigns the progognostic variables to noah's
!  model space. The state vector consists of veg
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
     NOAHMP36_struc(n)%noahmp36(t)%albd(:) = albd(t) 
     NOAHMP36_struc(n)%noahmp36(t)%albi(:) = albi(t) 
!     NOAHMP36_struc(n)%noahmp36(t)%albd(1) = albd(t) 
!     NOAHMP36_struc(n)%noahmp36(t)%albd(2) = albi(t) 
!     NOAHMP36_struc(n)%noahmp36(t)%albi(1) = albd(t) 
!     NOAHMP36_struc(n)%noahmp36(t)%albi(2) = albi(t) 
     if(albd(t).ne.-9999.0.and.albi(t).ne.-9999.0) then 
        NOAHMP36_struc(n)%noahmp36(t)%alb_upd_flag = .true.
     else
        NOAHMP36_struc(n)%noahmp36(t)%alb_upd_flag = .false. 
     endif
  enddo

    
end subroutine noahmp36_setalbedovars


