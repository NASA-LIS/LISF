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
! !ROUTINE: clsmf25_setsoilm
!  \label{clsmf25_setsoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine clsmf25_setsoilm(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only  : LIS_verify, LIS_logunit
  use clsmf25_lsmMod
  use clsmf25_model

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP

  type(ESMF_Field)       :: catdefField
  type(ESMF_Field)       :: rzexcField
  type(ESMF_Field)       :: srfexcField
  integer                :: t
  integer                :: status
  real, pointer          :: catdef(:)
  real, pointer          :: rzexc(:)
  real, pointer          :: srfexc(:)

  character*100          :: lsm_state_objs(4)

  call ESMF_StateGet(LSM_State,"Catchment deficit",catdefField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for catdef in clsmf25_getsoilm')
  call ESMF_StateGet(LSM_State,"Root zone excess",rzexcField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for rzexc in clsmf25_getsoilm')
  call ESMF_StateGet(LSM_State,"Surface excess",srfexcField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for srfexc in clsmf25_getsoilm')

  call ESMF_FieldGet(catdefField,localDE=0,farrayPtr=catdef,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for catdef in clsmf25_getsoilm')
  call ESMF_FieldGet(rzexcField,localDE=0,farrayPtr=rzexc,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for rzexc in clsmf25_getsoilm')
  call ESMF_FieldGet(srfexcField,localDE=0,farrayPtr=srfexc,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for srfexc in clsmf25_getsoilm')

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

!     if(catdef(t).gt.0.1.and.srfexc(t).gt.0.1) then 
     if(catdef(t).gt.0.1) then 
        clsmf25_struc(n)%cat_progn(t)%catdef = catdef(t)
        !     clsmf25_struc(n)%cat_progn(t)%rzexc  = rzexc(t)
        clsmf25_struc(n)%cat_progn(t)%srfexc = srfexc(t)
        !     if(t.eq.5641) print*, 'set ',catdef(t), rzexc(t), srfexc(t)
     endif
  enddo


  CALL CALC_SOIL_MOIST (                       &
       LIS_rc%npatch(n,LIS_rc%lsm_index),      &
       clsmf25_struc(n)%cat_param%vegcls,  &
       clsmf25_struc(n)%cat_param%dzsf,    &
       clsmf25_struc(n)%cat_param%vgwmax,  &
       clsmf25_struc(n)%cat_param%cdcr1,   &
       clsmf25_struc(n)%cat_param%cdcr2,   &
       clsmf25_struc(n)%cat_param%wpwet,   &
       clsmf25_struc(n)%cat_param%poros,   &
       clsmf25_struc(n)%cat_param%psis,    &
       clsmf25_struc(n)%cat_param%bee,     &
       clsmf25_struc(n)%cat_param%ars1,    & 
       clsmf25_struc(n)%cat_param%ars2,    & 
       clsmf25_struc(n)%cat_param%ars3,    & 
       clsmf25_struc(n)%cat_param%ara1,    & 
       clsmf25_struc(n)%cat_param%ara2,    & 
       clsmf25_struc(n)%cat_param%ara3,    & 
       clsmf25_struc(n)%cat_param%ara4,    & 
       clsmf25_struc(n)%cat_param%arw1,    & 
       clsmf25_struc(n)%cat_param%arw2,    & 
       clsmf25_struc(n)%cat_param%arw3,    & 
       clsmf25_struc(n)%cat_param%arw4,    &
       clsmf25_struc(n)%cat_progn%srfexc,  &
       clsmf25_struc(n)%cat_progn%rzexc,   &
       clsmf25_struc(n)%cat_progn%catdef,  &
       clsmf25_struc(n)%cat_diagn%sfmc,    & 
       clsmf25_struc(n)%cat_diagn%rzmc,    & 
       clsmf25_struc(n)%cat_diagn%prmc)

end subroutine clsmf25_setsoilm

