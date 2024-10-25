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
! !ROUTINE: clsmf25_updatesoilm
!  \label{clsmf25_updatesoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine clsmf25_updatesoilm(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use clsmf25_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP

  type(ESMF_Field)       :: catdefField,catdefIncrField
  type(ESMF_Field)       :: rzexcField,rzexcIncrField
  type(ESMF_Field)       :: srfexcField,srfexcIncrField
  integer                :: t
  integer                :: status
  real, pointer          :: catdef(:),catdefincr(:)
  real, pointer          :: rzexc(:),rzexcincr(:)
  real, pointer          :: srfexc(:),srfexcincr(:)

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


  call ESMF_StateGet(LSM_Incr_State,"Catchment deficit",&
       catdefincrField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for catdef in clsmf25_getsoilm')
  call ESMF_StateGet(LSM_Incr_State,"Root zone excess",&
       rzexcincrField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for rzexc in clsmf25_getsoilm')
  call ESMF_StateGet(LSM_Incr_State,"Surface excess",srfexcincrField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for srfexc in clsmf25_getsoilm')

  call ESMF_FieldGet(catdefincrField,localDE=0,farrayPtr=catdefincr,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for catdef in clsmf25_getsoilm')
  call ESMF_FieldGet(rzexcincrField,localDE=0,farrayPtr=rzexcincr,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for rzexc in clsmf25_getsoilm')
  call ESMF_FieldGet(srfexcincrField,localDE=0,farrayPtr=srfexcincr,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for srfexc in clsmf25_getsoilm')

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     if(t.eq.5641) print*, 'upd ',catdef(t), catdefincr(t)
     catdef(t)  = catdef(t) + catdefincr(t)
     rzexc(t)   = rzexc(t)  + rzexcincr(t)
     srfexc(t)  = srfexc(t) + srfexcincr(t)
  enddo
end subroutine clsmf25_updatesoilm

