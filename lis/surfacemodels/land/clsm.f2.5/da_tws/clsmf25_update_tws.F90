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
! !ROUTINE: clsmf25_update_tws
!  \label{clsmf25_update_tws}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 29Sep2011: Ben Zaitchik: Applied to GRACE
!
! !INTERFACE:
subroutine clsmf25_update_tws(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only  : LIS_verify, LIS_logunit
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

  type(ESMF_Field)       :: catdefField
  type(ESMF_Field)       :: rzexcField
  type(ESMF_Field)       :: srfexcField
  type(ESMF_Field)       :: wesn1Field
  type(ESMF_Field)       :: wesn2Field
  type(ESMF_Field)       :: wesn3Field
  type(ESMF_Field)       :: cdIncrField
  type(ESMF_Field)       :: rzIncrField
  type(ESMF_Field)       :: sfIncrField
  type(ESMF_Field)       :: ws1IncrField
  type(ESMF_Field)       :: ws2IncrField
  type(ESMF_Field)       :: ws3IncrField

  real, pointer          :: catdef(:)
  real, pointer          :: rzexc(:)
  real, pointer          :: srfexc(:)
  real, pointer          :: wesn1(:)
  real, pointer          :: wesn2(:)
  real, pointer          :: wesn3(:)
  real, pointer          :: cdIncr(:)
  real, pointer          :: rzIncr(:)
  real, pointer          :: sfIncr(:)
  real, pointer          :: ws1Incr(:)
  real, pointer          :: ws2Incr(:)
  real, pointer          :: ws3Incr(:)
  integer                :: t
  integer                :: status

  call ESMF_StateGet(LSM_State,"Catchment Deficit",catdefField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Root Zone Excess",rzexcField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Surface Excess",srfexcField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 1",wesn1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 2",wesn2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 3",wesn3Field,rc=status)
  call LIS_verify(status) 

  call ESMF_FieldGet(catdefField,localDE=0,farrayPtr=catdef,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(rzexcField,localDE=0,farrayPtr=rzexc,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(srfexcField,localDE=0,farrayPtr=srfexc,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(wesn1Field,localDE=0,farrayPtr=wesn1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(wesn2Field,localDE=0,farrayPtr=wesn2,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(wesn3Field,localDE=0,farrayPtr=wesn3,rc=status)
  call LIS_verify(status)  

  call ESMF_StateGet(LSM_Incr_State,"Catchment Deficit",cdIncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Root Zone Excess",rzIncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Surface Excess",sfIncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Water Equivalent Snow 1",ws1IncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Water Equivalent Snow 2",ws2IncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Water Equivalent Snow 3",ws3IncrField,rc=status)
  call LIS_verify(status)  

  call ESMF_FieldGet(cdIncrField,localDE=0,farrayPtr=cdIncr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(rzIncrField,localDE=0,farrayPtr=rzIncr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sfIncrField,localDE=0,farrayPtr=sfIncr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(ws1IncrField,localDE=0,farrayPtr=ws1Incr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(ws2IncrField,localDE=0,farrayPtr=ws2Incr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(ws3IncrField,localDE=0,farrayPtr=ws3Incr,rc=status)
  call LIS_verify(status)
  
! if(LIS_localPet.eq.208) then 
!     write(LIS_logunit,*) 'upd ',wesn1(16903),wesn2(16903),wesn3(16903),&
!          ws1Incr(16903),ws2Incr(160903),ws3Incr(16903), &
!          wesn1(16903) + ws1Incr(16903),&
!          wesn2(16903) + ws2Incr(16903),&
!          wesn3(16903) + ws3Incr(16903)
!  endif

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     if(t.ge.1.and.t.le.10) then 
!        write(LIS_logunit,*) 'upd ', cdIncr(t), rzIncr(t), sfIncr(t)
!     endif

     catdef(t) = catdef(t) + cdIncr(t) 
     rzexc(t) = rzexc(t) + rzIncr(t) 
     srfexc(t) = srfexc(t) + sfIncr(t)
     wesn1(t) = wesn1(t) + ws1Incr(t)
     wesn2(t) = wesn2(t) + ws2Incr(t)
     wesn3(t) = wesn3(t) + ws3Incr(t)

  enddo
end subroutine clsmf25_update_tws

