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
! !ROUTINE: noah33_getrunoffs_et
!  \label{noah33_getrunoffs_et}
!
! !REVISION HISTORY:
!  6 May 2011: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah33_getrunoffs_et(n)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_routingMod, only : LIS_runoff_state
  use LIS_logMod,     only : LIS_verify
  use LIS_historyMod
  use LIS_constantsMod
  use noah33_lsmMod, only : noah33_struc

  implicit none
! !ARGUMENTS: 
  integer,  intent(in)   :: n 
!
! !DESCRIPTION:
!  
!
! 
!EOP
  type(ESMF_Field)       :: sfrunoff_field
  type(ESMF_Field)       :: baseflow_field
  type(ESMF_Field)       :: et_field
  real, pointer          :: sfrunoff(:)
  real, pointer          :: baseflow(:)
  real, pointer          :: et(:)
  real, allocatable          :: gvar1(:)
  real, allocatable          :: gvar2(:)
  real, allocatable          :: gvar3(:)
  real, allocatable          :: runoff1(:)
  real, allocatable          :: runoff2(:)
  real, allocatable          :: runoff1_t(:)
  real, allocatable          :: runoff2_t(:)
  real, allocatable          :: et_t(:)
  integer                :: t
  integer                :: c,r
  integer                :: status

  allocate(runoff1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(runoff2(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(et(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(runoff1_t(LIS_rc%ntiles(n)))
  allocate(runoff2_t(LIS_rc%ntiles(n)))
  allocate(et_t(LIS_rc%ntiles(n)))

  if(LIS_masterproc) then 
     call ESMF_StateGet(LIS_runoff_state(n),"Surface Runoff",&
          sfrunoff_field,rc=status)
     call LIS_verify(status,'ESMF_StateGet failed for Surface Runoff')
     
     call ESMF_StateGet(LIS_runoff_state(n),"Subsurface Runoff",&
          baseflow_field,rc=status)
     call LIS_verify(status,'ESMF_StateGet failed for Subsurface Runoff')

     call ESMF_StateGet(LIS_runoff_state(n),"Total Evapotranspiration",&
          et_field,rc=status)
     call LIS_verify(status,'ESMF_StateGet failed for Total Evapotranspiration')

     call ESMF_FieldGet(sfrunoff_field,localDE=0,farrayPtr=sfrunoff,rc=status)
     call LIS_verify(status,'ESMF_FieldGet failed for Surface Runoff')
     
     call ESMF_FieldGet(baseflow_field,localDE=0,farrayPtr=baseflow,rc=status)
     call LIS_verify(status,'ESMF_FieldGet failed for Subsurface Runoff')

     call ESMF_FieldGet(et_field,localDE=0,farrayPtr=et,rc=status)
     call LIS_verify(status,'ESMF_FieldGet failed for Total Evapotranspiration')

  endif

  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)  !units?
     runoff1(t) = noah33_struc(n)%noah(t)%runoff1*LIS_CONST_RHOFW
     runoff2(t) = noah33_struc(n)%noah(t)%runoff2*LIS_CONST_RHOFW
     et(t) = noah33_struc(n)%noah(t)%evap !in kg/m2s
  enddo

  call LIS_patch2tile(n,1,runoff1_t, runoff1)
  call LIS_patch2tile(n,1,runoff2_t, runoff2)
  call LIS_patch2tile(n,1,et_t, et)

!gather the model tiles before assigning to the global data structure. 
!  call LIS_gather_gridded_output(n, 1, gvar1, runoff1)
!  call LIS_gather_gridded_output(n, 1, gvar2, runoff2)
  call LIS_gather_tiled_vector_withhalo_output(n, gvar1, runoff1_t)
  call LIS_gather_tiled_vector_withhalo_output(n, gvar2, runoff2_t)
  call LIS_gather_tiled_vector_withhalo_output(n, gvar3, et_t)

  if(LIS_masterproc) then
     
     sfrunoff = gvar1
     baseflow = gvar2
     et = gvar3

     deallocate(gvar1)
     deallocate(gvar2)
     deallocate(gvar3)
  endif
  

!  call gather_gridded_output(n, sfrunoff, runoff1)
!  call gather_gridded_output(n, baseflow, runoff2)

  deallocate(runoff1)
  deallocate(runoff2)
  deallocate(et)
  deallocate(runoff1_t)
  deallocate(runoff2_t)
  deallocate(et_t)
end subroutine noah33_getrunoffs_et


