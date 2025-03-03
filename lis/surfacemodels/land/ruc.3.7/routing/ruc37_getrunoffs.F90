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
! !ROUTINE: ruc37_getrunoffs
!  \label{ruc37_getrunoffs}
!
! !INTERFACE:
subroutine ruc37_getrunoffs(n)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_routingMod, only : LIS_runoff_state
  use LIS_logMod,     only : LIS_verify
  use LIS_constantsMod
  use LIS_historyMod
  use ruc37_lsmMod, only : ruc37_struc


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
  real, pointer          :: sfrunoff(:)
  real, pointer          :: baseflow(:)
  real, allocatable          :: gvar1(:)
  real, allocatable          :: gvar2(:)
  real, allocatable          :: runoff1(:)
  real, allocatable          :: runoff2(:)
  real, allocatable          :: runoff1_t(:)
  real, allocatable          :: runoff2_t(:)
  integer                :: t
  integer                :: c,r
  integer                :: status

  allocate(runoff1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(runoff2(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(runoff1_t(LIS_rc%ntiles(n)))
  allocate(runoff2_t(LIS_rc%ntiles(n)))

  if(LIS_masterproc) then 
     call ESMF_StateGet(LIS_runoff_state(n),"Surface Runoff",&
          sfrunoff_field,rc=status)
     call LIS_verify(status, 'ESMF_StateGet failed for Surface Runoff')
     
     call ESMF_StateGet(LIS_runoff_state(n),"Subsurface Runoff",&
          baseflow_field,rc=status)
     call LIS_verify(status, 'ESMF_StateGet failed for Subsurface Runoff')

     call ESMF_FieldGet(sfrunoff_field,localDE=0,farrayPtr=sfrunoff,rc=status)
     call LIS_verify(status,'ESMF_FieldGet failed for Surface Runoff')
     
     call ESMF_FieldGet(baseflow_field,localDE=0,farrayPtr=baseflow,rc=status)
     call LIS_verify(status,'ESMF_FieldGet failed for Subsurface Runoff')

  endif

  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)  !units - from mm s-1 to m/s
     runoff1(t) = RUC37_struc(n)%ruc37(t)%qs/LIS_CONST_RHOFW
     runoff2(t) = RUC37_struc(n)%ruc37(t)%qsb/LIS_CONST_RHOFW
  enddo

  call LIS_patch2tile(n,1,runoff1_t, runoff1)
  call LIS_patch2tile(n,1,runoff2_t, runoff2)

!gather the model tiles before assigning to the global data structure. 
  call LIS_gather_tiled_vector_withhalo_output(n, gvar1, runoff1_t)
  call LIS_gather_tiled_vector_withhalo_output(n, gvar2, runoff2_t)

  if(LIS_masterproc) then

     sfrunoff = gvar1
     baseflow = gvar2

     deallocate(gvar1)
     deallocate(gvar2)

  endif

  deallocate(runoff1)
  deallocate(runoff2)
  deallocate(runoff1_t)
  deallocate(runoff2_t)

end subroutine ruc37_getrunoffs
