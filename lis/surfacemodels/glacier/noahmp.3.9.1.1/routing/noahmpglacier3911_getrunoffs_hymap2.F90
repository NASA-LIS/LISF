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
! !ROUTINE: noahmpglacier3911_getrunoffs_hymap2
!  \label{noahmpglacier3911_getrunoffs_hymap2}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
! !INTERFACE:
subroutine noahmpglacier3911_getrunoffs_hymap2(n)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_routingMod, only : LIS_runoff_state
  use LIS_logMod,     only : LIS_verify
  use LIS_historyMod
  use noahmpglacier3911_Mod, only : noahmpgl3911_struc

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
  integer                :: t
  integer                :: c,r
  integer                :: status
  real, allocatable      :: runoff1(:)
  real, allocatable      :: runoff2(:)
  real, allocatable      :: runoff1_t(:)
  real, allocatable      :: runoff2_t(:)
  
  allocate(runoff1(LIS_rc%npatch(n,LIS_rc%glacier_index)))
  allocate(runoff2(LIS_rc%npatch(n,LIS_rc%glacier_index)))
  allocate(runoff1_t(LIS_rc%ntiles(n)))
  allocate(runoff2_t(LIS_rc%ntiles(n)))

  runoff1_t = -9999.0
  runoff2_t = -9999.0
  
  call ESMF_StateGet(LIS_runoff_state(n),"Surface Runoff",&
       sfrunoff_field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for Surface Runoff')
  
  call ESMF_StateGet(LIS_runoff_state(n),"Subsurface Runoff",&
       baseflow_field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for Subsurface Runoff')

  call ESMF_FieldGet(sfrunoff_field,localDE=0,farrayPtr=sfrunoff,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for Surface Runoff')
  
  call ESMF_FieldGet(baseflow_field,localDE=0,farrayPtr=baseflow,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for Subsurface Runoff')
  

  do t=1, LIS_rc%npatch(n,LIS_rc%glacier_index)  !units?
     runoff1(t) = noahmpgl3911_struc(n)%noahmpgl(t)%runsrf
     runoff2(t) = noahmpgl3911_struc(n)%noahmpgl(t)%runsub
  enddo

  runoff1_t = LIS_rc%udef
  runoff2_t = LIS_rc%udef

  call LIS_patch2tile(n,LIS_rc%glacier_index,runoff1_t, runoff1)
  call LIS_patch2tile(n,LIS_rc%glacier_index,runoff2_t, runoff2)

  do t=1,LIS_rc%ntiles(n)
!if the runoff fields are undefined from the LSM, modify them 
     if(sfrunoff(t).eq.-9999.0) then 
        sfrunoff(t) = runoff1_t(t)
     endif
     if(baseflow(t).eq.-9999.0) then 
        baseflow(t) = runoff2_t(t)
     endif
  enddo

  deallocate(runoff1)
  deallocate(runoff2)
  deallocate(runoff1_t)
  deallocate(runoff2_t)

end subroutine noahmpglacier3911_getrunoffs_hymap2
