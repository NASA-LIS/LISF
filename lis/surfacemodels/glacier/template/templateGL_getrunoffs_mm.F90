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
! !ROUTINE: templateGL_getrunoffs_mm
!  \label{templateGL_getrunoffs_mm}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
! !INTERFACE:
subroutine templateGL_getrunoffs_mm(n)

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
  real, allocatable      :: gvar1(:)
  real, allocatable      :: gvar2(:)
  integer                :: t
  integer                :: c,r
  integer                :: status
  real, allocatable      :: runoff1(:)
  real, allocatable      :: runoff2(:)
  real, allocatable      :: runoff1_t(:)
  real, allocatable      :: runoff2_t(:)
  
!  real                   :: temp(LIS_rc%lnc(n),LIS_rc%lnr(n))

  allocate(runoff1(LIS_rc%npatch(n,LIS_rc%glacier_index)))
  allocate(runoff2(LIS_rc%npatch(n,LIS_rc%glacier_index)))
  allocate(runoff1_t(LIS_rc%ntiles(n)))
  allocate(runoff2_t(LIS_rc%ntiles(n)))

  runoff1_t = -9999.0
  runoff2_t = -9999.0

  if(LIS_masterproc) then 
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

  endif

  do t=1, LIS_rc%npatch(n,LIS_rc%glacier_index)  !units?
     runoff1(t) = LIS_rc%udef
     runoff2(t) = LIS_rc%udef
  enddo

  call LIS_patch2tile(n,LIS_rc%glacier_index,runoff1_t, runoff1)
  call LIS_patch2tile(n,LIS_rc%glacier_index,runoff2_t, runoff2)

!gather the model tiles before assigning to the global data structure. 
  call LIS_gather_tiled_vector_withhalo_output(n, gvar1, runoff1_t)
  call LIS_gather_tiled_vector_withhalo_output(n, gvar2, runoff2_t)

!  open(100,file='test.bin',form='unformatted')
!  temp = -9999.0
!  do t=1,LIS_rc%ntiles(n)
!     c = LIS_domain(n)%tile(t)%col
!     r = LIS_domain(n)%tile(t)%row
!     temp(c,r) = gvar1(t)
!  enddo

!  write(100) temp
!  temp = -9999.0
!  do t=1,LIS_rc%ntiles(n)
!     c = LIS_domain(n)%tile(t)%col
!     r = LIS_domain(n)%tile(t)%row
!     temp(c,r) = sfrunoff(t)
!  enddo

!  write(100) temp


  if(LIS_masterproc) then

     do t=1,LIS_rc%ntiles(n)
!if the runoff fields are undefined from the LSM, modify them 
        if(sfrunoff(t).eq.-9999.0) then 
           sfrunoff(t) = gvar1(t)
        endif
        if(baseflow(t).eq.-9999.0) then 
           baseflow(t) = gvar2(t)
        endif
     enddo

     deallocate(gvar1)
     deallocate(gvar2)

  endif

!  temp = -9999.0
!  do t=1,LIS_rc%ntiles(n)
!     c = LIS_domain(n)%tile(t)%col
!     r = LIS_domain(n)%tile(t)%row
!     temp(c,r) = sfrunoff(t)
!  enddo

!  write(100) temp
!  stop

  deallocate(runoff1)
  deallocate(runoff2)
  deallocate(runoff1_t)
  deallocate(runoff2_t)
 
end subroutine templateGL_getrunoffs_mm
