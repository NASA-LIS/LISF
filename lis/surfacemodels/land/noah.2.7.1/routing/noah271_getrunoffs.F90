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
! !ROUTINE: noah271_getrunoffs
!  \label{noah271_getrunoffs}
!
! !REVISION HISTORY:
!  6 May 2011: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah271_getrunoffs(n)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_routingMod, only : LIS_runoff_state
  use LIS_logMod,     only : LIS_verify
  use LIS_historyMod, only : LIS_gather_gridded_output
  use noah271_lsmMod, only : noah271_struc

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
  type(ESMF_Field)       :: sfrunoff_count_field
  type(ESMF_Field)       :: baseflow_field
  real, pointer          :: sfrunoff(:,:)
  real, pointer          :: baseflow(:,:)
  real, pointer          :: runoff_count(:,:)
  real, allocatable          :: gvar1(:,:)
  real, allocatable          :: gvar2(:,:)
  real                   :: runoff1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                   :: runoff2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  integer                :: t
  integer                :: c,r
  integer                :: status

  if(LIS_masterproc) then 
     call ESMF_StateGet(LIS_runoff_state(n),"Surface Runoff",&
          sfrunoff_field,rc=status)
     call LIS_verify(status,'ESMF_StateGet failed for Surface Runoff')
     
     call ESMF_StateGet(LIS_runoff_state(n),"Subsurface Runoff",&
          baseflow_field,rc=status)
     call LIS_verify(status,'ESMF_StateGet failed for Subsurface Runoff')

     call ESMF_StateGet(LIS_runoff_state(n), "Surface Runoff count",&
          sfrunoff_count_field,rc=status)
     call LIS_verify(status, &
                     'ESMF_StateGet failed for Surface Runoff count')
     
     call ESMF_FieldGet(sfrunoff_field,localDE=0,farrayPtr=sfrunoff,rc=status)
     call LIS_verify(status,'ESMF_FieldGet failed for Surface Runoff')
     
     call ESMF_FieldGet(baseflow_field,localDE=0,farrayPtr=baseflow,rc=status)
     call LIS_verify(status,'ESMF_FieldGet failed for Subsurface Runoff')

     call ESMF_FieldGet(sfrunoff_count_field, localDE=0, farrayPtr=&
          runoff_count, rc=status)
     call LIS_verify(status, &
                   'ESMF_FieldGet failed for Subsurface Runoff count')
  endif

  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)  !units?
     runoff1(t) = noah271_struc(n)%noah(t)%runoff1
     runoff2(t) = noah271_struc(n)%noah(t)%runoff2 
  enddo
!gather the model tiles before assigning to the global data structure. 
  call LIS_gather_gridded_output(n, 1, gvar1, runoff1)
  call LIS_gather_gridded_output(n, 1, gvar2, runoff2)

  if(LIS_masterproc) then
     do r=1, LIS_rc%gnr(n)
        do c=1, LIS_rc%gnc(n)
           sfrunoff(c,r) = sfrunoff(c,r) + gvar1(c,r)
           baseflow(c,r) = baseflow(c,r) + gvar2(c,r)
           runoff_count(c,r) = runoff_count(c,r) + 1
        enddo
     enddo

     deallocate(gvar1)
     deallocate(gvar2)
     
  endif
  

!  call gather_gridded_output(n, sfrunoff, runoff1)
!  call gather_gridded_output(n, baseflow, runoff2)

end subroutine noah271_getrunoffs


