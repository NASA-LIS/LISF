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
! !ROUTINE: noahmp401_getrunoffs_hymap2
!  \label{noahmp401_getrunoffs_hymap2}
!
! !REVISION HISTORY:
!  6 May 2011: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp401_getrunoffs_hymap2(n)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_routingMod, only : LIS_runoff_state
  use LIS_logMod
  use LIS_historyMod
  use noahmp401_lsmMod, only : noahmp401_struc

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
  real, allocatable          :: runoff1(:)
  real, allocatable          :: runoff2(:)
  real, allocatable          :: runoff1_t(:)
  real, allocatable          :: runoff2_t(:)
  
  !ag (25Apr2017)
  type(ESMF_Field)       :: evapotranspiration_Field
  real,pointer           :: evapotranspiration(:)
  real, allocatable      :: evapotranspiration1(:)
  real, allocatable      :: evapotranspiration1_t(:)
  integer                :: evapflag

  allocate(runoff1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(runoff2(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(runoff1_t(LIS_rc%ntiles(n)))
  allocate(runoff2_t(LIS_rc%ntiles(n)))

  runoff1_t = -9999.0
  runoff2_t = -9999.0

  call ESMF_AttributeGet(LIS_runoff_state(n),"Routing model evaporation option",&
       evapflag, rc=status)
!if option is not defined, then assume that no evap calculations will be done
  if(status.ne.0)then 
     evapflag = 0
  endif
  
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

  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)  !units?
     runoff1(t) = NOAHMP401_struc(n)%noahmp401(t)%runsf
     runoff2(t) = NOAHMP401_struc(n)%noahmp401(t)%runsb
  enddo

  runoff1_t = LIS_rc%udef
  runoff2_t = LIS_rc%udef

  call LIS_patch2tile(n,1,runoff1_t, runoff1)
  call LIS_patch2tile(n,1,runoff2_t, runoff2)

  sfrunoff = runoff1_t
  baseflow = runoff2_t

  deallocate(runoff1)
  deallocate(runoff2)
  deallocate(runoff1_t)
  deallocate(runoff2_t)
 
  !ag (05Jun2017)
  !Including meteorological forcings + evapotranspiration for computing evaporation from open waters in HyMAP2)
  if(evapflag.ne.0)then
    allocate(evapotranspiration1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
    allocate(evapotranspiration1_t(LIS_rc%ntiles(n)))
    
    call ESMF_StateGet(LIS_runoff_state(n),"Total Evapotranspiration",&
         evapotranspiration_Field, rc=status)
    call LIS_verify(status, "noahmp401_getrunoffs_hymap2: ESMF_StateGet failed for Total Evapotranspiration")
        
    call ESMF_FieldGet(evapotranspiration_Field,localDE=0,&
         farrayPtr=evapotranspiration,rc=status)
    call LIS_verify(status, "noahmp401_getrunoffs_hymap2: ESMF_FieldGet failed for Total Evapotranspiration")
    
    do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)  
       evapotranspiration1(t)  =  &
            NOAHMP401_struc(n)%noahmp401(t)%ecan + &
            NOAHMP401_struc(n)%noahmp401(t)%etran + &
            NOAHMP401_struc(n)%noahmp401(t)%edir 
    enddo

    call LIS_patch2tile(n,1,evapotranspiration1_t, evapotranspiration1)
  
    evapotranspiration = evapotranspiration1_t

    deallocate(evapotranspiration1)
    deallocate(evapotranspiration1_t)
  endif  

end subroutine noahmp401_getrunoffs_hymap2
