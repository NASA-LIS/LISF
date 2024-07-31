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
! !ROUTINE: noah271_map_sm
! \label{noah271_map_sm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine noah271_map_sm(n,OBS_State,LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_logMod,  only : LIS_verify
  use noah271_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
  type(ESMF_State)         :: LSM_State
! !DESCRIPTION:
!
!  This subroutine directly maps the observation state to the corresponding 
!  variables in the LSM state for soil moisture data assimilation.
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF State for observations \newline
!  \item[LSM\_State] ESMF State for LSM state variables \newline
!  \end{description}
!
!EOP

  type(ESMF_Field)         :: sm1Field
  type(ESMF_Field)         :: obs_sm_field
  real, allocatable            :: soilm1(:)

  real, allocatable            :: smobs(:)
  integer                  :: t
  integer                  :: status
  integer                  :: obs_state_count
  character*100,allocatable    :: obs_state_objs(:)

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status)
  
!change this to ... retrieving data using "observation value"

  call ESMF_StateGet(OBS_State,itemCount=obs_state_count,rc=status)
  call LIS_verify(status)
  allocate(obs_state_objs(obs_state_count))
  
  call ESMF_StateGet(OBS_State,itemNameList=obs_state_objs,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(OBS_State,obs_state_objs(1),obs_sm_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(obs_sm_field,localDE=0,farrayPtr=smobs,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     if(smobs(LIS_domain(n)%tile(t)%index).gt.0) then 
        soilm1(t) = smobs(LIS_domain(n)%tile(t)%index)
     endif
  enddo
  deallocate(obs_state_objs)
end subroutine noah271_map_sm
   
