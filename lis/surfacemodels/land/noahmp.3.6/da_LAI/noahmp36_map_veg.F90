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
! !ROUTINE: noahmp36_map_veg
! \label{noahmp36_map_veg}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp36_map_veg(n,k,OBS_State,LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_surface
  use LIS_constantsMod, only  : LIS_CONST_TKFRZ
  use LIS_logMod,   only  : LIS_logunit, LIS_verify
  use noahmp36_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
  type(ESMF_State)         :: LSM_Incr_State
! !DESCRIPTION:
!
!  This subroutine directly maps the observation state to the corresponding 
!  variables in the LSM state for veg data assimilation.
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF State for observations \newline
!  \item[LSM\_State] ESMF State for LSM state variables \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: laiIncrField
  type(ESMF_Field)         :: obs_LAI_field
  real, pointer            :: laiincr(:)
  real, pointer            :: LAIobs(:)
  integer                  :: t
  integer                  :: status
  integer                  :: obs_state_count
  character*100,allocatable    :: obs_state_objs(:)
  real, allocatable            :: noahmp36_lai(:)

  allocate(noahmp36_lai(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  call ESMF_StateGet(LSM_Incr_State,"LAI",laiIncrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(laiIncrField,localDE=0,farrayPtr=laiincr,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(OBS_State,itemCount=obs_state_count,rc=status)
  call LIS_verify(status)
  allocate(obs_state_objs(obs_state_count))
  
  call ESMF_StateGet(OBS_State,itemNameList=obs_state_objs,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(OBS_State,obs_state_objs(1),obs_LAI_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(obs_LAI_field,localDE=0,farrayPtr=LAIobs,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     noahmp36_lai(t)  = noahmp36_struc(n)%noahmp36(t)%lai
  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(LAIobs(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index).ge.0) then 
        laiincr(t) = LAIobs(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index)&
             - noahmp36_lai(t)
     else
        laiincr(t) = 0 
     endif
  enddo
!  stop
  deallocate(obs_state_objs)
  deallocate(noahmp36_lai)

end subroutine noahmp36_map_veg
   
