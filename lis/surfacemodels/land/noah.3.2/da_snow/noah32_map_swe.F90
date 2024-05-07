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
! !ROUTINE: noah32_map_swe
! \label{noah32_map_swe}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!
! !INTERFACE:
subroutine noah32_map_swe(n,OBS_State,LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_constantsMod, only  : LIS_CONST_TKFRZ
  use LIS_logMod,   only  : LIS_logunit, LIS_verify
  use noah32_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
  type(ESMF_State)         :: LSM_State
! !DESCRIPTION:
!
!  This subroutine directly maps the observation state to the corresponding 
!  variables in the LSM state for SWE data assimilation.
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF State for observations \newline
!  \item[LSM\_State] ESMF State for LSM state variables \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: sweField
  type(ESMF_Field)         :: obs_swe_field
  real, allocatable            :: swe(:)
  type(ESMF_Field)         :: snodField
  real, allocatable            :: snod(:)

  real, allocatable            :: sweobs(:)
  integer                  :: t
  integer                  :: status
  integer                  :: obs_state_count
  character*100,allocatable    :: obs_state_objs(:)
  real                     :: newsn,newsnc,snowhc,tempc, dsnew, hnewc

  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(OBS_State,itemCount=obs_state_count,rc=status)
  call LIS_verify(status)
  allocate(obs_state_objs(obs_state_count))
  
  call ESMF_StateGet(OBS_State,itemNameList=obs_state_objs,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(OBS_State,obs_state_objs(1),obs_swe_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(obs_swe_field,localDE=0,farrayPtr=sweobs,rc=status)
  call LIS_verify(status)

! Based on SWE, we manually update snowdepth    

  do t=1,LIS_rc%ntiles(n)
     if(sweobs(LIS_domain(n)%tile(t)%index).gt.0) then 
        swe(t) = sweobs(LIS_domain(n)%tile(t)%index)
        newsn = swe(t)
        snowhc = noah32_struc(n)%noah(t)%snowh*100.0
        newsnc = newsn*100.0
        if(noah32_struc(n)%noah(t)%tair.lt.LIS_CONST_TKFRZ) then 
           tempc = noah32_struc(n)%noah(t)%tair-LIS_CONST_TKFRZ
        else
           tempc = 0.0
        endif
        if(tempc.le.-15.0) then 
           dsnew = 0.05
        else
           dsnew = 0.05+0.0017*(tempc+15)**1.5
        endif
        hnewc = newsnc/dsnew
        snowhc = snowhc+hnewc
        snod(t) = snowhc*0.01
     endif
  enddo
#if 0 
  do t=1,LIS_rc%ntiles(n)
     if(sweobs(LIS_domain(n)%tile(t)%index).gt.0) then 
        swe(t) = sweobs(LIS_domain(n)%tile(t)%index)
        
        rhos = 150.0
        if(swe(t).gt.0.013) then 
           snod_upd = rhos*swe(t)
        else
           snod_upd = 0.013*1000.0/rhos
        endif
        snod(t) = snod_upd
     endif
  enddo
#endif
  deallocate(obs_state_objs)
end subroutine noah32_map_swe
   
