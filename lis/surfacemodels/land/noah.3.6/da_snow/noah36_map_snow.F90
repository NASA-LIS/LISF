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
! !ROUTINE: noah36_map_snow
! \label{noah36_map_snow}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!
! !INTERFACE:
subroutine noah36_map_snow(n,OBS_State,LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_surface
  use LIS_constantsMod, only  : LIS_CONST_TKFRZ
  use LIS_logMod,   only  : LIS_logunit, LIS_verify
  use noah36_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
  type(ESMF_State)         :: LSM_Incr_State
! !DESCRIPTION:
!
!  This subroutine directly maps the observation state to the corresponding 
!  variables in the LSM state for SCA data assimilation.
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF State for observations \newline
!  \item[LSM\_State] ESMF State for LSM state variables \newline
!  \end{description}
!
!EOP
  real, parameter          :: nominal_swe = 10.0
  type(ESMF_Field)         :: sweincrField
  type(ESMF_Field)         :: obs_sca_field
  real, pointer            :: sweincr(:)
  type(ESMF_Field)         :: snodincrField
  real, pointer            :: snodincr(:)
  real, pointer            :: scaobs(:)
  integer                  :: t
  integer                  :: status
  integer                  :: obs_state_count
  character*100,allocatable    :: obs_state_objs(:)
  real                     :: tempc, dsnew, hnewc, snowhc
  real                     :: snowh, newsn, newsnc
  real, allocatable            :: swe(:)
  real, allocatable            :: snod(:)

  allocate(swe(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(snod(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  call ESMF_StateGet(LSM_Incr_State,"SWE",sweincrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweincrField,localDE=0,farrayPtr=sweincr,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(LSM_Incr_State,"Snowdepth",snodincrField,rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(snodincrField,localDE=0,farrayPtr=snodincr,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(OBS_State,itemCount=obs_state_count,rc=status)
  call LIS_verify(status)
  allocate(obs_state_objs(obs_state_count))
  
  call ESMF_StateGet(OBS_State,itemNameList=obs_state_objs,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(OBS_State,obs_state_objs(1),obs_sca_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(obs_sca_field,localDE=0,farrayPtr=scaobs,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     swe(t)  = noah36_struc(n)%noah(t)%sneqv
     snod(t) = noah36_struc(n)%noah(t)%snowh
  enddo
! Based on SCA, we update the SWE based on the rule

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
! Observation shows snow and model doesn't
     if(scaobs(LIS_surface(n,1)%tile(t)%index).ne.-9999.0) then 
        if((scaobs(LIS_surface(n,1)%tile(t)%index).gt.0.40).and.&
             (swe(t)*1000.0.le.1.0)) then 
!           print*, 'obs snow, no model snow',t,&
!                scaobs(LIS_surface(n,1)%tile(t)%index),&
!                swe(t)
           swe(t) = swe(t)+nominal_swe/1000.0
           snowh = snod(t)
           newsn = nominal_swe/1000.0
           snowhc = snowh*100.0
           newsnc = newsn*100.0
           if(noah36_struc(n)%noah(t)%tair.lt.LIS_CONST_TKFRZ) then 
              tempc = noah36_struc(n)%noah(t)%tair-LIS_CONST_TKFRZ
           else
              tempc = 0.0
           endif
           if(tempc .le. -15) then 
              dsnew = 0.05
           else
              dsnew = 0.05 + 0.0017 * (tempc + 15.0 )**1.5
           endif
           hnewc = newsnc/dsnew
           snowhc = snowhc + hnewc
           snowh = snowhc *0.01
           snod(t) = snowh
           if(swe(t).le.0.0) snod(t) = 0.0
           if(snod(t).lt.swe(t)) then 
              snod(t) = snod(t)
           endif
!Model shows snow and observation doesn't. 
! remove the snow only if the model has a thin snow pack. 
        elseif((scaobs(LIS_surface(n,1)%tile(t)%index).lt.0.10).and.&
          (swe(t)*1000.0.gt.0.0)) then 
!        write(LIS_logunit,*) 'model snow, no obs snow',t, &
!             scaobs(LIS_surface(n,1)%tile(t)%index),&
!             swe(t)
           if(swe(t)*1000.0.lt.20.0) then 
              swe(t) = 0.0
              snod(t) = 0.0
           endif
        endif
     endif
  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     sweincr(t)  = swe(t)  - noah36_struc(n)%noah(t)%sneqv
     snodincr(t) = snod(t) - noah36_struc(n)%noah(t)%snowh
  enddo 
  deallocate(obs_state_objs)
  deallocate(swe)
  deallocate(snod)

end subroutine noah36_map_snow
   
