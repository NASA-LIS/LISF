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
! !ROUTINE: noah36_qc_soilmobs
! \label{noah36_qc_soilmobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
! 28Aug2017: Mahdi Navari; Updated to take into account the latest developments in the SM DA 
!
! !INTERFACE:
subroutine noah36_qc_soilmobs(n,k,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use LIS_DAobservationsMod
  use noah36_lsmMod
  

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine performs any model-based QC of the observation 
!  prior to data assimilation. Here the soil moisture observations
!  are flagged when LSM indicates that (1) rain is falling (2)
!  soil is frozen or (3) ground is fully or partially covered 
!  with snow. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: obs_sm_field

  real, pointer            :: smobs(:)
  integer                  :: t
  integer                  :: gid
  integer                  :: status
  real                     :: lat,lon
  real                     :: smc1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: smc2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: smc3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: smc4(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sh2o1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sh2o2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sh2o3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sh2o4(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc4(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: vegt(LIS_rc%npatch(n,LIS_rc%lsm_index))

  real                     :: rainf_obs(LIS_rc%obs_ngrid(k))
  real                     :: sneqv_obs(LIS_rc%obs_ngrid(k))
  real                     :: sca_obs(LIS_rc%obs_ngrid(k))
  real                     :: shdfac_obs(LIS_rc%obs_ngrid(k))
  real                     :: t1_obs(LIS_rc%obs_ngrid(k))
  real                     :: smcwlt_obs(LIS_rc%obs_ngrid(k))
  real                     :: smcmax_obs(LIS_rc%obs_ngrid(k))
  real                     :: smc1_obs(LIS_rc%obs_ngrid(k))
  real                     :: smc2_obs(LIS_rc%obs_ngrid(k))
  real                     :: smc3_obs(LIS_rc%obs_ngrid(k))
  real                     :: smc4_obs(LIS_rc%obs_ngrid(k))
  real                     :: sh2o1_obs(LIS_rc%obs_ngrid(k))
  real                     :: sh2o2_obs(LIS_rc%obs_ngrid(k))
  real                     :: sh2o3_obs(LIS_rc%obs_ngrid(k))
  real                     :: sh2o4_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc1_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc2_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc3_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc4_obs(LIS_rc%obs_ngrid(k))
  real                     :: vegt_obs(LIS_rc%obs_ngrid(k))

  call ESMF_StateGet(OBS_State,"Observation01",obs_sm_field,&
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in noah36_qc_soilmobs")
  call ESMF_FieldGet(obs_sm_field,localDE=0,farrayPtr=smobs,rc=status)
  call LIS_verify(status,& 
       "ESMF_FieldGet failed in noah36_qc_soilmobs")

  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     smc1(t) = noah36_struc(n)%noah(t)%smc(1)
     smc2(t) = noah36_struc(n)%noah(t)%smc(2)
     smc3(t) = noah36_struc(n)%noah(t)%smc(3)
     smc4(t) = noah36_struc(n)%noah(t)%smc(4)

     sh2o1(t) = noah36_struc(n)%noah(t)%sh2o(1)
     sh2o2(t) = noah36_struc(n)%noah(t)%sh2o(2)
     sh2o3(t) = noah36_struc(n)%noah(t)%sh2o(3)
     sh2o4(t) = noah36_struc(n)%noah(t)%sh2o(4)

     stc1(t) = noah36_struc(n)%noah(t)%stc(1)
     stc2(t) = noah36_struc(n)%noah(t)%stc(2)
     stc3(t) = noah36_struc(n)%noah(t)%stc(3)
     stc4(t) = noah36_struc(n)%noah(t)%stc(4)

     vegt(t) = noah36_struc(n)%noah(t)%vegt
  enddo


  call LIS_convertPatchSpaceToObsSpace(n,k,&       
       LIS_rc%lsm_index, &
       noah36_struc(n)%noah(:)%rainf,&
       rainf_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       noah36_struc(n)%noah(:)%sneqv,&
       sneqv_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       noah36_struc(n)%noah(:)%sca,&
       sca_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       noah36_struc(n)%noah(:)%shdfac,&
       shdfac_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       noah36_struc(n)%noah(:)%t1,&
       t1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       noah36_struc(n)%noah(:)%smcmax,&
       smcmax_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       noah36_struc(n)%noah(:)%smcwlt,&
       smcwlt_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc1,&
       smc1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc2,&
       smc2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc3,&
       smc3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc4,&
       smc4_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sh2o1,&
       sh2o1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sh2o2,&
       sh2o2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sh2o3,&
       sh2o3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sh2o4,&
       sh2o4_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc1,&
       stc1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc2,&
       stc2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc3,&
       stc3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc4,&
       stc4_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       vegt,&
       vegt_obs)


  do t = 1,LIS_rc%obs_ngrid(k)
     if(smobs(t).ne.LIS_rc%udef) then 
        if(rainf_obs(t).gt.3E-6) then 
           smobs(t) = LIS_rc%udef
        elseif(abs(smc1_obs(t)- &
             sh2o1_obs(t)).gt.0.0001) then
           smobs(t) = LIS_rc%udef
        elseif(abs(smc2_obs(t)- &
             sh2o2_obs(t)).gt.0.0001) then
           smobs(t) = LIS_rc%udef
        elseif(abs(smc3_obs(t)- &
             sh2o3_obs(t)).gt.0.0001) then
           smobs(t) = LIS_rc%udef
        elseif(abs(smc4_obs(t)- &
             sh2o4_obs(t)).gt.0.0001) then
           smobs(t) = LIS_rc%udef
        elseif(stc1_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        elseif(stc2_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        elseif(stc3_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        elseif(stc4_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        elseif(t1_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        elseif(vegt_obs(t).le.4) then !forest types
           smobs(t) = LIS_rc%udef
        elseif(sneqv_obs(t).gt.0.001) then 
           smobs(t) = LIS_rc%udef
        elseif(sca_obs(t).gt.0.0001) then 
           smobs(t) = LIS_rc%udef
        elseif(shdfac_obs(t).gt.0.7) then 
           smobs(t) = LIS_rc%udef        
!too close to the tails, could be due to scaling, so reject. 
        elseif(smcmax_obs(t)-smobs(t).lt.0.02) then 
           smobs(t) = LIS_rc%udef
        elseif(smobs(t) - smcwlt_obs(t).lt.0.02) then 
           smobs(t) = LIS_rc%udef
        endif
     endif
  enddo

end subroutine noah36_qc_soilmobs

