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
! !ROUTINE: clsmf25_qc_soilmobs
! \label{clsmf25_qc_soilmobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
!
! !INTERFACE:
subroutine clsmf25_qc_soilmobs(n,k,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use LIS_vegDataMod, only : LIS_lai,LIS_gfrac
  use LIS_DAobservationsMod
  use clsmf25_lsmMod

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
  integer                  :: status
  real                     :: vegt(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: lai(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: swe(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: poros(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: rainf(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: tsurf(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: tp1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: tp2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: tp3(LIS_rc%npatch(n,LIS_rc%lsm_index))

  real                     :: rainf_obs(LIS_rc%obs_ngrid(k))
  real                     :: tp1_obs(LIS_rc%obs_ngrid(k))
  real                     :: tp2_obs(LIS_rc%obs_ngrid(k))
  real                     :: tp3_obs(LIS_rc%obs_ngrid(k))
  real                     :: tsurf_obs(LIS_rc%obs_ngrid(k))
  real                     :: vegt_obs(LIS_rc%obs_ngrid(k))
  real                     :: lai_obs(LIS_rc%obs_ngrid(k))
  real                     :: poros_obs(LIS_rc%obs_ngrid(k))
  real                     :: swe_obs(LIS_rc%obs_ngrid(k))

  
  call ESMF_StateGet(OBS_State,"Observation01",obs_sm_field,&
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in clsmf25_qc_soilmobs")
  call ESMF_FieldGet(obs_sm_field,localDE=0,farrayPtr=smobs,rc=status)
  call LIS_verify(status,& 
       "ESMF_FieldGet failed in clsmf25_qc_soilmobs")

 do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     vegt(t) = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
     swe(t) = clsmf25_struc(n)%cat_progn(t)%wesn(1) + &
          clsmf25_struc(n)%cat_progn(t)%wesn(2) + & 
          clsmf25_struc(n)%cat_progn(t)%wesn(3)

     lai(t) = LIS_lai(n)%tlai(t)
     poros(t) = clsmf25_struc(n)%cat_param(t)%poros
     rainf(t) = clsmf25_struc(n)%met_force(t)%rainf
     tsurf(t) = clsmf25_struc(n)%cat_diagn(t)%tsurf
     tp1(t) = clsmf25_struc(n)%cat_diagn(t)%tp(1)+273.16
     tp2(t) = clsmf25_struc(n)%cat_diagn(t)%tp(2)+273.16
     tp3(t) = clsmf25_struc(n)%cat_diagn(t)%tp(3)+273.16
  enddo
  
  call LIS_convertPatchSpaceToObsSpace(n,k,&       
       LIS_rc%lsm_index, &
       rainf,&
       rainf_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       tp1,&
       tp1_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       tp2,&
       tp2_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       tp3,&
       tp3_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       tsurf,&
       tsurf_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       vegt,&
       vegt_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       lai,&
       lai_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       swe,&
       swe_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       poros,&
       poros_obs)

  do t = 1,LIS_rc%obs_ngrid(k)

!     gid  = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
!     lat =  LIS_domain(n)%grid(gid)%lat
!     lon =  LIS_domain(n)%grid(gid)%lon

     if(smobs(t).ne.LIS_rc%udef) then 
        if(rainf_obs(t).gt.3E-6) then 
           smobs(t) = LIS_rc%udef
        endif
        if(tp1_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        endif
        if(tp2_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        endif
        if(tp3_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        endif
        if(tsurf_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        endif
        if(vegt_obs(t).le.3) then !forest types
           smobs(t) = LIS_rc%udef
        endif
        if(swe_obs(t).gt.0.0001) then 
           smobs(t) = LIS_rc%udef
        endif
        if(lai_obs(t).gt.2) then 
           smobs(t) = LIS_rc%udef        
        endif
        if(abs(smobs(t)- poros_obs(t)).lt.0.05) then 
           smobs(t) = LIS_rc%udef
        endif
     endif

  enddo

end subroutine clsmf25_qc_soilmobs
   
