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
! !ROUTINE: RUC37_qc_soilmobs
! \label{RUC37_qc_soilmobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
!
! !INTERFACE:
subroutine RUC37_qc_soilmobs(n,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use RUC37_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
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

  call ESMF_StateGet(OBS_State,"Observation01",obs_sm_field,&
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in RUC37_qc_soilmobs")
  call ESMF_FieldGet(obs_sm_field,localDE=0,farrayPtr=smobs,rc=status)
  call LIS_verify(status,& 
       "ESMF_FieldGet failed in RUC37_qc_soilmobs")

  do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid  = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
     lat =  LIS_domain(n)%grid(gid)%lat
     lon =  LIS_domain(n)%grid(gid)%lon

!     if(lat.gt.39.and.lon.gt.-110.and.lon.lt.-95) then 
!        smobs(gid) = LIS_rc%udef
!     endif

!     if(LIS_rc%mo.eq.9.or.LIS_rc%mo.eq.10.or.&
!          LIS_rc%mo.eq.11.or.LIS_rc%mo.eq.12 &
!          .or.LIS_rc%mo.eq.1.or.LIS_rc%mo.eq.2.or.LIS_rc%mo.eq.3) then 
!        smobs(gid) = LIS_rc%udef
!     endif
!     print*, t, RUC37_struc(n)%ruc37(t)%smc(1), RUC37_struc(n)%ruc37(t)%sho(1)
!     if(abs(RUC37_struc(n)%ruc37(t)%smc(1)-&
!          RUC37_struc(n)%ruc37(t)%sho(1)).gt.0.01) then 
!        smobs(gid) = LIS_rc%udef
!     endif
!     if(RUC37_struc(n)%ruc37(t)%smc(2).ne.&
!          RUC37_struc(n)%ruc37(t)%sho(2)) then 
!        smobs(gid) = LIS_rc%udef
!     endif   
     if(smobs(gid).ne.LIS_rc%udef) then 
        if(RUC37_struc(n)%ruc37(t)%rainf.gt.3E-6) then 
           smobs(gid) = LIS_rc%udef
!           print*, 'rainf ',gid,t,RUC37_struc(n)%ruc37(t)%rainf
        elseif(abs(RUC37_struc(n)%ruc37(t)%smc(1)- &
             RUC37_struc(n)%ruc37(t)%sho(1)).gt.0.0001) then
           smobs(gid) = LIS_rc%udef
!           print*, 'smc1.ne.sho1 ',gid,t,RUC37_struc(n)%ruc37(t)%smc(1),&
!                RUC37_struc(n)%ruc37(t)%sho(1)
        elseif(abs(RUC37_struc(n)%ruc37(t)%smc(2)- &
             RUC37_struc(n)%ruc37(t)%sho(2)).gt.0.0001) then
           smobs(gid) = LIS_rc%udef
!           print*, 'smc2.ne.sho2 ',gid,t,RUC37_struc(n)%ruc37(t)%smc(2),&
!                RUC37_struc(n)%ruc37(t)%sho(2)
        elseif(abs(RUC37_struc(n)%ruc37(t)%smc(3)-&
             RUC37_struc(n)%ruc37(t)%sho(3)).gt.0.0001) then
           smobs(gid) = LIS_rc%udef
!           print*, 'smc3.ne.sho3 ',gid,t,RUC37_struc(n)%ruc37(t)%smc(3),&
!                RUC37_struc(n)%ruc37(t)%sho(3)
        elseif(abs(RUC37_struc(n)%ruc37(t)%smc(4)-&
             RUC37_struc(n)%ruc37(t)%sho(4)).gt.0.0001) then
           smobs(gid) = LIS_rc%udef
        elseif(abs(RUC37_struc(n)%ruc37(t)%smc(5)-&
             RUC37_struc(n)%ruc37(t)%sho(5)).gt.0.0001) then
           smobs(gid) = LIS_rc%udef
        elseif(abs(RUC37_struc(n)%ruc37(t)%smc(6)-&
             RUC37_struc(n)%ruc37(t)%sho(6)).gt.0.0001) then
           smobs(gid) = LIS_rc%udef
        elseif(abs(RUC37_struc(n)%ruc37(t)%smc(7)-&
             RUC37_struc(n)%ruc37(t)%sho(7)).gt.0.0001) then
           smobs(gid) = LIS_rc%udef
        elseif(abs(RUC37_struc(n)%ruc37(t)%smc(8)-&
             RUC37_struc(n)%ruc37(t)%sho(8)).gt.0.0001) then
           smobs(gid) = LIS_rc%udef

!          print*, 'smc4.ne.sho4 ',gid,t,RUC37_struc(n)%ruc37(t)%smc(4),&
 !               RUC37_struc(n)%ruc37(t)%sho(4)
        elseif(RUC37_struc(n)%ruc37(t)%stc(1).le.LIS_CONST_TKFRZ) then
           smobs(gid) = LIS_rc%udef
 !          print*, 'stc1 < 0 ',gid,t, RUC37_struc(n)%ruc37(t)%stc(1)
        elseif(RUC37_struc(n)%ruc37(t)%stc(2).le.LIS_CONST_TKFRZ) then
           smobs(gid) = LIS_rc%udef
 !          print*, 'stc2 < 0 ',gid,t, RUC37_struc(n)%ruc37(t)%stc(2)
        elseif(RUC37_struc(n)%ruc37(t)%stc(3).le.LIS_CONST_TKFRZ) then
           smobs(gid) = LIS_rc%udef
 !          print*, 'stc2 < 0 ',gid,t, RUC37_struc(n)%ruc37(t)%stc(2)
        elseif(RUC37_struc(n)%ruc37(t)%stc(4).le.LIS_CONST_TKFRZ) then
           smobs(gid) = LIS_rc%udef
        elseif(RUC37_struc(n)%ruc37(t)%stc(5).le.LIS_CONST_TKFRZ) then
           smobs(gid) = LIS_rc%udef
        elseif(RUC37_struc(n)%ruc37(t)%stc(6).le.LIS_CONST_TKFRZ) then
           smobs(gid) = LIS_rc%udef
        elseif(RUC37_struc(n)%ruc37(t)%stc(7).le.LIS_CONST_TKFRZ) then
           smobs(gid) = LIS_rc%udef
        elseif(RUC37_struc(n)%ruc37(t)%stc(8).le.LIS_CONST_TKFRZ) then
           smobs(gid) = LIS_rc%udef
 !          print*, 'stc2 < 0 ',gid,t, RUC37_struc(n)%ruc37(t)%stc(2)
        elseif(RUC37_struc(n)%ruc37(t)%tskin.le.LIS_CONST_TKFRZ) then
           smobs(gid) = LIS_rc%udef
 !          print*, 't1 < 0 ',gid,t, RUC37_struc(n)%ruc37(t)%t1
        elseif(RUC37_struc(n)%ruc37(t)%vegetype.le.4) then !forest types
           smobs(gid) = LIS_rc%udef
 !          print*, 'dense veg ',gid,t, RUC37_struc(n)%ruc37(t)%vegt
        elseif(RUC37_struc(n)%ruc37(t)%sneqv.gt.0.001) then 
           smobs(gid) = LIS_rc%udef
 !          print*, 'sneqv>0 ',gid,t, RUC37_struc(n)%ruc37(t)%sneqv
        elseif(RUC37_struc(n)%ruc37(t)%snowc.gt.0.0001) then 
           smobs(gid) = LIS_rc%udef
 !          print*, 'sca>0 ',gid,t, RUC37_struc(n)%ruc37(t)%sca
        elseif(RUC37_struc(n)%ruc37(t)%shdfac.gt.0.5) then 
           smobs(gid) = LIS_rc%udef        
 !          print*, 'shdfac>0.7 ',gid,t, RUC37_struc(n)%ruc37(t)%shdfac
!too close to the tails, could be due to scaling, so reject. 
!        elseif(abs(smobs(gid)-RUC37_struc(n)%ruc37(t)%smcmax).lt.0.05) then 
!           smobs(gid) = LIS_rc%udef
!        elseif(abs(smobs(gid)-RUC37_struc(n)%ruc37(t)%smcwlt).lt.0.05) then 
!           smobs(gid) = LIS_rc%udef
        endif
     endif
  enddo

end subroutine RUC37_qc_soilmobs

