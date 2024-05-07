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
! !ROUTINE: noah271_qc_soilmobs
! \label{noah271_qc_soilmobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
!
! !INTERFACE:
subroutine noah271_qc_soilmobs(n,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use LIS_logMod,  only : LIS_verify
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use noah271_lsmMod

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

  real, allocatable            :: smobs(:)
  integer                  :: t
  integer                  :: N_obs_size
  integer                  :: status

  call ESMF_AttributeGet(OBS_State,name="Number Of Observations",&
       value=N_obs_size,rc=status)
  
  call ESMF_StateGet(OBS_State,"Observation01",obs_sm_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(obs_sm_field,localDE=0,farrayPtr=smobs,rc=status)
  call LIS_verify(status)
  
  do t=1,N_obs_size
     if(noah271_struc(n)%noah(t)%rainf.gt.3E-6.and.smobs(t).ne.LIS_rc%udef) then 
!        print*, 'flagging for rain ',t, smobs(t), noah271_struc(n)%noah(t)%rainf
        smobs(t) = LIS_rc%udef
     elseif(noah271_struc(n)%noah(t)%stc(1).le.LIS_CONST_TKFRZ.and.&
          smobs(t).ne.LIS_rc%udef) then 
!        print*, 'flagging for stc ',t, smobs(t), noah271_struc(n)%noah(t)%stc(1)
        smobs(t) = LIS_rc%udef
     elseif(noah271_struc(n)%noah(t)%sca.gt.0.001.and.&
          smobs(t).ne.LIS_rc%udef) then 
!        print*, 'flagging for sca ',t, smobs(t), noah271_struc(n)%noah(t)%sca
        smobs(t) = LIS_rc%udef
     endif
  enddo
  
end subroutine noah271_qc_soilmobs

