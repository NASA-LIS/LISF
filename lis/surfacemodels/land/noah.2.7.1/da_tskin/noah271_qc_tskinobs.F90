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
! !ROUTINE: noah271_qc_tskinobs
! \label{noah271_qc_tskinobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
!
! !INTERFACE:
subroutine noah271_qc_tskinobs(n,OBS_State)
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
!  prior to data assimilation. Here the tskin observations
!  are flagged when LSM indicates that (1) rain is falling (2)
!  ground is fully or partially covered with snow. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: obs_tskin_field

  real, allocatable            :: tskinobs(:)
  integer                  :: t,i
  integer                  :: N_obs_size
  integer                  :: status

  call ESMF_AttributeGet(OBS_State,name="Number Of Observations",&
       value=N_obs_size,rc=status)
  
  call ESMF_StateGet(OBS_State,"Observation01",obs_tskin_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(obs_tskin_field,localDE=0,farrayPtr=tskinobs,rc=status)
  call LIS_verify(status)
  
  do t=1, LIS_rc%ntiles(n)
     i = LIS_domain(n)%gindex(LIS_domain(n)%tile(t)%col, LIS_domain(n)%tile(t)%row)

     if(noah271_struc(n)%noah(t)%rainf.gt.5.78E-5.and.tskinobs(i).ne.LIS_rc%udef) then 
!        print*, 'flagging for rain ',t, smobs(t), noah271_struc(n)%noah(t)%rainf
        tskinobs(i) = LIS_rc%udef
     elseif(noah271_struc(n)%noah(t)%sca.gt.0.33.and.&
          tskinobs(i).ne.LIS_rc%udef) then 
!        print*, 'flagging for sca ',t, tskinobs(t), noah271_struc(n)%noah(t)%sca
        tskinobs(i) = LIS_rc%udef
     endif
  enddo
  
end subroutine noah271_qc_tskinobs

