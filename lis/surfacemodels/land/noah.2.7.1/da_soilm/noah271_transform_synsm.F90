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
! !ROUTINE: noah271_transform_synsm
! \label{noah271_transform_synsm}
!
! !REVISION HISTORY:
! 25Jun2006: Sujay Kumar: Initial Specification
!
! !INTERFACE:
subroutine noah271_transform_synsm(n,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use LIS_logMod,  only : LIS_verify
  use noah271_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine transforms the Synthetic observation state 
!  (mm) to a volumetric soil moisture state. 
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
     if(smobs(t).gt.0) then 
        smobs(t) = smobs(t)/&
             (1000.0*noah271_struc(n)%lyrthk(1))
     endif
  enddo
  
end subroutine noah271_transform_synsm

