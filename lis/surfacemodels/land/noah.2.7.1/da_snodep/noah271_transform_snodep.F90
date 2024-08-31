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
! !ROUTINE: noah271_transform_snodep
! \label{noah271_transform_snodep}
!
! !REVISION HISTORY:
! 25Jun2006: Sujay Kumar: Initial Specification
!
! !INTERFACE:
subroutine noah271_transform_snodep(n,OBS_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use LIS_logMod,  only : LIS_verify
  use noah271_lsmMod
  use SNODEPobs_Mod, only : SNODEP_obs_obj
!EOP
  implicit none

  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine transforms the SNODEP state
!  (mm) to the lsm state
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!EOP
  type(ESMF_Field)         :: obs_snodep_field
  real, pointer            :: snodepobs(:)
  integer                  :: t
  integer                  :: N_obs_size
  integer                  :: status
  
  call ESMF_AttributeGet(OBS_State,name="Number Of Observations",&
       value=N_obs_size,rc=status)
  call LIS_verify(status, 'attributeget error in noah271_transform_snodep')
  call ESMF_StateGet(OBS_State,"Observation01",obs_snodep_field,&
       rc=status)
  call LIS_verify(status,'stateget error in noah271_transform_snodep')
  call ESMF_FieldGet(obs_snodep_field,localDE=0,farrayPtr=snodepobs,rc=status)
  call LIS_verify(status,'fieldget error in noah271_transform_snodep')

  ! If using 8th mesh SNODEP data, convert it from inches to meters.
  ! 16th mesh SNODEP data are already in meters.
  if ( SNODEP_obs_obj(n)%mesh == 8 ) then
     do t=1,N_obs_size
        if(snodepobs(t).ne.LIS_rc%udef) then 
           if(snodepobs(t).gt.408.0) then 
              snodepobs(t) = 0.0        
           endif
           snodepobs(t) = snodepobs(t)*2.54E-2     !inches to meters
        endif
     enddo
  endif
end subroutine noah271_transform_snodep
