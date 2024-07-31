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
! !ROUTINE: noah271_transform_scf
! \label{noah271_transform_scf}
!
! !REVISION HISTORY:
! 25Jun2006: Sujay Kumar: Initial Specification
!
! !INTERFACE:
subroutine noah271_transform_scf(n,OBS_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use LIS_logMod,  only : LIS_verify
  use noah271_lsmMod
!EOP
  implicit none

  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine transforms the MODIS observation state 
!  (mm) to the lsm state
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!EOP
  type(ESMF_Field)         :: obs_sca_field

  real, allocatable            :: scaobs(:)
  integer                  :: t
  integer                  :: N_obs_size
  integer                  :: status

  call ESMF_AttributeGet(OBS_State,name="Number Of Observations",&
       value=N_obs_size,rc=status)
  
  call ESMF_StateGet(OBS_State,"Observation01",obs_sca_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(obs_sca_field,localDE=0,farrayPtr=scaobs,rc=status)
  call LIS_verify(status)
  
  do t=1,N_obs_size
     if(scaobs(t).gt.0) then 
        scaobs(t) = scaobs(t)/100.0
     endif
  enddo
end subroutine noah271_transform_scf

